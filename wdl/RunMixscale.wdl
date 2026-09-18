version 1.0

# Mixscale Perturb-seq, one perturbation per Terra job.
#
# CHANGES vs the previous version, all load-bearing:
#   1. nextflow_branch  main -> large-scale/stream   (main has none of the scaling work)
#   2. docker_image     1.1.8 -> 1.1.11              (HVG kernels + collapsed solver + gcsfs)
#   3. subsample        true -> false                (REQUIRED by the collapsed solver)
#   4. RETICULATE_PYTHON now exported                (step 1 dies without it)
#   5. memory_gb        64 -> 58                     (measured peak 43.94 GB x 1.3)
#   6. preflight checks moved before the h5ad download, and made fatal
#   7. h5ad_gcs_folder -> h5ad_gcs, and the gsutil cp is GONE  (see below)
#   8. disk_gb          100 -> 60                    (the h5ad no longer lands on disk)
#
# NO MORE PER-TARGET SHARDS.
#
# This job used to be handed pert_<gene>.h5ad: that target's cells plus the whole
# NT pool, pre-cut by an upstream sharding stage. At 3M cells / 305k NT that shard
# is ~99.8% NT cells, and building one per target means storing ~5,000 copies of
# the same NT pool -- ~16 TB, ~$320/month, and weeks of serial writing at anndata's
# ~12.8 MB/s gzip rate.
#
# Now the job is handed the WHOLE cell-line file (filtered_correct_pairs.h5ad) and
# step 1 selects its own cells: target_gene in {target, NT}, optionally AND-ed with
# a keep_cell allowlist column. It reads only those rows' nonzeros out of HDF5, one
# bounded slice per contiguous run of wanted cells. Gated bit-identical -- the
# Seurat object built this way is identical() to the one built from the old shard,
# counts and meta.data both.
#
# Two things make that cheap rather than catastrophic, and BOTH are required:
#
#   * The file must be CSR and sorted by target_gene (NT block first). Then every
#     target is exactly 2 contiguous runs. On the original csc_matrix, unsorted
#     file the same selection is 13,965 scattered runs touching every nonzero in
#     the file -- i.e. reading all 31 GB to get 3.2 GB of cells. Do the sort+CSR
#     conversion once per cell line: tools/repack_h5ad.py in this repo, or, better
#     at 3M cells, upstream in whatever writes filtered_correct_pairs.h5ad.
#
#   * The image must be 1.1.11+, which carries gcsfs. h5py opens the gs:// object
#     through it and issues ranged GETs for just those two runs, so the 31 GB file
#     is never localised -- transfer is ~3.2 GB, the same bytes the old shard cost,
#     and nothing touches local disk. On an older image a gs:// URI fails fast with
#     a message saying so; the preflight below catches it before anything expensive.
#
# This WDL runs Nextflow INSIDE one Terra VM, so all three steps share that VM and it
# must be sized to the worst of them (step 2). Splitting into three Cromwell tasks so
# each could be sized to its own peak saves only ~8% of GB-seconds -- step 2 dominates
# the wall clock anyway -- and two extra VM boots plus image pulls eat most of that
# back. Not worth the rewrite.

workflow RunMixscaleNextflow {
  input {
    String target_gene

    # The WHOLE cell-line h5ad -- one file shared by every target of this cell
    # line, NOT a per-target shard. Passed straight through to Nextflow and read
    # in place over gcsfs; it is never downloaded, so its size does not enter
    # disk_gb. Must be CSR + sorted by target_gene (see the header).
    # Example: gs://bucket/path/CM_rna_filtered_correct_pairs.repacked.h5ad
    String h5ad_gcs

    # Optional per-cell allowlist, AND-ed with the target/NT selection. Names a
    # boolean column in guide_csv; "" disables it. This replaces the
    # `annslicer filter --obs-column keep_cell` stage of the old sharding
    # pipeline -- a QC allowlist is now applied in the job instead of being
    # baked into a shard, so changing it costs nothing.
    String keep_cell_col = ""

    # Terra localizes this CSV as a WDL File.
    File guide_csv

    String target_gene_col = "target_gene"
    String cell_col = "cell"
    String guide_col = "pair_key"
    String nt_label = "ONE_INTERGENIC_SITE"

    # glm_gp option. MUST stay false whenever collapsed = true: the collapsed
    # solver is exact only against an unsubsampled fit, because glm_gp's
    # subsample=TRUE caps the overdispersion MLE at 1000 cells and that cap does
    # not survive the per-design-row rewrite. Setting this true silently changes
    # the answer -- no error, different DEGs.
    Boolean subsample = false

    # Fit step 3 with the collapsed Gamma-Poisson solver. Exact (same design,
    # same optimizer, same objective; only per-cell sums become per-design-row
    # sums) and it removes the dense n_genes x n_cells Mu, which at 305k NT is
    # 77 GB on its own -- the stock path does not fit at all at that scale.
    # This is already the default on large-scale/stream; passed explicitly so
    # the behaviour does not change silently if the branch is repointed.
    Boolean collapsed = true

    # CalcPerturbSig chunk size passed to Nextflow / 02_mixscale_preprocess.R
    Int chunk_cells = 3000

    # Final output GCS directory.
    String output_gcs_dir

    # Public GitHub repo containing mixscale_nextflow.
    # NOTE: large-scale/stream must be PUSHED before this will clone.
    String nextflow_repo = "https://github.com/xyscmbbb/mixscale_nextflow.git"
    String nextflow_branch = "large-scale/stream"

    # Contains R/Mixscale + Python/anndata + Java + Nextflow + gsutil + git.
    # 1.1.11 adds gcsfs, which is what lets h5ad_gcs be read in place instead of
    # downloaded -- REQUIRED for a gs:// h5ad. 1.1.10 compiles hvf_kernels.c AND
    # the collapsed solver at image-build time, exports MIXSCALE_HVF_SO /
    # MIXSCALE_CPP_CACHE / MIXSCALE_COLLAPSED_SRC, and puts the conda toolchain on
    # PATH. 1.1.9 compiles hvf_kernels.c and exports MIXSCALE_HVF_SO. On 1.1.8 the
    # run still works but rebuilds the kernels at job start and logs that it did.
    String docker_image = "docker.io/xyscmbbb/r-mixscale:1.1.11"

    # 8 vCPU is what every number below was measured at. Cores are nearly free
    # up to the GCE extended-memory line (6.5 GB RAM per vCPU = 52 GB at 8),
    # but do not buy cores to escape that line -- it does not work that way.
    Int cpu = 8

    # Measured step-2 container peak at 305,110 cells is 43.94 GB; 58 = that
    # x 1.3 headroom. Steps 1 and 3 peak lower (32.61 and 37.85 GB) but share
    # this VM. Raise to 64 if you want margin for the untested case below.
    #
    # CAVEAT: the 305k benchmark had only 4 DE genes, so CalcPerturbSig's DENSE
    # PRTB matrix (n_de_genes x n_cells x 8B) was 0.01 GB. A realistic target
    # with ~3,815 DE genes would allocate ~9.3 GB there. That stage measured
    # 20.03 GB, so it should land near 30 GB and stay under the 43.94 GB peak
    # that DE selection sets -- but that is arithmetic, not a measurement.
    Int memory_gb = 58

    # Disk is billed at the full rate with NO spot discount, so at these runtimes
    # it is a large fraction of the bill -- 250 GB pd-ssd is ~$0.058/hr, i.e.
    # ~$0.029 per 30-minute job, against ~$0.064 of spot compute.
    #
    # Measured footprint of one 305,110-cell run:
    #   input h5ad .................. 0 GB      (read over gcsfs, never localised)
    #   step-1 seurat_obj.rds .......  3.8 GB
    #   step-2 mixscale_obj.rds ..... 17.2 GB
    #   step-3 CSVs ................. ~5 MB
    #   total ....................... ~21 GB
    #
    # It does not double: publishDir mode "copy" is only on the step-3 process,
    # so the big intermediates live once in work/ and are never copied out, and
    # Nextflow's local executor symlinks staged inputs rather than copying them.
    # The Docker image goes on the boot disk, not on this local-disk.
    #
    # 60 GB is ~2.9x margin at 305k. Dropping the localised shard is what bought
    # the 100 -> 60: disk is billed at the full rate with NO spot discount, so at
    # a ~30 min job this is real money (~$0.017 saved per job, ~$85 across 5,000).
    # At 455k the step-1 RDS alone is 26 GB, so raise this for very deep targets.
    Int disk_gb = 60

    # Valid Terra/GCP local disk types are usually "SSD" or "HDD". HDD is much
    # cheaper but throttles the multi-GB saveRDS/readRDS between steps.
    String disk_type = "SSD"

    # Number of preemptible attempts before falling back to non-preemptible.
    # A ~30 min job has roughly 2-3% preemption exposure per attempt.
    Int preemptible = 2
  }

  call RunOneGene {
    input:
      target_gene = target_gene,
      h5ad_gcs = h5ad_gcs,
      keep_cell_col = keep_cell_col,
      guide_csv = guide_csv,
      target_gene_col = target_gene_col,
      cell_col = cell_col,
      guide_col = guide_col,
      nt_label = nt_label,
      subsample = subsample,
      collapsed = collapsed,
      chunk_cells = chunk_cells,
      output_gcs_dir = output_gcs_dir,
      nextflow_repo = nextflow_repo,
      nextflow_branch = nextflow_branch,
      docker_image = docker_image,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      disk_type = disk_type,
      preemptible = preemptible
  }

  output {
    String de_res_df_gcs = RunOneGene.de_res_df_gcs
    String meta_data_gcs = RunOneGene.meta_data_gcs
    String nextflow_trace_gcs = RunOneGene.nextflow_trace_gcs
  }
}


task RunOneGene {
  input {
    String target_gene
    String h5ad_gcs
    String keep_cell_col

    File guide_csv

    String target_gene_col
    String cell_col
    String guide_col
    String nt_label
    Boolean subsample
    Boolean collapsed
    Int chunk_cells

    String output_gcs_dir

    String nextflow_repo
    String nextflow_branch

    String docker_image

    Int cpu
    Int memory_gb
    Int disk_gb
    String disk_type
    Int preemptible
  }

  command <<<
    set -euo pipefail

    echo "============================================================"
    echo "[WDL] Running Mixscale Nextflow"
    echo "[WDL] target_gene: ~{target_gene}"
    echo "[WDL] h5ad_gcs: ~{h5ad_gcs}"
    echo "[WDL] keep_cell_col: ~{keep_cell_col}"
    echo "[WDL] guide_csv original localized path: ~{guide_csv}"
    echo "[WDL] target_gene_col: ~{target_gene_col}"
    echo "[WDL] cell_col: ~{cell_col}"
    echo "[WDL] guide_col: ~{guide_col}"
    echo "[WDL] nt_label: ~{nt_label}"
    echo "[WDL] subsample: ~{subsample}"
    echo "[WDL] collapsed: ~{collapsed}"
    echo "[WDL] chunk_cells: ~{chunk_cells}"
    echo "[WDL] cpu: ~{cpu}"
    echo "[WDL] memory_gb: ~{memory_gb}"
    echo "[WDL] disk_gb: ~{disk_gb}"
    echo "[WDL] disk_type: ~{disk_type}"
    echo "[WDL] preemptible: ~{preemptible}"
    echo "[WDL] output_gcs_dir: ~{output_gcs_dir}"
    echo "[WDL] docker_image: ~{docker_image}"
    echo "[WDL] nextflow_branch: ~{nextflow_branch}"
    echo "============================================================"

    # ---------------------------------------------------------------
    # reticulate does NOT find the right interpreter on its own.
    # py_discover_config() in this image resolves /opt/miniconda3/bin/python3,
    # which is the BASE conda env and has no anndata -- step 1 then dies with
    # ModuleNotFoundError after the h5ad has already been downloaded. The
    # mixscale env is the one with anndata/h5py/scipy.
    # ---------------------------------------------------------------
    export RETICULATE_PYTHON=/opt/miniconda3/envs/mixscale/bin/python
    echo "[WDL] RETICULATE_PYTHON: ${RETICULATE_PYTHON}"

    # ---------------------------------------------------------------
    # Deliberately NOT setting MIXSCALE_CPP_CACHE here.
    #
    # Step 3 compiles bin/collapsed/*.cpp at run time with Rcpp::sourceCpp.
    # Image 1.1.10 has already done that build and ships the cache at
    # /opt/mixscale_cpp, with the sources it was built from at
    # /opt/mixscale_collapsed (MIXSCALE_COLLAPSED_SRC). sourceCpp keys its cache
    # on the source PATH as well as its contents, which is why the image has to
    # ship the sources too -- load_collapsed.R compiles from that fixed path when
    # the files still md5-match the checked-out branch, and from the branch
    # itself when they do not. Pointing MIXSCALE_CPP_CACHE at the task directory
    # would miss the baked cache and pay the ~15 s build on every one of the
    # 5,000 jobs.
    #
    # On 1.1.9 and older there is no baked cache, and load_collapsed.R compiles
    # into a tempdir -- still correct, just ~15 s slower. It also puts the conda
    # toolchain on PATH itself: this image's conda R invokes its compiler as
    # x86_64-conda-linux-gnu-c++, which lives in /opt/miniconda3/envs/mixscale/bin,
    # a directory that is not on the default PATH. Before that fix step 3 died
    # with "x86_64-conda-linux-gnu-c++: not found" after steps 1 and 2 had
    # already run. It never showed up in benchmarking because every bench script
    # bind-mounts a prebuilt cache from the host; a Terra VM has no such cache.
    # ---------------------------------------------------------------
    echo "[WDL] MIXSCALE_CPP_CACHE: ${MIXSCALE_CPP_CACHE:-<unset, will compile into a tempdir>}"
    echo "[WDL] MIXSCALE_COLLAPSED_SRC: ${MIXSCALE_COLLAPSED_SRC:-<unset>}"

    echo "[WDL] Tool paths"
    which Rscript
    which python
    which java
    which nextflow
    which gsutil
    which git

    echo "[WDL] Versions"
    Rscript --version
    python --version
    java -version
    nextflow -version

    # ---------------------------------------------------------------
    # Preflight. Everything here is cheap and fails in seconds; it runs BEFORE
    # the multi-GB h5ad download so a misconfiguration does not cost a transfer.
    # ---------------------------------------------------------------
    echo "[WDL] Preflight: python env used by reticulate"
    "${RETICULATE_PYTHON}" -c "import anndata, h5py, scipy; print('anndata', anndata.__version__)"

    # A gs:// h5ad is read in place through gcsfs. Without it step 1 stops with a
    # clear message, but only after the repo clone and the solver build -- check it
    # here instead, where it costs seconds. Also confirms anndata.io.read_elem,
    # which is how step 1 pulls obs/var off the same handle.
    if [[ "~{h5ad_gcs}" == gs://* ]]; then
      echo "[WDL] Preflight: gcsfs (required, --h5ad is a gs:// URI)"
      "${RETICULATE_PYTHON}" -c "import gcsfs, fsspec; from anndata.io import read_elem; print('gcsfs', gcsfs.__version__, 'fsspec', fsspec.__version__)" || {
        echo "[ERROR] image ~{docker_image} has no gcsfs."
        echo "[ERROR] use docker.io/xyscmbbb/r-mixscale:1.1.11 or newer, or pass a local h5ad path."
        exit 1
      }
    fi

    echo "[WDL] Preflight: prebuilt HVG kernels in the image"
    if [[ -n "${MIXSCALE_HVF_SO:-}" && -f "${MIXSCALE_HVF_SO}" ]]; then
      echo "[WDL] MIXSCALE_HVF_SO: ${MIXSCALE_HVF_SO}"
      Rscript -e 'dyn.load(Sys.getenv("MIXSCALE_HVF_SO")); stopifnot(is.loaded("rowVar2_cm"), is.loaded("rowVarStd_cm")); cat("[WDL] HVG kernels load OK\n")'
    else
      echo "[WDL] WARN: MIXSCALE_HVF_SO not set or missing -- image is probably older than 1.1.9."
      echo "[WDL] WARN: the run will still work but recompiles the kernels at step-2 start."
    fi

    # -----------------------------
    # Normalize paths
    # -----------------------------
    OUTPUT_GCS_DIR="~{output_gcs_dir}"
    OUTPUT_GCS_DIR="${OUTPUT_GCS_DIR%/}"

    # Passed to Nextflow verbatim. Not localised, so there is no *_LOCAL twin and
    # no ../ prefix -- unlike the guide CSV below, this never becomes a file here.
    H5AD="~{h5ad_gcs}"

    GUIDE_CSV_INPUT="~{guide_csv}"
    GUIDE_CSV_LOCAL="guide_csv.csv"

    RESULTS_DIR="results/~{target_gene}"
    TRACE_FILE="${RESULTS_DIR}/pert_~{target_gene}_nextflow_trace.txt"

    echo "[WDL] H5AD: ${H5AD}"
    echo "[WDL] GUIDE_CSV_INPUT: ${GUIDE_CSV_INPUT}"
    echo "[WDL] GUIDE_CSV_LOCAL: ${GUIDE_CSV_LOCAL}"
    echo "[WDL] RESULTS_DIR: ${RESULTS_DIR}"
    echo "[WDL] TRACE_FILE: ${TRACE_FILE}"

    mkdir -p "${RESULTS_DIR}"

    # -----------------------------
    # Compute Nextflow memory
    # Use the smaller of:
    #   memory_gb * 0.96
    #   memory_gb - 4
    # Rounded down to integer GB.
    # -----------------------------
    MEMORY_GB="~{memory_gb}"

    NF_MEMORY_GB=$(awk -v mem="${MEMORY_GB}" 'BEGIN {
      a = mem * 0.96
      b = mem - 4
      c = (a < b) ? a : b
      if (c < 1) c = 1
      printf "%d", c
    }')

    echo "[WDL] MEMORY_GB: ${MEMORY_GB}"
    echo "[WDL] NF_MEMORY_GB: ${NF_MEMORY_GB}"

    # -----------------------------
    # Clone Nextflow repo FIRST -- a missing branch is the cheapest failure to
    # hit, and hitting it before the h5ad download saves a wasted transfer.
    # -----------------------------
    echo "[WDL] Clone Nextflow repo"
    git clone --branch "~{nextflow_branch}" --depth 1 "~{nextflow_repo}" mixscale_nextflow

    echo "[WDL] Cloned commit:"
    git -C mixscale_nextflow log --oneline -1

    # The steps source siblings out of bin/ via dirname(.self), and step 3
    # reaches into bin/collapsed/. A clone gets all of it; this asserts that a
    # repointed branch has not silently dropped the collapsed solver.
    for f in bin/01_convert_h5ad_to_obj.R \
             bin/02_mixscale_preprocess.R \
             bin/03_Run_wmvRegDE_scaled_debug.R \
             bin/hvf_lowmem.R bin/hvf_kernels.c bin/lognorm_blocked.R \
             bin/mixscale_de_genes_lean.R bin/collapsed/load_collapsed.R; do
      if [[ ! -s "mixscale_nextflow/${f}" ]]; then
        echo "[ERROR] branch ~{nextflow_branch} is missing ${f}"
        echo "[ERROR] this branch does not carry the scaling work; expected large-scale/stream"
        exit 1
      fi
    done
    echo "[WDL] Preflight: bin/ complete"

    # Load the collapsed solver now rather than discovering a broken toolchain
    # 20 minutes in, at step 3. On 1.1.10 this is a cache hit (~0.2 s); on older
    # images it compiles, and the log line below reports which happened.
    echo "[WDL] Preflight: compiling the collapsed Gamma-Poisson solver"
    Rscript -e 'source("mixscale_nextflow/bin/collapsed/load_collapsed.R"); load_collapsed_glm_gp("mixscale_nextflow/bin/collapsed"); stopifnot(exists("fit_gp_weighted_cpp"), exists("gp_se_weighted_cpp")); cat("[WDL] collapsed solver OK\n")'

    # The h5ad is deliberately NOT downloaded. Step 1 opens it where it lives and
    # reads only this target's row runs; localising it would move ~31 GB to read
    # ~3.2 GB of it, and would need the disk back.
    #
    # Cheap existence check so a typo'd URI fails here rather than 30 s into R.
    if [[ "${H5AD}" == gs://* ]]; then
      echo "[WDL] Checking the h5ad exists (not downloading it)"
      gsutil ls -l "${H5AD}"
    else
      ls -lh "${H5AD}"
    fi

    # -----------------------------
    # Explicitly copy Terra-localized guide_csv into the task working directory
    # -----------------------------
    echo "[WDL] Copying localized guide_csv into working directory"
    cp "${GUIDE_CSV_INPUT}" "${GUIDE_CSV_LOCAL}"

    echo "[WDL] Check localized input files"
    ls -lh "${GUIDE_CSV_LOCAL}"

    echo "[WDL] First few lines of guide_csv:"
    head -n 5 "${GUIDE_CSV_LOCAL}" || true

    echo "[WDL] Local files:"
    ls -lh

    cd mixscale_nextflow

    echo "[WDL] Repo contents:"
    ls -lh
    ls -lh bin

    # -----------------------------
    # Nextflow runs inside WDL Docker image r-mixscale.
    # Disable nested containers so it does not pull another Docker image.
    #
    # Consequence: main.nf's per-step env exports still apply (they are in the
    # process script bodies), but process.container is ignored, so the image
    # that matters is the WDL runtime one above.
    # -----------------------------
    cat > terra_override.config <<'EOF'
process {
  executor = 'local'
  container = null
}

docker {
  enabled = false
}

singularity {
  enabled = false
}

podman {
  enabled = false
}
EOF

    echo "[WDL] terra_override.config:"
    cat terra_override.config

    # -----------------------------
    # Run Nextflow
    #
    # Not passed, and deliberately left at the branch defaults:
    #   step2_threads  null -> task.cpus
    #   de_target_nnz  1.2e8   (wilcox row-block size; cutting it 3x moved the
    #                           peak 0.07 GB and cost 257 s -- not a useful knob)
    #   keep_cell_col  passed ONLY when non-empty. `--keep_cell_col ""` does NOT
    #                  give Nextflow an empty string -- it sees a flag with no
    #                  value and sets the param to the STRING "true", and step 1
    #                  then looks for a pair_csv column named `true` and dies.
    #                  Omitting the option leaves the branch default of "".
    #   de_threads     step 3 stays SERIAL. --threads > 1 there is a measured
    #                  regression: 1.8x slower, +44% memory, and numerically
    #                  divergent from the serial path.
    # -----------------------------
    echo "[WDL] Run Nextflow"

    nextflow run main.nf \
      -c terra_override.config \
      --h5ad "${H5AD}" \
      --pair_csv "../${GUIDE_CSV_LOCAL}" \
      --perturb_gene "~{target_gene}" \
      --target_gene_col "~{target_gene_col}" \
      --cell_col "~{cell_col}" \
      --guide_col "~{guide_col}" \
      --nt_label "~{nt_label}" \
      ~{if keep_cell_col != "" then "--keep_cell_col '" + keep_cell_col + "'" else ""} \
      --subsample "~{subsample}" \
      --collapsed "~{collapsed}" \
      --chunk_cells "~{chunk_cells}" \
      --cpus "~{cpu}" \
      --memory "${NF_MEMORY_GB} GB" \
      --outdir "../${RESULTS_DIR}" \
      -ansi-log false \
      -with-trace "../${TRACE_FILE}"

    cd ..

    echo "[WDL] Results generated by Nextflow:"
    find results -type f -print -exec ls -lh {} \;

    DE_FILE="results/~{target_gene}/pert_~{target_gene}_de_res_df.csv"
    META_FILE="results/~{target_gene}/pert_~{target_gene}_meta_data.csv"
    TRACE_FILE="results/~{target_gene}/pert_~{target_gene}_nextflow_trace.txt"

    if [[ ! -s "${DE_FILE}" ]]; then
      echo "[ERROR] Missing or empty DE file: ${DE_FILE}"
      exit 1
    fi

    if [[ ! -s "${META_FILE}" ]]; then
      echo "[ERROR] Missing or empty metadata file: ${META_FILE}"
      exit 1
    fi

    if [[ ! -s "${TRACE_FILE}" ]]; then
      echo "[ERROR] Missing or empty Nextflow trace file: ${TRACE_FILE}"
      echo "[WDL] Searching for any trace files:"
      find . -name "*trace*" -type f -print -exec ls -lh {} \; || true
      exit 1
    fi

    echo "[WDL] Preview Nextflow trace:"
    head -n 20 "${TRACE_FILE}" || true

    # -----------------------------
    # Upload outputs directly to final GCS directory
    # -----------------------------
    echo "[WDL] Upload outputs directly to GCS"

    gsutil -m cp "${DE_FILE}" "${OUTPUT_GCS_DIR}/pert_~{target_gene}_de_res_df.csv"
    gsutil -m cp "${META_FILE}" "${OUTPUT_GCS_DIR}/pert_~{target_gene}_meta_data.csv"
    gsutil -m cp "${TRACE_FILE}" "${OUTPUT_GCS_DIR}/pert_~{target_gene}_nextflow_trace.txt"

    echo "[WDL] Confirm uploaded files:"
    gsutil ls -lh "${OUTPUT_GCS_DIR}/pert_~{target_gene}_de_res_df.csv"
    gsutil ls -lh "${OUTPUT_GCS_DIR}/pert_~{target_gene}_meta_data.csv"
    gsutil ls -lh "${OUTPUT_GCS_DIR}/pert_~{target_gene}_nextflow_trace.txt"

    echo "[WDL] Done"
    echo "[WDL] DE result: ${OUTPUT_GCS_DIR}/pert_~{target_gene}_de_res_df.csv"
    echo "[WDL] metadata: ${OUTPUT_GCS_DIR}/pert_~{target_gene}_meta_data.csv"
    echo "[WDL] Nextflow trace: ${OUTPUT_GCS_DIR}/pert_~{target_gene}_nextflow_trace.txt"
  >>>

  output {
    String de_res_df_gcs = "~{output_gcs_dir}/pert_~{target_gene}_de_res_df.csv"
    String meta_data_gcs = "~{output_gcs_dir}/pert_~{target_gene}_meta_data.csv"
    String nextflow_trace_gcs = "~{output_gcs_dir}/pert_~{target_gene}_nextflow_trace.txt"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} ~{disk_type}"
    preemptible: preemptible
  }
}
