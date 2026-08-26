nextflow.enable.dsl = 2

// ============================================================
// Step 1: Convert h5ad to Seurat object
// ============================================================

process CONVERT_H5AD_TO_OBJ {
    tag "${perturb_gene}"
    debug true

    cpus   params.cpus_step1 ?: params.cpus
    memory params.mem_step1  ?: params.memory

    input:
    path h5ad
    path pair_csv
    val perturb_gene

    output:
    tuple val(perturb_gene), path("pert_${perturb_gene}.seurat_obj.rds")

    script:
    """
    set -euo pipefail

    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export BLIS_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}

    echo "============================================================"
    echo "[STEP 1/3] CONVERT_H5AD_TO_OBJ START: ${perturb_gene}"
    echo "[STEP 1/3] time: \$(date)"
    echo "[STEP 1/3] workdir: \$(pwd)"
    echo "[STEP 1/3] h5ad: ${h5ad}"
    echo "[STEP 1/3] pair_csv: ${pair_csv}"
    echo "[STEP 1/3] target_gene_col: ${params.target_gene_col}"
    echo "[STEP 1/3] cell_col: ${params.cell_col}"
    echo "[STEP 1/3] guide_col: ${params.guide_col}"
    echo "[STEP 1/3] nt_label: ${params.nt_label}"
    echo "[STEP 1/3] task.cpus: ${task.cpus}"
    echo "[STEP 1/3] task.memory: ${task.memory}"
    echo "[STEP 1/3] nproc: \$(nproc)"
    echo "============================================================"

    stdbuf -oL -eL Rscript ${projectDir}/bin/01_convert_h5ad_to_obj.R \\
      --h5ad ${h5ad} \\
      --pair_csv ${pair_csv} \\
      --perturb_gene ${perturb_gene} \\
      --target_gene_col ${params.target_gene_col} \\
      --cell_col ${params.cell_col} \\
      --guide_col ${params.guide_col} \\
      --nt_label ${params.nt_label} \\
      --max_features ${params.max_features} \\
      --subset_to_gene ${params.subset_to_gene} \\
      --out_rds pert_${perturb_gene}.seurat_obj.rds \\
      2>&1 | tee step1_convert_h5ad.log

    echo "============================================================"
    echo "[STEP 1/3] CONVERT_H5AD_TO_OBJ DONE: ${perturb_gene}"
    echo "[STEP 1/3] time: \$(date)"
    ls -lh pert_${perturb_gene}.seurat_obj.rds
    echo "============================================================"
    """
}

// ============================================================
// Step 2: Mixscale preprocessing
// ============================================================

process MIXSCALE_PREPROCESS {
    tag "${perturb_gene}"
    debug true

    cpus   params.cpus_step2 ?: params.cpus
    memory params.mem_step2  ?: params.memory

    input:
    tuple val(perturb_gene), path(seurat_obj)

    output:
    tuple val(perturb_gene), path("pert_${perturb_gene}.mixscale_obj.rds")

    script:
    """
    set -euo pipefail

    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export BLIS_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}

    echo "============================================================"
    echo "[STEP 2/3] MIXSCALE_PREPROCESS START: ${perturb_gene}"
    echo "[STEP 2/3] time: \$(date)"
    echo "[STEP 2/3] workdir: \$(pwd)"
    echo "[STEP 2/3] input: ${seurat_obj}"
    echo "[STEP 2/3] target_gene_col: ${params.target_gene_col}"
    echo "[STEP 2/3] guide_col: ${params.guide_col}"
    echo "[STEP 2/3] nt_label: ${params.nt_label}"
    echo "[STEP 2/3] ndims: ${params.ndims}"
    echo "[STEP 2/3] num_neighbors: ${params.num_neighbors}"
    echo "[STEP 2/3] chunk_cells: ${params.chunk_cells}"
    echo "[STEP 2/3] min_de_genes: ${params.min_de_genes}"
    echo "[STEP 2/3] max_de_genes: ${params.max_de_genes}"
    echo "[STEP 2/3] logfc_threshold: ${params.logfc_threshold}"
    echo "[STEP 2/3] task.cpus: ${task.cpus}"
    echo "[STEP 2/3] task.memory: ${task.memory}"
    echo "[STEP 2/3] nproc: \$(nproc)"
    echo "============================================================"

    stdbuf -oL -eL Rscript ${projectDir}/bin/02_mixscale_preprocess.R \\
      --obj_rds ${seurat_obj} \\
      --target_gene_col ${params.target_gene_col} \\
      --guide_col ${params.guide_col} \\
      --nt_label ${params.nt_label} \\
      --ndims ${params.ndims} \\
      --num_neighbors ${params.num_neighbors} \\
      --chunk_cells ${params.chunk_cells} \\
      --min_de_genes ${params.min_de_genes} \\
      --max_de_genes ${params.max_de_genes} \\
      --logfc_threshold ${params.logfc_threshold} \\
      --out_rds pert_${perturb_gene}.mixscale_obj.rds \\
      2>&1 | tee step2_mixscale_preprocess.log

    echo "============================================================"
    echo "[STEP 2/3] MIXSCALE_PREPROCESS DONE: ${perturb_gene}"
    echo "[STEP 2/3] time: \$(date)"
    ls -lh pert_${perturb_gene}.mixscale_obj.rds
    echo "============================================================"
    """
}

// ============================================================
// Step 3: Weighted Mixscale DE
// ============================================================

process RUN_WMVREGDE {
    tag "${perturb_gene}"
    debug true

    cpus   params.cpus_step3 ?: params.cpus
    memory params.mem_step3  ?: params.memory

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    tuple val(perturb_gene), path(mixscale_obj)

    output:
    path "pert_${perturb_gene}_de_res_df.csv"
    path "pert_${perturb_gene}_meta_data.csv"

    script:
    """
    set -euo pipefail

    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export BLIS_NUM_THREADS=${task.cpus}
    export VECLIB_MAXIMUM_THREADS=${task.cpus}
    export NUMEXPR_NUM_THREADS=${task.cpus}

    echo "============================================================"
    echo "[STEP 3/3] RUN_WMVREGDE START: ${perturb_gene}"
    echo "[STEP 3/3] time: \$(date)"
    echo "[STEP 3/3] workdir: \$(pwd)"
    echo "[STEP 3/3] input: ${mixscale_obj}"
    echo "[STEP 3/3] target_gene_col: ${params.target_gene_col}"
    echo "[STEP 3/3] nt_label: ${params.nt_label}"
    echo "[STEP 3/3] subsample: ${params.subsample}"
    echo "[STEP 3/3] min_pct: ${params.min_pct}"
    echo "[STEP 3/3] collapsed: ${params.collapsed}"
    echo "[STEP 3/3] min_cells_group: ${params.min_cells_group}"
    echo "[STEP 3/3] task.cpus: ${task.cpus}"
    echo "[STEP 3/3] task.memory: ${task.memory}"
    echo "[STEP 3/3] nproc: \$(nproc)"
    echo "[STEP 3/3] OMP_NUM_THREADS: \$OMP_NUM_THREADS"
    echo "[STEP 3/3] OPENBLAS_NUM_THREADS: \$OPENBLAS_NUM_THREADS"
    echo "[STEP 3/3] MKL_NUM_THREADS: \$MKL_NUM_THREADS"
    echo "============================================================"

    stdbuf -oL -eL Rscript ${projectDir}/bin/03_Run_wmvRegDE_scaled_debug.R \\
      --obj_rds ${mixscale_obj} \\
      --perturb_gene ${perturb_gene} \\
      --target_gene_col ${params.target_gene_col} \\
      --nt_label ${params.nt_label} \\
      --subsample ${params.subsample} \\
      --logfc_threshold ${params.logfc_threshold} \\
      --min_pct ${params.min_pct} \\
      --min_cells_group ${params.min_cells_group} \\
      --collapsed ${params.collapsed} \\
      2>&1 | tee step3_wmvregde.log

    echo "============================================================"
    echo "[STEP 3/3] RUN_WMVREGDE DONE: ${perturb_gene}"
    echo "[STEP 3/3] time: \$(date)"
    ls -lh pert_${perturb_gene}_de_res_df.csv pert_${perturb_gene}_meta_data.csv
    echo "============================================================"
    """
}

// ============================================================
// Workflow
// ============================================================

workflow {
    if (params.h5ad == null) {
        error "Please provide --h5ad"
    }

    if (params.pair_csv == null) {
        error "Please provide --pair_csv"
    }

    if (params.perturb_gene == null) {
        error "Please provide --perturb_gene"
    }

    h5ad_ch = Channel.fromPath(params.h5ad, checkIfExists: true)
    pair_csv_ch = Channel.fromPath(params.pair_csv, checkIfExists: true)

    convert_ch = CONVERT_H5AD_TO_OBJ(
        h5ad_ch,
        pair_csv_ch,
        params.perturb_gene
    )

    mixscale_ch = MIXSCALE_PREPROCESS(convert_ch)

    RUN_WMVREGDE(mixscale_ch)
}