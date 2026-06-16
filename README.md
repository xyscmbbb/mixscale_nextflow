# mixscale_nextflow

[Mixscale](https://github.com/satijalab/Mixscale) is an R package for analyzing single-cell Perturb-seq data. It quantifies cell-level perturbation strength and performs weighted differential expression analysis, which can improve statistical power in CRISPRi Perturb-seq datasets. The official Mixscale vignette describes two main steps: calculating Mixscale scores and running scoring-based weighted DE tests.

For large-scale Perturb-seq datasets, Mixscale can require substantial memory and long runtimes, especially for genome-wide screens with thousands of perturbed genes. This Nextflow workflow wraps the Mixscale analysis into a reproducible, cluster-friendly pipeline. The workflow subsets the full dataset to one perturbed target gene plus non-targeting controls, because each target gene can be analyzed independently. This substantially reduces memory requirements and makes it easy to parallelize across target genes.

## Workflow overview

The pipeline runs one target gene at a time in three steps:

1. `01_convert_h5ad_to_obj.R`  
   Converts a sliced `.h5ad` file into a Seurat object RDS.

2. `02_mixscale_preprocess.R`  
   Runs Mixscale preprocessing and computes Mixscale scores.

3. `03_Run_wmvRegDE_scaled_debug.R`  
   Runs weighted Mixscale differential expression analysis and writes final CSV outputs.

## Inputs

The main inputs are:

- `--h5ad`: sliced `.h5ad` file for one perturbation plus non-targeting controls.
- `--pair_csv`: per-cell guide assignment CSV.
- `--perturb_gene`: target gene to analyze.
- `--target_gene_col`: column containing the target gene label.
- `--cell_col`: column containing cell barcodes.
- `--guide_col`: column containing guide or guide-pair labels.
- `--nt_label`: non-targeting control label.
- `--outdir`: output directory.

The `.h5ad` input can be sliced from a larger `.h5ad` file, for example using [`annslicer`](https://github.com/cellarium-ai/annslicer), a low-memory utility for slicing and merging AnnData `.h5ad`/`.zarr` files.

The `pair_csv` should include at least:

| Column | Description |
|---|---|
| `cell_col` | Cell barcode / cell ID |
| `guide_col` | Guide or guide-pair assignment for each cell |
| `target_gene_col` | Target gene assignment for each cell |

## Run example

```bash
nextflow run main.nf -resume \
  --h5ad pert_RRP9.h5ad \
  --pair_csv NGN2-MN_per_cell_target.csv \
  --perturb_gene RRP9 \
  --target_gene_col target_gene \
  --cell_col cell \
  --guide_col pair_key \
  --nt_label ONE_INTERGENIC_SITE \
  --subsample false \
  --min_pct 0.02 \
  --logfc_threshold 0 \
  --min_cells_group 10 \
  --outdir results \
  --cpus 8 \
  --memory 120 \
  --chunk_cells 5000 \
  -ansi-log false \
  -with-trace
```

## Performance and bottlenecks

**Step 3 (weighted DE) dominates the pipeline's runtime, and is also a memory peak.**
It fits one negative-binomial GLM per gene with `glm_gp`, which is **single-threaded**
and, with
`on_disk = FALSE`, **materializes a dense genes × cells working matrix** (there is no
in-memory sparse IRLS path). Cost therefore scales with two axes:

- **Number of genes fitted** — the full transcriptome by default.
- **Number of cells** = perturbed + **all non-targeting controls**. The NT pool is
  usually the larger of the two, so it is a first-class driver of both time and memory.

For reference, on RRP9 the genome-wide step-3 fit takes ~10 min / ~30 GiB at 19k cells,
and ~46 min / ~100 GB at ~69k cells (mostly NT).

**Step 2 (`MIXSCALE_PREPROCESS`) is also memory-heavy**, though shorter. Its
`CalcPerturbSig` builds a perturbation-signature matrix whose peak RAM scales with the
cell count; on the ~69k-cell RRP9 run it reached **~64 GB**. This step is the reason
the pipeline still needs a large-memory machine even when step 3 is filtered. Its peak
is bounded by **`--chunk_cells`**, which processes cells in chunks: **lower
`chunk_cells` lowers peak RAM at the cost of runtime, higher speeds it up but raises
peak** (e.g. `--chunk_cells 5000` gave ~64 GB / ~6 min; a smaller chunk size trades
that for a longer wall time). Tune `chunk_cells` down if step 2 is the memory-limiting
step on your machine.

### What helps: pre-filter the gene set

The dominant lever is **fitting fewer genes**. `--min_pct` and `--logfc_threshold`
prune the gene set fed to `glm_gp`: a gene is fitted only if it is expressed in at
least `min_pct` of either the perturbed or NT cells, passes the `logfc_threshold`, and
is seen in at least `--min_cells_group` cells. On the RRP9 benchmark, raising
`min_pct` to `0.02` cut the fitted set from ~24k to ~900 genes, giving roughly a **12×
speedup and 4× lower memory** on the serial path. Defaults (`min_pct = 0`,
`logfc_threshold = 0`) reproduce the original unfiltered behaviour.

> **Caveat:** `p_adj_bh` is BH-corrected over the number of genes *actually fitted*, so
> adjusted p-values from a filtered run are **not** directly comparable to an unfiltered
> run (the smaller denominator inflates significance, and genes pruned by `min_pct` are
> never tested). For comparable statistics, fit the filtered set but apply BH correction
> over the full gene universe, or report the raw `p_weight` and adjust yourself.

### What does *not* help

- **Multi-threading `glm_gp`.** It has no internal parallelism, and fitting gene
  chunks in forked workers benchmarked as a **regression** — slower and more memory
  (each fork densifies its own chunk simultaneously, plus BLAS oversubscription) and
  numerically divergent from the serial fit. The pipeline is intentionally serial.
- **`--subsample true`.** This only subsamples `glm_gp`'s overdispersion estimation,
  which is not the cost driver, so it gives **no measurable speed or memory benefit**
  while perturbing the weighted p-values (it shifted ~⅓ of DEG calls on the benchmark).
  Leave `--subsample false`.

For very large NT pools, the most direct way to lower the step-3 memory ceiling is to
**cap the number of NT cells** before the fit (smaller dense matrix on both axes),
since neither threading nor subsampling addresses it.
