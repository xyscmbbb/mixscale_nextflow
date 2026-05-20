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
  --outdir results \
  --cpus 8 \
  --memory 120 \
  --chunk_cells 5000 \
  -ansi-log false \
  -with-trace
