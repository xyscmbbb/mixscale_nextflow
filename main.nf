nextflow.enable.dsl = 2

params.h5ad              = null
params.pair_csv          = null
params.perturb_gene      = null
params.outdir            = "results"
params.container         = "docker.io/xyscmbbb/r-mixscale:1.1.7"

params.target_gene_col   = "target_gene"
params.cell_col          = "cell"
params.guide_col         = "pair_key"
params.nt_label          = "ONE_INTERGENIC_SITE"

params.max_features      = 60609
params.subset_to_gene    = true
params.ndims             = 30
params.num_neighbors     = 20
params.min_de_genes      = 0
params.max_de_genes      = 100
params.logfc_threshold   = 0
params.min_pct           = 0
params.min_cells_group   = 10
params.seed              = 1

if (params.h5ad == null) {
    error "Missing required parameter: --h5ad /path/to/input.h5ad"
}
if (params.pair_csv == null) {
    error "Missing required parameter: --pair_csv /path/to/per_cell_target.csv"
}
if (params.perturb_gene == null) {
    error "Missing required parameter: --perturb_gene ZNHIT6"
}

process CONVERT_H5AD_TO_OBJ {
    tag "${perturb_gene}"
    container params.container
    debug true

    cpus 1
    memory '40 GB'
    time '12h'

    input:
    path h5ad
    path pair_csv
    val perturb_gene

    output:
    path "pert_${perturb_gene}.seurat_obj.rds", emit: seurat_obj

    script:
    """
    set -euo pipefail

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
    echo "============================================================"

    Rscript ${projectDir}/bin/01_convert_h5ad_to_obj.R \
      --h5ad ${h5ad} \
      --pair_csv ${pair_csv} \
      --perturb_gene ${perturb_gene} \
      --target_gene_col ${params.target_gene_col} \
      --cell_col ${params.cell_col} \
      --guide_col ${params.guide_col} \
      --nt_label ${params.nt_label} \
      --max_features ${params.max_features} \
      --subset_to_gene ${params.subset_to_gene} \
      --out_rds pert_${perturb_gene}.seurat_obj.rds \
      2>&1 | tee step1_convert_h5ad.log

    echo "============================================================"
    echo "[STEP 1/3] CONVERT_H5AD_TO_OBJ DONE: ${perturb_gene}"
    echo "[STEP 1/3] time: \$(date)"
    ls -lh pert_${perturb_gene}.seurat_obj.rds
    echo "============================================================"
    """
}

process MIXSCALE_PREPROCESS {
    tag "${perturb_gene}"
    container params.container
    debug true

    cpus 1
    memory '40 GB'
    time '24h'

    input:
    path seurat_obj
    val perturb_gene

    output:
    path "pert_${perturb_gene}.mixscale_obj.rds", emit: mixscale_obj

    script:
    """
    set -euo pipefail

    echo "============================================================"
    echo "[STEP 2/3] MIXSCALE_PREPROCESS START: ${perturb_gene}"
    echo "[STEP 2/3] time: \$(date)"
    echo "[STEP 2/3] workdir: \$(pwd)"
    echo "[STEP 2/3] input: ${seurat_obj}"
    echo "[STEP 2/3] target_gene_col: ${params.target_gene_col}"
    echo "[STEP 2/3] guide_col: ${params.guide_col}"
    echo "[STEP 2/3] nt_label: ${params.nt_label}"
    echo "============================================================"

    Rscript ${projectDir}/bin/02_mixscale_preprocess.R \
      --obj_rds ${seurat_obj} \
      --target_gene_col ${params.target_gene_col} \
      --nt_label ${params.nt_label} \
      --guide_col ${params.guide_col} \
      --ndims ${params.ndims} \
      --num_neighbors ${params.num_neighbors} \
      --min_de_genes ${params.min_de_genes} \
      --max_de_genes ${params.max_de_genes} \
      --logfc_threshold ${params.logfc_threshold} \
      --out_rds pert_${perturb_gene}.mixscale_obj.rds \
      2>&1 | tee step2_mixscale_preprocess.log

    echo "============================================================"
    echo "[STEP 2/3] MIXSCALE_PREPROCESS DONE: ${perturb_gene}"
    echo "[STEP 2/3] time: \$(date)"
    ls -lh pert_${perturb_gene}.mixscale_obj.rds
    echo "============================================================"
    """
}

process RUN_WMVREGDE_SCALED_DEBUG {
    tag "${perturb_gene}"
    publishDir params.outdir, mode: 'copy', overwrite: true
    container params.container
    debug true

    cpus 1
    memory '40 GB'
    time '24h'

    input:
    path mixscale_obj
    val perturb_gene

    output:
    path "pert_${perturb_gene}_de_res_df.csv"
    path "pert_${perturb_gene}_meta_data.csv"

    script:
    """
    set -euo pipefail

    echo "============================================================"
    echo "[STEP 3/3] RUN_WMVREGDE_SCALED_DEBUG START: ${perturb_gene}"
    echo "[STEP 3/3] time: \$(date)"
    echo "[STEP 3/3] workdir: \$(pwd)"
    echo "[STEP 3/3] input: ${mixscale_obj}"
    echo "============================================================"

    stdbuf -oL -eL Rscript --vanilla ${projectDir}/bin/03_Run_wmvRegDE_scaled_debug.R \
      --obj_rds ${mixscale_obj} \
      --perturb_gene ${perturb_gene} \
      --target_gene_col ${params.target_gene_col} \
      --nt_label ${params.nt_label} \
      --logfc_threshold ${params.logfc_threshold} \
      --min_pct ${params.min_pct} \
      --min_cells_group ${params.min_cells_group} \
      2>&1 | tee step3_wmvregde.log

    echo "============================================================"
    echo "[STEP 3/3] RUN_WMVREGDE_SCALED_DEBUG DONE: ${perturb_gene}"
    echo "[STEP 3/3] time: \$(date)"
    ls -lh pert_${perturb_gene}_de_res_df.csv pert_${perturb_gene}_meta_data.csv
    echo "============================================================"
    """
}

workflow {
    h5ad_ch = Channel.fromPath(params.h5ad, checkIfExists: true)
    pair_csv_ch = Channel.fromPath(params.pair_csv, checkIfExists: true)

    converted = CONVERT_H5AD_TO_OBJ(h5ad_ch, pair_csv_ch, params.perturb_gene)
    prepped   = MIXSCALE_PREPROCESS(converted.seurat_obj, params.perturb_gene)
    RUN_WMVREGDE_SCALED_DEBUG(prepped.mixscale_obj, params.perturb_gene)
}
