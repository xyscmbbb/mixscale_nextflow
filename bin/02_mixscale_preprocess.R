#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Mixscale)
})

option_list <- list(
  make_option("--obj_rds", type = "character"),
  make_option("--target_gene_col", type = "character", default = "target_gene"),
  make_option("--nt_label", type = "character", default = "ONE_INTERGENIC_SITE"),
  make_option("--guide_col", type = "character", default = "pair_key"),
  make_option("--ndims", type = "integer", default = 30),
  make_option("--num_neighbors", type = "integer", default = 20),
  make_option("--min_de_genes", type = "integer", default = 0),
  make_option("--max_de_genes", type = "integer", default = 100),
  make_option("--logfc_threshold", type = "double", default = 0),
  make_option("--out_rds", type = "character", default = "mixscale_obj.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$obj_rds), !is.null(opt$out_rds))

message("[02] Reading Seurat object: ", opt$obj_rds)
obj <- readRDS(opt$obj_rds)

if (!opt$target_gene_col %in% colnames(obj@meta.data)) {
  stop("target_gene_col '", opt$target_gene_col, "' not found in obj@meta.data.")
}

fine_mode <- FALSE
if (opt$guide_col %in% colnames(obj@meta.data)) {
  fine_mode <- TRUE
} else {
  message("[02][WARN] guide_col '", opt$guide_col, "' not found. RunMixscale will use fine.mode = FALSE.")
}

message("[02] NormalizeData / FindVariableFeatures / ScaleData / RunPCA")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

message("[02] CalcPerturbSig")
obj <- CalcPerturbSig(
  object = obj,
  assay = "RNA",
  slot = "data",
  gd.class = opt$target_gene_col,
  nt.cell.class = opt$nt_label,
  reduction = "pca",
  ndims = opt$ndims,
  num.neighbors = opt$num_neighbors,
  new.assay.name = "PRTB",
  split.by = NULL
)

message("[02] RunMixscale")
if (fine_mode) {
  obj <- RunMixscale(
    object = obj,
    assay = "PRTB",
    slot = "scale.data",
    labels = opt$target_gene_col,
    nt.class.name = opt$nt_label,
    min.de.genes = opt$min_de_genes,
    logfc.threshold = opt$logfc_threshold,
    de.assay = "RNA",
    max.de.genes = opt$max_de_genes,
    new.class.name = "mixscale_score",
    fine.mode = TRUE,
    fine.mode.labels = opt$guide_col,
    verbose = TRUE
  )
} else {
  obj <- RunMixscale(
    object = obj,
    assay = "PRTB",
    slot = "scale.data",
    labels = opt$target_gene_col,
    nt.class.name = opt$nt_label,
    min.de.genes = opt$min_de_genes,
    logfc.threshold = opt$logfc_threshold,
    de.assay = "RNA",
    max.de.genes = opt$max_de_genes,
    new.class.name = "mixscale_score",
    fine.mode = FALSE,
    verbose = TRUE
  )
}

saveRDS(obj, opt$out_rds)
message("[02] Wrote: ", opt$out_rds)
