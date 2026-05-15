#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(reticulate)
  library(Matrix)
  library(Seurat)
})

option_list <- list(
  make_option("--h5ad", type = "character"),
  make_option("--pair_csv", type = "character"),
  make_option("--perturb_gene", type = "character"),
  make_option("--target_gene_col", type = "character", default = "target_gene"),
  make_option("--cell_col", type = "character", default = "cell"),
  make_option("--guide_col", type = "character", default = "pair_key"),
  make_option("--nt_label", type = "character", default = "ONE_INTERGENIC_SITE"),
  make_option("--max_features", type = "integer", default = 60609),
  make_option("--subset_to_gene", type = "character", default = "true"),
  make_option("--out_rds", type = "character", default = "seurat_obj.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$h5ad), !is.null(opt$pair_csv), !is.null(opt$perturb_gene), !is.null(opt$out_rds))
subset_to_gene <- tolower(opt$subset_to_gene) %in% c("true", "t", "1", "yes", "y")

message("[01] Started at: ", Sys.time())
message("[01] h5ad: ", opt$h5ad)
message("[01] pair_csv: ", opt$pair_csv)
message("[01] perturb_gene: ", opt$perturb_gene)
message("[01] target_gene_col: ", opt$target_gene_col)
message("[01] cell_col: ", opt$cell_col)
message("[01] guide_col: ", opt$guide_col)
message("[01] nt_label: ", opt$nt_label)
message("[01] RETICULATE_PYTHON: ", Sys.getenv("RETICULATE_PYTHON"))

message("[01] reticulate Python config:")
print(reticulate::py_config())

message("[01] Reading h5ad: ", opt$h5ad)
anndata <- import("anndata", convert = FALSE)
scipy <- import("scipy.sparse", convert = FALSE)
ad <- anndata$read_h5ad(opt$h5ad)
X <- ad$X

if (py_to_r(scipy$issparse(X))) {
  message("[01] adata.X is sparse; converting scipy sparse matrix to R sparse Matrix.")
  X_coo <- X$tocoo()
  shape <- py_to_r(X_coo$shape)
  n_row <- as.integer(shape[[1]])
  n_col <- as.integer(shape[[2]])
  i <- as.integer(py_to_r(X_coo$row)) + 1L
  j <- as.integer(py_to_r(X_coo$col)) + 1L
  x <- as.numeric(py_to_r(X_coo$data))
  mat <- sparseMatrix(i = i, j = j, x = x, dims = c(n_row, n_col))
} else {
  message("[01][WARN] adata.X is dense; converting to sparse Matrix in R.")
  mat <- Matrix(py_to_r(X), sparse = TRUE)
}

cell_names <- py_to_r(ad$obs_names$to_list())
gene_names <- py_to_r(ad$var_names$to_list())
rownames(mat) <- cell_names
colnames(mat) <- gene_names
counts <- t(mat)

if (!is.na(opt$max_features) && opt$max_features > 0 && nrow(counts) > opt$max_features) {
  message("[01] Keeping first ", opt$max_features, " features out of ", nrow(counts), ".")
  counts <- counts[seq_len(opt$max_features), , drop = FALSE]
}

rownames(counts) <- gsub("_", "-", rownames(counts))
rownames(counts) <- make.unique(rownames(counts))
colnames(counts) <- make.unique(colnames(counts))

obs <- py_to_r(ad$obs)
rownames(obs) <- make.unique(cell_names)
obs <- obs[colnames(counts), , drop = FALSE]

if (!opt$target_gene_col %in% colnames(obs)) {
  stop("[01] target_gene_col '", opt$target_gene_col, "' not found in h5ad obs. Available columns: ", paste(colnames(obs), collapse = ", "))
}

if (!file.exists(opt$pair_csv)) {
  stop("[01] pair_csv does not exist: ", opt$pair_csv)
}

message("[01] Reading pair CSV: ", opt$pair_csv)
pair_df <- read.csv(opt$pair_csv, row.names = 1, check.names = FALSE)

if (!opt$cell_col %in% colnames(pair_df)) {
  stop("[01] cell_col '", opt$cell_col, "' not found in pair_csv. Available columns: ", paste(colnames(pair_df), collapse = ", "))
}
if (!opt$guide_col %in% colnames(pair_df)) {
  stop("[01] guide_col '", opt$guide_col, "' not found in pair_csv. Available columns: ", paste(colnames(pair_df), collapse = ", "))
}

# Add guide_col from external per-cell CSV into h5ad obs before creating Seurat object.
obs[[opt$guide_col]] <- pair_df[match(rownames(obs), pair_df[[opt$cell_col]]), opt$guide_col]
n_missing_guides <- sum(is.na(obs[[opt$guide_col]]))
message("[01] Added guide column to metadata: ", opt$guide_col)
message("[01] Missing ", opt$guide_col, ": ", n_missing_guides, " / ", nrow(obs))

if (n_missing_guides == nrow(obs)) {
  stop("[01] All guide_col values are NA. Cell names probably do not match between h5ad obs_names and pair_csv cell_col.")
}

if (subset_to_gene) {
  keep_cells <- obs[[opt$target_gene_col]] %in% c(opt$perturb_gene, opt$nt_label)
  message("[01] Subsetting cells to ", opt$perturb_gene, " + ", opt$nt_label, ": ", sum(keep_cells), " / ", nrow(obs))
  if (sum(keep_cells) == 0) stop("[01] No cells remain after subsetting. Check --perturb_gene, --target_gene_col, and --nt_label.")
  counts <- counts[, keep_cells, drop = FALSE]
  obs <- obs[keep_cells, , drop = FALSE]
}

obj <- CreateSeuratObject(counts = counts, meta.data = obs)
message("[01] Seurat object: ", nrow(obj), " genes x ", ncol(obj), " cells")
message("[01] Metadata columns include: ", paste(colnames(obj@meta.data), collapse = ", "))
message("[01] Saving: ", opt$out_rds)
saveRDS(obj, opt$out_rds)
message("[01] Wrote: ", opt$out_rds)
message("[01] Finished at: ", Sys.time())
