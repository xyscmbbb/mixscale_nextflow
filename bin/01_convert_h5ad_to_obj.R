#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(reticulate)
  library(Matrix)
  library(Seurat)
})

.self <- tryCatch(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                  error = function(e) "01_convert_h5ad_to_obj.R")
source(file.path(dirname(.self), "counts_stats.R"))
source(file.path(dirname(.self), "rds_io.R"))

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

# Per-stage wall clock. At 305,110 cells step 1 is 892 s end to end and the
# stages are wildly uneven (an HDF5 slice read, an O(nnz) build, and a gzip of a
# 17 GB object), so a single total says nothing about what to attack.
.t0 <- Sys.time(); .tprev <- .t0
tick01 <- function(what) {
  now <- Sys.time()
  message(sprintf("[01][time] %-28s %7.1f s  (total %7.1f s)", what,
                  as.numeric(difftime(now, .tprev, units = "secs")),
                  as.numeric(difftime(now, .t0, units = "secs"))))
  .tprev <<- now
}
message("[01] Reading h5ad: ", opt$h5ad)
anndata <- import("anndata", convert = FALSE)
np      <- import("numpy",   convert = FALSE)
bi      <- import_builtins(convert = FALSE)

# obs/var only -- backed="r" leaves X on disk so it is never materialised twice.
ad <- anndata$read_h5ad(opt$h5ad, backed = "r")
cell_names <- py_to_r(ad$obs_names$to_list())
gene_names <- py_to_r(ad$var_names$to_list())
obs        <- py_to_r(ad$obs)

# Read one HDF5 dataset into a preallocated R vector in slices, so python never
# holds the whole array alongside the R copy.
h5_read_vec <- function(ds, n, integer_out, chunk = 5e7) {
  out <- if (integer_out) integer(n) else numeric(n)
  a <- 0
  while (a < n) {
    b <- min(a + chunk, n)
    v <- py_to_r(np$asarray(ds$`__getitem__`(bi$slice(bi$int(a), bi$int(b)))))
    out[(a + 1):b] <- v
    a <- b
  }
  out
}

h5py <- import("h5py", convert = FALSE)
fh   <- h5py$File(opt$h5ad, "r")
Xg   <- fh[["X"]]
enc  <- tryCatch(as.character(py_to_r(Xg$attrs[["encoding-type"]])), error = function(e) "")

if (enc %in% c("csr_matrix", "csc_matrix")) {
  # An h5ad stores X as cells x genes. Either sparse encoding is just the triple
  # (indptr, indices, data), which maps straight onto a dgCMatrix's (p, i, x) --
  # so the matrix can be built from the HDF5 datasets with no scipy COO, no
  # dgTMatrix intermediate and no 1-based index vectors.
  #
  #   CSR(cells x genes) is BIT-IDENTICAL to CSC(genes x cells), i.e. exactly the
  #   counts matrix Seurat wants -- the transpose is free.
  #   CSC(cells x genes) gives the cells x genes dgCMatrix, and still needs one
  #   real t(). Exporting shards as CSR avoids that copy entirely.
  shape  <- py_to_r(Xg$attrs[["shape"]])
  n_cell <- as.integer(shape[[1]])
  n_gene <- as.integer(shape[[2]])
  nnz    <- as.numeric(py_to_r(Xg[["indices"]]$shape)[[1]])
  message(sprintf("[01] X is %s, %d cells x %d genes, nnz = %.0f (%.0f/cell)",
                  enc, n_cell, n_gene, nnz, nnz / n_cell))
  if (nnz >= 2147483647)
    stop("[01] nnz = ", format(nnz, scientific = FALSE), " exceeds the dgCMatrix ",
         "limit of 2^31-1. Shard the h5ad into fewer cells.")

  message("[01] reading indices ...")
  i_vec <- h5_read_vec(Xg[["indices"]], nnz, TRUE)
  p_vec <- as.integer(py_to_r(np$asarray(Xg[["indptr"]])))
  message("[01] reading data ...")
  x_vec <- h5_read_vec(Xg[["data"]], nnz, FALSE)
  fh$close()

  # new() shares these vectors rather than copying them; the validity check is
  # what guarantees the stored indices were sorted within each slice.
  dm <- if (enc == "csr_matrix") c(n_gene, n_cell) else c(n_cell, n_gene)
  dn <- if (enc == "csr_matrix") list(gene_names, cell_names)
        else                     list(cell_names, gene_names)
  counts <- new("dgCMatrix", i = i_vec, p = p_vec, x = x_vec, Dim = dm, Dimnames = dn)
  rm(i_vec, p_vec, x_vec); gc(FALSE)

  if (enc == "csc_matrix") {
    message("[01] transposing to genes x cells (one copy; export shards as CSR to skip this)")
    counts <- Matrix::t(counts)
    counts <- as(counts, "CsparseMatrix")
    gc(FALSE)
  }
} else {
  message("[01][WARN] X is not sparse (encoding-type '", enc, "'); ",
          "falling back to the scipy COO path (much heavier).")
  fh$close()
  scipy <- import("scipy.sparse", convert = FALSE)
  ad2 <- anndata$read_h5ad(opt$h5ad)
  X <- ad2$X
  if (py_to_r(scipy$issparse(X))) {
    X_coo <- X$tocoo()
    shape <- py_to_r(X_coo$shape)
    mat <- sparseMatrix(i = as.integer(py_to_r(X_coo$row)) + 1L,
                        j = as.integer(py_to_r(X_coo$col)) + 1L,
                        x = as.numeric(py_to_r(X_coo$data)),
                        dims = c(as.integer(shape[[1]]), as.integer(shape[[2]])))
  } else {
    mat <- Matrix(py_to_r(X), sparse = TRUE)
  }
  rownames(mat) <- cell_names
  colnames(mat) <- gene_names
  counts <- t(mat)
  rm(mat); gc(FALSE)
}

if (!is.na(opt$max_features) && opt$max_features > 0 && nrow(counts) > opt$max_features) {
  message("[01] Keeping first ", opt$max_features, " features out of ", nrow(counts), ".")
  counts <- counts[seq_len(opt$max_features), , drop = FALSE]
}

rownames(counts) <- gsub("_", "-", rownames(counts))
rownames(counts) <- make.unique(rownames(counts))
colnames(counts) <- make.unique(colnames(counts))

rownames(obs) <- make.unique(cell_names)
obs <- obs[colnames(counts), , drop = FALSE]

if (!opt$target_gene_col %in% colnames(obs)) {
  stop("[01] target_gene_col '", opt$target_gene_col, "' not found in h5ad obs. Available columns: ", paste(colnames(obs), collapse = ", "))
}

if (!file.exists(opt$pair_csv)) {
  stop("[01] pair_csv does not exist: ", opt$pair_csv)
}

tick01("read h5ad -> dgCMatrix")
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
  if (!all(keep_cells)) {
    counts <- counts[, keep_cells, drop = FALSE]
    obs <- obs[keep_cells, , drop = FALSE]
  } else {
    message("[01] all cells kept; skipping the subset copy")
  }
}

# CreateSeuratObject calls CalcN, which builds colSums(counts > 0) as a full
# lgCMatrix just to count nonzeros -- 28.77 of the 44.24 GiB the call costs at
# nnz = 1.43e9. Both statistics are computed off the (p, x) slots instead and
# handed in through meta.data, which is where CalcN would have put them. When
# the h5ad's obs already carries a column, the stock path overwrites CalcN's
# value with it, so it is left alone here too.
if (!"nCount_RNA" %in% colnames(obs))   obs$nCount_RNA   <- n_counts_per_cell(counts)
if (!"nFeature_RNA" %in% colnames(obs)) obs$nFeature_RNA <- n_features_per_cell(counts)

assay <- CreateAssay5Object(counts = counts)
rm(counts); gc(FALSE)

.op <- options(Seurat.object.assay.calcn = FALSE)
obj <- CreateSeuratObject(assay, assay = "RNA", meta.data = obs)
options(.op)
rm(assay); gc(FALSE)

# CalcN's two columns sit directly after orig.ident; keep that order so the
# metadata is column-for-column what the stock path produced.
.front <- intersect(c("orig.ident", "nCount_RNA", "nFeature_RNA"), colnames(obj@meta.data))
obj@meta.data <- obj@meta.data[, c(.front, setdiff(colnames(obj@meta.data), .front)),
                               drop = FALSE]
tick01("build Seurat object")
message("[01] Seurat object: ", nrow(obj), " genes x ", ncol(obj), " cells")
message("[01] Metadata columns include: ", paste(colnames(obj@meta.data), collapse = ", "))
message("[01] Saving: ", opt$out_rds)
save_rds_fast(obj, opt$out_rds)
tick01("saveRDS")
message("[01] Wrote: ", opt$out_rds)
message("[01] Finished at: ", Sys.time())
