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
  make_option("--keep_cell_col", type = "character", default = ""),
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
h5py    <- import("h5py",    convert = FALSE)
ad_io   <- import("anndata.io", convert = FALSE)

is_gcs <- grepl("^gs://", opt$h5ad)

# One handle for obs, var and X. h5py takes a Python file object as happily as a
# path, which is what lets --h5ad be a gs:// URI: gcsfs turns the object into a
# seekable file and h5py then issues ranged GETs for only the bytes it reads.
# Combined with the CSR row subset below, a job transfers its own cells and
# nothing else -- the whole h5ad never lands on the VM's disk.
gcs_file <- NULL
if (is_gcs) {
  gcsfs <- tryCatch(import("gcsfs", convert = FALSE), error = function(e)
    stop("[01] --h5ad is a gs:// URI but the python gcsfs package is missing. ",
         "Use image r-mixscale:1.1.11 or newer, or pass a local path."))
  # block_size is the ranged-GET size. Too small and a 5 GB slice becomes tens of
  # thousands of round trips; 16 MiB keeps the read sequential and the count sane.
  #
  # cache_type MATTERS MORE THAN block_size, and the default is the wrong one.
  # fsspec defaults to "readahead", which holds exactly ONE block and drops it on
  # the next miss. h5py does not read front to back -- it walks the superblock and
  # object headers, which live at scattered offsets -- so every metadata hop
  # evicts the block the previous hop just paid 16 MiB for. "blockcache" keeps an
  # LRU of blocks instead, so those hops hit. Measured on a 292 MB CSR h5ad in
  # GCS, reading every row:
  #
  #   readahead  16 MiB   405.3 MB fetched   1.40x amplification   4.80 s
  #   blockcache 16 MiB   302.0 MB fetched   1.05x                 2.80 s
  #
  # maxblocks caps the LRU, and hence the RAM this holds: 8 x 16 MiB = 134 MB.
  # Anything from 4 up measured the same wall clock, so this is the small end of
  # the flat region rather than a tuned value.
  fs_gcs   <- gcsfs$GCSFileSystem()
  gcs_file <- fs_gcs$open(sub("^gs://", "", opt$h5ad), "rb",
                          block_size  = bi$int(16777216),
                          cache_type  = "blockcache",
                          cache_options = bi$dict(maxblocks = bi$int(8)))
  fh <- h5py$File(gcs_file, "r")
  message("[01] reading over gcsfs (ranged GETs, 16 MiB blocks, blockcache); ",
          "the h5ad is NOT downloaded")
} else {
  fh <- h5py$File(opt$h5ad, "r")
}

# obs/var only -- X stays on disk (or in GCS) so it is never materialised twice.
obs_py     <- ad_io$read_elem(fh[["obs"]])
var_py     <- ad_io$read_elem(fh[["var"]])
cell_names <- py_to_r(obs_py$index$astype("str")$to_list())
gene_names <- py_to_r(var_py$index$astype("str")$to_list())
obs        <- py_to_r(obs_py)
rm(obs_py, var_py)

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

# One bounded HDF5 slice, ds[from, to), as a fresh R vector. Callers keep their
# own preallocated destination and assign into it at top level, where R can do
# the write in place -- passing the destination into a function would force a
# copy of the whole vector.
h5_slice <- function(ds, from, to) {
  py_to_r(np$asarray(ds$`__getitem__`(bi$slice(bi$int(from), bi$int(to)))))
}

# ---------------------------------------------------------------------------
# Which cells this job wants is decided HERE, before X is touched, because a
# CSR h5ad can then be read for just those rows. That is what lets a run start
# from the whole filtered_correct_pairs.h5ad instead of a per-target shard:
# the shard's only job was to pre-select these cells, and pre-selecting them is
# 99% NT cells copied once per target.
#
# Everything below needs obs and the pair CSV only -- both cheap.
# ---------------------------------------------------------------------------
cell_names_u <- make.unique(cell_names)
rownames(obs) <- cell_names_u

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

pair_row <- match(rownames(obs), pair_df[[opt$cell_col]])

# Add guide_col from external per-cell CSV into h5ad obs before creating Seurat object.
obs[[opt$guide_col]] <- pair_df[pair_row, opt$guide_col]
n_missing_guides <- sum(is.na(obs[[opt$guide_col]]))
message("[01] Added guide column to metadata: ", opt$guide_col)
message("[01] Missing ", opt$guide_col, ": ", n_missing_guides, " / ", nrow(obs))

if (n_missing_guides == nrow(obs)) {
  stop("[01] All guide_col values are NA. Cell names probably do not match between h5ad obs_names and pair_csv cell_col.")
}

keep_cells <- rep(TRUE, nrow(obs))

if (subset_to_gene) {
  keep_cells <- obs[[opt$target_gene_col]] %in% c(opt$perturb_gene, opt$nt_label)
  message("[01] Subsetting cells to ", opt$perturb_gene, " + ", opt$nt_label, ": ", sum(keep_cells), " / ", nrow(obs))
}

# Optional allowlist column on the pair CSV. This is the annslicer `filter`
# stage from the PotC shard pipeline (keep_cell = matched-pair AND QC-pass)
# moved into the job, so the big h5ad does not have to be pre-filtered either.
# A cell missing from the CSV is not in the allowlist and is dropped.
if (nzchar(opt$keep_cell_col)) {
  if (!opt$keep_cell_col %in% colnames(pair_df)) {
    stop("[01] keep_cell_col '", opt$keep_cell_col, "' not found in pair_csv. Available columns: ", paste(colnames(pair_df), collapse = ", "))
  }
  kc <- pair_df[pair_row, opt$keep_cell_col]
  kc <- if (is.logical(kc)) !is.na(kc) & kc
        else tolower(as.character(kc)) %in% c("true", "t", "1", "yes", "y")
  message("[01] Applying keep_cell_col '", opt$keep_cell_col, "': ", sum(kc), " / ", length(kc), " cells allowed")
  keep_cells <- keep_cells & kc
}

keep_idx <- which(keep_cells)
if (length(keep_idx) == 0) {
  stop("[01] No cells remain after subsetting. Check --perturb_gene, --target_gene_col, --nt_label and --keep_cell_col.")
}
subset_rows <- length(keep_idx) < length(keep_cells)
message("[01] Cells selected: ", length(keep_idx), " / ", length(keep_cells))

rows_already_subset <- FALSE

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
  # The 2^31-1 cap is on the matrix we BUILD, not on the file. Reading a row
  # subset out of a file whose total nnz is over the cap is fine and is the
  # whole point of the row path -- so this only guards the whole-read path
  # (the row path checks nnz_keep in its own branch).
  if (nnz >= 2147483647 && !(enc == "csr_matrix" && subset_rows))
    stop("[01] nnz = ", format(nnz, scientific = FALSE), " exceeds the dgCMatrix ",
         "limit of 2^31-1. Repack the h5ad as CSR so the wanted cells can be ",
         "read directly, or shard it into fewer cells.")

  if (enc == "csr_matrix" && subset_rows) {
    # CSR indptr is over CELLS, so the wanted cells are a set of row slices and
    # only their nonzeros have to be read. This is what makes running from the
    # whole h5ad cost the same as running from a shard: a shard is 99% NT cells,
    # and so is the slice we read here.
    #
    # Cost is one HDF5 read per RUN of consecutive wanted cells. Repacking the
    # h5ad sorted by target_gene (tools/repack_h5ad.py) makes that two runs --
    # the NT block and the target's block -- i.e. two sequential reads. Without
    # the sort it is still correct and still moves only nnz_kept bytes, just in
    # many small reads; the log line below is how you tell which you are in.
    ip_ds  <- Xg[["indptr"]]
    p_all  <- as.numeric(py_to_r(np$asarray(ip_ds)))
    row_nnz  <- p_all[keep_idx + 1L] - p_all[keep_idx]
    nnz_keep <- sum(row_nnz)

    brk       <- which(diff(keep_idx) != 1L)
    run_first <- c(1L, brk + 1L)
    run_last  <- c(brk, length(keep_idx))
    message(sprintf("[01] CSR row subset: %d cells in %d run(s), nnz = %.0f of %.0f (%.1f%%)",
                    length(keep_idx), length(run_first), nnz_keep, nnz,
                    100 * nnz_keep / nnz))
    if (nnz_keep >= 2147483647)
      stop("[01] selected nnz = ", format(nnz_keep, scientific = FALSE),
           " exceeds the dgCMatrix limit of 2^31-1.")

    idx_ds <- Xg[["indices"]]
    dat_ds <- Xg[["data"]]
    i_vec <- integer(nnz_keep)
    x_vec <- numeric(nnz_keep)
    off <- 0
    chunk <- 5e7
    for (r in seq_along(run_first)) {
      a <- p_all[keep_idx[run_first[r]]]
      b <- p_all[keep_idx[run_last[r]] + 1L]
      # Bounded slices, and assigned at top level so R writes in place rather
      # than copying the destination vector.
      while (a < b) {
        e <- min(a + chunk, b)
        i_vec[(off + 1):(off + (e - a))] <- h5_slice(idx_ds, a, e)
        x_vec[(off + 1):(off + (e - a))] <- h5_slice(dat_ds, a, e)
        off <- off + (e - a)
        a <- e
      }
    }
    fh$close()
    if (!is.null(gcs_file)) gcs_file$close()

    p_vec <- as.integer(c(0, cumsum(row_nnz)))
    counts <- new("dgCMatrix", i = i_vec, p = p_vec, x = x_vec,
                  Dim = c(n_gene, length(keep_idx)),
                  Dimnames = list(gene_names, cell_names_u[keep_idx]))
    rm(i_vec, p_vec, x_vec, p_all, row_nnz); gc(FALSE)
    rows_already_subset <- TRUE
  } else {
    if (subset_rows && enc == "csc_matrix")
      message("[01][WARN] X is csc_matrix, so a cell subset cannot be read directly ",
              "(indptr is over genes); reading the whole matrix. Repack the h5ad as ",
              "CSR with tools/repack_h5ad.py to read only the wanted cells.")
    message("[01] reading indices ...")
    i_vec <- h5_read_vec(Xg[["indices"]], nnz, TRUE)
    p_vec <- as.integer(py_to_r(np$asarray(Xg[["indptr"]])))
    message("[01] reading data ...")
    x_vec <- h5_read_vec(Xg[["data"]], nnz, FALSE)
    fh$close()
    if (!is.null(gcs_file)) gcs_file$close()

    # new() shares these vectors rather than copying them; the validity check is
    # what guarantees the stored indices were sorted within each slice.
    dm <- if (enc == "csr_matrix") c(n_gene, n_cell) else c(n_cell, n_gene)
    dn <- if (enc == "csr_matrix") list(gene_names, cell_names)
          else                     list(cell_names, gene_names)
    counts <- new("dgCMatrix", i = i_vec, p = p_vec, x = x_vec, Dim = dm, Dimnames = dn)
    rm(i_vec, p_vec, x_vec); gc(FALSE)
  }

  if (enc == "csc_matrix" && !rows_already_subset) {
    message("[01] transposing to genes x cells (one copy; export shards as CSR to skip this)")
    counts <- Matrix::t(counts)
    counts <- as(counts, "CsparseMatrix")
    gc(FALSE)
  }
} else {
  message("[01][WARN] X is not sparse (encoding-type '", enc, "'); ",
          "falling back to the scipy COO path (much heavier).")
  fh$close()
  if (!is.null(gcs_file)) gcs_file$close()
  if (is_gcs)
    stop("[01] X is not sparse and --h5ad is a gs:// URI. The dense fallback ",
         "re-reads the file by path; localise it first, or store X sparse.")
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

# The row path already named its columns cell_names_u[keep_idx]; the whole-read
# paths named them from the raw obs_names, and stock applied make.unique here.
if (!rows_already_subset) colnames(counts) <- cell_names_u

if (!is.na(opt$max_features) && opt$max_features > 0 && nrow(counts) > opt$max_features) {
  message("[01] Keeping first ", opt$max_features, " features out of ", nrow(counts), ".")
  counts <- counts[seq_len(opt$max_features), , drop = FALSE]
}

rownames(counts) <- gsub("_", "-", rownames(counts))
rownames(counts) <- make.unique(rownames(counts))

# make.unique is applied to the FULL cell list above (cell_names_u) and only
# then subset, so a name stays whatever it was called in the unsubset file --
# which cells happen to be selected must not change any cell's name.
obs <- obs[keep_idx, , drop = FALSE]

if (subset_rows) {
  # anndata drops unused categories when it subsets, so a per-target shard
  # arrived with target_gene carrying exactly the two levels in it. Subsetting
  # in-job has to do the same or the factor keeps every level in the whole file
  # -- at 5,000 targets that is 4,999 empty classes handed to step 2, which is
  # not a cosmetic difference. This is what makes running from the whole h5ad
  # produce the identical object to running from the shard.
  obs <- droplevels(obs)
}

if (!rows_already_subset && subset_rows) {
  # Whole matrix was read (CSC source, or no subsetting possible at read time).
  counts <- counts[, keep_idx, drop = FALSE]
  gc(FALSE)
}
stopifnot(identical(colnames(counts), rownames(obs)))

tick01("read h5ad -> dgCMatrix")

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
