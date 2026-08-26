#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Mixscale)
  library(Matrix)
  library(RANN)
})

# Resolve the helper next to this script regardless of how it was invoked.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(.self), "mixscale_de_genes.R"))
source(file.path(dirname(.self), "mixscale_loo_patch.R"))

# stage timer -- the log needs to show where step 2's time goes so the cost can
# be projected to a larger NT pool.
.stage_t0 <- Sys.time()
.stage_nm <- NULL
stage <- function(name) {
  if (!is.null(.stage_nm))
    message(sprintf("[02][time] %-34s %7.1f s", .stage_nm,
                    as.numeric(difftime(Sys.time(), .stage_t0, units = "secs"))))
  .stage_nm <<- name
  .stage_t0 <<- Sys.time()
}


option_list <- list(
  make_option("--obj_rds", type = "character"),
  make_option("--target_gene_col", type = "character", default = "target_gene"),
  make_option("--nt_label", type = "character", default = "ONE_INTERGENIC_SITE"),
  make_option("--guide_col", type = "character", default = "pair_key"),

  make_option("--ndims", type = "integer", default = 30),
  make_option("--num_neighbors", type = "integer", default = 20),
  make_option("--chunk_cells", type = "integer", default = 3000),

  make_option("--min_de_genes", type = "integer", default = 0),
  make_option("--max_de_genes", type = "integer", default = 100),
  make_option("--logfc_threshold", type = "double", default = 0),

  make_option("--out_rds", type = "character", default = "mixscale_obj.rds")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$obj_rds), !is.null(opt$out_rds))

# ------------------------------------------------------------
# Chunked CalcPerturbSig for Mixscale
# ------------------------------------------------------------

sparse_expm1 <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    x <- as(x, "dgCMatrix")
    x@x <- expm1(x@x)
    return(x)
  } else {
    return(expm1(x))
  }
}

CalcPerturbSig_chunked_for_Mixscale <- function(
    object,
    assay = "RNA",
    layer = "data",
    gd.class = "target_gene",
    nt.cell.class = "NT",
    reduction = "pca",
    ndims = 30,
    num.neighbors = 20,
    features = NULL,
    new.assay.name = "PRTB",
    chunk.cells = 3000,
    split.by = NULL,
    store.neighbors = FALSE,
    verbose = TRUE
) {
  DefaultAssay(object) <- assay

  if (is.null(features)) {
    features <- rownames(object[[assay]])
  }

  features <- intersect(features, rownames(object[[assay]]))

  if (length(features) == 0) {
    stop("No valid features found.")
  }

  if (verbose) {
    message("[CalcPerturbSig_chunked] features: ", length(features))
    message("[CalcPerturbSig_chunked] cells: ", ncol(object))
    message("[CalcPerturbSig_chunked] chunk.cells: ", chunk.cells)
  }

  expr <- GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )[features, , drop = FALSE]

  emb <- Embeddings(object, reduction = reduction)

  if (ndims > ncol(emb)) {
    stop("ndims is larger than available dimensions in reduction: ", reduction)
  }

  emb <- emb[, seq_len(ndims), drop = FALSE]

  meta <- object[[]]

  if (!gd.class %in% colnames(meta)) {
    stop(gd.class, " not found in object metadata.")
  }

  if (is.null(split.by)) {
    split.list <- list(rep1 = colnames(object))
  } else {
    if (!split.by %in% colnames(meta)) {
      stop(split.by, " not found in object metadata.")
    }
    split.list <- split(rownames(meta), meta[[split.by]])
  }

  one_matrix_gb <- length(features) * ncol(object) * 8 / 1024^3
  message("[CalcPerturbSig_chunked] Preallocating PRTB matrix: ~", round(one_matrix_gb, 2), " GB")

  prtb_mat <- matrix(
    NA_real_,
    nrow = length(features),
    ncol = ncol(object),
    dimnames = list(features, colnames(object))
  )

  neighbor_list <- list()

  for (rep_name in names(split.list)) {
    message("[CalcPerturbSig_chunked] Processing split: ", rep_name)

    cells_use <- intersect(split.list[[rep_name]], colnames(object))

    nt_cells <- cells_use[
      meta[cells_use, gd.class, drop = TRUE] == nt.cell.class
    ]

    if (length(nt_cells) < num.neighbors) {
      stop(
        "Split ", rep_name, " has only ", length(nt_cells),
        " NT cells, fewer than num.neighbors = ", num.neighbors
      )
    }

    message("[CalcPerturbSig_chunked] cells in split: ", length(cells_use))
    message("[CalcPerturbSig_chunked] NT cells in split: ", length(nt_cells))

    nt_emb <- emb[nt_cells, , drop = FALSE]
    all_emb <- emb[cells_use, , drop = FALSE]

    nt_expr_linear <- sparse_expm1(expr[, nt_cells, drop = FALSE])

    cell_chunks <- split(
      cells_use,
      ceiling(seq_along(cells_use) / chunk.cells)
    )

    if (store.neighbors) {
      neighbor_chunks <- vector("list", length(cell_chunks))
    }

    # RANN::nn2 rebuilds a kd-tree over `data` on every call, so querying once
    # per chunk rebuilt the tree over all NT cells length(cell_chunks) times
    # (102 times at 305k cells / chunk.cells = 3000). One call for the whole
    # split returns the same neighbours -- same data, same query order, same k
    # -- and the index matrix is only n_cells x k integers.
    message("[CalcPerturbSig_chunked] kNN for all ", length(cells_use),
            " cells in one nn2 call")
    nn_all <- RANN::nn2(
      data = nt_emb,
      query = all_emb,
      k = num.neighbors
    )$nn.idx
    rownames(nn_all) <- cells_use

    for (i in seq_along(cell_chunks)) {
      query_cells <- cell_chunks[[i]]

      message(
        "[CalcPerturbSig_chunked] chunk ", i, "/", length(cell_chunks),
        " | cells = ", length(query_cells)
      )

      nn <- nn_all[query_cells, , drop = FALSE]

      W <- Matrix::sparseMatrix(
        i = as.vector(nn),
        j = rep(seq_len(nrow(nn)), times = ncol(nn)),
        x = 1 / num.neighbors,
        dims = c(length(nt_cells), length(query_cells))
      )

      avg_nt <- nt_expr_linear %*% W
      avg_nt <- log1p(avg_nt)

      query_expr <- expr[, query_cells, drop = FALSE]

      prtb_chunk <- as.matrix(query_expr - avg_nt)

      rownames(prtb_chunk) <- features
      colnames(prtb_chunk) <- query_cells

      prtb_mat[, query_cells] <- prtb_chunk

      if (store.neighbors) {
        rownames(nn) <- query_cells
        neighbor_chunks[[i]] <- nn
      }

      rm(nn, W, avg_nt, query_expr, prtb_chunk)
      gc()
    }

    if (store.neighbors) {
      neighbor_list[[make.names(paste0(new.assay.name, "_", rep_name))]] <- do.call(
        rbind,
        neighbor_chunks
      )
      rm(neighbor_chunks)
    }

    rm(nt_expr_linear, nt_emb, all_emb, cell_chunks, nn_all)
    gc()
  }

  message("[CalcPerturbSig_chunked] Creating assay: ", new.assay.name)

  object[[new.assay.name]] <- CreateAssayObject(
    data = prtb_mat,
    check.matrix = FALSE
  )

  DefaultAssay(object) <- new.assay.name

  if (store.neighbors) {
    object@tools[[paste0("CalcPerturbSig_chunked_for_Mixscale.", assay, ".", reduction)]] <- neighbor_list
  }

  rm(prtb_mat)
  gc()

  return(object)
}

# ------------------------------------------------------------
# Main script
# ------------------------------------------------------------

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

message("[02] Parameters")
message("[02] target_gene_col: ", opt$target_gene_col)
message("[02] guide_col: ", opt$guide_col)
message("[02] nt_label: ", opt$nt_label)
message("[02] ndims: ", opt$ndims)
message("[02] num_neighbors: ", opt$num_neighbors)
message("[02] chunk_cells: ", opt$chunk_cells)
message("[02] min_de_genes: ", opt$min_de_genes)
message("[02] max_de_genes: ", opt$max_de_genes)
message("[02] logfc_threshold: ", opt$logfc_threshold)

DefaultAssay(obj) <- "RNA"

stage("normalize/HVG/scale/PCA")
message("[02] NormalizeData / FindVariableFeatures / ScaleData / RunPCA")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = max(50, opt$ndims), verbose = FALSE)

message("[02] Slim object before CalcPerturbSig_chunked")
obj <- DietSeurat(
  obj,
  assays = "RNA",
  dimreducs = "pca",
  graphs = NULL,
  misc = FALSE
)

DefaultAssay(obj) <- "RNA"

# scale.data fed RunPCA and nothing after it. It is dense -- 2000 variable
# features x n_cells x 8 B, i.e. ~4.6 GB at 305k cells -- so drop it before
# CalcPerturbSig allocates the PRTB matrix.
try({
  sd_dim <- dim(SeuratObject::LayerData(obj[["RNA"]], layer = "scale.data"))
  if (!is.null(sd_dim) && prod(sd_dim) > 0) {
    message(sprintf("[02] dropping scale.data (%d x %d, %.2f GB)",
                    sd_dim[1], sd_dim[2], prod(sd_dim) * 8 / 1024^3))
    SeuratObject::LayerData(obj[["RNA"]], layer = "scale.data") <- NULL
  }
}, silent = TRUE)
gc()

# RunMixscale reads the PRTB assay only at its body line 259, immediately subset
# to de.genes, and de.genes is computed from the RNA assay -- it never touches
# PRTB. Selecting those genes up front lets CalcPerturbSig build a PRTB matrix of
# ~max_de_genes rows instead of all ~60k, which is where step 2's time and memory
# were going. The list is handed back to RunMixscale so it does not recompute it.
stage("DE gene selection (wilcox)")
message("[02] Precomputing DE genes (RNA assay) to restrict CalcPerturbSig")
de_genes_list <- mixscale_de_genes(
  object = obj,
  labels = opt$target_gene_col,
  nt.class.name = opt$nt_label,
  de.assay = "RNA",
  logfc.threshold = opt$logfc_threshold,
  fine.mode = fine_mode,
  fine.mode.labels = if (fine_mode) opt$guide_col else NULL,
  verbose = TRUE
)
prtb_features <- unique(unlist(de_genes_list, use.names = FALSE))
message("[02] DE genes across all targets: ", length(prtb_features),
        " (of ", nrow(obj[["RNA"]]), " -> ",
        round(nrow(obj[["RNA"]]) / max(1, length(prtb_features))), "x fewer PRTB rows)")

if (length(prtb_features) == 0) {
  stop("[02] No DE genes selected for any target; cannot build a PRTB assay.")
}

stage("CalcPerturbSig (kNN + PRTB)")
message("[02] CalcPerturbSig_chunked_for_Mixscale")
obj <- CalcPerturbSig_chunked_for_Mixscale(
  object = obj,
  assay = "RNA",
  layer = "data",
  gd.class = opt$target_gene_col,
  nt.cell.class = opt$nt_label,
  reduction = "pca",
  ndims = opt$ndims,
  num.neighbors = opt$num_neighbors,
  features = prtb_features,
  new.assay.name = "PRTB",
  chunk.cells = opt$chunk_cells,
  split.by = NULL,
  store.neighbors = FALSE,
  verbose = TRUE
)

DefaultAssay(obj) <- "PRTB"

stage("RunMixscale (scores + LOO)")
patch_runmixscale_loo()
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
    DE.gene = de_genes_list,
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
    DE.gene = de_genes_list,
    new.class.name = "mixscale_score",
    fine.mode = FALSE,
    verbose = TRUE
  )
}

message("[02] Slim object after RunMixscale for DE")

mixscale_meta <- obj@meta.data

obj_slim <- obj

DefaultAssay(obj_slim) <- "RNA"

if ("PRTB" %in% Assays(obj_slim)) {
  obj_slim[["PRTB"]] <- NULL
}

obj_slim <- DietSeurat(
  obj_slim,
  assays = "RNA",
  layers = "counts",
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)

obj_slim@meta.data <- mixscale_meta
DefaultAssay(obj_slim) <- "RNA"

gc()

stage("saveRDS")
message("[02] Saving slim object: ", opt$out_rds)
saveRDS(
  obj_slim,
  opt$out_rds,
  compress = FALSE
)

message("[02] Wrote slim object: ", opt$out_rds)
message("[02] Final assays: ", paste(Assays(obj_slim), collapse = ", "))
message("[02] Final RNA layers: ", paste(Layers(obj_slim[["RNA"]]), collapse = ", "))
stage(NULL)
