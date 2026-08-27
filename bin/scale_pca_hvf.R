# ScaleData + RunPCA without paying for the data layer
# ----------------------------------------------------
# Both functions read ONLY the variable features -- 2,000 rows of the data layer --
# but they run through the Seurat object, and ScaleData stores its result by writing
# the assay back. That write-back hands the assay out through a temporary, so R
# duplicates it, and duplicating an Assay5 deep-copies @layers: the whole 17 GB data
# layer is copied to store a 4.6 GB scale.data. Measured at 305,110 cells the stage
# peaks at 51.78 GB while its live data is 16.5 GB.
#
# The two are run on a temporary object holding only the variable-feature rows, and
# the reduction is transplanted back. ScaleData centres and scales each feature
# independently and RunPCA reads nothing but scale.data, so with the same rows in the
# same order the scale.data, embeddings, loadings and stdev are identical.
#
# scale.data is deliberately NOT transplanted. Nothing between RunPCA and the
# DietSeurat that follows reads it -- the caller drops it explicitly -- so putting it
# back would allocate 4.6 GB to delete it two lines later.

scale_pca_hvf <- function(object, assay = "RNA", npcs = 50, features = NULL,
                          reduction.name = "pca", verbose = FALSE) {
  vf <- features
  if (is.null(vf)) vf <- SeuratObject::VariableFeatures(object)
  stopifnot(length(vf) > 0)

  a <- object[[assay]]
  # Positional access to the stored layer: LayerData() re-attaches dimnames on
  # fetch, which duplicates the whole matrix before a handful of rows are taken.
  d <- a@layers[["data"]]
  if (is.null(d)) d <- SeuratObject::LayerData(a, layer = "data")
  ridx <- match(vf, rownames(a))
  stopifnot(!anyNA(ridx))
  sub <- d[ridx, , drop = FALSE]
  dimnames(sub) <- list(vf, colnames(object))
  rm(d, a, ridx)
  gc(FALSE)

  tmp <- SeuratObject::CreateSeuratObject(
    SeuratObject::CreateAssay5Object(data = sub), assay = assay)
  rm(sub)
  SeuratObject::VariableFeatures(tmp) <- vf
  tmp <- Seurat::ScaleData(tmp, verbose = verbose)
  tmp <- Seurat::RunPCA(tmp, npcs = npcs, verbose = verbose,
                        reduction.name = reduction.name)

  red <- tmp[[reduction.name]]
  # RunPCA stamps the reduction with the assay it ran on; that is the same assay
  # name here, but set it explicitly so the transplant does not depend on it.
  red@assay.used <- assay
  rm(tmp)
  gc(FALSE)
  object@reductions[[reduction.name]] <- red
  object
}
