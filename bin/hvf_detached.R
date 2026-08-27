# FindVariableFeatures computes a 2000-element ranking and writes it back into
# the object -- and the write-back, not the computation, is what costs the
# memory. `obj[["RNA"]]@meta.data <- v` hands the assay out through a temporary,
# so when the slot is set R sees two references and duplicates the assay; an S4
# duplicate deep-copies the @layers list, counts matrix and all.
#
# Measured at CM scale and projected to 305,110 cells (nnz = 1.430e9):
#   VST kernels alone (no write-back)          15.3 B/nnz  ->  20.4 GB
# CORRECTION: that 20.4 GB came from R's gc(), which cannot see the Eigen copies
# Seurat's kernels make. The cgroup trace at 305k shows the kernels actually cost
# ~48 GB -- they, not the write-back, set the container peak. See hvf_kernels.c.
#   obj <- FindVariableFeatures(obj)           36.6 B/nnz  ->  48.7 GB
#   assay <- FindVariableFeatures(assay)       30.1 B/nnz  ->  40.1 GB
#   detach assay, modify, re-attach            15.7 B/nnz  ->  20.9 GB
#
# So the fix is not to make VST cheaper -- it is already free -- but to hold the
# sole reference to the assay while writing, which lets R write in place.
#
# The body below is FindVariableFeatures.StdAssay's own tail, verbatim: the same
# `vf_vst_<layer>_` column naming, the same var.features reset, and the same
# two-step VariableFeatures() re-derivation. The only change is that it runs on
# a detached assay and assigns through slots rather than through `[[<-`, which
# would re-introduce the temporary.
find_variable_features_detached <- function(obj, assay = "RNA",
                                            nfeatures = 2000L, span = 0.3,
                                            clip = NULL, verbose = FALSE) {
  a <- obj@assays[[assay]]
  # Drop the object's reference before touching the assay, so `a` is the only
  # thing pointing at the counts matrix while it is modified.
  obj@assays[[assay]] <- NULL

  layer <- SeuratObject::Layers(object = a, search = "counts")
  # fast = TRUE skips re-attaching dimnames, which would copy the whole matrix.
  d <- SeuratObject::LayerData(object = a, layer = layer, fast = TRUE)
  # StdAssay picks the method by class: a dgCMatrix inherits V3Matrix, so it
  # dispatches to FindVariableFeatures.default (the `method=` form), NOT the
  # generic, whose V3Matrix method takes selection.method/clip.max instead.
  # vst_lowmem is VST.dgCMatrix with Seurat's two Eigen-transposing kernels
  # replaced by transpose-free ones; gated identical() on the full VST
  # data.frame. Seurat's own path costs three extra copies of the counts matrix
  # per call, twice -- 48 GB at 305k, and invisible to R's gc().
  load_hvf_kernels(dirname(.self))
  hvf.info <- vst_lowmem(data = d, nselect = nfeatures,
                         span = span, clip = clip, verbose = verbose)
  rm(d)

  colnames(hvf.info) <- paste("vf", "vst", layer, colnames(hvf.info), sep = "_")
  rownames(hvf.info) <- SeuratObject::Features(x = a, layer = layer)
  md <- a@meta.data
  md[["var.features"]] <- NULL
  md[["var.features.rank"]] <- NULL
  for (cn in colnames(hvf.info)) md[[cn]] <- NULL
  for (cn in colnames(hvf.info)) md[[cn]] <- hvf.info[[cn]]
  a@meta.data <- md
  SeuratObject::VariableFeatures(a) <-
    SeuratObject::VariableFeatures(object = a, nfeatures = nfeatures, method = "vst")

  obj@assays[[assay]] <- a
  obj
}
