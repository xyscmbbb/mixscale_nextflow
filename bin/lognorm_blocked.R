# Blockwise LogNormalize, bit-identical to NormalizeData(method = "LogNormalize").
#
# Seurat routes a dgCMatrix through Seurat:::LogNorm, an Rcpp function that takes
#   Eigen::SparseMatrix<double> data
# BY VALUE. So the counts matrix exists three times at once -- the R original,
# Rcpp's Eigen copy, and the R matrix built from the return -- which is why
# NormalizeData peaks at 59.1 GB on a 16 GB matrix at 305,110 cells.
#
# The kernel it runs is elementwise given the column sums:
#   colSum = sum of the column's stored values, in storage order
#   x     -> log1p(x / colSum * scale_factor)
# Reproducing that expression in the same association order gives bit-identical
# output while touching one column block at a time. Verified with identical()
# against NormalizeData in validate_lognorm.R.
#
# `dimnames` is taken separately because Seurat v5 stores layers WITHOUT them:
# attaching them to the input first would copy the whole 17 GB matrix just to
# label it, whereas the result is built here and can be labelled for free. The
# new matrix also SHARES @i and @p with the input, so the data layer costs only
# its own @x (8 B per nonzero), not a second full matrix.
lognorm_blocked <- function(cnt, scale.factor = 1e4, block = 20000L,
                            dimnames = NULL) {
  stopifnot(inherits(cnt, "dgCMatrix"))
  p  <- cnt@p
  cs <- Matrix::colSums(cnt)
  xn <- numeric(length(cnt@x))
  n  <- ncol(cnt)
  a  <- 1L
  while (a <= n) {
    b  <- min(a + block - 1L, n)
    k1 <- p[a] + 1L
    k2 <- p[b + 1L]
    if (k2 >= k1) {
      reps <- diff(p[a:(b + 1L)])
      xn[k1:k2] <- log1p(cnt@x[k1:k2] / rep.int(cs[a:b], reps) * scale.factor)
    }
    a <- b + 1L
  }
  new("dgCMatrix", i = cnt@i, p = cnt@p, x = xn,
      Dim = cnt@Dim,
      Dimnames = if (is.null(dimnames)) cnt@Dimnames else dimnames)
}
