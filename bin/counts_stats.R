# nCount_RNA / nFeature_RNA without materialising anything nonzero-sized.
#
# CreateSeuratObject gets these from SeuratObject:::CalcN, whose nFeature term
# is colSums(counts > 0) -- a whole lgCMatrix (8 B per nonzero) built only to
# count nonzeros. Measured at CM scale and projected to nnz = 1.430e9
# (305,110 cells), CalcN is 28.77 GiB of the 44.24 GiB that
# CreateSeuratObject(counts = ...) costs; with it switched off the same call
# costs 1.79 GiB.
#
# A dgCMatrix column j holds p[j+1] - p[j] stored entries, so nFeature is
# diff(@p) minus any stored entry that is not positive. Explicit zeros are rare
# but are not assumed away: @x is scanned for them in slices, so the largest
# allocation is one slice.
n_features_per_cell <- function(m, chunk = 2e7) {
  x   <- m@x
  n   <- length(x)
  bad <- integer(0)
  a   <- 1
  while (a <= n) {
    b <- min(a + chunk - 1, n)
    w <- which(x[a:b] <= 0)
    if (length(w)) bad <- c(bad, as.integer(w + a - 1))
    a <- b + 1
  }
  nf <- diff(m@p)
  if (length(bad)) {
    # bad - 1L is the 0-based offset into @x; findInterval against @p maps it
    # to the column that holds it (empty columns cannot, so ties are safe).
    nf <- nf - tabulate(findInterval(bad - 1L, m@p), nbins = ncol(m))
  }
  # colSums() on the logical matrix CalcN builds returns integer, and that type
  # reaches meta.data; diff(@p) is integer too, so it is left alone.
  nf
}

# Matrix::colSums on a dgCMatrix allocates only its result; it is what CalcN
# uses for nCount and is kept as-is. unname() because CalcN's value reaches
# meta.data unnamed.
n_counts_per_cell <- function(m) unname(Matrix::colSums(m))
