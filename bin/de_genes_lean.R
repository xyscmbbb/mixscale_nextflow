# Memory-lean, hoisted replacement for the FoldChange half of
# Seurat:::TopDEGenesMixscape -> FindMarkers, for the "data" layer only.
#
# Measured at 305,110 cells (38,292 genes, nnz = 1.43e9): stock FoldChange adds
# 26.2 GB on top of a 48.2 GB baseline, and does it once per guide. Three
# separate wastes, all removed here:
#
#   1. cells.2 is the whole NT pool -- the SAME cells for every guide and every
#      target -- yet pct.2 / data.2 are recomputed per call. Hoisted to one
#      pass, reused for every guide.
#   2. object[features, cells.2] materialises a full-size dgCMatrix (15.8 GB)
#      four times per call. rowSums(x > 0) and rowSums(expm1(x)) are both
#      additive over disjoint column blocks, so one block at a time is enough.
#   3. pct is rowSums(x > 0), which builds a 10.4 GB lgCMatrix just to count
#      stored nonzeros. tabulate(x@i + 1L) returns the identical integer counts
#      with no intermediate matrix.
#
# Exactness: (2) changes the association order of the expm1 row sums, so data.2
# can differ from the stock value in the last ulp. That value never leaves
# TopDEGenesMixscape -- it returns rownames() only, the row order is
# order(p_val, -abs(pct.1 - pct.2)), and fc is read solely by the
# abs(fc) >= logfc.threshold gate. pct.1/pct.2 (integer counts) and the wilcox
# p-values are bit-identical. Verified against the stock path in
# validate_de_genes.R.

# rowSums(x > 0) and rowSums(expm1(x)) over a column set, one block at a time.
# cols is an integer index vector into d's columns.
nt_row_stats_blocked <- function(d, cols, block = 20000L) {
  ng  <- nrow(d)
  nz  <- integer(ng)
  ex  <- numeric(ng)
  n   <- length(cols)
  a   <- 1L
  while (a <= n) {
    b <- min(a + block - 1L, n)
    s <- d[, cols[a:b], drop = FALSE]
    nz <- nz + tabulate(s@i[s@x > 0] + 1L, nbins = ng)
    ex <- ex + rowSums(expm1(s))
    rm(s)
    a <- b + 1L
  }
  list(nz = nz, ex = ex, n = n)
}

# Same as the block above but for a small cell set, in one shot.
row_stats_direct <- function(d, cols) {
  s  <- d[, cols, drop = FALSE]
  ng <- nrow(d)
  list(nz = tabulate(s@i[s@x > 0] + 1L, nbins = ng),
       ex = rowSums(expm1(s)),
       n  = length(cols))
}

# Reproduces FoldChange.StdAssay(slot = "data", pseudocount.use = 1, base = 2):
#   log((rowSums(expm1(x)) + 1) / NCOL(x), 2)
fc_from_stats <- function(s1, s2) {
  d1 <- log((s1$ex + 1) / s1$n, base = 2)
  d2 <- log((s2$ex + 1) / s2$n, base = 2)
  data.frame(avg_log2FC = d1 - d2,
             pct.1      = round(s1$nz / s1$n, 3),
             pct.2      = round(s2$nz / s2$n, 3))
}
