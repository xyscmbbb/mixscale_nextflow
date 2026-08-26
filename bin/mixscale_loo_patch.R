# Exact rewrite of RunMixscale's leave-one-out PRTB score.
#
# Stock body (Mixscale::RunMixscale):
#   omit_mat <- outer(de.genes, de.genes, `!=`)
#   calc_pvec2 <- function(include_gene)
#     matrixStats::rowSums2(pvec_mat[, include_gene, drop = F]) / sum(vec_mat[include_gene])
#   gv.list[[gene]][[s]][, de.genes] <- apply(omit_mat, 2, calc_pvec2)
#
# For each of the D de.genes it re-sums the other D-1 columns of pvec_mat, and
# pvec_mat[, include_gene] materialises an (n_cells x D-1) copy every time. That
# is O(n * D^2) work and D allocations of ~n*D doubles. With D = 100 and 600k
# NT cells each copy is ~475 MB and the block dominates step 2.
#
# But leaving one term out of a sum is just the total minus that term, in both
# the numerator and the denominator. The rewrite below is the same arithmetic in
# one pass -- O(n * D), no per-gene copies -- not an approximation.
patch_runmixscale_loo <- function(verbose = TRUE) {
  stopifnot(requireNamespace("Mixscale", quietly = TRUE))
  src <- deparse(body(Mixscale::RunMixscale), width.cutoff = 500L)

  old <- c(
    "omit_mat <- outer(de.genes, de.genes, `!=`)",
    "calc_pvec2 <- function(include_gene) {",
    "pvec2 <- matrixStats::rowSums2(pvec_mat[, include_gene, drop = F])/sum(vec_mat[include_gene])",
    "return(pvec2)",
    "}",
    "gv.list[[gene]][[s]][, de.genes] <- apply(omit_mat, 2, calc_pvec2)")
  norm <- function(x) gsub("[[:space:]]+", " ", trimws(x))
  hit <- which(norm(src) == norm(old[1]))
  if (length(hit) != 1L)
    stop("[loo] expected exactly one LOO block in RunMixscale, found ", length(hit),
         " -- Mixscale version changed, refusing to patch")
  span <- hit:(hit + length(old) - 1L)
  if (!identical(norm(src[span]), norm(old)))
    stop("[loo] LOO block does not match the expected source; refusing to patch")

  new <- c(
    "pvec_tot <- matrixStats::rowSums2(pvec_mat)",
    "vec_den <- sum(vec_mat) - vec_mat",
    "loo <- sweep(pvec_tot - pvec_mat, 2, vec_den, `/`)",
    "colnames(loo) <- de.genes",
    "gv.list[[gene]][[s]][, de.genes] <- loo")

  src <- append(src[-span], new, after = hit - 1L)
  f <- Mixscale::RunMixscale
  body(f) <- parse(text = paste(src, collapse = "\n"))[[1]]
  environment(f) <- environment(Mixscale::RunMixscale)
  assignInNamespace("RunMixscale", f, ns = "Mixscale")
  if (verbose) message("[02] RunMixscale LOO block patched (O(n*D^2) -> O(n*D), exact)")
  invisible(TRUE)
}
