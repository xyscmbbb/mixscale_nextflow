# Collapsed drop-in for glmGamPoi::glm_gp on the mixscale step-3 design.
#
# Scope: size_factors = FALSE (offset identically zero), overdispersion = TRUE,
# overdispersion_shrinkage = TRUE, subsample = FALSE, no ridge penalty, no groups.
# That is exactly what wmvRegDE calls.
#
# ---------------------------------------------------------------------------
# 1. The collapse
#
# mixscale's design is ~ 1 + weight + log_ct (or, with several cell lines,
# ~ 0 + cell_type + weight:cell_type + log_ct). NT cells carry weight == 0 in
# every weight column, so an NT cell's design row is determined by
# (cell_type, log_ct) alone, and log_ct = log1p(integer nCount) takes only a few
# tens of thousands of distinct values no matter how many cells there are.
# Cells sharing a design row share mu, so every per-cell sum in glm_gp reduces to
# a per-group sum. The reduction is exact, not an approximation, and it applies to
# all three expensive blocks: the beta IRLS, the overdispersion MLE, and the
# standard-error pass.
#
# What has to be carried per gene, per group:
#   S_r = sum of counts        (beta IRLS, overdispersion objective)
#   L_r = sum of log(y + 1)    (reproduces estimate_betas_roughly exactly)
#   f_r = group size           (shared across genes)
# plus, per gene, a count-value histogram (which upstream builds itself, see
# make_table_if_small) and sum(y), sum(y^2), which give the moment dispersion
# estimate and nlminb's start value.
#
# ---------------------------------------------------------------------------
# 2. Why the dense matrices go away
#
# Mu is what makes stock glm_gp blow up: calculate_mu() materialises a dense
# n_genes x n_cells matrix, twice. Mu has no zeros -- every cell has a positive
# fitted mean -- so sparse storage cannot help it. What it has instead is massive
# column DUPLICATION, and that is what the collapse removes: at 3M cells there are
# only ~31k distinct design rows, so 23000 x 3e6 (552 GB) becomes 23000 x 31000
# (5.7 GB).
#
# That still leaves four n_genes x n_groups blocks (S, L, Ybar, Mu ~ 23 GB at
# that size), so two further things are done here:
#
#   * Ybar and Mu are never formed. fit_gp_weighted_s_cpp takes the raw group
#     sums S and divides by f per gene, and returns dot(mu_g, f)/n directly --
#     the only thing the full Mu was ever needed for was rowMeans2(Mu) feeding
#     the shrinkage trend.
#   * Genes are processed in chunks. The fit is independent per gene, so nothing
#     but per-gene scalars (disp_est, gene_means) and the p-column Beta/se have
#     to span the whole gene axis. Peak memory becomes O(gene_chunk x n_groups)
#     instead of O(n_genes x n_groups).
#
# Keeping S sparse was measured and rejected: the collapse itself fills S in
# (15.8% dense at 1.8 cells/group, 44% at 97 cells/group), so at the sizes that
# matter sparse storage saves ~1.5x and costs a great deal of complexity.
#
# The cost of chunking is one extra pass of row-subsetting over Y (the shrinkage
# is global, so beta pass 2 cannot start until every gene's disp_est exists).
# No IRLS arithmetic is repeated, so results are unchanged.

collapsed_glm_gp <- function(Y, design, gene_idx = NULL, gnames = NULL,
                             verbose = TRUE,
                             do_cox_reid_adjustment = TRUE,
                             max_iter_od = 200,
                             gene_chunk = NULL, chunk_budget_gb = 0.5,
                             threads = 1L, solver = c("qr", "chol")) {
  solver_id <- if (match.arg(solver) == "chol") 1L else 0L
  stopifnot(inherits(Y, "dgCMatrix"))
  # `gene_idx` selects which rows of Y to fit. Taking an index rather than a
  # pre-subset matrix keeps the caller from having to materialise one: at
  # min_pct = 0 the selection is ~99.7% of rows, so the subset would be a second
  # near-complete copy of the counts. Chunks are slices of `gsel`, so nothing
  # downstream changes -- output row k is gene gsel[k] either way.
  gsel <- if (is.null(gene_idx)) seq_len(nrow(Y)) else as.integer(gene_idx)
  n <- ncol(Y); ng <- length(gsel); p <- ncol(design)
  stopifnot(nrow(design) == n)
  tick <- function(msg, t0) if (verbose)
    message(sprintf("[collapsed] %-22s %6.1f s", msg, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # ---- 1. group identical design rows -------------------------------------
  # This block runs once per fit, and step 3 does ~100 leave-one-out fits after
  # the genome-wide one. Each LOO fit touches a single gene but re-groups every
  # cell, so at 455k cells the grouping dominated the whole step: 10.2 s per LOO
  # fit against 0.7 s of solver work.
  #
  # The old route pasted each design row into a string and matched on it, which
  # is ~n string allocations plus two hash passes. Sorting the rows numerically
  # and marking the boundaries is the same partition for a fraction of the cost.
  # Group ids are renumbered back to first-appearance order afterwards, so `grp`
  # is identical to what match(key, unique(key)) produced -- group order sets the
  # summation order in the kernels, and keeping it identical keeps the floating
  # point identical too.
  t0 <- Sys.time()
  if (n == 1L) {
    grp <- 1L
  } else {
    o <- do.call(order, c(asplit(design, 2), list(method = "radix")))
    brk <- c(TRUE, rowSums(design[o[-1L], , drop = FALSE] !=
                           design[o[-n], , drop = FALSE]) > 0)
    gs <- integer(n); gs[o] <- cumsum(brk)          # ids in sorted order
    # renumber sorted ids by where each group first appears in the original order
    remap <- integer(max(gs))
    remap[order(match(seq_len(max(gs)), gs))] <- seq_len(max(gs))
    grp <- remap[gs]
  }
  R <- max(grp)
  Xc <- design[match(seq_len(R), grp), , drop = FALSE]
  # group sizes: tabulate is O(n) into an R-vector; the old route built an
  # n-nonzero sparse indicator matrix and row-summed it.
  f <- as.numeric(tabulate(grp, R))
  grp0 <- as.integer(grp - 1L)
  # gene-major copy: chunking is then a contiguous column slice (cheap) and the
  # stats kernel gets one gene per thread with no write races. Y belongs to the
  # caller now, so this is a genuine doubling of the counts while both live --
  # the largest remaining term in step 3's peak.
  # Taking the gene names as an argument lets the caller hand in an UNNAMED
  # counts matrix. Seurat v5 stores the layer without dimnames, so naming it
  # duplicates it -- 17 GB at 305k cells -- purely to label rows that only
  # this line and the output ever read.
  gnames <- if (is.null(gnames)) rownames(Y)[gsel] else gnames[gsel]
  Yt <- Matrix::t(Y); gc(FALSE)
  if (verbose) message(sprintf("[collapsed] %d cells -> %d design rows (%.1fx)", n, R, n / R))

  if (is.null(gene_chunk))
    gene_chunk <- max(1L, min(ng, as.integer(chunk_budget_gb * 1024^3 / (R * 8))))
  chunks <- split(seq_len(ng), ceiling(seq_len(ng) / gene_chunk))
  if (verbose) message(sprintf("[collapsed] %d genes in %d chunk(s) of <= %d (%.2f GB/chunk block)",
                               ng, length(chunks), gene_chunk, gene_chunk * R * 8 / 1024^3))

  # per-gene scalars and the p-column results are all that span the full gene axis
  Beta      <- matrix(NA_real_, ng, p)
  disp_est  <- numeric(ng)
  gene_mean <- numeric(ng)

  # ---- 2-5. pass A, per gene chunk ----------------------------------------
  t0 <- Sys.time(); t_stats <- 0; t_beta1 <- 0; t_od <- 0
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    ta <- Sys.time()
    Ytc <- Yt[, gsel[idx], drop = FALSE]
    st <- chunk_stats_cpp(Ytc@p, Ytc@i, Ytc@x, grp0, length(idx), R, threads)
    t_stats <- t_stats + as.numeric(difftime(Sys.time(), ta, units = "secs"))

    # rough dispersion (moment) and exact rough beta. estimate_betas_roughly
    # solves the LS problem X b = log(y + 1) over all cells; its normal equations
    # X'X b = X'z are reproduced exactly by the group sums, because
    # X'X = Xc' diag(f) Xc and X'z = Xc' L.
    bm <- st$sumy / n
    bv <- (st$sumy2 - n * bm^2) / (n - 1)                    # rowVars, n-1 denom
    disp_init <- (bv - bm) / bm^2                            # xim == 1 at zero offset
    disp_init[is.na(disp_init) | disp_init < 0] <- 0
    beta_init <- t(solve(crossprod(Xc * f, Xc), crossprod(Xc, t(st$L))))
    st$L <- NULL

    ta <- Sys.time()
    r1 <- fit_gp_weighted_s_cpp(st$S, Xc, f, disp_init,
                                dev_const_hist_cpp(st$Hp, st$Hi, st$Hx, disp_init, threads), beta_init, n,
                                nthreads = threads, solver = solver_id)
    Beta[idx, ] <- r1$Beta
    gene_mean[idx] <- r1$mu_mean                             # == rowMeans2(dense Mu)
    t_beta1 <- t_beta1 + as.numeric(difftime(Sys.time(), ta, units = "secs"))

    # The overdispersion MLE is the one block that cannot go under OpenMP: it calls
    # R's nlminb per gene, so it has to stay on the R side. Forking is safe here in a
    # way the old --threads regression was not -- that forked around glm_gp, and each
    # worker densified its own n_genes x n_cells Mu, so peak RAM scaled with workers.
    # Here a worker only reads the chunk's S (already allocated, shared copy-on-write)
    # and returns one double per gene, so peak RAM is flat in `threads`. No RNG is
    # involved and results are collected by index, so the output is deterministic.
    ta <- Sys.time()
    one_gene <- function(k) {
      mu <- pmax(as.numeric(exp(Xc %*% Beta[idx[k], ])), 1e-6)
      od_mle_one(st$S[k, ], f, mu, Xc, hist_row(st, k, n, st$nz[k]),
                 st$sumy[k], st$sumy2[k], n, do_cox_reid_adjustment, max_iter_od)
    }
    disp_est[idx] <- if (threads > 1L) {
      unlist(parallel::mclapply(seq_along(idx), one_gene, mc.cores = threads,
                                mc.preschedule = TRUE), use.names = FALSE)
    } else {
      vapply(seq_along(idx), one_gene, numeric(1))
    }
    t_od <- t_od + as.numeric(difftime(Sys.time(), ta, units = "secs"))
    rm(Ytc, st, r1); gc(FALSE)
    if (verbose) message(sprintf("[collapsed]   pass A chunk %d/%d", ci, length(chunks)))
  }
  if (verbose) message(sprintf("[collapsed] pass A: stats %.1fs  beta1 %.1fs  overdisp %.1fs",
                               t_stats, t_beta1, t_od))

  # ---- 6. shrinkage (unchanged upstream code, gene-level so already cheap) ---
  t0 <- Sys.time()
  shr <- glmGamPoi:::overdispersion_shrinkage(
    disp_est, gene_means = gene_mean, df = n - p,
    ql_disp_trend = length(disp_est) >= 100,
    npoints = max(0.1 * length(disp_est), 100), verbose = FALSE)
  disp_latest <- shr$dispersion_trend
  tick("shrinkage", t0)

  # ---- 7. pass B: beta pass 2, then standard errors -------------------------
  # predict.glmGamPoi computes the se from fit$overdispersions, which is the
  # PRE-shrinkage disp_est -- not the dispersion_trend the beta refit uses.
  # Reproduce that split rather than the (more self-consistent) single-theta se,
  # otherwise the standard errors come out too small and every marginal gene
  # gains significance.
  t0 <- Sys.time()
  se  <- matrix(NA_real_, ng, p)
  dev <- numeric(ng); iters <- numeric(ng)
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    Ytc <- Yt[, gsel[idx], drop = FALSE]
    stb <- chunk_stats_cpp(Ytc@p, Ytc@i, Ytc@x, grp0, length(idx), R, threads)
    S  <- stb$S
    r2 <- fit_gp_weighted_s_cpp(S, Xc, f, disp_latest[idx],
                                dev_const_hist_cpp(stb$Hp, stb$Hi, stb$Hx, disp_latest[idx], threads), Beta[idx, , drop = FALSE], n,
                                nthreads = threads, solver = solver_id)
    Beta[idx, ] <- r2$Beta
    dev[idx] <- r2$deviance; iters[idx] <- r2$iter
    se[idx, ] <- gp_se_weighted_omp_cpp(r2$Beta, Xc, f, disp_est[idx], nthreads = threads)
    rm(Ytc, stb, S, r2); gc(FALSE)
  }
  tick("pass B (beta2 + se)", t0)

  dimnames(Beta) <- list(gnames, colnames(design))
  dimnames(se) <- dimnames(Beta)
  list(Beta = Beta, se = se, deviances = dev, iterations = iters,
       overdispersions = disp_est, overdispersion_shrinkage_list = shr,
       gene_means = gene_mean, design_collapsed = Xc, group_sizes = f,
       gene_chunk = gene_chunk, n_groups = R, threads = threads, solver = match.arg(solver))
}

# Per-gene, per-group sufficient statistics for one chunk of genes.
chunk_stats <- function(Yc, Mt, f, n) {
  ngc <- nrow(Yc)
  S <- as.matrix(Yc %*% Mt)                                  # sum(y)       per group
  Ylog <- Yc; Ylog@x <- log1p(Ylog@x)
  L <- as.matrix(Ylog %*% Mt); rm(Ylog)                      # sum(log1p y) per group
  sumy  <- Matrix::rowSums(Yc)
  sumy2 <- Matrix::rowSums({ Y2 <- Yc; Y2@x <- Y2@x^2; Y2 })
  # count-value histogram per gene: H[g, v] = #cells with count exactly v.
  # sparseMatrix sums duplicate (i, j), so this is one pass over the nonzeros.
  # Transposed (counts x genes) so a gene is one dgCMatrix column -> O(nnz_g).
  maxc <- if (length(Yc@x)) max(Yc@x) else 1
  Ht <- Matrix::sparseMatrix(j = Yc@i + 1L, i = as.integer(Yc@x), x = 1,
                             dims = c(maxc, ngc))
  list(S = S, L = L, sumy = sumy, sumy2 = sumy2, Ht = Ht,
       nz = tabulate(Yc@i + 1L, nbins = ngc))
}

# mu-independent part of the collapsed Gamma-Poisson deviance, per gene.
dev_const_w <- function(Y, theta) {
  A <- Matrix::rowSums({ Z <- Y; Z@x <- 2 * Z@x * log(Z@x); Z })
  th <- theta[Y@i + 1L]
  B <- Matrix::rowSums({ Z <- Y; Z@x <- 2 * (Z@x + 1 / th) * log1p(Z@x * th); Z })
  out <- A - B
  pois <- theta < 1e-6
  if (any(pois)) out[pois] <- A[pois] - 2 * Matrix::rowSums(Y[pois, , drop = FALSE])
  out
}

# count-value table for one gene, in the (unique_counts, frequencies) form that
# glmGamPoi's loglik/score/hessian already accept. The zero bin is included.
hist_row <- function(st, g, n, nz) {
  if (st$Hp[g + 1L] == st$Hp[g]) return(list(counts = 0, freq = n))
  k <- (st$Hp[g] + 1L):st$Hp[g + 1L]
  list(counts = c(0, st$Hi[k] + 1), freq = c(n - nz, st$Hx[k]))
}

# Mirrors glmGamPoi:::conventional_overdispersion_mle exactly -- same guards, same
# start value, same nlminb calls and fallbacks -- but evaluated on group sums.
od_mle_one <- function(S, f, mu, X, tab, sumy, sumy2, n, do_cr, max_iter) {
  if (sumy == 0) return(0)
  uc <- tab$counts; cf <- tab$freq
  if (gp_od_score_w(S, f, mu, log(1e-8), X, do_cr, uc, cf, n) < 0) return(0)
  m <- sumy / n
  sv <- ((sumy2 - n * m^2) / (n - 1) - m) / m^2
  if (is.na(sv) || sv <= 0) sv <- 0.5
  obj  <- function(lt, cr = do_cr) -gp_od_loglik_w(S, f, mu, lt, X, cr, uc, cf, n)
  grad <- function(lt, cr = do_cr) -gp_od_score_w(S, f, mu, lt, X, cr, uc, cf, n)
  hess <- function(lt) matrix(-gp_od_deriv_score_w(S, f, mu, lt, X, do_cr, uc, cf, n), 1, 1)
  ctl <- list(iter.max = max_iter); lo <- log(1e-16); hi <- log(1e16)
  res <- nlminb(log(sv), obj, grad, hess, lower = lo, upper = hi, control = ctl)
  if (res$convergence != 0)
    res <- nlminb(log(sv), obj, grad, lower = lo, upper = hi, control = ctl)
  if (res$convergence != 0)
    res <- nlminb(log(sv), function(lt) obj(lt, FALSE), function(lt) grad(lt, FALSE),
                  lower = lo, upper = hi, control = ctl)
  exp(res$par)
}
