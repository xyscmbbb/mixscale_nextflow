// Frequency-weighted Gamma-Poisson IRLS.
//
// This is glmGamPoi's fitBeta_fisher_scoring specialised for the mixscale step-3
// design (~ 1 + weight + log_ct, size_factors = FALSE so the offset is identically
// zero) and extended with per-column frequency weights f.
//
// Why frequency weights: NT cells all have weight == 0 exactly, so they differ
// only in log_ct = log1p(nCount). Cells sharing a design row share mu, so every
// per-cell sum in the IRLS reduces to a per-group statistic. At 3M NT cells the
// number of distinct design rows saturates near 31k, a ~96x column reduction.
//
// The collapse is exact, not an approximation:
//   sum_{i in g} dev(y_i, mu, theta)
//     = -2*S_g*log(mu) + 2*(S_g + f_g/theta)*log1p(mu*theta) + D_g
// with S_g = sum of counts in the group and D_g = sum_i [2*y_i*log(y_i)
//   - 2*(y_i + 1/theta)*log1p(y_i*theta)] independent of mu, so it is computed
// once per gene and passed in as dev_const.
//
// Everything else -- the QR step, the step-halving line search, the convergence
// test, the tolerances -- is kept identical to upstream so the fits agree.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Weighted Fisher scoring step. Identical to fisher_scoring_qr_step except that
// the working weight carries the group size f and counts are group means.
inline arma::vec fs_qr_step_w(const arma::mat& X, const arma::vec& ybar,
                              const arma::vec& mu, double theta, const arma::vec& f) {
  arma::mat q, r;
  arma::vec w_vec = f % (mu / (1.0 + theta * mu));
  arma::vec w_sqrt = arma::sqrt(w_vec);
  arma::mat Xw = X.each_col() % w_sqrt;
  arma::qr_econ(q, r, Xw);
  arma::vec score = (q.each_col() % w_sqrt).t() * ((ybar - mu) / mu);
  return arma::solve(arma::trimatu(r), score);
}

// Same step, via the normal equations instead of a QR.
//
// fs_qr_step_w allocates six 31000-row temporaries per IRLS iteration (w, sqrt(w),
// Xw, Q, R, the score vector). With ~8 iterations x 38000 genes that allocation
// churn, not the arithmetic, is what the profile is made of.
//
// The algebra is the same, not an approximation. qr_econ gives Xw = QR with
// Xw = diag(sqrt w) X, so R'R = X'WX, and
//   solve(trimatu(R), (diag(sqrt w) Q)' v) = R^-1 R^-T X'W v = (X'WX)^-1 X'W v
// which is exactly what solving the normal equations returns. QR is the more
// stable route in general; at p = 3 with columns (1, weight in [0,1], log_ct ~ 8)
// the Gram matrix is nowhere near ill-conditioned, and the parity run is the check.
//
// Xt is the design TRANSPOSED (p x n_groups) so one design row is contiguous --
// the accumulation is a single streaming pass with no temporaries at all.
inline void accum_normal_eq(const arma::mat& Xt, const arma::vec& ybar, const arma::vec& mu,
                            double theta, const arma::vec& f, arma::mat& A, arma::vec& s) {
  const arma::uword n = Xt.n_cols, p = Xt.n_rows;
  A.zeros(); s.zeros();
  const double* xp = Xt.memptr();
  for (arma::uword i = 0; i < n; i++, xp += p) {
    const double m = mu(i);
    const double w = f(i) * m / (1.0 + theta * m);
    const double r = w * (ybar(i) - m) / m;
    for (arma::uword a = 0; a < p; a++) {
      s(a) += r * xp[a];
      const double wxa = w * xp[a];
      for (arma::uword b = a; b < p; b++) A(a, b) += wxa * xp[b];
    }
  }
  for (arma::uword a = 0; a < p; a++)
    for (arma::uword b = 0; b < a; b++) A(a, b) = A(b, a);
}

inline arma::vec fs_chol_step_w(const arma::mat& Xt, const arma::vec& ybar,
                                const arma::vec& mu, double theta, const arma::vec& f,
                                arma::mat& A, arma::vec& s) {
  accum_normal_eq(Xt, ybar, mu, theta, f, A, s);
  arma::vec step;
  if (!arma::solve(step, A, s, arma::solve_opts::likely_sympd)) step.zeros(A.n_rows);
  return step;
}

// Collapsed Gamma-Poisson deviance. dev_const carries the mu-independent part.
inline double gp_dev_w(const arma::vec& ybar, const arma::vec& f, const arma::vec& mu,
                       double theta, double dev_const) {
  double dev = 0.0;
  if (theta < 1e-6) {
    // Poisson limit: sum 2*(y*log(y/mu) - (y - mu)); the y*log(y) part is in dev_const.
    for (arma::uword i = 0; i < mu.n_elem; i++)
      dev += 2.0 * f(i) * (mu(i) - ybar(i) * std::log(mu(i)));
  } else {
    for (arma::uword i = 0; i < mu.n_elem; i++)
      dev += 2.0 * f(i) * (-ybar(i) * std::log(mu(i))
                           + (ybar(i) + 1.0 / theta) * std::log1p(mu(i) * theta));
  }
  return dev + dev_const;
}

// Step-halving line search. Mirrors decrease_deviance in beta_estimation.cpp.
inline double decrease_deviance_w(arma::vec& beta, arma::vec& mu, const arma::vec& step,
                                  const arma::mat& X, const arma::vec& ybar, const arma::vec& f,
                                  double theta, double dev_const, double dev_old,
                                  double tolerance, double max_rel_mu_change) {
  arma::vec mu_old = mu;
  double speed = 1.0;
  beta = beta + step;
  int line_iter = 0;
  double dev = 0.0;
  while (true) {
    mu = arma::exp(X * beta);
    dev = gp_dev_w(ybar, f, mu, theta, dev_const);
    double conv = std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1);
    double mu_rel_change = arma::max(mu / mu_old);
    if ((dev < dev_old && mu_rel_change < max_rel_mu_change) || conv < tolerance) break;
    if (line_iter >= 100) { dev = std::numeric_limits<double>::quiet_NaN(); break; }
    speed = speed / 2.0;
    beta = beta - step * speed;
    line_iter++;
  }
  return dev;
}

// [[Rcpp::export]]
List fit_gp_weighted_cpp(const arma::mat& Ybar, const arma::mat& X, const arma::vec& f,
                         const arma::vec& thetas, const arma::vec& dev_const,
                         const arma::mat& beta_init,
                         double tolerance = 1e-8, double max_rel_mu_change = 1e5,
                         int max_iter = 1000) {
  int n_genes = Ybar.n_rows;
  int p = X.n_cols;
  if ((int)Ybar.n_cols != (int)X.n_rows) stop("Ybar columns must match X rows");
  if ((int)f.n_elem != (int)X.n_rows) stop("f must have one entry per design row");

  arma::mat beta_mat = beta_init;
  NumericVector iterations(n_genes), deviance(n_genes);
  // sqrt(diag(solve(X'WX))) at the optimum -- this is exactly what
  // predict(se.fit = TRUE, newdata = diag(p)) reduces to when the offset is zero.
  arma::mat se_mat(n_genes, p, arma::fill::value(NA_REAL));

  for (int g = 0; g < n_genes; g++) {
    if (g % 100 == 0) checkUserInterrupt();
    arma::vec ybar = Ybar.row(g).t();
    arma::vec beta = beta_mat.row(g).t();
    double theta = thetas(g);

    if (beta.has_nan() || !arma::is_finite(beta) || ISNAN(theta)) {
      beta.fill(NA_REAL); beta_mat.row(g) = beta.t();
      iterations(g) = 0; deviance(g) = NA_REAL; continue;
    }
    arma::vec mu = arma::exp(X * beta);
    double dev_old = gp_dev_w(ybar, f, mu, theta, dev_const(g));

    for (int t = 0; t < max_iter; t++) {
      iterations(g)++;
      arma::vec step = fs_qr_step_w(X, ybar, mu, theta, f);
      double dev = decrease_deviance_w(beta, mu, step, X, ybar, f, theta,
                                       dev_const(g), dev_old, tolerance, max_rel_mu_change);
      double conv = std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1);
      dev_old = dev;
      if (std::isnan(conv)) { beta.fill(NA_REAL); iterations(g) = max_iter; break; }
      if (conv < tolerance) break;
    }
    beta_mat.row(g) = beta.t();
    deviance(g) = dev_old;

    if (beta.is_finite()) {
      arma::vec w_vec = f % (mu / (1.0 + theta * mu));
      arma::mat XtWX = X.t() * (X.each_col() % w_vec);
      arma::mat inv;
      if (arma::inv_sympd(inv, XtWX)) se_mat.row(g) = arma::sqrt(inv.diag()).t();
    }
  }
  return List::create(Named("Beta", beta_mat), Named("se", se_mat),
                      Named("iter", iterations), Named("deviance", deviance));
}

// Standard errors, as a separate pass.
//
// This has to be its own function because glm_gp and predict.glmGamPoi disagree
// about which overdispersion to use: the second beta pass is fitted with the
// SHRUNKEN dispersion_trend (glm_gp_impl.R:117-139), while predict() computes
// the standard error with the PRE-shrinkage disp_est (predict.R:205,
// `disp <- object$overdispersions[gene_idx]`). Folding the se into the fit would
// silently use the shrunken value for both and report smaller standard errors
// than stock -- so theta_se is passed in separately.
//
// The formula is predict.glmGamPoi's: with newdata = diag(p) it reduces to
// sqrt(diag(solve(X'WX))), W = f * mu/(1 + mu*theta_se), at zero offset.
// [[Rcpp::export]]
arma::mat gp_se_weighted_cpp(const arma::mat& Beta, const arma::mat& X, const arma::vec& f,
                             const arma::vec& theta_se) {
  int n_genes = Beta.n_rows, p = X.n_cols;
  arma::mat se_mat(n_genes, p, arma::fill::value(NA_REAL));
  for (int g = 0; g < n_genes; g++) {
    if (g % 100 == 0) checkUserInterrupt();
    arma::vec beta = Beta.row(g).t();
    if (!beta.is_finite() || ISNAN(theta_se(g))) continue;
    arma::vec mu = arma::exp(X * beta);
    arma::vec w = f % (mu / (1.0 + mu * theta_se(g)));
    arma::mat XtWX = X.t() * (X.each_col() % w);
    arma::mat inv;
    if (arma::inv_sympd(inv, XtWX)) se_mat.row(g) = arma::sqrt(inv.diag()).t();
  }
  return se_mat;
}

// Chunk-friendly variant of the beta fit.
//
// Two allocations are removed relative to fit_gp_weighted_cpp:
//   * Ybar (n_genes x n_groups) -- the caller passes the raw group sums S and
//     the division by f happens per gene, inside the loop.
//   * Mu (n_genes x n_groups) -- the only thing the caller needed the full Mu
//     for was rowMeans2(Mu) to feed the shrinkage trend, so that reduction is
//     returned directly as mu_mean(g) = dot(mu_g, f) / n_total.
// At 23k genes x 31k groups each of those is 5.7 GB, so dropping both is 11.4 GB
// off the peak before any chunking is even applied.
//
// The arithmetic is otherwise identical to fit_gp_weighted_cpp -- same step,
// same line search, same convergence test -- so results are unchanged.
//
// nthreads > 1 runs the gene loop under OpenMP. Genes are fully independent (each
// touches only its own row of S/beta_mat and its own scalars), there is no reduction
// and no shared mutable state, so the results are bit-identical to nthreads = 1 --
// unlike the old R-level --threads fork, which diverged. Nothing in the loop calls
// the R API: checkUserInterrupt is gone and the output vectors are std::vector.
// Pin the BLAS to one thread (OPENBLAS_NUM_THREADS=1) or the LAPACK call inside
// qr_econ will oversubscribe against this loop.
// [[Rcpp::export]]
List fit_gp_weighted_s_cpp(const arma::mat& S, const arma::mat& X, const arma::vec& f,
                           const arma::vec& thetas, const arma::vec& dev_const,
                           const arma::mat& beta_init, double n_total,
                           double tolerance = 1e-8, double max_rel_mu_change = 1e5,
                           int max_iter = 1000, int nthreads = 1, int solver = 0) {
  int n_genes = S.n_rows;
  int p = X.n_cols;
  if ((int)S.n_cols != (int)X.n_rows) stop("S columns must match X rows");
  if ((int)f.n_elem != (int)X.n_rows) stop("f must have one entry per design row");
  // Transposed once for the whole chunk: the normal-equations step wants design
  // rows contiguous, and this is the only place the layout is decided.
  const arma::mat Xt = X.t();

  arma::mat beta_mat = beta_init;
  // Plain std::vector, not Rcpp::NumericVector: nothing in the parallel region may
  // touch the R API, and that includes Rcpp proxy assignment.
  std::vector<double> iterations(n_genes, 0.0), deviance(n_genes, 0.0), mu_mean(n_genes, 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(nthreads)
#endif
  for (int g = 0; g < n_genes; g++) {
    // Per-iteration workspace for the normal-equations solver, allocated once per
    // gene (per thread) rather than once per IRLS iteration.
    arma::mat A(p, p); arma::vec sc(p);
    arma::vec ybar = S.row(g).t() / f;
    arma::vec beta = beta_mat.row(g).t();
    double theta = thetas(g);

    if (beta.has_nan() || !arma::is_finite(beta) || ISNAN(theta)) {
      beta.fill(NA_REAL); beta_mat.row(g) = beta.t();
      iterations[g] = 0; deviance[g] = NA_REAL; mu_mean[g] = NA_REAL; continue;
    }
    arma::vec mu = arma::exp(X * beta);
    double dev_old = gp_dev_w(ybar, f, mu, theta, dev_const(g));

    for (int t = 0; t < max_iter; t++) {
      iterations[g]++;
      arma::vec step = (solver == 1) ? fs_chol_step_w(Xt, ybar, mu, theta, f, A, sc)
                                     : fs_qr_step_w(X, ybar, mu, theta, f);
      double dev = decrease_deviance_w(beta, mu, step, X, ybar, f, theta,
                                       dev_const(g), dev_old, tolerance, max_rel_mu_change);
      double conv = std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1);
      dev_old = dev;
      if (std::isnan(conv)) { beta.fill(NA_REAL); iterations[g] = max_iter; break; }
      if (conv < tolerance) break;
    }
    beta_mat.row(g) = beta.t();
    deviance[g] = dev_old;
    mu_mean[g] = beta.is_finite() ? arma::dot(mu, f) / n_total : NA_REAL;
  }
  return List::create(Named("Beta", beta_mat), Named("iter", wrap(iterations)),
                      Named("deviance", wrap(deviance)), Named("mu_mean", wrap(mu_mean)));
}

// [[Rcpp::export]]
arma::mat gp_se_weighted_omp_cpp(const arma::mat& Beta, const arma::mat& X, const arma::vec& f,
                                 const arma::vec& theta_se, int nthreads = 1) {
  int n_genes = Beta.n_rows, p = X.n_cols;
  arma::mat se_mat(n_genes, p, arma::fill::value(NA_REAL));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(nthreads)
#endif
  for (int g = 0; g < n_genes; g++) {
    arma::vec beta = Beta.row(g).t();
    if (!beta.is_finite() || ISNAN(theta_se(g))) continue;
    arma::vec mu = arma::exp(X * beta);
    arma::vec w = f % (mu / (1.0 + mu * theta_se(g)));
    arma::mat XtWX = X.t() * (X.each_col() % w);
    arma::mat inv;
    if (arma::inv_sympd(inv, XtWX)) se_mat.row(g) = arma::sqrt(inv.diag()).t();
  }
  return se_mat;
}

// ---------------------------------------------------------------------------
// Fused sufficient statistics.
//
// The R version made five separate passes over the chunk (S, L, sumy, sumy2,
// histogram), two of which copied the whole @x vector.  This does all of them
// in one pass, parallel over genes.
//
// The chunk arrives gene-major (CSC, columns = genes) so each thread owns a
// disjoint set of columns and there are no races.  Per-group accumulators are
// kept in thread-local contiguous buffers and only the touched entries are
// written out / reset, so the cost is O(nnz) not O(n_genes * n_groups).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List chunk_stats_cpp(IntegerVector Yp, IntegerVector Yi, NumericVector Yx,
                     IntegerVector grp, int ngc, int R, int nthreads = 1) {
  arma::mat S(ngc, R, arma::fill::zeros);
  arma::mat L(ngc, R, arma::fill::zeros);
  arma::vec sumy(ngc, arma::fill::zeros), sumy2(ngc, arma::fill::zeros);
  IntegerVector nz(ngc);

  // per-gene (value, count) histogram of the nonzeros, assembled into CSC below
  std::vector< std::vector<int> >    hval(ngc);
  std::vector< std::vector<double> > hcnt(ngc);

  const int* pp = Yp.begin();
  const int* ii = Yi.begin();
  const double* xx = Yx.begin();
  const int* gg = grp.begin();

  int maxc = 1;
  for (R_xlen_t k = 0; k < Yx.size(); k++) {
    int v = (int) xx[k];
    if (v > maxc) maxc = v;
  }

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    std::vector<double> accS(R, 0.0), accL(R, 0.0);
    std::vector<int> touched;
    std::vector<int> cnt(maxc + 1, 0);
    std::vector<int> seen;

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
    for (int g = 0; g < ngc; g++) {
      const int lo = pp[g], hi = pp[g + 1];
      nz[g] = hi - lo;
      touched.clear();
      seen.clear();
      double sy = 0.0, sy2 = 0.0;
      for (int k = lo; k < hi; k++) {
        const double y = xx[k];
        const int r = gg[ii[k]];
        if (accS[r] == 0.0 && accL[r] == 0.0) touched.push_back(r);
        accS[r] += y;
        accL[r] += std::log1p(y);
        sy  += y;
        sy2 += y * y;
        const int v = (int) y;
        if (cnt[v] == 0) seen.push_back(v);
        cnt[v]++;
      }
      sumy(g) = sy; sumy2(g) = sy2;
      for (size_t t = 0; t < touched.size(); t++) {
        const int r = touched[t];
        S(g, r) = accS[r]; L(g, r) = accL[r];
        accS[r] = 0.0; accL[r] = 0.0;
      }
      std::sort(seen.begin(), seen.end());
      hval[g].reserve(seen.size());
      hcnt[g].reserve(seen.size());
      for (size_t t = 0; t < seen.size(); t++) {
        hval[g].push_back(seen[t]);
        hcnt[g].push_back((double) cnt[seen[t]]);
        cnt[seen[t]] = 0;
      }
    }
  }

  // assemble the histogram as a (maxc x ngc) dgCMatrix: one gene = one column
  IntegerVector Hp(ngc + 1);
  R_xlen_t tot = 0;
  for (int g = 0; g < ngc; g++) { Hp[g] = (int) tot; tot += hval[g].size(); }
  Hp[ngc] = (int) tot;
  IntegerVector Hi(tot);
  NumericVector Hx(tot);
  R_xlen_t o = 0;
  for (int g = 0; g < ngc; g++)
    for (size_t t = 0; t < hval[g].size(); t++, o++) {
      Hi[o] = hval[g][t] - 1;      // value v -> row v-1
      Hx[o] = hcnt[g][t];
    }

  return List::create(Named("S", S), Named("L", L),
                      Named("sumy", sumy), Named("sumy2", sumy2),
                      Named("nz", nz), Named("maxc", maxc),
                      Named("Hp", Hp), Named("Hi", Hi), Named("Hx", Hx));
}

// ---------------------------------------------------------------------------
// Deviance constant from the count histogram.
//
// dev_const is a sum over cells of a function of y alone (given the gene's
// theta), and cells with y == 0 contribute nothing, so the histogram of the
// nonzeros determines it exactly -- no second pass over the data.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec dev_const_hist_cpp(IntegerVector Hp, IntegerVector Hi, NumericVector Hx,
                             const arma::vec& theta, int nthreads = 1) {
  const int ngc = Hp.size() - 1;
  arma::vec out(ngc, arma::fill::zeros);
  const int* pp = Hp.begin();
  const int* ii = Hi.begin();
  const double* xx = Hx.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for (int g = 0; g < ngc; g++) {
    const double th = theta(g);
    const bool pois = th < 1e-6;
    double acc = 0.0;
    for (int k = pp[g]; k < pp[g + 1]; k++) {
      const double y = (double)(ii[k] + 1), fr = xx[k];
      const double A = 2.0 * y * std::log(y);
      const double B = pois ? 2.0 * y
                            : 2.0 * (y + 1.0 / th) * std::log1p(y * th);
      acc += fr * (A - B);
    }
    out(g) = acc;
  }
  return out;
}
