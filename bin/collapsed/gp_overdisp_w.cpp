// Frequency-weighted Gamma-Poisson overdispersion objective.
//
// These are glmGamPoi's conventional_loglikelihood_fast / _score_function_fast /
// _deriv_score_function_fast (src/overdispersion.cpp), rewritten so the per-cell
// sums run over design-row GROUPS instead of cells.
//
// Every per-cell loop in the originals has the form
//     sum_i g(mu_i) + sum_i y_i * h(mu_i)
// and cells sharing a design row share mu, so it collapses exactly to
//     sum_r f_r * g(mu_r) + sum_r S_r * h(mu_r)
// with f_r the group size and S_r the group's count sum. The lgamma/digamma/
// trigamma terms already depend only on the count-value table, which upstream
// builds itself (make_table_if_small), so those are passed through unchanged.
//
// Branch structure, tolerances, the cr_correction_factor and the small-mu_theta
// Taylor guards are kept byte-for-byte identical to upstream so the objective
// surface -- and therefore nlminb's path over it -- is the same.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static const double cr_correction_factor = 0.99;

// Cox-Reid: -0.5 * log det(X' diag(f*w) X) * cr_correction_factor
static inline double cr_logdet_w(const arma::mat& X, const arma::vec& f,
                                 const arma::vec& mu, double theta) {
  arma::vec w = f / (1.0 / mu + theta);
  arma::mat b = X.t() * (X.each_col() % w);
  arma::mat L, U, P;
  arma::lu(L, U, P, b);
  double ld = arma::sum(arma::log(arma::diagvec(L)));
  arma::vec u_diag = arma::diagvec(U);
  for (double e : u_diag) ld += e < 1e-50 ? std::log(1e-50) : std::log(e);
  return -0.5 * ld * cr_correction_factor;
}

// [[Rcpp::export]]
double gp_od_loglik_w(const arma::vec& S, const arma::vec& f, const arma::vec& mu,
                      double log_theta, const arma::mat& X, bool do_cr_adj,
                      const arma::vec& uniq_counts, const arma::vec& count_freq,
                      double n_total) {
  double theta = std::exp(log_theta);
  double cr_term = do_cr_adj ? cr_logdet_w(X, f, mu, theta) : 0.0;
  double tn1 = 1.0 / theta;

  double lgamma_term = 0.0;
  for (arma::uword k = 0; k < uniq_counts.n_elem; k++)
    lgamma_term += count_freq(k) * std::lgamma(uniq_counts(k) + tn1);
  lgamma_term -= n_total * std::lgamma(tn1);

  double ll_part = 0.0;
  for (arma::uword r = 0; r < mu.n_elem; r++)
    ll_part += (-S(r) - f(r) * tn1) * std::log(mu(r) + tn1);
  ll_part -= n_total * tn1 * log_theta;

  return lgamma_term + ll_part + cr_term;
}

// [[Rcpp::export]]
double gp_od_score_w(const arma::vec& S, const arma::vec& f, const arma::vec& mu,
                     double log_theta, const arma::mat& X, bool do_cr_adj,
                     const arma::vec& uniq_counts, const arma::vec& count_freq,
                     double n_total) {
  double theta = std::exp(log_theta);
  double tn1 = 1.0 / theta;

  double cr_term = 0.0;
  if (do_cr_adj) {
    arma::vec w = f / (1.0 / mu + theta);
    arma::vec w1 = 1.0 / (1.0 / mu + theta);   // dw = -w1^2, carried with f below
    arma::vec dw = -f % w1 % w1;
    arma::mat b = X.t() * (X.each_col() % w);
    arma::mat db = X.t() * (X.each_col() % dw);
    arma::mat b_inv = arma::inv_sympd(b + arma::eye(b.n_rows, b.n_cols) * 1e-6);
    cr_term = -0.5 * arma::trace(b_inv * db) * cr_correction_factor;
  }

  // digamma block: identical to upstream, it only ever used the count table
  double digamma_term = 0.0, max_y = 0.0, sum_y = 0.0, sum_prod_y = 0.0;
  for (arma::uword k = 0; k < uniq_counts.n_elem; k++) {
    digamma_term += count_freq(k) * Rf_digamma(uniq_counts(k) + tn1);
    sum_y += count_freq(k) * uniq_counts(k);
    sum_prod_y += count_freq(k) * (uniq_counts(k) - 1) * uniq_counts(k);
    max_y = std::max(max_y, uniq_counts(k));
  }
  double corr = tn1 > 1e5 ? sum_prod_y / (2 * tn1) : 0.0;
  if (max_y * 1e6 < tn1) {
    digamma_term = sum_y - corr;
  } else {
    digamma_term -= n_total * Rf_digamma(tn1);
    digamma_term *= tn1;
    digamma_term = std::min(digamma_term, sum_y - corr);
  }

  double ll_part = 0.0;
  for (arma::uword r = 0; r < mu.n_elem; r++) {
    double mu_theta = mu(r) * theta;
    if (mu_theta < 1e-10) {
      ll_part += f(r) * (mu_theta * mu_theta * (1.0 / (1.0 + mu_theta) - 0.5));
    } else if (mu_theta < 1e-4) {
      double inv = 1.0 / (1.0 + mu_theta);
      double upper_bound = mu_theta * mu_theta * inv;
      double lower_bound = mu_theta * mu_theta * (inv - 0.5);
      double suggest = std::log1p(mu_theta) - mu(r) / (mu(r) + tn1);
      ll_part += f(r) * std::max(std::min(suggest, upper_bound), lower_bound);
    } else {
      ll_part += f(r) * (std::log(1.0 + mu_theta) - mu(r) / (mu(r) + tn1));
    }
    ll_part += S(r) / (mu(r) + tn1);
  }
  ll_part *= tn1;
  return ll_part - digamma_term + cr_term * theta;
}

// [[Rcpp::export]]
double gp_od_deriv_score_w(const arma::vec& S, const arma::vec& f, const arma::vec& mu,
                           double log_theta, const arma::mat& X, bool do_cr_adj,
                           const arma::vec& uniq_counts, const arma::vec& count_freq,
                           double n_total) {
  double theta = std::exp(log_theta);
  double cr_term = 0.0, cr_term2 = 0.0;
  if (do_cr_adj) {
    arma::vec w1 = 1.0 / (1.0 / mu + theta);
    arma::vec w = f % w1;
    arma::vec dw = -f % w1 % w1;
    arma::vec d2w = 2.0 * f % w1 % w1 % w1;
    arma::mat b = X.t() * (X.each_col() % w);
    arma::mat db = X.t() * (X.each_col() % dw);
    arma::mat d2b = X.t() * (X.each_col() % d2w);
    arma::mat b_inv = arma::inv_sympd(b + arma::eye(b.n_rows, b.n_cols) * 1e-6);
    arma::mat d_i_db = b_inv * db;
    double ddetb = arma::trace(d_i_db);
    double d2detb = ddetb * ddetb - arma::trace(d_i_db * d_i_db) + arma::trace(b_inv * d2b);
    cr_term = (0.5 * ddetb * ddetb - 0.5 * d2detb) * cr_correction_factor;
    cr_term2 = -0.5 * ddetb * cr_correction_factor;
  }

  double tn1 = 1.0 / theta;
  double tn2 = tn1 * tn1;
  double digamma_term = 0.0, trigamma_term = 0.0;
  for (arma::uword k = 0; k < uniq_counts.n_elem; k++) {
    digamma_term += count_freq(k) * Rf_digamma(uniq_counts(k) + tn1);
    trigamma_term += count_freq(k) * Rf_trigamma(uniq_counts(k) + tn1);
  }
  trigamma_term *= tn2;
  digamma_term -= n_total * Rf_digamma(tn1);
  trigamma_term -= tn2 * n_total * Rf_trigamma(tn1);

  double ll_part_1 = 0.0, ll_part_2 = 0.0;
  for (arma::uword r = 0; r < mu.n_elem; r++) {
    double omt = 1.0 + mu(r) * theta;
    ll_part_1 += f(r) * (std::log(omt) - mu(r) / (mu(r) + tn1)) + S(r) / (mu(r) + tn1);
    ll_part_2 += (f(r) * mu(r) * mu(r) * theta + S(r)) / (omt * omt);
  }
  double ll_part = -2.0 * tn1 * (ll_part_1 - digamma_term) + (ll_part_2 + trigamma_term);
  return ll_part + cr_term * theta * theta + (ll_part_1 - digamma_term) * tn1 + cr_term2 * theta;
}
