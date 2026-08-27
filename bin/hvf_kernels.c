/* Transpose-free replacements for Seurat's SparseRowVar2 / SparseRowVarStd.
 *
 * Seurat takes `Eigen::SparseMatrix<double> mat` BY VALUE and then runs
 * `mat = mat.transpose()`, because it wants per-gene statistics out of a
 * cell-major matrix. At 305k cells that is three full-size copies of the counts
 * matrix live at once (Rcpp's copy into Eigen, the transpose temporary, and the
 * assignment result), none of which show up in R's gc() accounting. Measured on
 * the cgroup: the R matrix sits at 16.63 GB and each of the two calls climbs
 * 16.63 -> 32.6 -> 48.6 -> 64.6 GB and back down. That excursion alone set step
 * 2's container peak; every other stage stays under 43.8 GB.
 *
 * One column-major pass needs no transpose and no copy -- it reads the
 * dgCMatrix slots in place. Bit-identity is the whole point, so the
 * accumulation order matches Seurat's: Seurat walks column k of the TRANSPOSED
 * matrix, i.e. gene k's nonzeros in ascending original-column order, and
 * visiting original columns j = 0..n-1 in order hands each gene its terms in
 * that same ascending-j sequence. Every partial sum is therefore the identical
 * double at every step. The tail terms (zero contribution, then the divide) are
 * applied in Seurat's order too, and pow() is used where Seurat uses pow().
 *
 * Plain C, not C++: R's Makeconf in this image points at a conda toolchain that
 * is absent, and building with the system g++ would link a different libstdc++
 * than R itself uses. C against R's own API has no such ABI surface.
 */
#include <R.h>
#include <Rinternals.h>
#include <math.h>

SEXP rowVar2_cm(SEXP i_, SEXP p_, SEXP x_, SEXP ng_, SEXP nc_, SEXP mu_) {
  const int *i = INTEGER(i_), *p = INTEGER(p_);
  const double *x = REAL(x_), *mu = REAL(mu_);
  const int ng = asInteger(ng_), nc = asInteger(nc_);
  double *acc   = (double *) R_alloc((size_t) ng, sizeof(double));
  int    *nZero = (int *)    R_alloc((size_t) ng, sizeof(int));
  for (int g = 0; g < ng; ++g) { acc[g] = 0.0; nZero[g] = nc; }
  for (int j = 0; j < nc; ++j) {
    const R_xlen_t lo = p[j], hi = p[j + 1];
    for (R_xlen_t k = lo; k < hi; ++k) {
      const int g = i[k];
      nZero[g] -= 1;
      acc[g] += pow(x[k] - mu[g], 2);
    }
  }
  SEXP out = PROTECT(allocVector(REALSXP, ng));
  double *o = REAL(out);
  for (int g = 0; g < ng; ++g) {
    double colSum = acc[g] + pow(mu[g], 2) * nZero[g];
    o[g] = colSum / (nc - 1);
  }
  UNPROTECT(1);
  return out;
}

SEXP rowVarStd_cm(SEXP i_, SEXP p_, SEXP x_, SEXP ng_, SEXP nc_,
                  SEXP mu_, SEXP sd_, SEXP vmax_) {
  const int *i = INTEGER(i_), *p = INTEGER(p_);
  const double *x = REAL(x_), *mu = REAL(mu_), *sd = REAL(sd_);
  const int ng = asInteger(ng_), nc = asInteger(nc_);
  const double vmax = asReal(vmax_);
  double *acc   = (double *) R_alloc((size_t) ng, sizeof(double));
  int    *nZero = (int *)    R_alloc((size_t) ng, sizeof(int));
  for (int g = 0; g < ng; ++g) { acc[g] = 0.0; nZero[g] = nc; }
  for (int j = 0; j < nc; ++j) {
    const R_xlen_t lo = p[j], hi = p[j + 1];
    for (R_xlen_t k = lo; k < hi; ++k) {
      const int g = i[k];
      nZero[g] -= 1;                 /* Seurat decrements before the sd check */
      if (sd[g] == 0) continue;
      double v = (x[k] - mu[g]) / sd[g];
      if (vmax < v) v = vmax;        /* std::min(vmax, v) */
      acc[g] += pow(v, 2);
    }
  }
  SEXP out = PROTECT(allocVector(REALSXP, ng));
  double *o = REAL(out);
  for (int g = 0; g < ng; ++g) {
    if (sd[g] == 0) { o[g] = 0.0; continue; }   /* Seurat leaves the 0-init */
    double colSum = acc[g] + pow((0 - mu[g]) / sd[g], 2) * nZero[g];
    o[g] = colSum / (nc - 1);
  }
  UNPROTECT(1);
  return out;
}
