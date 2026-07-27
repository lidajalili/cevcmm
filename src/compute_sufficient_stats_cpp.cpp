// =============================================================================
// Day 16: C++ port of compute_sufficient_stats() core math.
//
// Computes the six per-node cross-products
//   a   = y' y                          (scalar)
//   b   = X' y                          (p x 1 matrix)
//   C   = X' X                          (p x p)
//   ZtZ = Z' Z                          (q x q)
//   Zty = Z' y                          (q x 1 matrix)
//   XtZ = X' Z                          (p x q)
// plus the observation count n_obs, packaged as an Rcpp::List with the
// EXACT shapes (matrix not vec for b and Zty) that the R-side downstream
// code in fit_ss.R / fit_csl.R / etc. expects from crossprod().
//
// The R wrapper compute_sufficient_stats(..., use_cpp = TRUE) attaches
// the "vcmm_ss" class after calling this function.
//
// Bit-equivalence to the R path is verified by
// inst/validation/day16_rcpp_suffstats_validation.R to a tolerance of
// 1e-10 (typical observed difference is ~1e-13 from BLAS summation order).
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Sufficient-statistics core in C++
//'
//' Internal RcppArmadillo implementation of the six cross-products
//' that make up a VCMM per-node summary. Called by
//' \code{compute_sufficient_stats()} when \code{use_cpp = TRUE}
//' (the default since Day 16).
//'
//' @param y Numeric vector, length n.
//' @param X Numeric matrix, n by p.
//' @param Z Numeric matrix, n by q.
//' @return A plain \code{list} (no class) with elements
//'   \code{a, b, C, ZtZ, Zty, XtZ, n_obs}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List compute_sufficient_stats_cpp(
    const arma::vec& y,
    const arma::mat& X,
    const arma::mat& Z
) {
  const arma::uword n = X.n_rows;

  // R's crossprod(X, y) returns a p x 1 *matrix*. The downstream
  // arithmetic in fit_ss / fit_csl treats b and Zty as matrices (matrix
  // - matrix product compatibility). Force the same shape by copying
  // y (an n-vec, i.e. arma::Col<double>) into an n x 1 arma::mat.
  //
  // We use Armadillo's own copy constructor rather than std::memcpy
  // on y.memptr(). When n == 0 the vector's memptr() is NULL, and
  // std::memcpy(NULL, NULL, 0) is undefined behavior per the C
  // standard even though nothing is copied. gcc-UBSAN on the CRAN
  // Fedora tests-gcc-SAN farm reports this as
  //   runtime error: null pointer passed as argument 1,
  //   which is declared to never be null
  // even when the destination copy is zero bytes. Using the
  // arma::mat(arma::vec) constructor handles the n == 0 case
  // internally without dereferencing any pointer, and produces the
  // same n x 1 matrix for all n > 0.
  const arma::mat y_mat(y);

  // Six cross-products. Armadillo dispatches each to BLAS:
  //   X.t() * y_mat  -> dgemv
  //   X.t() * X      -> dgemm (Armadillo may use dsyrk if it detects
  //                     the X.t() * X symmetric pattern, depending on
  //                     the compile-time settings of the linked BLAS).
  //   X.t() * Z      -> dgemm
  const double    a   = arma::dot(y, y);
  const arma::mat b   = X.t() * y_mat;   // p x 1 matrix
  const arma::mat C   = X.t() * X;       // p x p
  const arma::mat ZtZ = Z.t() * Z;       // q x q
  const arma::mat Zty = Z.t() * y_mat;   // q x 1 matrix
  const arma::mat XtZ = X.t() * Z;       // p x q

  return Rcpp::List::create(
    Rcpp::Named("a")     = a,
    Rcpp::Named("b")     = b,
    Rcpp::Named("C")     = C,
    Rcpp::Named("ZtZ")   = ZtZ,
    Rcpp::Named("Zty")   = Zty,
    Rcpp::Named("XtZ")   = XtZ,
    Rcpp::Named("n_obs") = static_cast<int>(n)
  );
}
