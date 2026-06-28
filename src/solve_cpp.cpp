// =============================================================================
// Day 17: RcppArmadillo backends for invert_matrix() dispatch.
//
// Three functions cover the inversion methods used in VCMM:
//
//   invert_spd_cpp(A):
//     Cholesky-based inverse for symmetric positive-definite matrices.
//     LAPACK dpotrf + dpotri. ~2x faster than LU for SPD inputs.
//     This is the *expected* fast path for VCMM K matrices (which are
//     prior-augmented Gram matrices, always symmetric PD by construction).
//
//   invert_general_cpp(A):
//     LU-based inverse for general square matrices.
//     LAPACK dgetrf + dgetri. Used as a fallback if Cholesky fails (i.e.
//     if the matrix turns out not to be PD).
//
//   pinv_cpp(A, tol):
//     SVD-based Moore-Penrose pseudo-inverse.
//     LAPACK dgesvd. Used as the final fallback if both Cholesky and LU
//     fail, and for the small-matrix branch of the SVD path. The
//     split-merge SVD variant (paper Algorithm 2) stays in R because
//     it's an algorithmic choice driven by numerical stability, not a
//     perf hotspot.
//
// All three return arma::mat. On Armadillo internal failure (bool false
// from the safe-form inv/pinv), they call Rcpp::stop() so the R-side
// tryCatch wrapper can fall through to the next method.
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Cholesky-based inverse for symmetric positive-definite matrices
//'
//' Internal RcppArmadillo backend for the fast path of
//' \code{\link{invert_matrix}}. Uses LAPACK \code{dpotrf} + \code{dpotri}.
//'
//' Throws an R error if \code{A} is not symmetric positive-definite; the
//' R wrapper catches this and falls back to \code{invert_general_cpp}.
//'
//' @param A Symmetric positive-definite square matrix.
//' @return Inverse of A.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat invert_spd_cpp(const arma::mat& A) {
  // Defensively symmetrize. The K matrices passed in from fit_ss /
  // fit_csl are mathematically symmetric (X'X plus a symmetric prior /
  // penalty) but the floating-point BLAS summation order leaves
  // |A[i,j] - A[j,i]| at ~1e-14 in some entries. Armadillo's
  // inv_sympd() strictly checks symmetry and logs a "given matrix is
  // not symmetric" warning when the asymmetry exceeds its tolerance --
  // it still produces a correct inverse using the upper triangle, but
  // the warning is noise users shouldn't see. Average upper and lower
  // triangles to make the input bit-exactly symmetric.
  const arma::mat A_sym = 0.5 * (A + A.t());

  arma::mat A_inv;
  const bool ok = arma::inv_sympd(A_inv, A_sym);
  if (!ok) {
    Rcpp::stop("invert_spd_cpp: matrix is not symmetric positive-definite");
  }
  return A_inv;
}

//' General matrix inverse via LU factorisation
//'
//' Internal RcppArmadillo backend used as a fallback by
//' \code{\link{invert_matrix}} when Cholesky fails. Uses LAPACK
//' \code{dgetrf} + \code{dgetri}.
//'
//' @param A Square invertible matrix.
//' @return Inverse of A.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat invert_general_cpp(const arma::mat& A) {
  arma::mat A_inv;
  const bool ok = arma::inv(A_inv, A);
  if (!ok) {
    Rcpp::stop("invert_general_cpp: matrix is singular or near-singular");
  }
  return A_inv;
}

//' SVD-based Moore-Penrose pseudo-inverse
//'
//' Internal RcppArmadillo backend used as the final fallback by
//' \code{\link{invert_matrix}}. Uses LAPACK \code{dgesvd} via
//' \code{arma::pinv}.
//'
//' @param A Numeric matrix.
//' @param tol Tolerance below which singular values are treated as zero.
//'   If \code{<= 0}, use Armadillo's default
//'   (\code{max(dim(A)) * eps_machine * max(svd(A))}).
//' @return Moore-Penrose pseudo-inverse of A.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat pinv_cpp(const arma::mat& A, double tol = -1.0) {
  arma::mat A_pinv;
  bool ok;
  if (tol <= 0.0) {
    ok = arma::pinv(A_pinv, A);
  } else {
    ok = arma::pinv(A_pinv, A, tol);
  }
  if (!ok) {
    Rcpp::stop("pinv_cpp: SVD did not converge");
  }
  return A_pinv;
}
