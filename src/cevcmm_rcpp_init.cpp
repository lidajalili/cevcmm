// Day 15: minimal RcppArmadillo scaffolding for cevcmm.
//
// This file's only purpose is to verify the C++ toolchain end-to-end:
//   1. The compiler invokes correctly (-std=c++14 or higher).
//   2. Rcpp dispatch generates the right .Call glue.
//   3. RcppArmadillo's headers and library link.
//
// Days 16-18 replace the placeholders with real ports:
//   - accumulate_ss_cpp()  (Day 16)
//   - ss_solve_cpp()       (Day 17)
//   - truncated SVD path   (Day 18)
//
// Do not delete cevcmm_rcpp_check() until Days 16-18 are done; it is
// the canary the day-15 validation script checks.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <sstream>
#include <string>

//' Verify the Rcpp / RcppArmadillo toolchain is wired up
//'
//' Internal stub used by Day 15's profiling script to confirm
//' the C++ compiler, Rcpp dispatch, and Armadillo linking all work
//' before the real ports on Days 16-18. Performs a trivial 3x3
//' identity-matrix trace to exercise the linker.
//'
//' @return A character string of the form
//'   "OK (Armadillo X.Y.Z; trace(I_3) = 3)".
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
std::string cevcmm_rcpp_check() {
  arma::mat M = arma::eye<arma::mat>(3, 3);
  double tr = arma::trace(M);
  std::ostringstream out;
  out << "OK (Armadillo "
      << ARMA_VERSION_MAJOR << "."
      << ARMA_VERSION_MINOR << "."
      << ARMA_VERSION_PATCH
      << "; trace(I_3) = " << tr << ")";
  return out.str();
}
