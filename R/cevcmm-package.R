#' cevcmm: Communication-Efficient Varying Coefficient Mixed-Effects Models
#'
#' Scalable inference for Varying Coefficient Mixed-Effects Models
#' (VCMMs) with large, correlated random effects. The package implements:
#'
#' * Sufficient-statistics (SS) iterative estimator (Algorithm 1).
#' * One-step communication-efficient surrogate likelihood (CSL) estimator.
#' * SVD-stabilized variants for ill-conditioned random-effect Gram matrices.
#' * Kronecker and separable covariance structures for origin-destination
#'   and group-shared random effects.
#'
#' The package is in early development. See the project ROADMAP for the
#' current status.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed Effect Models: Methodology, Theory, and
#' Applications.
#'
#' @keywords internal
#' @importFrom stats coef cov nobs pnorm printCoefmat vcov
#' @importFrom Rcpp sourceCpp
#' @useDynLib cevcmm, .registration = TRUE
"_PACKAGE"
