#===============================================================================
# Kronecker random-effects covariance for VCMMs
#
# Implements the structured covariance Sigma_alpha = Sigma_2x2 ⊗ Sigma_spatial
# used in origin-destination (OD) settings, where each of G groups has a
# 2-dimensional random effect (origin / destination) and the q = 2 * G
# random-effect vector alpha is alpha = vec( G x 2 matrix ).
#
# Per Algorithm 1 of Lin and Jalili (2026), the nuisance parameters
# eta = (Sigma_alpha, sigma_eps^2) are updated at each iteration. For
# Sigma_alpha = Sigma_2x2 ⊗ Sigma_spatial, the moment-based M_eta rule
# is implemented here as estimate_kronecker_components(); fit_ss() and
# fit_csl() invoke it once per iteration when control$update_variance is TRUE.
#===============================================================================

#-------------------------------------------------------------------------------
# Build the random-effect precision matrix from Kronecker components.
#
# Returns sigma_eps^2 * (Sigma_2x2_inv ⊗ Sigma_spatial_inv), which goes into
# the alpha-update equation in place of the diagonal ridge (sigma_eps^2 /
# sigma_alpha^2) I_q used for re_cov = "diag".
#-------------------------------------------------------------------------------

#' Build the Kronecker-structured random-effect precision matrix
#'
#' For Sigma_alpha = Sigma_2x2 \%x\% Sigma_spatial (the OD setting), returns
#' \deqn{
#'   \sigma_\varepsilon^2 \cdot
#'   (\Sigma_{2 \times 2}^{-1} \otimes \Sigma_{\text{spatial}}^{-1})
#' }
#' which is added to \code{crossprod(Z)} inside the random-effect block of the
#' VCMM Hessian. This is the structured analogue of
#' \code{(sigma_eps^2 / sigma_alpha^2) * I_q} used under
#' \code{re_cov = "diag"}.
#'
#' @param Sigma_2x2 A 2 by 2 positive-definite numeric matrix.
#' @param Sigma_spatial A G by G positive-definite numeric matrix
#'   describing spatial correlation across G groups.
#' @param sigma_eps Positive numeric. Residual standard deviation
#'   (used to scale the precision into the implicit-loss parametrisation
#'   that fit_ss and fit_csl share).
#'
#' @return A 2G by 2G numeric matrix.
#'
#' @references
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
build_kronecker_precision <- function(Sigma_2x2, Sigma_spatial, sigma_eps) {
  if (!is.matrix(Sigma_2x2) || !identical(dim(Sigma_2x2), c(2L, 2L))) {
    stop("`Sigma_2x2` must be a 2 by 2 matrix.", call. = FALSE)
  }
  if (!is.matrix(Sigma_spatial) ||
      nrow(Sigma_spatial) != ncol(Sigma_spatial)) {
    stop("`Sigma_spatial` must be a square matrix.", call. = FALSE)
  }
  if (!is.numeric(sigma_eps) || length(sigma_eps) != 1L ||
      sigma_eps <= 0 || !is.finite(sigma_eps)) {
    stop("`sigma_eps` must be a single positive finite numeric.",
         call. = FALSE)
  }

  Sigma_2x2_inv     <- solve(Sigma_2x2)
  Sigma_spatial_inv <- solve(Sigma_spatial)
  (sigma_eps^2) * kronecker(Sigma_2x2_inv, Sigma_spatial_inv)
}

#-------------------------------------------------------------------------------
# Moment-based M_eta rule for Sigma_2x2 and Sigma_spatial from alpha_hat.
#
# Given alpha_hat = vec(G x 2 matrix) where column 1 is origin effects and
# column 2 is destination effects:
#   Sigma_2x2_hat     = cov(alpha_mat)                    (2 x 2)
#   Sigma_spatial_hat = average of normalized tcrossprod  (G x G)
#
# Both are projected to be symmetric positive-definite via a minimal ridge.
# Sigma_spatial_hat is rescaled to a correlation matrix (unit diagonal) for
# identifiability of the Kronecker decomposition.
#
# Note on identifiability: with one alpha_hat sample of length 2G, the row
# covariance is a rank-1 quantity. The result captures the dominant spatial
# pattern of alpha_hat but is not a full G x G covariance in the strict
# statistical sense. Users with a parametric spatial kernel (e.g.
# exponential decay exp(-d/phi)) should pass it via Sigma_spatial_init
# and set update_variance = FALSE to keep it fixed.
#-------------------------------------------------------------------------------

#' Moment-based estimator for Kronecker covariance components
#'
#' Given a length-\eqn{2G} random-effects estimate
#' \eqn{\hat\alpha = \mathrm{vec}(G \times 2\ \mathrm{matrix})} (column 1
#' = origin effects, column 2 = destination effects), returns
#' moment-based estimates of \code{Sigma_2x2} and \code{Sigma_spatial}
#' that together specify \eqn{\hat\Sigma_\alpha = \Sigma_{2\times 2}
#' \otimes \Sigma_{\text{spatial}}}.
#'
#' This is the \eqn{M_\eta} update rule the paper permits (one of several
#' allowed choices listed alongside ML, REML, quasi-likelihood, and
#' Fisher-scoring). It is the form actually used in the simulation code's
#' post-hoc covariance estimator, lifted into the iterative loop here so
#' that the package works for real OD data where \eqn{\Sigma_\alpha} is
#' unknown.
#'
#' \strong{Identifiability caveat.} With only one \eqn{\hat\alpha}
#' sample, the spatial component is reconstructed from a rank-1
#' quantity. The output captures the dominant spatial pattern in
#' \eqn{\hat\alpha} but is a low-rank approximation rather than a full
#' rank-G covariance. For applications where a parametric spatial
#' kernel is available, supply it via \code{Sigma_spatial_init} in
#' \code{vcmm()} and fit with \code{update_variance = FALSE}.
#'
#' @param alpha Numeric vector of length \eqn{2G}.
#' @param n_groups Integer \eqn{G}, the number of groups. Must satisfy
#'   \code{length(alpha) == 2 * n_groups}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{Sigma_2x2}: 2 by 2 symmetric positive-definite matrix.
#'   \item \code{Sigma_spatial}: G by G symmetric positive-definite
#'     correlation matrix (unit diagonal).
#' }
#'
#' @references
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
estimate_kronecker_components <- function(alpha, n_groups) {
  if (!is.numeric(alpha)) {
    stop("`alpha` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(n_groups) || length(n_groups) != 1L ||
      n_groups < 1 || n_groups != as.integer(n_groups)) {
    stop("`n_groups` must be a single positive integer.", call. = FALSE)
  }
  n_groups <- as.integer(n_groups)
  if (length(alpha) != 2L * n_groups) {
    stop(sprintf(
      "length(alpha) = %d does not match 2 * n_groups = %d.",
      length(alpha), 2L * n_groups), call. = FALSE)
  }

  # Reshape: G x 2  (column 1 = origin, column 2 = destination)
  alpha_mat <- matrix(alpha, nrow = n_groups, ncol = 2L)

  # ----- Sigma_2x2 = column covariance, projected to PD --------------------
  Sigma_2x2 <- stats::cov(alpha_mat)
  Sigma_2x2 <- .project_pd(Sigma_2x2)

  # ----- Sigma_spatial = average normalised row covariance ----------------
  alpha_O <- alpha_mat[, 1L]
  alpha_D <- alpha_mat[, 2L]

  norm_O <- max(sum(alpha_O^2), 1e-8)
  norm_D <- max(sum(alpha_D^2), 1e-8)
  Sigma_spatial_O <- tcrossprod(alpha_O) / norm_O
  Sigma_spatial_D <- tcrossprod(alpha_D) / norm_D
  Sigma_spatial   <- (Sigma_spatial_O + Sigma_spatial_D) / 2

  # Rescale to unit diagonal (correlation matrix) for Kronecker identifiability
  diag_vals <- diag(Sigma_spatial)
  diag_vals[diag_vals < 1e-8] <- 1
  D_inv <- diag(1 / sqrt(diag_vals), n_groups)
  Sigma_spatial <- D_inv %*% Sigma_spatial %*% D_inv
  Sigma_spatial <- (Sigma_spatial + t(Sigma_spatial)) / 2   # symmetry

  # Project to PD with minimal ridge
  Sigma_spatial <- .project_pd(Sigma_spatial)

  list(Sigma_2x2 = Sigma_2x2, Sigma_spatial = Sigma_spatial)
}

#-------------------------------------------------------------------------------
# Internal: project a symmetric matrix to positive-definite with minimal ridge
#-------------------------------------------------------------------------------
.project_pd <- function(M, ridge_floor = 1e-8) {
  M <- (M + t(M)) / 2
  eig_min <- min(eigen(M, only.values = TRUE, symmetric = TRUE)$values)
  if (eig_min < ridge_floor) {
    ridge <- max(ridge_floor, -eig_min + ridge_floor)
    M <- M + diag(ridge, nrow(M))
  }
  M
}

#-------------------------------------------------------------------------------
# Internal: build a fresh re_cov_state object (passed through fit_ss/fit_csl)
#-------------------------------------------------------------------------------
.new_re_cov_state <- function(type,
                              Sigma_2x2     = NULL,
                              Sigma_spatial = NULL,
                              n_groups      = NULL) {
  list(
    type          = type,
    Sigma_2x2     = Sigma_2x2,
    Sigma_spatial = Sigma_spatial,
    n_groups      = n_groups
  )
}

#-------------------------------------------------------------------------------
# Internal: assemble the prior precision matrix for the alpha-update.
#
# For diag:      sigma_eps^2 / sigma_alpha^2 * I_q
# For kronecker: sigma_eps^2 * (Sigma_2x2_inv ⊗ Sigma_spatial_inv)
#-------------------------------------------------------------------------------
.build_prior_precision <- function(re_cov_state, sigma_eps, sigma_alpha, q) {
  type <- if (is.null(re_cov_state)) "diag" else re_cov_state$type
  switch(
    type,
    diag = {
      ridge <- (sigma_eps^2) / (sigma_alpha^2)
      diag(ridge, q)
    },
    kronecker = build_kronecker_precision(re_cov_state$Sigma_2x2,
                                          re_cov_state$Sigma_spatial,
                                          sigma_eps),
    stop(sprintf("Unknown re_cov_state type: '%s'.", type), call. = FALSE)
  )
}

#-------------------------------------------------------------------------------
# Internal: update re_cov_state in-place for one M_eta step.
#
# Design choice (v0.1):
#   For re_cov = "kronecker" we update Sigma_2x2 only via column covariance
#   of the (G x 2) reshape of alpha_hat. Sigma_spatial is HELD FIXED at the
#   user-supplied initial value because the moment-based estimator from
#   one alpha_hat is rank-2 at most, which is not identifiable and causes
#   iterative reinforcement of noise. For real OD data, Sigma_spatial
#   should be a parametric spatial kernel supplied by the user (e.g.
#   exp(-d/phi) with a known distance matrix). A future release will add
#   parametric spatial estimation.
#-------------------------------------------------------------------------------
.update_re_cov_state <- function(re_cov_state, alpha) {
  if (is.null(re_cov_state) || re_cov_state$type == "diag") {
    return(re_cov_state)
  }
  if (re_cov_state$type == "kronecker") {
    # Update Sigma_2x2 only (low-dim, identifiable from G samples)
    alpha_mat <- matrix(alpha, nrow = re_cov_state$n_groups, ncol = 2L)
    S2 <- stats::cov(alpha_mat)
    S2 <- .project_pd(S2)
    re_cov_state$Sigma_2x2 <- S2
    # Sigma_spatial is intentionally NOT updated -- held fixed at the
    # user-supplied (or identity-default) value.
    return(re_cov_state)
  }
  re_cov_state
}
