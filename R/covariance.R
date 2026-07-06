#===============================================================================
# Kronecker-structured random-effects covariance for VCMMs
#
# Implements the structured covariance Sigma_alpha = Sigma_left ⊗ Sigma_right
# in the column-stacking convention:
#
#   alpha = vec_col(M),   M is (G x q_left),
#   Var(alpha) = Sigma_left (q_left x q_left)  ⊗  Sigma_right (G x G).
#
# Two user-facing names route through the same internals:
#
#   re_cov = "kronecker": OD-style; q_left defaults to 2 (origin / dest);
#                         user-facing matrices are Sigma_2x2 and Sigma_spatial.
#   re_cov = "separable": group-shared dense; q_left is required and is
#                         the per-group random-effect dimension;
#                         user-facing matrices are Sigma_q and Omega_G.
#
# Mathematically identical -- only the parameter names and defaults differ.
# Per Algorithm 1 of Jalili and Lin (2025) the M_eta rule updates the
# random-effects covariance nuisance once per iteration; here we implement
# the weighted moment estimator + EM-style posterior-variance correction
# (Theorem 1 permits any fixed map of the summaries).
#===============================================================================

#-------------------------------------------------------------------------------
# Build the prior precision matrix sigma_eps^2 * (Sigma_left^{-1} ⊗ Sigma_right^{-1})
# which augments crossprod(Z) in the alpha-update equation.
#-------------------------------------------------------------------------------

#' Build the Kronecker-structured random-effect precision matrix
#'
#' For \eqn{\Sigma_\alpha = \Sigma_{\mathrm{left}} \otimes \Sigma_{\mathrm{right}}}
#' (column-stacking convention), returns
#' \deqn{
#'   \sigma_\varepsilon^2 \cdot
#'   (\Sigma_{\mathrm{left}}^{-1} \otimes \Sigma_{\mathrm{right}}^{-1})
#' }
#' which is added to \code{crossprod(Z)} inside the random-effect block of
#' the VCMM Hessian. This is the structured analogue of
#' \code{(sigma_eps^2 / sigma_alpha^2) * I_q} used under
#' \code{re_cov = "diag"}.
#'
#' @param Sigma_left A \eqn{k \times k} positive-definite numeric matrix
#'   (\eqn{k = 2} for OD-style; arbitrary \eqn{k} for general separable).
#' @param Sigma_right A \eqn{G \times G} positive-definite numeric matrix.
#' @param sigma_eps Positive numeric. Residual standard deviation.
#'
#' @return A \eqn{kG \times kG} numeric matrix.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
build_kronecker_precision <- function(Sigma_left, Sigma_right, sigma_eps) {
  if (!is.matrix(Sigma_left) || nrow(Sigma_left) != ncol(Sigma_left)) {
    stop("`Sigma_left` must be a square matrix.", call. = FALSE)
  }
  if (!is.matrix(Sigma_right) || nrow(Sigma_right) != ncol(Sigma_right)) {
    stop("`Sigma_right` must be a square matrix.", call. = FALSE)
  }
  if (!is.numeric(sigma_eps) || length(sigma_eps) != 1L ||
      sigma_eps <= 0 || !is.finite(sigma_eps)) {
    stop("`sigma_eps` must be a single positive finite numeric.",
         call. = FALSE)
  }

  Sigma_left_inv  <- solve(Sigma_left)
  Sigma_right_inv <- solve(Sigma_right)
  (sigma_eps^2) * kronecker(Sigma_left_inv, Sigma_right_inv)
}

#-------------------------------------------------------------------------------
# Moment-based M_eta rule for Sigma_left from alpha_hat with Sigma_right fixed.
#
# Estimator (weighted, GLS-style):
#   Sigma_left_hat = (1/G) alpha_mat' Sigma_right^{-1} alpha_mat
#
# where alpha_mat = matrix(alpha, nrow = G, ncol = q_left) under
# column-stacking. Unbiased for the TRUE alpha when Sigma_right is correct;
# combine with EM correction below when plugging in the BLUP.
#-------------------------------------------------------------------------------

#' Moment-based estimator for the Kronecker left component
#'
#' Given a length-\eqn{kG} random-effects estimate
#' \eqn{\hat\alpha = \mathrm{vec}_{\mathrm{col}}(M)} (column-stacked,
#' \eqn{M \in \mathbb R^{G \times k}}) and the right-side covariance
#' \code{Sigma_right}, returns the weighted moment estimate
#' \deqn{\hat\Sigma_{\mathrm{left}} =
#'   \tfrac{1}{G}\, M^{\top}\, \Sigma_{\mathrm{right}}^{-1}\, M.}
#'
#' This is unbiased under the Kronecker model
#' \eqn{\alpha \sim N(0, \Sigma_{\mathrm{left}} \otimes
#' \Sigma_{\mathrm{right}})} when \code{Sigma_right} is correct.
#' When \eqn{\hat\alpha} is the BLUP rather than the true \eqn{\alpha},
#' apply the EM-style correction in \code{\link{vcmm}} (handled
#' automatically by \code{fit_ss} / \code{fit_csl}).
#'
#' Backwards-compatible alias: if \code{q_left = 2} (the default), this is
#' the same estimator as the previous \code{estimate_kronecker_components}
#' for the OD setting.
#'
#' @param alpha Numeric vector of length \eqn{kG}.
#' @param n_groups Integer \eqn{G}.
#' @param q_left Integer \eqn{k}, the left (within) dimension. Defaults to
#'   2 for backward compatibility with the OD setting.
#' @param Sigma_right Optional \eqn{G \times G} positive-definite covariance.
#'   If \code{NULL}, the unweighted sample covariance \code{cov(alpha_mat)}
#'   is returned (unbiased only when rows are uncorrelated).
#'
#' @return A \eqn{k \times k} symmetric positive-definite matrix.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
estimate_kronecker_components <- function(alpha,
                                          n_groups,
                                          q_left      = 2L,
                                          Sigma_right = NULL) {
  if (!is.numeric(alpha)) {
    stop("`alpha` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(n_groups) || length(n_groups) != 1L ||
      n_groups < 1 || n_groups != as.integer(n_groups)) {
    stop("`n_groups` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(q_left) || length(q_left) != 1L ||
      q_left < 1 || q_left != as.integer(q_left)) {
    stop("`q_left` must be a single positive integer.", call. = FALSE)
  }
  G <- as.integer(n_groups)
  k <- as.integer(q_left)

  if (length(alpha) != k * G) {
    stop(sprintf(
      "length(alpha) = %d does not match q_left * n_groups = %d * %d = %d.",
      length(alpha), k, G, k * G), call. = FALSE)
  }
  if (!is.null(Sigma_right)) {
    if (!is.matrix(Sigma_right) ||
        !isTRUE(all.equal(dim(Sigma_right), c(G, G)))) {
      stop(sprintf("`Sigma_right` must be a %d by %d matrix.", G, G),
           call. = FALSE)
    }
  }

  # Column-stacking: alpha_mat[g, k] = alpha[(k-1)*G + g]
  alpha_mat <- matrix(alpha, nrow = G, ncol = k)

  if (is.null(Sigma_right)) {
    S <- stats::cov(alpha_mat)
  } else {
    Sinv_alpha <- solve(Sigma_right, alpha_mat)
    S <- crossprod(alpha_mat, Sinv_alpha) / G
  }
  S <- (S + t(S)) / 2
  .project_pd(S)
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
# Internal: build a fresh re_cov_state object.
#
# Canonical field names: Sigma_left, Sigma_right, n_groups, q_left.
# Legacy aliases populated when type == "kronecker" with q_left = 2
# (Sigma_2x2 / Sigma_spatial), and when type == "separable"
# (Sigma_q / Omega_G), so existing scripts and the print method work.
#-------------------------------------------------------------------------------
.new_re_cov_state <- function(type,
                              Sigma_left  = NULL,
                              Sigma_right = NULL,
                              n_groups    = NULL,
                              q_left      = NULL) {
  state <- list(
    type        = type,
    Sigma_left  = Sigma_left,
    Sigma_right = Sigma_right,
    n_groups    = n_groups,
    q_left      = q_left
  )
  .sync_legacy_aliases(state)
}

#-------------------------------------------------------------------------------
# Internal: keep legacy field aliases in sync with the canonical fields.
# Called after every update to re_cov_state$Sigma_left or $Sigma_right.
#-------------------------------------------------------------------------------
.sync_legacy_aliases <- function(state) {
  if (is.null(state) || identical(state$type, "diag")) return(state)

  if (identical(state$type, "kronecker") && identical(state$q_left, 2L)) {
    state$Sigma_2x2     <- state$Sigma_left
    state$Sigma_spatial <- state$Sigma_right
  }
  if (identical(state$type, "separable")) {
    state$Sigma_q  <- state$Sigma_left
    state$Omega_G  <- state$Sigma_right
  }
  state
}

#-------------------------------------------------------------------------------
# Internal: assemble the prior precision matrix for the alpha-update.
#
#   diag                -> sigma_eps^2 / sigma_alpha^2 * I_q
#   kronecker/separable -> sigma_eps^2 * (Sigma_left_inv ⊗ Sigma_right_inv)
#-------------------------------------------------------------------------------
.build_prior_precision <- function(re_cov_state, sigma_eps, sigma_alpha, q) {
  type <- if (is.null(re_cov_state)) "diag" else re_cov_state$type
  if (identical(type, "diag")) {
    ridge <- (sigma_eps^2) / (sigma_alpha^2)
    return(diag(ridge, q))
  }
  if (identical(type, "kronecker") || identical(type, "separable")) {
    return(build_kronecker_precision(re_cov_state$Sigma_left,
                                     re_cov_state$Sigma_right,
                                     sigma_eps))
  }
  stop(sprintf("Unknown re_cov_state type: '%s'.", type), call. = FALSE)
}

#-------------------------------------------------------------------------------
# Internal: weighted-moment update (no EM correction).
# Used inside the iterative loop; EM correction is applied once at convergence
# via .apply_em_correction_if_kronecker().
#-------------------------------------------------------------------------------
.update_re_cov_state <- function(re_cov_state, alpha) {
  if (is.null(re_cov_state) || identical(re_cov_state$type, "diag")) {
    return(re_cov_state)
  }
  if (identical(re_cov_state$type, "kronecker") ||
      identical(re_cov_state$type, "separable")) {

    G <- re_cov_state$n_groups
    k <- re_cov_state$q_left
    alpha_mat <- matrix(alpha, nrow = G, ncol = k)

    Sinv_alpha <- solve(re_cov_state$Sigma_right, alpha_mat)
    S <- crossprod(alpha_mat, Sinv_alpha) / G
    S <- (S + t(S)) / 2
    S <- .project_pd(S)

    re_cov_state$Sigma_left <- S
    return(.sync_legacy_aliases(re_cov_state))
  }
  re_cov_state
}

#-------------------------------------------------------------------------------
# Internal: EM-style correction for Sigma_left.
#
# Plugging the BLUP into the moment estimator gives
#   E[(1/G) M' Sinv M] = Sigma_left - C_corr
# where C_corr is the posterior-variance partial trace:
#   C_corr[i, j] = (1/G) sum_{g, h} Sinv[g, h] *
#                  V_alpha_alpha[(i-1)G + g, (j-1)G + h]
# Adding C_corr removes the BLUP shrinkage bias exactly under the Gaussian
# model when Sigma_right is correct. This is the M_eta rule used by fit_ss
# and fit_csl when re_cov is "kronecker" or "separable".
#-------------------------------------------------------------------------------
.estimate_sigma_left_em <- function(alpha,
                                    V_aa_post,
                                    Sigma_right,
                                    n_groups,
                                    q_left) {
  G <- as.integer(n_groups)
  k <- as.integer(q_left)
  stopifnot(length(alpha) == k * G,
            identical(dim(V_aa_post), c(k * G, k * G)),
            identical(dim(Sigma_right), c(G, G)))

  alpha_mat <- matrix(alpha, nrow = G, ncol = k)
  Sinv      <- solve(Sigma_right)

  # Moment term: (1/G) M' Sinv M
  S_moment <- crossprod(alpha_mat, Sinv %*% alpha_mat) / G

  # Posterior-variance correction (partial trace).
  C_corr <- matrix(0, k, k)
  for (i in seq_len(k)) {
    rows_i <- ((i - 1L) * G + 1L):(i * G)
    for (j in seq_len(k)) {
      cols_j       <- ((j - 1L) * G + 1L):(j * G)
      V_block      <- V_aa_post[rows_i, cols_j, drop = FALSE]
      C_corr[i, j] <- sum(Sinv * V_block)
    }
  }
  C_corr <- C_corr / G

  S <- S_moment + C_corr
  S <- (S + t(S)) / 2
  .project_pd(S)
}

#-------------------------------------------------------------------------------
# Internal: apply the EM correction once at convergence. Called from both
# fit_ss (after K_inv is built) and fit_csl (after the Newton step).
# No-op for re_cov = "diag" or when update_variance = FALSE.
#-------------------------------------------------------------------------------
.apply_em_correction_if_kronecker <- function(re_cov_state,
                                              alpha,
                                              sigma_eps,
                                              K_inv,
                                              p,
                                              q,
                                              control) {
  if (!isTRUE(control$update_variance) || is.null(re_cov_state)) {
    return(re_cov_state)
  }
  if (!(identical(re_cov_state$type, "kronecker") ||
        identical(re_cov_state$type, "separable"))) {
    return(re_cov_state)
  }

  idx_a     <- seq.int(p + 1L, p + q)
  V_aa_post <- (sigma_eps^2) * K_inv[idx_a, idx_a, drop = FALSE]

  re_cov_state$Sigma_left <- .estimate_sigma_left_em(
    alpha       = alpha,
    V_aa_post   = V_aa_post,
    Sigma_right = re_cov_state$Sigma_right,
    n_groups    = re_cov_state$n_groups,
    q_left      = re_cov_state$q_left
  )
  .sync_legacy_aliases(re_cov_state)
}
