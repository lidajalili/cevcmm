#===============================================================================
# predict.vcmm_fit and logLik.vcmm_fit
#
# Implements the paper's subject-specific prediction formula (Section 5):
#   hat_y_i = X_design_i . beta_hat + Z_i . alpha_hat
# for groups seen during training (include_random = TRUE, the default).
# Provides a clean opt-out (include_random = FALSE) for the
# "new groups not seen in training" case, which is the marginal-mean
# predictor hat_y_i = X_design_i . beta_hat used by the package_comparison
# code under use_re = FALSE.
#
# logLik.vcmm_fit returns the (cached) marginal log-likelihood at
# convergence with attributes df and nobs, so AIC / BIC just work.
#===============================================================================

# ---------------------------------------------------------------------------
# Internal: marginal log-likelihood of the linear normal VCMM, evaluated at
# the converged (beta, alpha, sigma_eps, Sigma_alpha).
#
# Marginal model: y ~ N(X beta, Sigma_y), Sigma_y = sigma_eps^2 I + Z Sigma_alpha Z'.
#
# Using the matrix-determinant identity
#     |Sigma_y| = sigma_eps^{2(n - q)} * |Sigma_alpha| * |V_aa|,
# where V_aa = Z'Z + sigma_eps^2 * Sigma_alpha^{-1} (the alpha block of K
# we already assembled at convergence), and the alpha-update identity
#     Z'(y - X beta) = V_aa alpha,
# the quadratic form reduces to
#     (y - X beta)' Sigma_y^{-1} (y - X beta)
#       = (1 / sigma_eps^2) * [ ||y - X beta||^2 - alpha' V_aa alpha ].
# Both pieces are computable from the sufficient statistics + (beta, alpha),
# so no raw data is required.
# ---------------------------------------------------------------------------
.compute_marginal_loglik <- function(stats, beta, alpha,
                                     sigma_eps, sigma_alpha,
                                     re_cov_state, V_aa,
                                     n_obs, q) {

  # ||y - X beta||^2 from stats
  y_minus_Xb_norm2 <- as.numeric(
    stats$a -
    2 * crossprod(beta, stats$b) +
        crossprod(beta, stats$C %*% beta)
  )

  # alpha' V_aa alpha
  alpha_Vaa_alpha <- as.numeric(crossprod(alpha, V_aa %*% alpha))

  # Quadratic form
  quad <- (y_minus_Xb_norm2 - alpha_Vaa_alpha) / (sigma_eps^2)

  # log|Sigma_alpha|
  log_det_Sa <- if (is.null(re_cov_state) ||
                   identical(re_cov_state$type, "diag")) {
    q * log(sigma_alpha^2)
  } else {
    G  <- as.integer(re_cov_state$n_groups)
    kL <- as.integer(re_cov_state$q_left)
    # |Sigma_left ⊗ Sigma_right| = |Sigma_left|^G * |Sigma_right|^kL
    ld_left  <- as.numeric(determinant(re_cov_state$Sigma_left,
                                       logarithm = TRUE)$modulus)
    ld_right <- as.numeric(determinant(re_cov_state$Sigma_right,
                                       logarithm = TRUE)$modulus)
    G * ld_left + kL * ld_right
  }

  # log|V_aa|
  log_det_Vaa <- as.numeric(determinant(V_aa, logarithm = TRUE)$modulus)

  # log|Sigma_y|
  log_det_Sy <- 2 * (n_obs - q) * log(sigma_eps) + log_det_Sa + log_det_Vaa

  ll <- -0.5 * (n_obs * log(2 * pi) + log_det_Sy + quad)
  as.numeric(ll)
}

# ---------------------------------------------------------------------------
# logLik.vcmm_fit
# ---------------------------------------------------------------------------

#' Log-likelihood of a vcmm fit
#'
#' Returns the marginal log-likelihood
#' \deqn{\ell(\hat\beta, \hat\sigma_\varepsilon, \hat\Sigma_\alpha)
#'   = -\tfrac{n}{2}\log(2\pi)
#'     - \tfrac{1}{2}\log|\Sigma_y|
#'     - \tfrac{1}{2}(y - X\hat\beta)^{\top} \Sigma_y^{-1}(y - X\hat\beta),}
#' evaluated at the fitted parameter values, with
#' \eqn{\Sigma_y = \sigma_\varepsilon^2 I + Z\,\Sigma_\alpha\,Z^{\top}}.
#' The value is computed once at convergence and cached on the fit
#' object as \code{object$marginal_loglik}; this method simply retrieves
#' it and attaches \code{df} and \code{nobs} attributes so that
#' \code{AIC()} and \code{BIC()} work out of the box.
#'
#' Degrees of freedom counted are \eqn{p} (fixed-effects, including all
#' spline basis coefficients) plus the number of free variance-component
#' parameters:
#' \itemize{
#'   \item \code{re_cov = "diag"}: 2 (\eqn{\sigma_\varepsilon},
#'     \eqn{\sigma_\alpha}).
#'   \item \code{re_cov = "kronecker"} / \code{"separable"}:
#'     \eqn{1 + q_{\mathrm{left}}(q_{\mathrm{left}} + 1)/2}
#'     (\eqn{\sigma_\varepsilon} plus the free entries of
#'     \eqn{\Sigma_{\mathrm{left}}}). \eqn{\Sigma_{\mathrm{right}}} is held
#'     fixed at its user-supplied value and contributes \code{0} df.
#' }
#'
#' @param object A \code{vcmm_fit} object.
#' @param ... Unused.
#'
#' @return An object of class \code{"logLik"}; numeric scalar with
#'   \code{df} and \code{nobs} attributes.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' t <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 3), n, 3)
#' y <- 2 + sin(2 * pi * t) * x +
#'      as.vector(Z %*% rnorm(3, sd = 0.5)) + rnorm(n, sd = 0.5)
#' fit <- vcmm(y, X = x, Z = Z, t = t,
#'             control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
logLik.vcmm_fit <- function(object, ...) {
  ll <- object$marginal_loglik
  if (is.null(ll)) {
    stop("`object$marginal_loglik` is missing; refit with this version of cevcmm.",
         call. = FALSE)
  }

  # Variance-component df
  state <- object$re_cov_state
  n_vc  <- if (is.null(state) || identical(state$type, "diag")) {
    2L  # sigma_eps + sigma_alpha
  } else {
    kL <- as.integer(state$q_left)
    1L + as.integer(kL * (kL + 1L) / 2)
  }

  df_total <- as.integer(object$p) + n_vc

  out <- as.numeric(ll)
  attr(out, "df")   <- df_total
  attr(out, "nobs") <- as.integer(object$n_obs)
  class(out) <- "logLik"
  out
}

# ---------------------------------------------------------------------------
# predict.vcmm_fit
# ---------------------------------------------------------------------------

#' Predictions from a fitted VCMM
#'
#' Produces predicted responses at new \code{(X, t)} (and optionally
#' \code{Z}) values, using the paper's subject-specific predictor
#' \deqn{\hat y_i = \tilde{\mathbf x}_i^{\top} \hat{\tilde\beta}
#'                  + \mathbf z_i^{\top} \hat\alpha}
#' by default (Jalili and Lin, 2025, Section 5). For the alternative
#' "new groups not seen in training" scenario, pass
#' \code{include_random = FALSE} to use the marginal predictor
#' \eqn{\hat y_i = \tilde{\mathbf x}_i^{\top} \hat{\tilde\beta}}.
#'
#' \strong{newdata format.} A named list (or data frame whose columns
#' match these names) containing:
#' \itemize{
#'   \item \code{t}: numeric vector of length \eqn{N_{\mathrm{new}}}.
#'   \item \code{X}: numeric matrix \eqn{N_{\mathrm{new}} \times K} (or
#'     length-\eqn{N_{\mathrm{new}}} vector when \eqn{K = 1}) of varying-
#'     coefficient covariates, in the same column order used at fit time.
#'   \item \code{Z}: optional \eqn{N_{\mathrm{new}} \times q} random-
#'     effects design matrix. Must follow the same column-stacking
#'     convention as the training \code{Z} so that \code{Z \%*\% alpha}
#'     references the appropriate random-effect slots.
#' }
#'
#' \strong{Standard errors.} With \code{se.fit = TRUE} the per-prediction
#' standard error is the square root of
#' \eqn{[W, Z]\,\hat\sigma_\varepsilon^2 K^{-1}\,[W, Z]^{\top}} when
#' \code{include_random = TRUE} (joint uncertainty of fixed and random
#' effects), where \eqn{W} is the spline-expanded fixed-effects row;
#' the random-effect block is omitted when \code{include_random = FALSE}.
#' These are confidence-interval SEs on the mean; for prediction intervals
#' add \eqn{\hat\sigma_\varepsilon^2} to the variance.
#'
#' @param object A \code{vcmm_fit} object.
#' @param newdata A named list or data frame. See Details.
#' @param include_random Logical. If \code{TRUE} (default) and
#'   \code{newdata\$Z} is supplied, adds \eqn{Z\hat\alpha} to the
#'   prediction. Set \code{FALSE} for the marginal predictor (new groups
#'   scenario).
#' @param se.fit Logical. If \code{TRUE}, also returns per-prediction
#'   standard errors.
#' @param ... Unused.
#'
#' @return Either a numeric vector of length \eqn{N_{\mathrm{new}}}, or
#'   (when \code{se.fit = TRUE}) a list with components \code{fit} and
#'   \code{se.fit}.
#'
#' @seealso \code{\link{varying_coef}}, \code{\link{vcov.vcmm_fit}}.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' t <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 3), n, 3)
#' y <- 2 + sin(2 * pi * t) * x +
#'      as.vector(Z %*% rnorm(3, sd = 0.5)) + rnorm(n, sd = 0.5)
#'
#' fit <- vcmm(y, X = x, Z = Z, t = t,
#'             control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#'
#' # Default: subject-specific predictor (training-group prediction)
#' yhat_train <- predict(fit, newdata = list(t = t, X = x, Z = Z))
#' mean((y - yhat_train)^2)  # ~ sigma_eps^2 = 0.25
#'
#' # Marginal predictor (new groups scenario)
#' yhat_marg <- predict(fit, newdata = list(t = t, X = x, Z = Z),
#'                       include_random = FALSE)
predict.vcmm_fit <- function(object, newdata,
                              include_random = TRUE,
                              se.fit         = FALSE,
                              ...) {
  if (!inherits(object, "vcmm_fit")) {
    stop("`object` must be a 'vcmm_fit'.", call. = FALSE)
  }
  if (missing(newdata) || is.null(newdata)) {
    stop("predict.vcmm_fit requires `newdata`.", call. = FALSE)
  }
  if (!is.list(newdata)) {
    stop("`newdata` must be a named list or data frame.", call. = FALSE)
  }

  t_new <- newdata$t
  X_new <- newdata$X
  Z_new <- newdata$Z   # may be NULL

  if (is.null(t_new)) {
    stop("`newdata` must contain a numeric `t`.", call. = FALSE)
  }
  if (is.null(X_new)) {
    stop("`newdata` must contain `X` (the varying-coefficient covariates).",
         call. = FALSE)
  }
  if (is.vector(X_new) && is.numeric(X_new)) {
    X_new <- matrix(X_new, ncol = 1L)
  }
  if (!is.matrix(X_new) || !is.numeric(X_new)) {
    stop("`newdata$X` must be a numeric matrix (or vector when K = 1).",
         call. = FALSE)
  }
  n_new <- length(t_new)
  if (nrow(X_new) != n_new) {
    stop(sprintf("nrow(newdata$X) = %d does not match length(newdata$t) = %d.",
                 nrow(X_new), n_new),
         call. = FALSE)
  }
  if (anyNA(t_new) || anyNA(X_new)) {
    stop("NA values are not allowed in `newdata$t` or `newdata$X`.",
         call. = FALSE)
  }

  ds <- object$design
  if (is.null(ds)) {
    stop("`object$design` is missing; was the fit produced by vcmm()?",
         call. = FALSE)
  }

  if (ncol(X_new) != ds$K) {
    stop(sprintf("ncol(newdata$X) = %d does not match the fit's K = %d.",
                 ncol(X_new), ds$K),
         call. = FALSE)
  }

  # -------- Build X_design at newdata using the SAME basis as training --
  t_use <- if (isTRUE(ds$normalize_t)) {
    if (ds$t_max <= ds$t_min) {
      stop("Stored t-range has zero width.", call. = FALSE)
    }
    (t_new - ds$t_min) / (ds$t_max - ds$t_min)
  } else {
    t_new
  }

  B_new <- withCallingHandlers(
    splines::bs(
      t_use,
      degree         = as.integer(ds$degree),
      knots          = ds$internal_knots,
      Boundary.knots = ds$boundary_knots,
      intercept      = FALSE
    ),
    warning = function(w) {
      if (grepl("beyond boundary knots", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  B_new <- unclass(B_new); attributes(B_new) <- list(dim = dim(B_new))

  K <- as.integer(ds$K)
  blocks <- vector("list", K)
  for (k in seq_len(K)) {
    blocks[[k]] <- B_new * X_new[, k]
  }
  X_design_new <- cbind(1, do.call(cbind, blocks))   # n_new x (1 + K*m)

  if (ncol(X_design_new) != object$p) {
    stop(sprintf(
      "Internal: reconstructed design has %d columns, fit expects %d. ",
      ncol(X_design_new), object$p),
      "The basis spec stored on the fit may not match the original design.",
      call. = FALSE)
  }

  # -------- Fixed-effects contribution --------------------------------
  fit_part <- as.vector(X_design_new %*% object$beta)

  # -------- Random-effects contribution -------------------------------
  use_re <- isTRUE(include_random) && !is.null(Z_new)
  if (use_re) {
    if (!is.matrix(Z_new) || !is.numeric(Z_new)) {
      stop("`newdata$Z` must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(Z_new) != n_new) {
      stop(sprintf("nrow(newdata$Z) = %d does not match length(newdata$t) = %d.",
                   nrow(Z_new), n_new),
           call. = FALSE)
    }
    if (ncol(Z_new) != object$q) {
      stop(sprintf("ncol(newdata$Z) = %d does not match the fit's q = %d.",
                   ncol(Z_new), object$q),
           call. = FALSE)
    }
    if (anyNA(Z_new)) {
      stop("NA values are not allowed in `newdata$Z`.", call. = FALSE)
    }
    fit_part <- fit_part + as.vector(Z_new %*% object$alpha)
  }

  if (!isTRUE(se.fit)) {
    return(fit_part)
  }

  # -------- Standard errors -------------------------------------------
  # Per-row variance: w' (sigma_eps^2 * K_inv) w where w is the joint
  # row [X_design_new_i, Z_new_i] (or just X_design_new_i when no RE).
  if (is.null(object$K_inv)) {
    stop("`object$K_inv` is missing; cannot compute se.fit.", call. = FALSE)
  }
  se2  <- object$sigma_eps^2
  if (use_re) {
    W <- cbind(X_design_new, Z_new)   # n_new x (p + q)
    # var per row = rowSums((W %*% K_inv) * W)
    var_row <- rowSums((W %*% object$K_inv) * W) * se2
  } else {
    V_beta <- se2 * object$K_inv[seq_len(object$p),
                                 seq_len(object$p), drop = FALSE]
    var_row <- rowSums((X_design_new %*% V_beta) * X_design_new)
  }
  se_row <- sqrt(pmax(var_row, 0))

  list(fit = fit_part, se.fit = se_row)
}
