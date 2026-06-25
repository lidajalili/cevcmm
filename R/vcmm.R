#===============================================================================
# vcmm: main user-facing fit function for the Communication-Efficient
# Varying Coefficient Mixed-Effects Model
#
# Wraps:
#   build_vcmm_design()   -> X_design, penalty
#   compute_sufficient_stats() -> aggregated SS
#   fit_ss() or fit_csl() -> the chosen estimator
#
# Returns a vcmm_fit object enriched with design metadata so that
# downstream predict() / coef() / plot() (coming Days 11-13) can
# reconstruct the basis at new t values.
#===============================================================================

#' Fit a Varying Coefficient Mixed-Effects Model
#'
#' The main user-facing fit function. Builds the B-spline design and
#' penalty for one or more varying coefficients in \code{t}, computes
#' the aggregated sufficient statistics, and fits the model using either
#' the one-step CSL estimator (default) or the iterative SS estimator.
#'
#' \strong{Model.} For observations \eqn{i = 1, \ldots, n} the fitted
#' model is
#' \deqn{
#'   y_i
#'   = \beta_0(t_i) + \sum_{k=1}^{K} x_{ik}\, \beta_k(t_i)
#'   + z_i^\top \alpha + \varepsilon_i,
#'   \quad
#'   \alpha \sim N(0, \sigma_\alpha^2 I_q),
#'   \quad
#'   \varepsilon_i \sim N(0, \sigma_\varepsilon^2),
#' }
#' where each \eqn{\beta_k(t)} is a cubic B-spline with \code{n_basis}
#' basis functions and a second-order difference penalty.
#'
#' \strong{Method selection.} The default is \code{method = "auto"},
#' which picks \code{"csl"} when \eqn{N \cdot q > 10^5} or \eqn{q > 50}
#' and \code{"ss"} otherwise. The CSL estimator is first-order
#' equivalent to SS (Theorem 3.1 of the paper) while needing only one
#' Newton step from a \eqn{\sqrt{N}}-consistent pilot; it dominates at
#' scale. SS is the exact-MLE benchmark and the cleaner choice for
#' small problems. Use \code{method = "csl"} or \code{method = "ss"}
#' to override.
#'
#' \strong{Random-effects covariance.} Two structures are supported:
#' \itemize{
#'   \item \code{re_cov = "diag"} (default): \eqn{\alpha \sim N(0,
#'     \sigma_\alpha^2 I_q)}. Smallest, fastest.
#'   \item \code{re_cov = "kronecker"}: \eqn{\alpha \sim N(0, \Sigma_{2\times 2}
#'     \otimes \Sigma_{\text{spatial}})} for origin-destination data.
#'     Each of \code{n_groups} districts has a 2-dimensional random
#'     effect (origin, destination), so \code{ncol(Z)} must equal
#'     \code{2 * n_groups}. When \code{control$update_variance = TRUE},
#'     the 2-by-2 origin-destination cross-effect covariance
#'     \code{Sigma_2x2} is re-estimated at every iteration from the
#'     column covariance of \code{matrix(alpha_hat, G, 2)}. The G-by-G
#'     spatial covariance \code{Sigma_spatial} is \strong{held fixed}
#'     at the user-supplied initial value (default \code{I_G}). The
#'     reason is identifiability: with a single \eqn{\hat\alpha}, the
#'     moment-based spatial estimator is rank-2 at most and iterating
#'     it is unstable. For real OD data, supply a parametric spatial
#'     kernel via \code{Sigma_spatial_init} (e.g.
#'     \code{exp(-D / phi)} for a known distance matrix \code{D}).
#' }
#' (A third option, \code{"separable"} for group-shared dense random
#' effects, is planned for a later release.)
#'
#' @param y Numeric response vector of length \eqn{n}.
#' @param X Numeric \eqn{n \times K} matrix (or length-\eqn{n} vector if
#'   \eqn{K = 1}) of covariates that get varying coefficients in
#'   \code{t}.
#' @param Z Numeric \eqn{n \times q} random-effects design matrix.
#'   For \code{re_cov = "kronecker"}, \eqn{q} must equal
#'   \code{2 * n_groups}.
#' @param t Numeric vector of length \eqn{n} in which the coefficients
#'   vary smoothly.
#' @param method Character: \code{"auto"} (default), \code{"csl"}, or
#'   \code{"ss"}. See Details.
#' @param re_cov Character: random-effects covariance structure;
#'   \code{"diag"} (default) or \code{"kronecker"}. See Details.
#' @param n_groups Integer. Number of groups for
#'   \code{re_cov = "kronecker"} only. Must satisfy
#'   \code{ncol(Z) == 2 * n_groups}. Ignored otherwise.
#' @param Sigma_2x2_init Optional 2 by 2 initial Sigma_2x2 (only used
#'   for \code{re_cov = "kronecker"}). Defaults to
#'   \code{sigma_alpha^2 * I_2}. With \code{update_variance = TRUE} the
#'   final estimate is data-driven; the initial value only affects
#'   iteration 1.
#' @param Sigma_spatial_init Optional G by G initial Sigma_spatial
#'   (only used for \code{re_cov = "kronecker"}). Defaults to
#'   \code{I_G}. With \code{update_variance = TRUE} the final estimate
#'   is data-driven; the initial value only affects iteration 1.
#' @param n_basis Integer or \code{NULL}. Number of B-spline basis
#'   functions per varying coefficient. \code{NULL} (default)
#'   auto-picks \code{max(floor(n^(1/3)) + 4, 10)}.
#' @param degree Integer. B-spline degree (default 3 = cubic).
#' @param lambda Non-negative numeric. Smoothing parameter for the
#'   penalty (default 1).
#' @param control A \code{vcmm_control} object with fitting options.
#'   Pass \code{vcmm_control()} for defaults.
#' @param normalize_t Logical. If \code{TRUE} (default), \code{t} is
#'   linearly mapped to \code{[0, 1]} before building the basis.
#' @param ... Further arguments passed to \code{fit_csl()} (e.g.
#'   \code{pilot_max_iter}) when CSL is used. Ignored under SS.
#'
#' @return A \code{vcmm_fit} object (same class as the output of
#'   \code{fit_ss()} / \code{fit_csl()}) with an additional
#'   \code{design} element carrying the spline knots, degree, and
#'   \code{K} so that downstream prediction can reconstruct the basis,
#'   and a \code{K_inv} element holding the cached inverse Hessian
#'   used by \code{vcov()} and \code{summary()}.
#'
#' @references
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' # ---- Single varying coefficient ---------------------------------------
#' set.seed(1)
#' n <- 500
#' t  <- runif(n)
#' x  <- runif(n)
#' Z  <- matrix(rnorm(n * 3), n, 3)
#' alpha_true <- rnorm(3, sd = 0.5)
#' y  <- 2 + sin(2 * pi * t) * x +
#'       as.vector(Z %*% alpha_true) + rnorm(n, sd = 0.5)
#'
#' fit <- vcmm(y, X = x, Z = Z, t = t,
#'             control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' fit
#'
#' # ---- Force the SS estimator -------------------------------------------
#' fit_ss_ref <- vcmm(y, X = x, Z = Z, t = t, method = "ss",
#'                    control = vcmm_control(sigma_eps = 0.5,
#'                                           sigma_alpha = 0.5))
vcmm <- function(y,
                 X,
                 Z,
                 t,
                 method             = c("auto", "csl", "ss"),
                 re_cov             = c("diag", "kronecker", "separable"),
                 n_groups           = NULL,
                 Sigma_2x2_init     = NULL,
                 Sigma_spatial_init = NULL,
                 n_basis            = NULL,
                 degree             = 3L,
                 lambda             = 1,
                 control            = vcmm_control(),
                 normalize_t        = TRUE,
                 ...) {

  # ----- Argument matching and basic validation ---------------------------
  method <- match.arg(method)
  re_cov <- match.arg(re_cov)

  if (re_cov == "separable") {
    stop("re_cov = 'separable' is not yet implemented; planned for a later release.",
         call. = FALSE)
  }

  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("NA values are not allowed in `y`.", call. = FALSE)
  }
  n <- length(y)

  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(Z) != n) {
    stop(sprintf("nrow(Z) = %d does not match length(y) = %d.",
                 nrow(Z), n),
         call. = FALSE)
  }
  if (anyNA(Z)) {
    stop("NA values are not allowed in `Z`.", call. = FALSE)
  }
  q_dim <- ncol(Z)

  # ----- re_cov-specific validation and initial state --------------------
  re_cov_state <- NULL
  if (re_cov == "kronecker") {
    if (is.null(n_groups) ||
        !is.numeric(n_groups) || length(n_groups) != 1L ||
        n_groups < 1L || n_groups != as.integer(n_groups)) {
      stop("re_cov = 'kronecker' requires `n_groups` (a single positive integer).",
           call. = FALSE)
    }
    n_groups <- as.integer(n_groups)
    if (q_dim != 2L * n_groups) {
      stop(sprintf(
        "For re_cov = 'kronecker', ncol(Z) must equal 2 * n_groups; got ncol(Z) = %d and n_groups = %d.",
        q_dim, n_groups),
        call. = FALSE)
    }

    # Default initial Sigma matrices: identity for both
    # (steady-state value will be estimated from data if update_variance = TRUE)
    if (is.null(Sigma_2x2_init)) {
      Sigma_2x2_init <- (control$sigma_alpha^2) * diag(2L)
    }
    if (is.null(Sigma_spatial_init)) {
      Sigma_spatial_init <- diag(n_groups)
    }
    if (!isTRUE(all.equal(dim(Sigma_2x2_init), c(2L, 2L)))) {
      stop("`Sigma_2x2_init` must be a 2 by 2 matrix.", call. = FALSE)
    }
    if (!isTRUE(all.equal(dim(Sigma_spatial_init),
                          c(n_groups, n_groups)))) {
      stop(sprintf("`Sigma_spatial_init` must be %d by %d.",
                   n_groups, n_groups),
           call. = FALSE)
    }

    re_cov_state <- .new_re_cov_state(
      type          = "kronecker",
      Sigma_2x2     = Sigma_2x2_init,
      Sigma_spatial = Sigma_spatial_init,
      n_groups      = n_groups
    )

    if (!isTRUE(control$update_variance)) {
      message(
        "Note: re_cov = 'kronecker' is typically run with control$update_variance = TRUE ",
        "so that Sigma_2x2 and Sigma_spatial are estimated from the data. ",
        "Pass vcmm_control(..., update_variance = TRUE) to enable iterative estimation."
      )
    }
  }

  # ----- Resolve method = "auto" -----------------------------------------
  if (method == "auto") {
    if (n * q_dim > 1e5 || q_dim > 50L) {
      method_used <- "csl"
    } else {
      method_used <- "ss"
    }
    if (isTRUE(control$verbose)) {
      cat(sprintf("  [vcmm] method='auto' resolved to '%s' (N=%d, q=%d)\n",
                  method_used, n, q_dim))
    }
  } else {
    method_used <- method
  }

  # ----- Build design and penalty -----------------------------------------
  design <- build_vcmm_design(X           = X,
                              t           = t,
                              n_basis     = n_basis,
                              degree      = degree,
                              lambda      = lambda,
                              normalize_t = normalize_t)

  if (nrow(design$X_design) != n) {
    stop(sprintf("Internal: design rows (%d) != length(y) (%d).",
                 nrow(design$X_design), n),
         call. = FALSE)
  }

  # ----- Aggregated sufficient statistics ---------------------------------
  stats_obj <- compute_sufficient_stats(y, design$X_design, Z)

  # ----- Dispatch to the chosen estimator ---------------------------------
  fit <- switch(
    method_used,
    ss  = fit_ss (stats_obj, design$penalty, control,
                  re_cov_state = re_cov_state),
    csl = fit_csl(stats_obj, design$penalty, control,
                  re_cov_state = re_cov_state, ...)
  )

  # ----- Attach design metadata for downstream predict/plot ---------------
  fit$design <- design[c("internal_knots", "boundary_knots",
                         "degree", "n_basis", "K", "lambda",
                         "normalize_t", "t_min", "t_max")]
  fit$re_cov <- re_cov
  fit$call <- match.call()
  fit
}
