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
#' \strong{Default method.} The default is \code{method = "csl"},
#' matching the headline algorithm of Lin and Jalili (2026). By
#' Theorem 3.1, the CSL estimator is first-order equivalent to the
#' full-data SS estimator while needing only one Newton step from a
#' \eqn{\sqrt{N}}-consistent pilot. Use \code{method = "ss"} to compute
#' the exact-MLE benchmark; for large \eqn{N \cdot q} this can be much
#' slower without changing the answer beyond \eqn{O(1/\sqrt N)}.
#'
#' @param y Numeric response vector of length \eqn{n}.
#' @param X Numeric \eqn{n \times K} matrix (or length-\eqn{n} vector if
#'   \eqn{K = 1}) of covariates that get varying coefficients in
#'   \code{t}.
#' @param Z Numeric \eqn{n \times q} random-effects design matrix.
#' @param t Numeric vector of length \eqn{n} in which the coefficients
#'   vary smoothly.
#' @param method Character: \code{"csl"} (default) for the one-step
#'   communication-efficient estimator, or \code{"ss"} for the iterative
#'   sufficient-statistics estimator (exact MLE benchmark).
#' @param n_basis Integer or \code{NULL}. Number of B-spline basis
#'   functions per varying coefficient. \code{NULL} (default) auto-picks
#'   \code{max(floor(n^(1/3)) + 4, 10)}.
#' @param degree Integer. B-spline degree (default 3 = cubic).
#' @param lambda Non-negative numeric. Smoothing parameter for the
#'   penalty (default 1).
#' @param control A \code{vcmm_control} object with fitting options.
#'   Pass \code{vcmm_control()} for defaults.
#' @param normalize_t Logical. If \code{TRUE} (default), \code{t} is
#'   linearly mapped to \code{[0, 1]} before building the basis.
#' @param ... Further arguments passed to \code{fit_csl()} (e.g.
#'   \code{pilot_max_iter}) when \code{method = "csl"}. Ignored when
#'   \code{method = "ss"}.
#'
#' @return A \code{vcmm_fit} object (same class as the output of
#'   \code{fit_ss()} / \code{fit_csl()}) with an additional
#'   \code{design} element carrying the spline knots, degree, and
#'   \code{K} so that downstream prediction can reconstruct the basis.
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
#' # ---- Two varying coefficients -----------------------------------------
#' x1 <- runif(n); x2 <- runif(n)
#' y2 <- 2 + sin(2 * pi * t) * x1 + cos(2 * pi * t) * x2 +
#'       as.vector(Z %*% alpha_true) + rnorm(n, sd = 0.5)
#'
#' fit2 <- vcmm(y2, X = cbind(x1, x2), Z = Z, t = t,
#'              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' fit2
#'
#' # ---- Exact-MLE benchmark via SS ---------------------------------------
#' fit_ss_ref <- vcmm(y, X = x, Z = Z, t = t, method = "ss",
#'                    control = vcmm_control(sigma_eps = 0.5,
#'                                           sigma_alpha = 0.5))
vcmm <- function(y,
                 X,
                 Z,
                 t,
                 method      = c("csl", "ss"),
                 n_basis     = NULL,
                 degree      = 3L,
                 lambda      = 1,
                 control     = vcmm_control(),
                 normalize_t = TRUE,
                 ...) {

  # ----- Argument matching and basic validation ---------------------------
  method <- match.arg(method)

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

  # ----- Build design and penalty -----------------------------------------
  design <- build_vcmm_design(X           = X,
                              t           = t,
                              n_basis     = n_basis,
                              degree      = degree,
                              lambda      = lambda,
                              normalize_t = normalize_t)

  # nrow check after design (design also validates length(t) == nrow(X))
  if (nrow(design$X_design) != n) {
    stop(sprintf("Internal: design rows (%d) != length(y) (%d).",
                 nrow(design$X_design), n),
         call. = FALSE)
  }

  # ----- Aggregated sufficient statistics ---------------------------------
  stats <- compute_sufficient_stats(y, design$X_design, Z)

  # ----- Dispatch to the chosen estimator ---------------------------------
  fit <- switch(
    method,
    ss  = fit_ss (stats, design$penalty, control),
    csl = fit_csl(stats, design$penalty, control, ...)
  )

  # ----- Attach design metadata for downstream predict/plot ---------------
  fit$design <- design[c("internal_knots", "boundary_knots",
                         "degree", "n_basis", "K", "lambda",
                         "normalize_t", "t_min", "t_max")]
  fit$call <- match.call()
  fit
}
