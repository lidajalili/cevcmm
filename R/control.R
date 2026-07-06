#===============================================================================
# vcmm_control: validated fitting options for VCMM estimators
#
# A single options object passed to fit_ss() (and later fit_csl(), fit_svd()).
# Centralises every fitting knob so user-facing functions take few arguments.
#===============================================================================

#' Control parameters for VCMM fitting
#'
#' Builds a validated options list controlling the iterative SS estimator
#' (and, later, the CSL and SVD-stabilised estimators). Pass the returned
#' object as the \code{control} argument of \code{fit_ss()}.
#'
#' @param max_iter Integer. Maximum number of iterations (default 200).
#' @param tol_beta Positive numeric. Relative-change convergence tolerance
#'   for the fixed-effects coefficient vector \code{beta} (default 1e-6).
#' @param tol_alpha Positive numeric. Relative-change convergence tolerance
#'   for the random-effects vector \code{alpha} (default 1e-6).
#' @param sigma_eps Positive numeric. Initial residual standard deviation.
#'   If \code{update_variance = FALSE}, this value is held fixed throughout
#'   fitting (default 1).
#' @param sigma_alpha Positive numeric. Initial random-effect standard
#'   deviation. If \code{update_variance = FALSE}, this value is held fixed
#'   throughout fitting (default 1).
#' @param update_variance Logical. If \code{FALSE} (default), \code{sigma_eps}
#'   and \code{sigma_alpha} are held fixed at the supplied values --
#'   matching Algorithm 1 of Jalili and Lin (2025) as written. If
#'   \code{TRUE}, both are re-estimated at every iteration using the
#'   residual sum of squares formula (for \code{sigma_eps}) and a
#'   method-of-moments update (for \code{sigma_alpha}).
#' @param verbose Logical. If \code{TRUE}, print progress every 20
#'   iterations (default \code{FALSE}).
#'
#' @return A list of class \code{"vcmm_control"} containing the validated
#'   options.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' # Defaults: fix variances at 1, iterate up to 200 times.
#' ctrl <- vcmm_control()
#' ctrl
#'
#' # Fix variances at user-supplied values.
#' ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5)
#'
#' # Re-estimate variances each iteration.
#' ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5,
#'                      update_variance = TRUE)
vcmm_control <- function(max_iter        = 200L,
                         tol_beta        = 1e-6,
                         tol_alpha       = 1e-6,
                         sigma_eps       = 1.0,
                         sigma_alpha     = 1.0,
                         update_variance = FALSE,
                         verbose         = FALSE) {

  # --- Validation -----------------------------------------------------------
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop("`max_iter` must be a single positive integer.", call. = FALSE)
  }
  .check_pos <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop(sprintf("`%s` must be a single positive finite number.", name),
           call. = FALSE)
    }
  }
  .check_pos(tol_beta,    "tol_beta")
  .check_pos(tol_alpha,   "tol_alpha")
  .check_pos(sigma_eps,   "sigma_eps")
  .check_pos(sigma_alpha, "sigma_alpha")

  if (!is.logical(update_variance) || length(update_variance) != 1L ||
      is.na(update_variance)) {
    stop("`update_variance` must be a single TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE or FALSE.", call. = FALSE)
  }

  out <- list(
    max_iter        = as.integer(max_iter),
    tol_beta        = tol_beta,
    tol_alpha       = tol_alpha,
    sigma_eps       = sigma_eps,
    sigma_alpha     = sigma_alpha,
    update_variance = update_variance,
    verbose         = verbose
  )
  class(out) <- c("vcmm_control", "list")
  out
}

#' @rdname vcmm_control
#' @param x A \code{vcmm_control} object.
#' @param ... Unused.
#' @export
print.vcmm_control <- function(x, ...) {
  cat("<vcmm_control>  fitting options\n")
  cat(sprintf("  max_iter        : %d\n",    x$max_iter))
  cat(sprintf("  tol_beta        : %.2e\n",  x$tol_beta))
  cat(sprintf("  tol_alpha       : %.2e\n",  x$tol_alpha))
  cat(sprintf("  sigma_eps       : %.4f\n",  x$sigma_eps))
  cat(sprintf("  sigma_alpha     : %.4f\n",  x$sigma_alpha))
  cat(sprintf("  update_variance : %s\n",    x$update_variance))
  cat(sprintf("  verbose         : %s\n",    x$verbose))
  invisible(x)
}
