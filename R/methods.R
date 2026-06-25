#===============================================================================
# S3 methods for vcmm_fit: vcov, summary, print.summary
#
# vcov.vcmm_fit  -> variance-covariance matrix of the fixed-effects beta
# summary.vcmm_fit -> a vcmm_summary object holding the coefficient table
# print.vcmm_summary -> pretty printing of summary
#===============================================================================

#' Variance-covariance matrix of the fixed-effects from a vcmm fit
#'
#' Returns the asymptotic variance-covariance matrix of the
#' fixed-effects coefficient vector \eqn{\hat\beta} from a fitted
#' \code{vcmm_fit}. The matrix is computed as
#' \deqn{
#'   \widehat{\mathrm{Var}}(\hat\beta)
#'   = \hat\sigma_\varepsilon^2 \cdot [K^{-1}]_{1:p,\, 1:p},
#' }
#' where \eqn{K} is the prior-augmented Hessian assembled at
#' convergence and cached in \code{object$K_inv}. This is the standard
#' plug-in asymptotic-normal variance estimator for the linear normal
#' VCMM with fixed variance components.
#'
#' Pass \code{which = "alpha"} for the random-effect block,
#' \code{which = "both"} for the full \eqn{(p+q) \times (p+q)} joint
#' matrix.
#'
#' @param object A \code{vcmm_fit} object.
#' @param which Character: \code{"beta"} (default), \code{"alpha"}, or
#'   \code{"both"}.
#' @param ... Unused.
#'
#' @return A numeric matrix:
#' \itemize{
#'   \item \code{"beta"}: p by p.
#'   \item \code{"alpha"}: q by q.
#'   \item \code{"both"}: (p+q) by (p+q), joint.
#' }
#'
#' @references
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' t <- runif(n); x <- runif(n)
#' Z <- matrix(rnorm(n * 3), n, 3)
#' y <- 2 + sin(2 * pi * t) * x +
#'      as.vector(Z %*% rnorm(3, sd = 0.5)) + rnorm(n, sd = 0.5)
#'
#' fit <- vcmm(y, X = x, Z = Z, t = t,
#'             control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#'
#' V_beta <- vcov(fit)
#' dim(V_beta)
vcov.vcmm_fit <- function(object, which = c("beta", "alpha", "both"), ...) {
  which <- match.arg(which)
  if (is.null(object$K_inv)) {
    stop("`object$K_inv` is missing; was the fit produced by this version of cevcmm?",
         call. = FALSE)
  }
  p   <- object$p
  q   <- object$q
  se2 <- object$sigma_eps^2

  switch(
    which,
    beta  = se2 * object$K_inv[seq_len(p), seq_len(p), drop = FALSE],
    alpha = se2 * object$K_inv[seq.int(p + 1L, p + q),
                               seq.int(p + 1L, p + q), drop = FALSE],
    both  = se2 * object$K_inv
  )
}

#-------------------------------------------------------------------------------
# summary
#-------------------------------------------------------------------------------

#' Summarise a vcmm fit
#'
#' Produces a \code{vcmm_summary} object with a coefficient table for
#' the fixed-effects, the variance components, and basic fit
#' diagnostics. \code{print()} renders it in the style of \code{lm} /
#' \code{lmer}.
#'
#' The standard errors are the square roots of the diagonal of
#' \code{vcov(object, which = "beta")}, treating \eqn{\sigma_\varepsilon}
#' and \eqn{\sigma_\alpha} as fixed at their fitted values. They are
#' first-order valid (Theorem 3.1) but do not adjust for variance
#' component uncertainty.
#'
#' @param object A \code{vcmm_fit} object.
#' @param ... Unused.
#'
#' @return A list of class \code{"vcmm_summary"}.
#'
#' @export
summary.vcmm_fit <- function(object, ...) {
  V_beta <- vcov.vcmm_fit(object, which = "beta")
  se     <- sqrt(pmax(diag(V_beta), 0))
  beta   <- object$beta

  # Coefficient table; under H_0: beta_j = 0 we use a Wald z-statistic
  z      <- ifelse(se > 0, beta / se, NA_real_)
  pval   <- 2 * stats::pnorm(-abs(z))

  tab <- data.frame(
    Estimate   = beta,
    `Std. Error` = se,
    `z value`  = z,
    `Pr(>|z|)` = pval,
    check.names = FALSE
  )

  # Friendly row names: intercept, then per-covariate basis index
  design <- object$design
  if (!is.null(design)) {
    K       <- design$K
    n_basis <- design$n_basis
    if (K * n_basis + 1L == length(beta)) {
      rn <- character(length(beta))
      rn[1] <- "(Intercept)"
      for (k in seq_len(K)) {
        idx <- (1L + (k - 1L) * n_basis + 1L):(1L + k * n_basis)
        rn[idx] <- sprintf("beta_%d(t)[%d]", k, seq_len(n_basis))
      }
      rownames(tab) <- rn
    }
  }

  out <- list(
    call        = object$call,
    method      = object$method,
    coefficients = tab,
    sigma_eps   = object$sigma_eps,
    sigma_alpha = object$sigma_alpha,
    n_obs       = object$n_obs,
    p           = object$p,
    q           = object$q,
    iterations  = object$iterations,
    converged   = object$converged,
    re_cov      = if (!is.null(object$re_cov)) object$re_cov else "diag"
  )
  class(out) <- c("vcmm_summary", "list")
  out
}

#' @rdname summary.vcmm_fit
#' @param x A \code{vcmm_summary} object.
#' @param digits Integer. Number of significant digits to display.
#' @param signif.stars Logical. If \code{TRUE}, append significance
#'   stars to the coefficient table (\code{R}'s default convention).
#' @export
print.vcmm_summary <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  cat("VCMM fit\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  cat(sprintf("\nMethod:  %s   |  RE covariance: %s   |  N = %d   |  p = %d   |  q = %d\n",
              x$method, x$re_cov, x$n_obs, x$p, x$q))
  cat(sprintf("Converged: %s   |  iterations: %d\n",
              x$converged, x$iterations))

  cat("\nVariance components:\n")
  cat(sprintf("  sigma_eps   = %.4f\n", x$sigma_eps))
  cat(sprintf("  sigma_alpha = %.4f\n", x$sigma_alpha))

  cat("\nFixed-effects coefficients:\n")
  stats::printCoefmat(x$coefficients,
                      digits = digits,
                      signif.stars = signif.stars,
                      P.values = TRUE,
                      has.Pvalue = TRUE,
                      na.print = "NA")
  invisible(x)
}
