#===============================================================================
# S3 methods for vcmm_fit
#
# Existing (Day 4-era):
#   vcov.vcmm_fit         -- joint or marginal asymptotic covariance
#   summary.vcmm_fit      -- coefficient table + diagnostics object
#   print.vcmm_summary    -- pretty printer for the above
#
# Added on Day 11:
#   nobs.vcmm_fit         -- standard generic
#   coef.vcmm_fit         -- named numeric vector of beta
#   fixef / fixef.vcmm_fit -- intercept + reshaped basis-coef matrix
#   ranef / ranef.vcmm_fit -- alpha reshaped to match re_cov structure
#   varying_coef          -- evaluate beta_k(t) at new t values
#
# Refreshed on Day 11:
#   summary.vcmm_fit     -- now propagates re_cov_state
#   print.vcmm_summary   -- displays Sigma_left for kronecker / separable
#                          rather than the stale sigma_alpha
#
# Notes
# -----
# `fixef` and `ranef` exist as generics in nlme and lme4. To avoid a
# dependency, we define them here. Users with nlme/lme4 loaded after
# cevcmm should call cevcmm::fixef(fit) explicitly to avoid masking.
#===============================================================================

# ----------------------------------------------------------------------------
# vcov
# ----------------------------------------------------------------------------

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
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
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

# ----------------------------------------------------------------------------
# nobs
# ----------------------------------------------------------------------------

#' Number of observations from a vcmm fit
#'
#' Standard \code{stats::nobs} generic. Returns \eqn{N}, the total number
#' of observations used to compute the fit (or the sum of node sample
#' sizes for fits produced via \code{\link{fit_from_summaries}}).
#'
#' @param object A \code{vcmm_fit} object.
#' @param ... Unused.
#'
#' @return Integer.
#' @export
nobs.vcmm_fit <- function(object, ...) {
  if (is.null(object$n_obs)) return(NA_integer_)
  as.integer(object$n_obs)
}

# ----------------------------------------------------------------------------
# coef
# ----------------------------------------------------------------------------

#' Fixed-effects coefficient vector from a vcmm fit
#'
#' Returns \eqn{\hat\beta} as a named numeric vector. Under the package's
#' design (\code{X_design = cbind(1, B*x_1, ..., B*x_K)}, see
#' \code{\link{build_vcmm_design}}), entry 1 is the constant intercept
#' and entries \code{(1 + (k-1)*m + 1):(1 + k*m)} are the spline basis
#' coefficients for the \eqn{k}-th varying coefficient
#' \eqn{\beta_k(t)}.
#'
#' To evaluate \eqn{\beta_k(t)} at user-supplied \code{t} values, use
#' \code{\link{varying_coef}}. To get the same vector reshaped into a
#' (basis x covariate) matrix with the intercept split out, use
#' \code{\link{fixef.vcmm_fit}}.
#'
#' @param object A \code{vcmm_fit} object.
#' @param ... Unused.
#'
#' @return Named numeric vector of length \eqn{p = 1 + K \cdot m}.
#' @export
coef.vcmm_fit <- function(object, ...) {
  beta <- as.numeric(object$beta)
  names(beta) <- .build_beta_names(object$design, length(beta))
  beta
}

# ----------------------------------------------------------------------------
# fixef -- generic + method
# ----------------------------------------------------------------------------

#' Extract fixed-effects from a fitted model object
#'
#' Generic in the style of \code{nlme::fixef} / \code{lme4::fixef},
#' redefined here so \pkg{cevcmm} avoids a hard dependency on either.
#' If you also have \pkg{nlme} or \pkg{lme4} loaded, call
#' \code{cevcmm::fixef(fit)} explicitly.
#'
#' @param object A model object.
#' @param ... Method-specific arguments.
#'
#' @export
fixef <- function(object, ...) {
  UseMethod("fixef")
}

#' Fixed effects of a VCMM, reshaped by varying-coefficient block
#'
#' Splits the coefficient vector returned by \code{\link{coef.vcmm_fit}}
#' into:
#' \itemize{
#'   \item \code{intercept}: the constant scalar \eqn{\hat\beta_0}.
#'   \item \code{varying}: an \eqn{m \times K} matrix of B-spline basis
#'     coefficients, with row names \code{"basis1", ..., "basisM"} and
#'     column names \code{"X1", ..., "XK"}.
#' }
#' For \eqn{K = 0} (no varying covariate; intercept-only model), the
#' \code{varying} slot is \code{NULL}.
#'
#' @param object A \code{vcmm_fit} object.
#' @param ... Unused.
#'
#' @return A two-element list.
#'
#' @seealso \code{\link{coef.vcmm_fit}}, \code{\link{varying_coef}},
#'   \code{\link{ranef.vcmm_fit}}.
#'
#' @export
fixef.vcmm_fit <- function(object, ...) {
  beta <- coef(object)
  ds   <- object$design
  if (is.null(ds) || is.null(ds$K) || ds$K == 0L) {
    return(list(intercept = unname(beta[1L]), varying = NULL))
  }
  K <- as.integer(ds$K)
  m <- as.integer(ds$n_basis)
  if (length(beta) != 1L + K * m) {
    # Layout doesn't match expected; bail out gracefully.
    return(list(intercept = unname(beta[1L]), varying = NULL))
  }

  varying <- matrix(beta[-1L], nrow = m, ncol = K)
  rownames(varying) <- paste0("basis", seq_len(m))
  colnames(varying) <- paste0("X",     seq_len(K))

  list(intercept = unname(beta[1L]), varying = varying)
}

# ----------------------------------------------------------------------------
# ranef -- generic + method
# ----------------------------------------------------------------------------

#' Extract random effects from a fitted model object
#'
#' Generic in the style of \code{nlme::ranef} / \code{lme4::ranef},
#' redefined here so \pkg{cevcmm} avoids a hard dependency on either.
#' If you also have \pkg{nlme} or \pkg{lme4} loaded, call
#' \code{cevcmm::ranef(fit)} explicitly.
#'
#' @param object A model object.
#' @param ... Method-specific arguments.
#'
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' Random effects of a VCMM, reshaped by re_cov structure
#'
#' For \code{re_cov = "diag"}, returns a named numeric vector of length
#' \eqn{q}. For \code{"kronecker"} and \code{"separable"}, reshapes
#' \eqn{\hat\alpha} into a \eqn{G \times q_{\mathrm{left}}} matrix
#' (column-stacking convention used throughout the package). Row names
#' are \code{"g1", ..., "gG"}; column names are \code{c("origin",
#' "dest")} when \code{re_cov = "kronecker"} and \code{q_left = 2}, and
#' \code{"k1", ..., "kK"} otherwise.
#'
#' @param object A \code{vcmm_fit} object.
#' @param ... Unused.
#'
#' @return Numeric vector (\code{"diag"}) or numeric matrix
#'   (\code{"kronecker"} / \code{"separable"}).
#'
#' @seealso \code{\link{coef.vcmm_fit}}, \code{\link{fixef.vcmm_fit}}.
#'
#' @export
ranef.vcmm_fit <- function(object, ...) {
  alpha <- as.numeric(object$alpha)
  state <- object$re_cov_state

  if (is.null(state) || identical(state$type, "diag")) {
    names(alpha) <- paste0("a", seq_along(alpha))
    return(alpha)
  }

  G <- as.integer(state$n_groups)
  k <- as.integer(state$q_left)

  if (length(alpha) != k * G) {
    # Unexpected layout: return flat vector with a warning.
    warning("length(alpha) does not match q_left * n_groups; returning flat vector.",
            call. = FALSE)
    names(alpha) <- paste0("a", seq_along(alpha))
    return(alpha)
  }

  M <- matrix(alpha, nrow = G, ncol = k)
  rownames(M) <- paste0("g", seq_len(G))
  colnames(M) <- if (identical(state$type, "kronecker") && identical(k, 2L)) {
    c("origin", "dest")
  } else {
    paste0("k", seq_len(k))
  }
  M
}

# ----------------------------------------------------------------------------
# varying_coef -- evaluate beta_k(t) at user-supplied t values
# ----------------------------------------------------------------------------

#' Evaluate the varying coefficient(s) at new t values
#'
#' For a VCMM with varying coefficients \eqn{\beta_k(t), k = 1, ..., K},
#' this function evaluates \eqn{\hat\beta_k(t)} at any vector of
#' \code{t_new} values, using the same B-spline basis the fit was built
#' on. Useful for diagnostic plots, predictions, and reporting; the
#' Day-13 plot method uses it internally.
#'
#' The constant intercept \eqn{\hat\beta_0} is \emph{not} returned (it
#' does not vary in \code{t}); use \code{\link{fixef}(fit)\$intercept}
#' for that.
#'
#' Pointwise standard errors are available via \code{se.fit = TRUE},
#' computed as
#' \deqn{
#'   \mathrm{SE}(\hat\beta_k(t)) =
#'   \sqrt{B(t)^\top \widehat{\mathrm{Var}}(\beta_{(k)})\, B(t)}
#' }
#' where \eqn{\beta_{(k)}} is the basis-coefficient sub-vector for
#' coefficient \eqn{k} and the covariance comes from
#' \code{vcov(object, which = "beta")}.
#'
#' @param object A \code{vcmm_fit} object.
#' @param t_new Numeric vector at which to evaluate. Same scale as the
#'   original \code{t}; the package's normalisation is applied
#'   internally.
#' @param k Integer vector of which varying coefficients to evaluate
#'   (1-based). Default \code{NULL} = all \code{K}.
#' @param se.fit Logical. If \code{TRUE}, also returns pointwise
#'   standard errors.
#' @param ... Unused.
#'
#' @return Either a numeric matrix (\code{length(t_new)} by
#'   \code{length(k)}, default), or a list with components \code{fit}
#'   and \code{se.fit} when \code{se.fit = TRUE}.
#'
#' @seealso \code{\link{fixef.vcmm_fit}}, \code{\link{coef.vcmm_fit}}.
#'
#' @export
varying_coef <- function(object, t_new, k = NULL, se.fit = FALSE, ...) {
  if (!inherits(object, "vcmm_fit")) {
    stop("`object` must be a 'vcmm_fit'.", call. = FALSE)
  }
  ds <- object$design
  if (is.null(ds)) {
    stop("`object` has no `design` metadata; was it produced by vcmm()?",
         call. = FALSE)
  }
  if (!is.numeric(t_new)) {
    stop("`t_new` must be numeric.", call. = FALSE)
  }
  if (is.null(k)) {
    k <- seq_len(as.integer(ds$K))
  } else {
    k <- as.integer(k)
    if (any(k < 1L) || any(k > ds$K)) {
      stop(sprintf("`k` must lie in [1, %d].", ds$K), call. = FALSE)
    }
  }

  # Normalize t to the same scale used at fit time
  t_use <- if (isTRUE(ds$normalize_t)) {
    if (ds$t_max <= ds$t_min) {
      stop("Stored t-range has zero width.", call. = FALSE)
    }
    (t_new - ds$t_min) / (ds$t_max - ds$t_min)
  } else {
    t_new
  }

  # Re-evaluate the same B-spline basis using the stored knots.
  # Silence the harmless "x values beyond boundary knots" warning
  # that fires when the user passes a t_new spanning slightly past the
  # stored training boundary (very common with seq(0, 1, ...) when
  # min(t_train) is 0.001-ish). The basis evaluates correctly via
  # polynomial extrapolation; the warning is purely informational.
  B <- withCallingHandlers(
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
  B <- unclass(B); attributes(B) <- list(dim = dim(B))

  m    <- as.integer(ds$n_basis)
  beta <- as.numeric(object$beta)

  # Point estimates
  out <- matrix(0, nrow = length(t_new), ncol = length(k))
  for (i in seq_along(k)) {
    ki         <- k[i]
    coef_block <- beta[(1L + (ki - 1L) * m + 1L):(1L + ki * m)]
    out[, i]   <- as.vector(B %*% coef_block)
  }
  colnames(out) <- sprintf("beta_%d(t)", k)

  if (!isTRUE(se.fit)) return(out)

  # Standard errors: SE(B(t)' beta_k) = sqrt(B(t)' V_k B(t)) per row
  V_beta <- vcov.vcmm_fit(object, which = "beta")
  se_mat <- matrix(0, nrow = length(t_new), ncol = length(k))
  for (i in seq_along(k)) {
    ki  <- k[i]
    idx <- (1L + (ki - 1L) * m + 1L):(1L + ki * m)
    V_k <- V_beta[idx, idx, drop = FALSE]
    se_mat[, i] <- sqrt(pmax(rowSums((B %*% V_k) * B), 0))
  }
  colnames(se_mat) <- sprintf("se_beta_%d(t)", k)

  list(fit = out, se.fit = se_mat)
}

# ----------------------------------------------------------------------------
# summary  --  refreshed to carry re_cov_state through
# ----------------------------------------------------------------------------

#' Summarise a vcmm fit
#'
#' Produces a \code{vcmm_summary} object with a coefficient table for
#' the fixed-effects, the variance components, and basic fit
#' diagnostics. \code{print()} renders it in the style of \code{lm} /
#' \code{lmer}.
#'
#' The standard errors are the square roots of the diagonal of
#' \code{vcov(object, which = "beta")}, treating \eqn{\sigma_\varepsilon}
#' and \eqn{\sigma_\alpha} (or \eqn{\Sigma_\alpha}) as fixed at their
#' fitted values. They are first-order valid (Theorem 3.1) but do not
#' adjust for variance-component uncertainty.
#'
#' For \code{re_cov \%in\% c("kronecker", "separable")}, the variance-
#' components block of the printed summary displays the estimated
#' \eqn{\Sigma_{\mathrm{left}}} (\eqn{\Sigma_{2\times 2}} for OD or
#' \eqn{\Sigma_q} for group-shared dense) instead of a scalar
#' \eqn{\sigma_\alpha}.
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

  z      <- ifelse(se > 0, beta / se, NA_real_)
  pval   <- 2 * stats::pnorm(-abs(z))

  tab <- data.frame(
    Estimate     = beta,
    `Std. Error` = se,
    `z value`    = z,
    `Pr(>|z|)`   = pval,
    check.names  = FALSE
  )

  # Friendly row names
  rn <- .build_beta_names(object$design, length(beta))
  rownames(tab) <- rn

  out <- list(
    call         = object$call,
    method       = object$method,
    coefficients = tab,
    sigma_eps    = object$sigma_eps,
    sigma_alpha  = object$sigma_alpha,
    re_cov       = if (!is.null(object$re_cov)) object$re_cov else "diag",
    re_cov_state = object$re_cov_state,
    n_obs        = object$n_obs,
    p            = object$p,
    q            = object$q,
    iterations   = object$iterations,
    converged    = object$converged
  )
  class(out) <- c("vcmm_summary", "list")
  out
}

#' @rdname summary.vcmm_fit
#' @param x A \code{vcmm_summary} object.
#' @param digits Integer. Number of significant digits to display.
#' @param signif.stars Logical. If \code{TRUE}, append significance
#'   stars to the coefficient table.
#' @export
print.vcmm_summary <- function(x,
                               digits       = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  cat("VCMM fit\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  cat(sprintf(
    "\nMethod:  %s   |  RE covariance: %s   |  N = %s   |  p = %s   |  q = %s\n",
    .nv(x$method), .nv(x$re_cov), .nv(x$n_obs), .nv(x$p), .nv(x$q)))
  cat(sprintf("Converged: %s   |  iterations: %s\n",
              .nv(x$converged), .nv(x$iterations)))

  cat("\nVariance components:\n")
  cat(sprintf("  sigma_eps    = %.4f\n", x$sigma_eps))

  state <- x$re_cov_state
  is_kron <- !is.null(state) &&
             (identical(state$type, "kronecker") ||
              identical(state$type, "separable"))

  if (is_kron) {
    k  <- as.integer(state$q_left)
    G  <- as.integer(state$n_groups)
    SL <- state$Sigma_left

    if (identical(x$re_cov, "kronecker") && identical(k, 2L)) {
      rho <- if (SL[1, 1] > 0 && SL[2, 2] > 0) {
        SL[1, 2] / sqrt(SL[1, 1] * SL[2, 2])
      } else NA_real_
      cat("  Sigma_2x2    :\n")
      cat(sprintf("    [%.4f  %.4f]\n", SL[1, 1], SL[1, 2]))
      cat(sprintf("    [%.4f  %.4f]\n", SL[2, 1], SL[2, 2]))
      cat(sprintf("  OD corr      : %.4f\n", rho))
      cat(sprintf("  Sigma_spatial: %d x %d (G = %d groups, held fixed)\n",
                  G, G, G))
    } else {
      d  <- diag(SL)
      ds <- sqrt(pmax(d, 0))
      if (all(ds > 0)) {
        C   <- SL / outer(ds, ds)
        mac <- mean(abs(C[lower.tri(C)]))
      } else mac <- NA_real_
      left_name  <- if (identical(x$re_cov, "separable")) "Sigma_q" else "Sigma_left"
      right_name <- if (identical(x$re_cov, "separable")) "Omega_G" else "Sigma_right"
      cat(sprintf("  %-12s : %d x %d  (diag range %.4f .. %.4f, |corr| ~ %.3f)\n",
                  left_name, k, k, min(d), max(d), mac))
      cat(sprintf("  %-12s : %d x %d (G = %d groups, held fixed)\n",
                  right_name, G, G, G))
    }
  } else if (!is.null(x$sigma_alpha)) {
    cat(sprintf("  sigma_alpha  = %.4f\n", x$sigma_alpha))
  }

  cat("\nFixed-effects coefficients:\n")
  stats::printCoefmat(x$coefficients,
                      digits       = digits,
                      signif.stars = signif.stars,
                      P.values     = TRUE,
                      has.Pvalue   = TRUE,
                      na.print     = "NA")
  invisible(x)
}

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

# Build pretty names for the beta vector: "(Intercept)" then per-covariate
# basis labels "X1.basis1", ..., "XK.basisM". Falls back to "beta_1, ..."
# when the design metadata is missing or inconsistent.
.build_beta_names <- function(design, p) {
  if (is.null(design) ||
      is.null(design$K) || is.null(design$n_basis)) {
    return(paste0("beta_", seq_len(p)))
  }
  K <- as.integer(design$K)
  m <- as.integer(design$n_basis)
  if (1L + K * m != p) {
    return(paste0("beta_", seq_len(p)))
  }
  nm <- character(p)
  nm[1L] <- "(Intercept)"
  for (k in seq_len(K)) {
    idx <- (1L + (k - 1L) * m + 1L):(1L + k * m)
    nm[idx] <- sprintf("X%d.basis%d", k, seq_len(m))
  }
  nm
}

# "Value or NA-string" helper for print method robustness
.nv <- function(x) if (is.null(x)) "NA" else format(x)
