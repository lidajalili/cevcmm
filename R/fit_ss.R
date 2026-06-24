#===============================================================================
# fit_ss: iterative SS estimator (Algorithm 1 of Lin and Jalili, 2026)
#
# Solves alternately:
#   beta-update:  (C + P) beta = b - XtZ alpha
#   alpha-update: (ZtZ + (sigma_eps^2 / sigma_alpha^2) I) alpha
#                       = Zty - t(XtZ) beta
# All matrices are read off the aggregated sufficient-statistics object;
# the raw response and design matrices are not needed.
#===============================================================================

#' Fit a VCMM via the sufficient-statistics estimator
#'
#' Iteratively solves for the fixed-effects coefficient vector
#' \code{beta} and the random-effects vector \code{alpha} using only the
#' aggregated sufficient summary, implementing Algorithm 1 of Lin and
#' Jalili (2026) for the normal linear VCMM. The raw response and design
#' matrices are not needed: everything reads off \code{stats}, which may
#' come from a single call to \code{compute_sufficient_stats()} or from a
#' streaming accumulator (\code{init_accumulator()} plus repeated
#' \code{accumulate_stats()}).
#'
#' At each iteration the algorithm solves:
#' \itemize{
#'   \item beta-update: \code{(C + P) beta = b - XtZ \%*\% alpha}
#'   \item alpha-update:
#'     \code{(ZtZ + (sigma_eps^2 / sigma_alpha^2) * I) alpha = Zty - t(XtZ) \%*\% beta}
#' }
#' Both linear systems are solved via \code{invert_matrix()}, which
#' automatically dispatches between \code{solve()} and SVD pseudo-inverse
#' depending on dimension and condition number. Convergence is declared
#' when the relative change in both \code{beta} and \code{alpha} falls
#' below their respective tolerances.
#'
#' If \code{control$update_variance} is \code{TRUE}, the residual variance
#' is re-estimated each iteration from the SS residual sum of squares,
#' and the random-effect variance is re-estimated as
#' \code{mean(alpha^2)}.
#'
#' @param stats A \code{vcmm_ss} or \code{vcmm_accumulator} object
#'   containing the aggregated sufficient statistics.
#' @param penalty A symmetric \eqn{p \times p} penalty matrix from
#'   \code{build_penalty_matrix()}.
#' @param control A \code{vcmm_control} object with fitting options. Pass
#'   \code{vcmm_control()} to use defaults.
#'
#' @return A list of class \code{"vcmm_fit"} with elements:
#' \itemize{
#'   \item \code{beta}: fitted fixed-effects vector, length p.
#'   \item \code{alpha}: fitted random-effects vector, length q.
#'   \item \code{sigma_eps}: final residual standard deviation.
#'   \item \code{sigma_alpha}: final random-effect standard deviation.
#'   \item \code{iterations}: number of iterations performed.
#'   \item \code{converged}: \code{TRUE} if convergence tolerances were met.
#'   \item \code{elapsed_sec}: wall-clock fitting time in seconds.
#'   \item \code{n_obs}, \code{p}, \code{q}: data and design dimensions.
#'   \item \code{method}: character, \code{"SS"}.
#'   \item \code{control}: the \code{vcmm_control} object used.
#'   \item \code{call}: the matched call.
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
#' n <- 200; p <- 4; q <- 3
#' X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
#' Z <- matrix(rnorm(n * q), n, q)
#' alpha_true <- rnorm(q, sd = 0.5)
#' y <- as.vector(
#'   X %*% c(2, 0.5, -0.3, 0.8) + Z %*% alpha_true + rnorm(n, sd = 0.5)
#' )
#'
#' # Build sufficient statistics and penalty
#' stats <- compute_sufficient_stats(y, X, Z)
#' P     <- build_penalty_matrix(n_basis = p - 1, lambda = 0.1)
#'
#' # Fit with fixed variances
#' fit <- fit_ss(stats, P,
#'               vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' fit
#' coef(fit)
fit_ss <- function(stats, penalty, control = vcmm_control()) {

  # --- Input validation ----------------------------------------------------
  if (!(inherits(stats, "vcmm_ss") || inherits(stats, "vcmm_accumulator"))) {
    stop("`stats` must be a 'vcmm_ss' or 'vcmm_accumulator' object.",
         call. = FALSE)
  }
  if (!inherits(control, "vcmm_control")) {
    stop("`control` must be a 'vcmm_control' object (see vcmm_control()).",
         call. = FALSE)
  }
  if (!is.matrix(penalty) || !is.numeric(penalty)) {
    stop("`penalty` must be a numeric matrix.", call. = FALSE)
  }

  p     <- nrow(stats$C)
  q     <- nrow(stats$ZtZ)
  n_obs <- stats$n_obs

  if (!identical(dim(penalty), c(p, p))) {
    stop(sprintf(
      "`penalty` must be %d x %d to match `stats$C`; got %d x %d.",
      p, p, nrow(penalty), ncol(penalty)
    ), call. = FALSE)
  }
  if (n_obs <= 0L) {
    stop("`stats$n_obs` is 0; nothing to fit.", call. = FALSE)
  }

  # --- Initialize ----------------------------------------------------------
  beta        <- rep(0, p)
  alpha       <- rep(0, q)
  sigma_eps   <- control$sigma_eps
  sigma_alpha <- control$sigma_alpha

  converged  <- FALSE
  iter       <- 0L
  start_time <- proc.time()

  # --- Iterate --------------------------------------------------------------
  while (iter < control$max_iter) {
    iter      <- iter + 1L
    beta_old  <- beta
    alpha_old <- alpha

    ## Beta update: (C + P) beta = b - XtZ alpha
    rhs_beta <- as.vector(stats$b) - as.vector(stats$XtZ %*% alpha)
    lhs_beta <- stats$C + penalty
    beta     <- as.vector(invert_matrix(lhs_beta, q = p) %*% rhs_beta)

    ## Alpha update: (ZtZ + (s2_eps / s2_alpha) I) alpha = Zty - t(XtZ) beta
    ridge     <- (sigma_eps^2) / (sigma_alpha^2)
    rhs_alpha <- as.vector(stats$Zty) - as.vector(crossprod(stats$XtZ, beta))
    lhs_alpha <- stats$ZtZ + diag(ridge, q)
    alpha     <- as.vector(invert_matrix(lhs_alpha, q = q) %*% rhs_alpha)

    ## Optional variance updates
    if (control$update_variance) {
      # sigma_eps^2 = RSS / n, where RSS is written in terms of SS
      rss <- stats$a -
        2 * as.numeric(crossprod(beta,  stats$b)) +
            as.numeric(crossprod(beta,  stats$C   %*% beta)) -
        2 * as.numeric(crossprod(alpha, stats$Zty)) +
        2 * as.numeric(crossprod(beta,  stats$XtZ %*% alpha)) +
            as.numeric(crossprod(alpha, stats$ZtZ %*% alpha))
      sigma_eps   <- sqrt(max(rss / n_obs, 1e-8))

      # sigma_alpha^2 = mean(alpha^2)  (method of moments)
      sigma_alpha <- sqrt(max(mean(alpha^2), 1e-8))
    }

    ## Convergence check (relative change)
    d_beta  <- sqrt(sum((beta  - beta_old )^2)) / max(1, sqrt(sum(beta_old^2)))
    d_alpha <- sqrt(sum((alpha - alpha_old)^2)) / max(1, sqrt(sum(alpha_old^2)))

    if (control$verbose && (iter %% 20L == 0L)) {
      cat(sprintf("    [SS] iter %3d  d_beta=%.2e  d_alpha=%.2e\n",
                  iter, d_beta, d_alpha))
    }

    if (d_beta < control$tol_beta && d_alpha < control$tol_alpha) {
      converged <- TRUE
      break
    }
  }

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (control$verbose) {
    cat(sprintf("  [SS] converged=%s | iter=%d | elapsed=%.4fs\n",
                converged, iter, elapsed))
  }

  out <- list(
    beta        = beta,
    alpha       = alpha,
    sigma_eps   = sigma_eps,
    sigma_alpha = sigma_alpha,
    iterations  = iter,
    converged   = converged,
    elapsed_sec = elapsed,
    n_obs       = n_obs,
    p           = p,
    q           = q,
    method      = "SS",
    control     = control,
    call        = match.call()
  )
  class(out) <- c("vcmm_fit", "list")
  out
}

#' @rdname fit_ss
#' @param x A \code{vcmm_fit} object.
#' @export
print.vcmm_fit <- function(x, ...) {
  cat("<vcmm_fit>  Varying Coefficient Mixed-Effects Model fit\n")
  cat(sprintf("  method      : %s\n", x$method))
  cat(sprintf("  n_obs       : %d\n", x$n_obs))
  cat(sprintf("  p (fixed)   : %d\n", x$p))
  cat(sprintf("  q (random)  : %d\n", x$q))
  cat(sprintf("  iterations  : %d %s\n", x$iterations,
              if (x$converged) "(converged)" else "(NOT converged)"))
  cat(sprintf("  sigma_eps   : %.4f\n", x$sigma_eps))
  cat(sprintf("  sigma_alpha : %.4f\n", x$sigma_alpha))
  cat(sprintf("  elapsed     : %.4f sec\n", x$elapsed_sec))
  invisible(x)
}

#' @rdname fit_ss
#' @param object A \code{vcmm_fit} object.
#' @export
coef.vcmm_fit <- function(object, ...) {
  object$beta
}
