#===============================================================================
# fit_ss: iterative SS estimator (Algorithm 1 of Jalili and Lin, 2025)
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
#' aggregated sufficient summary, implementing Algorithm 1 of Jalili and
#' Lin (2025) for the normal linear VCMM. The raw response and design
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
#' @param re_cov_state Optional. An internal random-effects covariance
#'   state object (NULL for diagonal, or constructed for kronecker via
#'   \code{vcmm()} with \code{re_cov = "kronecker"}). When NULL, the
#'   prior precision is \code{(sigma_eps^2 / sigma_alpha^2) * I_q}
#'   matching the diagonal case. Advanced users typically reach
#'   \code{re_cov_state} via \code{vcmm()} rather than calling
#'   \code{fit_ss()} directly.
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
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
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
fit_ss <- function(stats, penalty, control = vcmm_control(),
                   re_cov_state = NULL) {

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

    ## Alpha update: (ZtZ + Sigma_alpha^{-1} prior precision) alpha = Zty - t(XtZ) beta
    ## The precision matrix is built from re_cov_state (diag or kronecker).
    prior_precision <- .build_prior_precision(re_cov_state, sigma_eps,
                                              sigma_alpha, q)
    rhs_alpha <- as.vector(stats$Zty) - as.vector(crossprod(stats$XtZ, beta))
    lhs_alpha <- stats$ZtZ + prior_precision
    alpha     <- as.vector(invert_matrix(lhs_alpha, q = q) %*% rhs_alpha)

    ## Optional variance / covariance updates (M_eta rule)
    if (control$update_variance) {
      # sigma_eps^2 = RSS / n  -- same formula for all re_cov types
      rss <- stats$a -
        2 * as.numeric(crossprod(beta,  stats$b)) +
            as.numeric(crossprod(beta,  stats$C   %*% beta)) -
        2 * as.numeric(crossprod(alpha, stats$Zty)) +
        2 * as.numeric(crossprod(beta,  stats$XtZ %*% alpha)) +
            as.numeric(crossprod(alpha, stats$ZtZ %*% alpha))
      sigma_eps <- sqrt(max(rss / n_obs, 1e-8))

      # Re-cov update:
      #   diag      -> sigma_alpha via method of moments
      #   kronecker -> update re_cov_state via estimate_kronecker_components()
      if (is.null(re_cov_state) || identical(re_cov_state$type, "diag")) {
        sigma_alpha <- sqrt(max(mean(alpha^2), 1e-8))
      } else {
        re_cov_state <- .update_re_cov_state(re_cov_state, alpha)
      }
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

  ## Assemble the joint Hessian once at convergence so summary() / vcov()
  ## don't have to recompute. K is in the same parametrisation as fit_csl.
  prior_precision_final <- .build_prior_precision(re_cov_state, sigma_eps,
                                                  sigma_alpha, q)
  V_bb_final <- stats$C   + penalty
  V_ba_final <- stats$XtZ
  V_aa_final <- stats$ZtZ + prior_precision_final
  K_final <- rbind(cbind(V_bb_final,   V_ba_final),
                   cbind(t(V_ba_final), V_aa_final))
  K_final <- (K_final + t(K_final)) / 2   # numerical symmetry
  K_inv   <- invert_matrix(K_final, q = p + q)

  ## EM-style correction for the Kronecker / separable prior covariance.
  ## Removes the BLUP shrinkage bias in Sigma_left by adding the
  ## posterior-variance partial-trace term, using K_inv we just built.
  ## No-op for re_cov = "diag" or when update_variance = FALSE.
  re_cov_state <- .apply_em_correction_if_kronecker(
    re_cov_state, alpha, sigma_eps, K_inv, p, q, control
  )

  ## If the EM correction mutated re_cov_state (only happens for
  ## kronecker/separable with update_variance = TRUE), the prior precision
  ## changed, so the K_inv we cached above is no longer consistent with
  ## the stored re_cov_state$Sigma_left. Rebuild K_inv from the post-EM
  ## precision so downstream vcov() / summary() report Wald SEs that
  ## match the final Sigma_left. Cheap: one extra (p+q) x (p+q)
  ## inversion in the kron/separable branch; no-op for diag.
  if (!is.null(re_cov_state) && !identical(re_cov_state$type, "diag") &&
      isTRUE(control$update_variance)) {
    prior_precision_final <- .build_prior_precision(re_cov_state, sigma_eps,
                                                    sigma_alpha, q)
    V_aa_final <- stats$ZtZ + prior_precision_final
    K_final    <- rbind(cbind(V_bb_final,   V_ba_final),
                        cbind(t(V_ba_final), V_aa_final))
    K_final    <- (K_final + t(K_final)) / 2
    K_inv      <- invert_matrix(K_final, q = p + q)
  }

  ## Marginal log-likelihood at convergence (used by logLik.vcmm_fit, AIC, BIC).
  marginal_loglik <- .compute_marginal_loglik(
    stats = stats, beta = beta, alpha = alpha,
    sigma_eps = sigma_eps, sigma_alpha = sigma_alpha,
    re_cov_state = re_cov_state,
    V_aa = V_aa_final, n_obs = n_obs, q = q
  )

  if (control$verbose) {
    cat(sprintf("  [SS] converged=%s | iter=%d | elapsed=%.4fs\n",
                converged, iter, elapsed))
  }

  out <- list(
    beta            = beta,
    alpha           = alpha,
    sigma_eps       = sigma_eps,
    sigma_alpha     = sigma_alpha,
    iterations      = iter,
    converged       = converged,
    elapsed_sec     = elapsed,
    n_obs           = n_obs,
    p               = p,
    q               = q,
    method          = "SS",
    control         = control,
    K_inv           = K_inv,
    re_cov_state    = re_cov_state,
    marginal_loglik = marginal_loglik,
    call            = match.call()
  )
  class(out) <- c("vcmm_fit", "list")
  out
}

## Format an elapsed time in seconds for printing. Sub-millisecond
## proc.time() values underflow to "0.0000" with the naive %.4f format
## (typical on Apple Silicon for small N); show them as "<0.001" so the
## user sees that the fit was very fast rather than thinking the timer
## is broken. Returns "NA" for NULL / NA inputs.
.fmt_elapsed <- function(s) {
  if (is.null(s) || is.na(s)) return("NA")
  if (s < 0.0005) return("<0.001")
  sprintf("%.4f", s)
}

#' @rdname fit_ss
#' @param x A \code{vcmm_fit} object.
#' @param ... Unused; present for S3 method consistency.
#' @export
print.vcmm_fit <- function(x, ...) {
  cat("<vcmm_fit>  Varying Coefficient Mixed-Effects Model fit\n")
  cat(sprintf("  method      : %s\n", x$method))
  cat(sprintf("  n_obs       : %d\n", x$n_obs))
  cat(sprintf("  p (fixed)   : %d\n", x$p))
  cat(sprintf("  q (random)  : %d\n", x$q))
  cat(sprintf("  RE cov      : %s\n",
              if (!is.null(x$re_cov)) x$re_cov else "diag"))

  if (identical(x$method, "CSL") && !is.null(x$pilot)) {
    cat(sprintf("  pilot iter  : %d (%s)\n",
                x$pilot$iterations,
                if (x$pilot$converged) "converged" else "loose"))
    cat(sprintf("  newton step : 1\n"))
  } else {
    cat(sprintf("  iterations  : %d %s\n", x$iterations,
                if (x$converged) "(converged)" else "(NOT converged)"))
  }

  cat(sprintf("  sigma_eps   : %.4f\n", x$sigma_eps))

  is_kron_state <- !is.null(x$re_cov_state) &&
    (identical(x$re_cov_state$type, "kronecker") ||
     identical(x$re_cov_state$type, "separable"))

  if (is_kron_state) {
    k <- x$re_cov_state$q_left
    G <- x$re_cov_state$n_groups
    SL <- x$re_cov_state$Sigma_left

    if (identical(x$re_cov, "kronecker") && identical(as.integer(k), 2L)) {
      ## OD-style display: full 2x2 plus the correlation.
      rho <- if (SL[1, 1] > 0 && SL[2, 2] > 0) {
        SL[1, 2] / sqrt(SL[1, 1] * SL[2, 2])
      } else NA_real_
      cat("  Sigma_2x2   :\n")
      cat(sprintf("    [%.4f  %.4f]\n", SL[1, 1], SL[1, 2]))
      cat(sprintf("    [%.4f  %.4f]\n", SL[2, 1], SL[2, 2]))
      cat(sprintf("  OD corr     : %.4f\n", rho))
      cat(sprintf("  Sigma_spatial: %d x %d (G = %d groups)\n", G, G, G))
    } else {
      ## General display for q_left > 2 or re_cov = "separable":
      ## show dimensions, diagonal range, mean abs off-diag correlation.
      sl_diag <- diag(SL)
      d_sd <- sqrt(pmax(sl_diag, 0))
      if (all(d_sd > 0)) {
        Corr <- SL / outer(d_sd, d_sd)
        off  <- Corr[lower.tri(Corr)]
        mean_abs_corr <- mean(abs(off))
      } else {
        mean_abs_corr <- NA_real_
      }
      left_name  <- if (identical(x$re_cov, "separable")) "Sigma_q" else "Sigma_left"
      right_name <- if (identical(x$re_cov, "separable")) "Omega_G" else "Sigma_right"
      cat(sprintf("  %-12s: %d x %d  (diag range %.4f .. %.4f, |corr| ~ %.3f)\n",
                  left_name, k, k, min(sl_diag), max(sl_diag), mean_abs_corr))
      cat(sprintf("  %-12s: %d x %d (G = %d groups, held fixed)\n",
                  right_name, G, G, G))
    }
  } else {
    cat(sprintf("  sigma_alpha : %.4f\n", x$sigma_alpha))
  }

  if (identical(x$method, "CSL") &&
      !is.null(x$pilot_elapsed_sec) && !is.null(x$step_elapsed_sec)) {
    cat(sprintf("  elapsed     : %s sec (pilot %ss + newton %ss)\n",
                .fmt_elapsed(x$elapsed_sec),
                .fmt_elapsed(x$pilot_elapsed_sec),
                .fmt_elapsed(x$step_elapsed_sec)))
  } else {
    cat(sprintf("  elapsed     : %s sec\n", .fmt_elapsed(x$elapsed_sec)))
  }

  invisible(x)
}
