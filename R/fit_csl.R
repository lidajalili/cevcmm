#===============================================================================
# fit_csl: one-step communication-efficient surrogate-likelihood estimator
#
# Implements the one-step refinement of Jalili and Lin (2025, Sec. 3):
#   theta_csl = theta_pilot - K^{-1} g(theta_pilot)
# where g and K use the FULL aggregated sufficient statistics. The pilot
# estimator is either user-supplied or produced internally by running
# fit_ss() to a loose tolerance.
#===============================================================================

#' Fit a VCMM via the one-step CSL estimator
#'
#' Performs a single Newton refinement of a pilot estimator using the
#' aggregated sufficient statistics, implementing the one-step
#' communication-efficient surrogate-likelihood (CSL) estimator of
#' Jalili and Lin (2025, Section 3). For the normal linear VCMM the
#' update is
#' \deqn{
#'   \widehat{\theta}_{CSL}
#'   = \widehat{\theta}_0
#'   - K^{-1} \, g(\widehat{\theta}_0),
#' }
#' where \eqn{\theta = (\beta, \alpha)}, the gradient \eqn{g} uses the
#' full aggregated stats, and the Hessian \eqn{K} also uses the full
#' aggregated stats with prior augmentation. Theorem 3.1 of the paper
#' shows that whenever the pilot estimator is \eqn{\sqrt{N}}-consistent,
#' the one-step CSL estimator is first-order equivalent to the full SS
#' estimator.
#'
#' \strong{Pilot estimator.} If \code{pilot = NULL} (default), an
#' internal SS pilot is run with loose tolerances and a small number of
#' iterations (\code{pilot_max_iter = 5}, \code{pilot_tol_beta = 1e-3},
#' \code{pilot_tol_alpha = 1e-3}); these defaults match the dense
#' simulation study of the paper. Setting \code{pilot_max_iter = 1L}
#' gives the cheapest possible pilot (the OLS-like single-step pilot
#' used in the origin-destination simulation), still
#' \eqn{\sqrt{N}}-consistent in the normal linear case. Advanced users
#' may pass any \code{vcmm_fit} object via the \code{pilot} argument --
#' e.g. a previously fitted \code{fit_ss()} result.
#'
#' \strong{Hessian.} This implementation uses the full aggregated
#' Hessian, not the reference-node curvature approximation \eqn{\tilde K}
#' from the paper. The two are first-order equivalent under the
#' conditions of Theorem 3.1; the full-aggregated form is the most
#' accurate and is the natural default for a single-server fit, which
#' is the typical use case of this package.
#'
#' @param stats A \code{vcmm_ss} or \code{vcmm_accumulator} object
#'   containing the aggregated sufficient statistics.
#' @param penalty A symmetric \eqn{p \times p} penalty matrix from
#'   \code{build_penalty_matrix()}.
#' @param control A \code{vcmm_control} object with fitting options.
#'   The \code{sigma_eps} and \code{sigma_alpha} entries provide the
#'   initial variance values used by the internal pilot (when
#'   \code{pilot = NULL}).
#' @param pilot Optional \code{vcmm_fit} object to use as the pilot
#'   estimator. If \code{NULL} (default), an internal SS pilot is run.
#' @param pilot_max_iter Integer. Maximum iterations for the internal
#'   SS pilot (default 5). Ignored if \code{pilot} is supplied.
#' @param pilot_tol_beta Positive numeric. Loose tolerance for the
#'   internal SS pilot (default 1e-3). Ignored if \code{pilot} is
#'   supplied.
#' @param pilot_tol_alpha Positive numeric. Loose tolerance for the
#'   internal SS pilot (default 1e-3). Ignored if \code{pilot} is
#'   supplied.
#' @param re_cov_state Optional. An internal random-effects covariance
#'   state object (NULL for diagonal, or constructed for kronecker via
#'   \code{vcmm()} with \code{re_cov = "kronecker"}). Passed through to
#'   the internal SS pilot and reused for the Newton-step Hessian.
#'   Advanced users typically reach \code{re_cov_state} via
#'   \code{vcmm()} rather than calling \code{fit_csl()} directly.
#'
#' @return A list of class \code{"vcmm_fit"} with the same fields as
#'   \code{fit_ss()}, plus:
#' \itemize{
#'   \item \code{pilot}: the \code{vcmm_fit} pilot used.
#'   \item \code{pilot_elapsed_sec}: pilot fitting time.
#'   \item \code{step_elapsed_sec}: Newton step time.
#' }
#'   The \code{method} field is \code{"CSL"}.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500; p <- 4; q <- 3
#' X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
#' Z <- matrix(rnorm(n * q), n, q)
#' alpha_true <- rnorm(q, sd = 0.5)
#' y <- as.vector(
#'   X %*% c(2, 0.5, -0.3, 0.8) + Z %*% alpha_true + rnorm(n, sd = 0.5)
#' )
#'
#' stats <- compute_sufficient_stats(y, X, Z)
#' P     <- build_penalty_matrix(n_basis = p - 1, lambda = 0.1)
#'
#' # Default: internal SS pilot (5 loose iterations) then one Newton step
#' fit <- fit_csl(stats, P,
#'                vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
#' fit
#'
#' # Cheapest pilot: single SS step
#' fit_one <- fit_csl(stats, P,
#'                    vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5),
#'                    pilot_max_iter = 1L)
#'
#' # Advanced: user-supplied pilot
#' my_pilot <- fit_ss(stats, P,
#'                    vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5,
#'                                 max_iter = 3))
#' fit_user <- fit_csl(stats, P,
#'                     vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5),
#'                     pilot = my_pilot)
fit_csl <- function(stats,
                    penalty,
                    control         = vcmm_control(),
                    pilot           = NULL,
                    pilot_max_iter  = 5L,
                    pilot_tol_beta  = 1e-3,
                    pilot_tol_alpha = 1e-3,
                    re_cov_state    = NULL) {

  # --- Input validation ---------------------------------------------------
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

  # --- Phase 1: pilot estimator ------------------------------------------
  start_total <- proc.time()

  if (is.null(pilot)) {
    # Build a loose control for the internal SS pilot
    pilot_ctrl <- vcmm_control(
      max_iter        = pilot_max_iter,
      tol_beta        = pilot_tol_beta,
      tol_alpha       = pilot_tol_alpha,
      sigma_eps       = control$sigma_eps,
      sigma_alpha     = control$sigma_alpha,
      update_variance = control$update_variance,
      verbose         = FALSE
    )
    if (control$verbose) {
      cat(sprintf("  [CSL] Phase 1: internal SS pilot (max_iter=%d, tol=%.0e)\n",
                  pilot_max_iter, pilot_tol_beta))
    }
    pilot <- fit_ss(stats, penalty, pilot_ctrl, re_cov_state = re_cov_state)

    # If the pilot updated the re_cov_state (via update_variance), inherit it
    if (!is.null(pilot$re_cov_state)) {
      re_cov_state <- pilot$re_cov_state
    }
  } else {
    if (!inherits(pilot, "vcmm_fit")) {
      stop("`pilot` must be a 'vcmm_fit' object (e.g. from fit_ss()).",
           call. = FALSE)
    }
    if (length(pilot$beta) != p) {
      stop(sprintf(
        "pilot$beta has length %d but stats$C is %d x %d.",
        length(pilot$beta), p, p
      ), call. = FALSE)
    }
    if (length(pilot$alpha) != q) {
      stop(sprintf(
        "pilot$alpha has length %d but stats$ZtZ is %d x %d.",
        length(pilot$alpha), q, q
      ), call. = FALSE)
    }
    if (control$verbose) {
      cat("  [CSL] Phase 1: using user-supplied pilot\n")
    }
  }

  pilot_elapsed <- pilot$elapsed_sec

  # --- Phase 2: one Newton step using full aggregated Hessian ------------
  if (control$verbose) {
    cat("  [CSL] Phase 2: one Newton step (full aggregated Hessian)\n")
  }
  start_step <- proc.time()

  beta_0      <- pilot$beta
  alpha_0     <- pilot$alpha
  sigma_eps   <- pilot$sigma_eps
  sigma_alpha <- pilot$sigma_alpha
  se2 <- sigma_eps^2
  sa2 <- sigma_alpha^2

  ## Prior precision for alpha (diag or kronecker, depending on re_cov_state)
  prior_precision <- .build_prior_precision(re_cov_state, sigma_eps,
                                            sigma_alpha, q)

  ## Gradient at (beta_0, alpha_0) using FULL stats.
  ## Form matches the implicit loss minimised by fit_ss():
  ##   L = -(1/2 se2) ( ||y - X beta - Z alpha||^2 + beta' P beta )
  ##       - (1/2) alpha' Sigma_alpha^{-1} alpha
  ## After multiplying through by se2 (the Newton direction is invariant),
  ## the working gradient and Hessian are:
  g_beta  <- as.vector((stats$C + penalty) %*% beta_0) +
             as.vector(stats$XtZ %*% alpha_0) -
             as.vector(stats$b)
  g_alpha <- as.vector(crossprod(stats$XtZ, beta_0)) +
             as.vector(stats$ZtZ %*% alpha_0) -
             as.vector(stats$Zty) +
             as.vector(prior_precision %*% alpha_0)
  g_vec <- c(g_beta, g_alpha)

  ## Hessian (full aggregated, prior-augmented), in the same parametrisation
  V_bb <- stats$C   + penalty
  V_ba <- stats$XtZ
  V_aa <- stats$ZtZ + prior_precision
  K <- rbind(cbind(V_bb,   V_ba),
             cbind(t(V_ba), V_aa))
  K <- (K + t(K)) / 2   # enforce symmetry numerically

  ## Newton step: theta_csl = theta_0 - K^{-1} g
  K_inv <- invert_matrix(K, q = p + q)
  step  <- as.vector(K_inv %*% g_vec)

  beta_csl  <- beta_0  - step[seq_len(p)]
  alpha_csl <- alpha_0 - step[seq.int(p + 1L, p + q)]

  ## Optional post-step variance / covariance update (same M_eta rule as fit_ss)
  if (control$update_variance) {
    rss <- stats$a -
      2 * as.numeric(crossprod(beta_csl,  stats$b)) +
          as.numeric(crossprod(beta_csl,  stats$C   %*% beta_csl)) -
      2 * as.numeric(crossprod(alpha_csl, stats$Zty)) +
      2 * as.numeric(crossprod(beta_csl,  stats$XtZ %*% alpha_csl)) +
          as.numeric(crossprod(alpha_csl, stats$ZtZ %*% alpha_csl))
    sigma_eps <- sqrt(max(rss / n_obs, 1e-8))

    if (is.null(re_cov_state) || identical(re_cov_state$type, "diag")) {
      sigma_alpha <- sqrt(max(mean(alpha_csl^2), 1e-8))
    } else {
      re_cov_state <- .update_re_cov_state(re_cov_state, alpha_csl)
    }
  }

  ## EM-style correction for the Kronecker / separable prior covariance.
  ## Overrides the moment-based update above with the unbiased
  ## posterior-variance-corrected estimator, using the K_inv we just built.
  ## No-op for re_cov = "diag" or when update_variance = FALSE.
  re_cov_state <- .apply_em_correction_if_kronecker(
    re_cov_state, alpha_csl, sigma_eps, K_inv, p, q, control
  )

  ## Marginal log-likelihood at the final (β_csl, α_csl, σ_ε, Σ_α). Rebuild
  ## V_aa from the FINAL prior precision so log|Σ_y| reflects the post-EM
  ## Sigma_left; this is the value summary() / logLik() / AIC() will report.
  prior_precision_final <- .build_prior_precision(re_cov_state, sigma_eps,
                                                  sigma_alpha, q)
  V_aa_final <- stats$ZtZ + prior_precision_final

  ## Rebuild K_inv from the post-EM precision so vcov() / summary() use
  ## Wald SEs consistent with the stored re_cov_state$Sigma_left. The
  ## one-step Newton update above used the pre-EM K_inv by design (the
  ## EM correction is computed from it), so this rebuild is purely for
  ## downstream inference, not the point estimate. No-op cost when
  ## re_cov = "diag" or update_variance = FALSE; one extra (p+q)x(p+q)
  ## inversion in the kron/separable + update_variance branch.
  if (!is.null(re_cov_state) && !identical(re_cov_state$type, "diag") &&
      isTRUE(control$update_variance)) {
    K_post <- rbind(cbind(V_bb,   V_ba),
                    cbind(t(V_ba), V_aa_final))
    K_post <- (K_post + t(K_post)) / 2
    K_inv  <- invert_matrix(K_post, q = p + q)
  }

  marginal_loglik <- .compute_marginal_loglik(
    stats = stats, beta = beta_csl, alpha = alpha_csl,
    sigma_eps = sigma_eps, sigma_alpha = sigma_alpha,
    re_cov_state = re_cov_state,
    V_aa = V_aa_final, n_obs = n_obs, q = q
  )

  step_elapsed  <- as.numeric((proc.time() - start_step)["elapsed"])
  total_elapsed <- as.numeric((proc.time() - start_total)["elapsed"])

  if (control$verbose) {
    cat(sprintf("  [CSL] pilot=%.4fs | newton=%.4fs | total=%.4fs\n",
                pilot_elapsed, step_elapsed, total_elapsed))
  }

  out <- list(
    beta              = beta_csl,
    alpha             = alpha_csl,
    sigma_eps         = sigma_eps,
    sigma_alpha       = sigma_alpha,
    iterations        = 1L,
    converged         = TRUE,
    elapsed_sec       = total_elapsed,
    pilot_elapsed_sec = pilot_elapsed,
    step_elapsed_sec  = step_elapsed,
    n_obs             = n_obs,
    p                 = p,
    q                 = q,
    method            = "CSL",
    pilot             = pilot,
    control           = control,
    K_inv             = K_inv,
    re_cov_state      = re_cov_state,
    marginal_loglik   = marginal_loglik,
    call              = match.call()
  )
  class(out) <- c("vcmm_fit", "list")
  out
}
