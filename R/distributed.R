#===============================================================================
# Distributed sufficient-statistics API for VCMM fitting
#
# Implements the two-stage workflow in Algorithm 1 of Jalili and Lin (2025):
#
#   STAGE 1 (each node, in parallel):
#     gamma_s = node_summary(y_s, X_s, Z_s)
#
#   STAGE 2 (central server):
#     gamma = Reduce(`+`, list(gamma_1, ..., gamma_K))
#     fit   = fit_from_summaries(gamma, penalty, control, ...)
#
# Communication cost per node: one (a, b, C, ZtZ, Zty, XtZ) tuple of
# constant size in N (depends only on p and q). Raw y, X, Z never leave
# the node. The aggregated summary is mathematically identical to the
# pooled-data sufficient statistic, so the resulting fit is bit-equivalent
# (up to floating-point summation noise of ~1e-12) to a single-node
# vcmm() fit on the pooled data.
#===============================================================================

#' Compute one node's sufficient-statistics summary
#'
#' Convenience alias for \code{\link{compute_sufficient_stats}}, intended
#' for distributed-computing workflows. Each compute node calls this on
#' its local data and transmits the small returned summary; the central
#' server aggregates summaries with \code{+} (or
#' \code{Reduce("+", summaries)}) and fits via
#' \code{\link{fit_from_summaries}}.
#'
#' \strong{What's transmitted.} The returned object contains six fixed-size
#' arrays whose dimensions depend only on \eqn{p = \mathrm{ncol}(X)} and
#' \eqn{q = \mathrm{ncol}(Z)}, never on the node-local sample size. For a
#' typical VCMM with \eqn{p \approx 15} and \eqn{q \approx 50}, one
#' summary is a few kilobytes regardless of \eqn{N_s}.
#'
#' \strong{Basis alignment.} All nodes must build their local
#' \code{X} (the spline-expanded fixed-effects design) using the
#' \emph{same} basis specification. The recommended workflow is to call
#' \code{\link{build_vcmm_design}} once with the full \code{t}-range
#' (or pre-agreed knots), broadcast the resulting basis, and have each
#' node evaluate it on local \code{t}. The package itself does not enforce
#' this; mis-aligned bases will silently give incorrect fits.
#'
#' @inheritParams compute_sufficient_stats
#'
#' @return A \code{vcmm_ss} object (additive via \code{+.vcmm_ss}).
#'
#' @seealso \code{\link{fit_from_summaries}},
#'   \code{\link{compute_sufficient_stats}},
#'   \code{\link{build_vcmm_design}}
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n_per_node <- 100; p <- 5; q <- 3
#' X <- cbind(1, matrix(rnorm(n_per_node * (p - 1)), n_per_node, p - 1))
#' Z <- matrix(rnorm(n_per_node * q), n_per_node, q)
#' y <- rnorm(n_per_node)
#'
#' gamma_1 <- node_summary(y, X, Z)
#' gamma_2 <- node_summary(y, X, Z)
#' gamma_pooled <- gamma_1 + gamma_2
#' gamma_pooled
node_summary <- function(y, X, Z) {
  compute_sufficient_stats(y, X, Z)
}

#' Additive aggregation of sufficient-statistics summaries
#'
#' Defines the \code{+} method for \code{vcmm_ss} objects, so that
#' \code{Reduce("+", summaries)} aggregates a list of node summaries
#' element-wise across the six sufficient-statistics components.
#'
#' Dimensions (\eqn{p}, \eqn{q}) must match between operands; an
#' informative error is thrown otherwise.
#'
#' @param e1,e2 \code{vcmm_ss} objects.
#'
#' @return A \code{vcmm_ss} object holding the component-wise sums.
#'
#' @export
"+.vcmm_ss" <- function(e1, e2) {
  if (!inherits(e2, "vcmm_ss")) {
    stop("Both operands of `+` must be 'vcmm_ss' objects.", call. = FALSE)
  }
  if (!identical(dim(e1$C), dim(e2$C))) {
    stop(sprintf(
      "Dimension mismatch in `+.vcmm_ss`: C is %s vs %s.",
      paste(dim(e1$C), collapse = "x"),
      paste(dim(e2$C), collapse = "x")
    ), call. = FALSE)
  }
  if (!identical(dim(e1$ZtZ), dim(e2$ZtZ))) {
    stop(sprintf(
      "Dimension mismatch in `+.vcmm_ss`: ZtZ is %s vs %s.",
      paste(dim(e1$ZtZ), collapse = "x"),
      paste(dim(e2$ZtZ), collapse = "x")
    ), call. = FALSE)
  }

  out <- list(
    a     = e1$a     + e2$a,
    b     = e1$b     + e2$b,
    C     = e1$C     + e2$C,
    ZtZ   = e1$ZtZ   + e2$ZtZ,
    Zty   = e1$Zty   + e2$Zty,
    XtZ   = e1$XtZ   + e2$XtZ,
    n_obs = e1$n_obs + e2$n_obs
  )
  class(out) <- c("vcmm_ss", "list")
  out
}

#' Fit a VCMM from aggregated sufficient-statistics summaries
#'
#' Server-side fit using only the small per-node summaries; no raw data
#' required. The mathematical guarantee (Theorem 1 of Jalili and Lin,
#' 2025) is that for a fixed partition of the data into nodes, the fit
#' obtained here is identical to the one a centralised
#' \code{\link{vcmm}()} call would produce on the pooled data, up to
#' floating-point summation noise.
#'
#' \strong{Identifiability re-centering.} \code{\link{vcmm}()} applies a
#' post-fit re-centering when \code{rowSums(Z)} is constant (indicator-Z
#' OD designs and similar), absorbing \code{rowsum_constant *
#' mean(alpha_hat)} into \code{beta_0}. The server has no access to
#' \code{Z}, so the caller must signal that the original \code{Z} had
#' constant row sums by passing \code{rowsum_constant}. Default
#' \code{NULL} means "no re-centering" (appropriate for dense-Z and
#' block-Z designs).
#'
#' @param summaries Either a list of \code{vcmm_ss} objects (one per
#'   node, the typical case) OR a single \code{vcmm_ss} that has already
#'   been aggregated (e.g., via \code{Reduce("+", ...)}) OR a
#'   \code{vcmm_accumulator}.
#' @param penalty The \eqn{p \times p} smoothing-penalty matrix used in
#'   the spline basis the nodes shared. Build it via the package's
#'   \code{build_vcmm_design()} (returned as \code{design$penalty}) or
#'   directly via \code{build_penalty_matrix()}.
#' @param control A \code{\link{vcmm_control}} object.
#' @param method Either \code{"csl"} (default) or \code{"ss"}.
#' @param re_cov Either \code{"diag"}, \code{"kronecker"}, or
#'   \code{"separable"}. See \code{\link{vcmm}} for full semantics.
#' @param n_groups Required if \code{re_cov} is \code{"kronecker"} or
#'   \code{"separable"}.
#' @param q_left Required for \code{"separable"}; defaults to 2 for
#'   \code{"kronecker"}.
#' @param Sigma_left_init,Sigma_right_init,Sigma_2x2_init,Sigma_spatial_init,Sigma_q_init,Omega_G_init
#'   Same aliases accepted as in \code{\link{vcmm}()}.
#' @param rowsum_constant Optional numeric. If supplied and non-zero,
#'   apply the same identifiability re-centering that
#'   \code{\link{vcmm}()} applies post-fit. See Details.
#' @param ... Passed to \code{\link{fit_csl}}.
#'
#' @return A \code{vcmm_fit} object.
#'
#' @seealso \code{\link{node_summary}}, \code{\link{vcmm}}.
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 600
#' t  <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 4), n, 4)
#' a_true <- rnorm(4, sd = 0.3)
#' y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% a_true) + rnorm(n, sd = 0.5)
#'
#' design <- build_vcmm_design(X = x, t = t)
#' Xd     <- design$X_design
#' ctrl   <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.3)
#'
#' # Split into 3 nodes and aggregate
#' idx_node <- sample.int(3, n, replace = TRUE)
#' summaries <- lapply(seq_len(3), function(s) {
#'   ii <- which(idx_node == s)
#'   node_summary(y[ii], Xd[ii, , drop = FALSE], Z[ii, , drop = FALSE])
#' })
#' fit <- fit_from_summaries(summaries,
#'                            penalty = design$penalty,
#'                            control = ctrl)
fit_from_summaries <- function(summaries,
                                penalty,
                                control            = vcmm_control(),
                                method             = c("csl", "ss"),
                                re_cov             = c("diag", "kronecker",
                                                       "separable"),
                                n_groups           = NULL,
                                q_left             = NULL,
                                Sigma_left_init    = NULL,
                                Sigma_right_init   = NULL,
                                Sigma_2x2_init     = NULL,
                                Sigma_spatial_init = NULL,
                                Sigma_q_init       = NULL,
                                Omega_G_init       = NULL,
                                rowsum_constant    = NULL,
                                ...) {
  method <- match.arg(method)
  re_cov <- match.arg(re_cov)

  # ----- Aggregate ---------------------------------------------------------
  if (inherits(summaries, "vcmm_ss") ||
      inherits(summaries, "vcmm_accumulator")) {
    agg <- summaries
  } else if (is.list(summaries)) {
    if (length(summaries) == 0L) {
      stop("`summaries` is an empty list.", call. = FALSE)
    }
    is_ss <- vapply(summaries, inherits, logical(1L), "vcmm_ss")
    if (!all(is_ss)) {
      stop("Every element of `summaries` must be a 'vcmm_ss' object ",
           "(from node_summary() or compute_sufficient_stats()).",
           call. = FALSE)
    }
    agg <- Reduce(`+`, summaries)
  } else {
    stop("`summaries` must be a 'vcmm_ss', a 'vcmm_accumulator', or a ",
         "list of 'vcmm_ss' objects.", call. = FALSE)
  }

  p_dim <- nrow(agg$C)
  q_dim <- nrow(agg$ZtZ)

  # ----- Penalty validation ------------------------------------------------
  if (!is.matrix(penalty) || !is.numeric(penalty)) {
    stop("`penalty` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(penalty) != p_dim || ncol(penalty) != p_dim) {
    stop(sprintf("`penalty` must be %d by %d (inferred from summaries).",
                 p_dim, p_dim), call. = FALSE)
  }

  # ----- Build re_cov_state (logic mirrors vcmm(); KEEP IN SYNC) ----------
  re_cov_state <- .build_re_cov_state_for_fit(
    re_cov             = re_cov,
    q_dim              = q_dim,
    n_groups           = n_groups,
    q_left             = q_left,
    control            = control,
    Sigma_left_init    = Sigma_left_init,
    Sigma_right_init   = Sigma_right_init,
    Sigma_2x2_init     = Sigma_2x2_init,
    Sigma_spatial_init = Sigma_spatial_init,
    Sigma_q_init       = Sigma_q_init,
    Omega_G_init       = Omega_G_init
  )

  # ----- Dispatch ----------------------------------------------------------
  fit <- switch(
    method,
    ss  = fit_ss (agg, penalty, control, re_cov_state = re_cov_state),
    csl = fit_csl(agg, penalty, control, re_cov_state = re_cov_state, ...)
  )

  fit$re_cov <- re_cov

  # ----- Optional identifiability re-centering ----------------------------
  # vcmm() does this internally using rowSums(Z); the server can't, so the
  # caller supplies `rowsum_constant` explicitly. NULL means "skip".
  if (!is.null(rowsum_constant)) {
    if (!is.numeric(rowsum_constant) || length(rowsum_constant) != 1L) {
      stop("`rowsum_constant` must be a single numeric value or NULL.",
           call. = FALSE)
    }
    if (abs(rowsum_constant) > 1e-10) {
      alpha_shift   <- mean(fit$alpha)
      fit$alpha     <- fit$alpha - alpha_shift
      fit$beta[1L]  <- fit$beta[1L] + rowsum_constant * alpha_shift
    }
  }

  fit$call <- match.call()
  fit
}

#-------------------------------------------------------------------------------
# Internal: build a re_cov_state from the same set of args vcmm() accepts.
# Extracted so vcmm() and fit_from_summaries() share one code path. Lives
# here (and not in covariance.R) because it is purely an arg-resolution
# helper, not a covariance algorithm.
#-------------------------------------------------------------------------------
.build_re_cov_state_for_fit <- function(re_cov,
                                        q_dim,
                                        n_groups,
                                        q_left,
                                        control,
                                        Sigma_left_init,
                                        Sigma_right_init,
                                        Sigma_2x2_init,
                                        Sigma_spatial_init,
                                        Sigma_q_init,
                                        Omega_G_init) {

  if (identical(re_cov, "diag")) {
    return(NULL)
  }
  if (!(identical(re_cov, "kronecker") || identical(re_cov, "separable"))) {
    stop(sprintf("Unknown re_cov '%s'.", re_cov), call. = FALSE)
  }

  # n_groups
  if (is.null(n_groups) ||
      !is.numeric(n_groups) || length(n_groups) != 1L ||
      n_groups < 1L || n_groups != as.integer(n_groups)) {
    stop(sprintf(
      "re_cov = '%s' requires `n_groups` (a single positive integer).",
      re_cov), call. = FALSE)
  }
  n_groups <- as.integer(n_groups)

  # q_left
  if (identical(re_cov, "kronecker")) {
    if (is.null(q_left)) q_left <- 2L
  } else { # separable
    if (is.null(q_left)) {
      stop("re_cov = 'separable' requires `q_left` (per-group dim).",
           call. = FALSE)
    }
  }
  if (!is.numeric(q_left) || length(q_left) != 1L ||
      q_left < 1L || q_left != as.integer(q_left)) {
    stop("`q_left` must be a single positive integer.", call. = FALSE)
  }
  q_left <- as.integer(q_left)

  if (q_dim != q_left * n_groups) {
    stop(sprintf(
      "Summary q (%d) != q_left * n_groups (%d * %d = %d).",
      q_dim, q_left, n_groups, q_left * n_groups), call. = FALSE)
  }

  # Resolve Sigma_left_init from aliases
  left_aliases <- list(Sigma_left_init = Sigma_left_init,
                       Sigma_2x2_init  = Sigma_2x2_init,
                       Sigma_q_init    = Sigma_q_init)
  supplied_left <- names(left_aliases)[!vapply(left_aliases, is.null, logical(1L))]
  if (length(supplied_left) > 1L) {
    stop(sprintf("Pass only ONE of: %s.",
                 paste(supplied_left, collapse = ", ")), call. = FALSE)
  }
  if (length(supplied_left) == 1L) {
    Sigma_left_init <- left_aliases[[supplied_left]]
  }
  if (is.null(Sigma_left_init)) {
    Sigma_left_init <- (control$sigma_alpha^2) * diag(q_left)
  }
  if (!is.matrix(Sigma_left_init) ||
      !isTRUE(all.equal(dim(Sigma_left_init), c(q_left, q_left)))) {
    stop(sprintf("`Sigma_left_init` (or alias) must be %d by %d.",
                 q_left, q_left), call. = FALSE)
  }

  # Resolve Sigma_right_init from aliases
  right_aliases <- list(Sigma_right_init   = Sigma_right_init,
                        Sigma_spatial_init = Sigma_spatial_init,
                        Omega_G_init       = Omega_G_init)
  supplied_right <- names(right_aliases)[!vapply(right_aliases, is.null, logical(1L))]
  if (length(supplied_right) > 1L) {
    stop(sprintf("Pass only ONE of: %s.",
                 paste(supplied_right, collapse = ", ")), call. = FALSE)
  }
  if (length(supplied_right) == 1L) {
    Sigma_right_init <- right_aliases[[supplied_right]]
  }
  if (is.null(Sigma_right_init)) {
    Sigma_right_init <- diag(n_groups)
  }
  if (!is.matrix(Sigma_right_init) ||
      !isTRUE(all.equal(dim(Sigma_right_init), c(n_groups, n_groups)))) {
    stop(sprintf("`Sigma_right_init` (or alias) must be %d by %d.",
                 n_groups, n_groups), call. = FALSE)
  }

  .new_re_cov_state(
    type        = re_cov,
    Sigma_left  = Sigma_left_init,
    Sigma_right = Sigma_right_init,
    n_groups    = n_groups,
    q_left      = q_left
  )
}
