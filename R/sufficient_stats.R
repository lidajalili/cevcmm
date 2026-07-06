#===============================================================================
# Sufficient statistics for VCMM estimation
#
# Implements the per-node summary Gamma_s = (a, b, C, ZtZ, Zty, XtZ) from
# Algorithm 1 of Jalili and Lin (2025). For the normal linear VCMM these
# summaries are computed once per data batch and aggregated additively.
#===============================================================================

#-------------------------------------------------------------------------------
# Internal: validate (y, X, Z) inputs to compute_sufficient_stats()
#-------------------------------------------------------------------------------
.check_yxz <- function(y, X, Z) {
  if (!is.numeric(y)) {
    stop("`y` must be numeric.", call. = FALSE)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix.", call. = FALSE)
  }

  n <- length(y)
  if (nrow(X) != n) {
    stop(sprintf("nrow(X) = %d does not match length(y) = %d.", nrow(X), n),
         call. = FALSE)
  }
  if (nrow(Z) != n) {
    stop(sprintf("nrow(Z) = %d does not match length(y) = %d.", nrow(Z), n),
         call. = FALSE)
  }
  if (anyNA(y) || anyNA(X) || anyNA(Z)) {
    stop("NA values are not allowed in `y`, `X`, or `Z`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Compute one batch or node sufficient statistics for a normal linear VCMM
#'
#' Computes the six-component summary that fully encodes one node's
#' contribution to the normal linear VCMM likelihood. These statistics are
#' additive across nodes and batches, so a central server can recover the
#' full-data estimator by summing per-node summaries -- without ever seeing
#' the raw response or design matrices.
#'
#' The returned components are:
#' \itemize{
#'   \item \code{a}: \code{sum(y^2)}, a scalar.
#'   \item \code{b}: \code{crossprod(X, y)}, dimension p by 1.
#'   \item \code{C}: \code{crossprod(X)}, dimension p by p.
#'   \item \code{ZtZ}: \code{crossprod(Z)}, dimension q by q.
#'   \item \code{Zty}: \code{crossprod(Z, y)}, dimension q by 1.
#'   \item \code{XtZ}: \code{crossprod(X, Z)}, dimension p by q.
#' }
#'
#' Since Day 16 the default backend is a RcppArmadillo implementation
#' (\code{use_cpp = TRUE}). Pass \code{use_cpp = FALSE} to fall back to
#' the pure-R reference path; this is used by the Day-16 bit-equivalence
#' validation and is also useful for debugging.
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric n by p fixed-effects design matrix (intercept plus
#'   spline basis columns).
#' @param Z Numeric n by q random-effects design matrix.
#' @param use_cpp Logical. If \code{TRUE} (default), dispatch to the
#'   RcppArmadillo backend \code{compute_sufficient_stats_cpp()}.
#'   If \code{FALSE}, use the pure-R \code{crossprod()} reference path.
#'   The two paths agree to within floating-point summation order
#'   (typically below 1e-13 relative).
#'
#' @return A list of class \code{"vcmm_ss"} with elements \code{a}, \code{b},
#'   \code{C}, \code{ZtZ}, \code{Zty}, \code{XtZ}, and \code{n_obs} (the
#'   number of observations summarized).
#'
#' @references
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @family sufficient statistics
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 3; q <- 2
#' X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
#' Z <- matrix(rnorm(n * q), n, q)
#' y <- rnorm(n)
#'
#' ss <- compute_sufficient_stats(y, X, Z)
#' str(ss)
compute_sufficient_stats <- function(y, X, Z, use_cpp = TRUE) {
  .check_yxz(y, X, Z)

  if (isTRUE(use_cpp)) {
    # RcppArmadillo backend. as.numeric() flattens y if it arrived as a
    # 1-column matrix; as.matrix() guarantees X and Z aren't data.frames.
    out <- compute_sufficient_stats_cpp(
      as.numeric(y),
      as.matrix(X),
      as.matrix(Z)
    )
  } else {
    # Pure-R reference path. Kept for the Day-16 bit-equivalence check
    # and as a fallback if the compiled .so is unavailable (e.g. when
    # users source individual R files directly without re-installing).
    out <- list(
      a     = sum(y^2),
      b     = crossprod(X, y),
      C     = crossprod(X),
      ZtZ   = crossprod(Z),
      Zty   = crossprod(Z, y),
      XtZ   = crossprod(X, Z),
      n_obs = length(y)
    )
  }
  class(out) <- c("vcmm_ss", "list")
  out
}

#' Initialize an empty sufficient-statistics accumulator
#'
#' Allocates a zero-filled accumulator with the correct dimensions to receive
#' batched calls to \code{accumulate_stats()}. Use this once before a
#' streaming loop over data batches or nodes.
#'
#' @param p Integer. Number of fixed-effects columns (intercept plus spline
#'   basis).
#' @param q Integer. Number of random-effects columns (length of alpha).
#'
#' @return A list of class \code{"vcmm_accumulator"} with the same six
#'   matrix slots as \code{compute_sufficient_stats()}, all initialised to
#'   zero, plus \code{n_obs = 0L}.
#'
#' @family sufficient statistics
#' @export
#'
#' @examples
#' acc <- init_accumulator(p = 5, q = 3)
#' dim(acc$C)   # 5 x 5
#' dim(acc$ZtZ) # 3 x 3
#' acc$n_obs    # 0
init_accumulator <- function(p, q) {
  if (!is.numeric(p) || length(p) != 1L || p < 1 || p != as.integer(p)) {
    stop("`p` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || q < 1 || q != as.integer(q)) {
    stop("`q` must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(p); q <- as.integer(q)

  out <- list(
    a     = 0,
    b     = matrix(0, p, 1),
    C     = matrix(0, p, p),
    ZtZ   = matrix(0, q, q),
    Zty   = matrix(0, q, 1),
    XtZ   = matrix(0, p, q),
    n_obs = 0L
  )
  class(out) <- c("vcmm_accumulator", "list")
  out
}

#' Add one batch or node statistics to a running accumulator
#'
#' Performs the additive aggregation
#' \code{acc <- acc + stats} component by component. After processing all
#' batches, the accumulator holds the full-data sufficient summary.
#'
#' @param acc A \code{vcmm_accumulator} object from
#'   \code{init_accumulator()}.
#' @param stats A \code{vcmm_ss} object from
#'   \code{compute_sufficient_stats()}.
#'
#' @return The updated accumulator (class \code{"vcmm_accumulator"}).
#'
#' @family sufficient statistics
#' @export
#'
#' @examples
#' set.seed(1)
#' n_batch <- 50; p <- 3; q <- 2
#' acc <- init_accumulator(p, q)
#'
#' for (b in 1:3) {
#'   X <- cbind(1, matrix(rnorm(n_batch * (p - 1)), n_batch, p - 1))
#'   Z <- matrix(rnorm(n_batch * q), n_batch, q)
#'   y <- rnorm(n_batch)
#'   ss <- compute_sufficient_stats(y, X, Z)
#'   acc <- accumulate_stats(acc, ss)
#' }
#'
#' acc$n_obs  # 150
accumulate_stats <- function(acc, stats) {
  if (!inherits(acc, "vcmm_accumulator")) {
    stop("`acc` must be a 'vcmm_accumulator' (from init_accumulator()).",
         call. = FALSE)
  }
  if (!inherits(stats, "vcmm_ss")) {
    stop("`stats` must be a 'vcmm_ss' object (from compute_sufficient_stats()).",
         call. = FALSE)
  }

  # Dimension checks
  if (!identical(dim(acc$C), dim(stats$C))) {
    stop(sprintf(
      "Dimension mismatch: acc$C is %s, stats$C is %s.",
      paste(dim(acc$C), collapse = "x"),
      paste(dim(stats$C), collapse = "x")
    ), call. = FALSE)
  }
  if (!identical(dim(acc$ZtZ), dim(stats$ZtZ))) {
    stop(sprintf(
      "Dimension mismatch: acc$ZtZ is %s, stats$ZtZ is %s.",
      paste(dim(acc$ZtZ), collapse = "x"),
      paste(dim(stats$ZtZ), collapse = "x")
    ), call. = FALSE)
  }

  acc$a     <- acc$a     + stats$a
  acc$b     <- acc$b     + stats$b
  acc$C     <- acc$C     + stats$C
  acc$ZtZ   <- acc$ZtZ   + stats$ZtZ
  acc$Zty   <- acc$Zty   + stats$Zty
  acc$XtZ   <- acc$XtZ   + stats$XtZ
  acc$n_obs <- acc$n_obs + stats$n_obs
  acc
}

#' @export
print.vcmm_ss <- function(x, ...) {
  cat("<vcmm_ss>  one-batch sufficient statistics\n")
  cat(sprintf("  n_obs : %d\n", x$n_obs))
  cat(sprintf("  p     : %d   (fixed-effects columns)\n", nrow(x$C)))
  cat(sprintf("  q     : %d   (random-effects columns)\n", nrow(x$ZtZ)))
  cat(sprintf("  a     : %.4g\n", x$a))
  invisible(x)
}

#' @export
print.vcmm_accumulator <- function(x, ...) {
  cat("<vcmm_accumulator>  accumulated sufficient statistics\n")
  cat(sprintf("  n_obs : %d\n", x$n_obs))
  cat(sprintf("  p     : %d   (fixed-effects columns)\n", nrow(x$C)))
  cat(sprintf("  q     : %d   (random-effects columns)\n", nrow(x$ZtZ)))
  cat(sprintf("  a     : %.4g\n", x$a))
  invisible(x)
}
