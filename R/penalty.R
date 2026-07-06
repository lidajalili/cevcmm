#===============================================================================
# B-spline second-order difference penalty matrix
#
# Builds P_lambda used by SS and CSL estimators for penalised B-spline
# coefficients. The intercept column is left unpenalised.
#===============================================================================

#' Build a B-spline second-order difference penalty matrix
#'
#' Constructs a symmetric \eqn{p \times p} penalty matrix
#' \eqn{\mathbf{P}_\lambda} for use with penalised B-spline varying
#' coefficients. For a single varying coefficient
#' (\code{n_blocks = 1}, default), \eqn{p = n\_basis + 1}; the intercept
#' (row/column 1) is unpenalised, and the remaining block is
#' \eqn{\lambda \mathbf{D}_2^\top \mathbf{D}_2} where \eqn{\mathbf{D}_2}
#' is the second-order difference operator on \code{n_basis}
#' coefficients.
#'
#' For \code{n_blocks > 1} (multiple covariates each with their own
#' varying coefficient), the result is block-diagonal: one zero entry
#' for the intercept, followed by \code{n_blocks} copies of the same
#' \code{n_basis x n_basis} penalty block. The total dimension is
#' \code{1 + n_blocks * n_basis}.
#'
#' A small ridge is added if any block is not positive-definite, to
#' guarantee numerical stability in downstream solves.
#'
#' @param n_basis Integer (\eqn{\geq 3}). Number of B-spline basis columns
#'   per varying coefficient.
#' @param lambda Non-negative numeric. Smoothing parameter \eqn{\lambda}.
#' @param n_blocks Integer (\eqn{\geq 1}). Number of varying coefficients
#'   (covariates) sharing the same penalty structure. Default 1.
#'
#' @return A symmetric \code{(1 + n_blocks * n_basis) x (1 + n_blocks * n_basis)}
#'   numeric matrix. Row/column 1 is zero (intercept is unpenalised).
#'
#' @references
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
#' B-splines and penalties. \emph{Statistical Science}, 11(2), 89--121.
#'
#' Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' # Single varying coefficient
#' P1 <- build_penalty_matrix(n_basis = 10, lambda = 1)
#' dim(P1)            # 11 x 11
#' P1[1, 1]           # 0 -- intercept is not penalised
#'
#' # Three varying coefficients sharing the same penalty
#' P3 <- build_penalty_matrix(n_basis = 10, lambda = 1, n_blocks = 3)
#' dim(P3)            # 31 x 31  (1 + 3*10)
#' isSymmetric(P3)
build_penalty_matrix <- function(n_basis, lambda, n_blocks = 1L) {
  if (!is.numeric(n_basis) || length(n_basis) != 1L ||
      n_basis < 3 || n_basis != as.integer(n_basis)) {
    stop("`n_basis` must be a single integer >= 3.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0 ||
      !is.finite(lambda)) {
    stop("`lambda` must be a single non-negative finite numeric.",
         call. = FALSE)
  }
  if (!is.numeric(n_blocks) || length(n_blocks) != 1L ||
      n_blocks < 1 || n_blocks != as.integer(n_blocks)) {
    stop("`n_blocks` must be a single integer >= 1.", call. = FALSE)
  }
  n_basis  <- as.integer(n_basis)
  n_blocks <- as.integer(n_blocks)

  # Second-order difference penalty block (shared across covariates)
  D       <- diff(diag(n_basis), differences = 2)
  P_block <- lambda * crossprod(D)
  P_block <- (P_block + t(P_block)) / 2   # enforce symmetry

  # Minimal ridge so the block stays positive-definite
  min_eig <- min(eigen(P_block, only.values = TRUE, symmetric = TRUE)$values)
  if (min_eig < 1e-10) {
    ridge   <- max(1e-8, -min_eig + 1e-8)
    P_block <- P_block + diag(ridge, n_basis)
  }

  # Assemble block-diagonal: intercept (0) + n_blocks copies of P_block
  p      <- 1L + n_blocks * n_basis
  P_full <- matrix(0, p, p)
  for (k in seq_len(n_blocks)) {
    start <- 2L + (k - 1L) * n_basis
    end   <- 1L + k * n_basis
    P_full[start:end, start:end] <- P_block
  }

  P_full
}
