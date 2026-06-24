#===============================================================================
# B-spline second-order difference penalty matrix
#
# Builds P_lambda used by SS and CSL estimators for penalised B-spline
# coefficients. The intercept column is left unpenalised.
#===============================================================================

#' Build a B-spline second-order difference penalty matrix
#'
#' Constructs a symmetric \eqn{p \times p} penalty matrix
#' \eqn{\mathbf{P}_\lambda} for use with cubic B-spline varying coefficients,
#' where \eqn{p = n_{\text{basis}} + 1} (an intercept column is prepended).
#' The penalty matrix is
#' \deqn{
#'   \mathbf{P}_\lambda
#'   = \mathrm{diag}\left(0, \; \lambda\, \mathbf{D}_2^\top \mathbf{D}_2 \right)
#' }
#' where \eqn{\mathbf{D}_2} is the second-order difference operator on
#' \code{n_basis} spline coefficients. The intercept (row/column 1) is not
#' penalised. A small ridge is added if the spline block is not
#' positive-definite, to guarantee numerical stability in downstream
#' solves.
#'
#' @param n_basis Integer (\eqn{\geq 3}). Number of B-spline basis columns.
#' @param lambda Non-negative numeric. Smoothing parameter \eqn{\lambda}.
#'
#' @return A symmetric \eqn{(n\_basis + 1) \times (n\_basis + 1)} numeric
#'   matrix. The first row and column are zero (intercept is unpenalised).
#'
#' @references
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
#' B-splines and penalties. \emph{Statistical Science}, 11(2), 89--121.
#'
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' P <- build_penalty_matrix(n_basis = 10, lambda = 1)
#' dim(P)       # 11 x 11 (intercept + 10 basis columns)
#' P[1, 1]      # 0 -- intercept is not penalised
#' isSymmetric(P)
build_penalty_matrix <- function(n_basis, lambda) {
  if (!is.numeric(n_basis) || length(n_basis) != 1L ||
      n_basis < 3 || n_basis != as.integer(n_basis)) {
    stop("`n_basis` must be a single integer >= 3.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0 ||
      !is.finite(lambda)) {
    stop("`lambda` must be a single non-negative finite numeric.",
         call. = FALSE)
  }
  n_basis <- as.integer(n_basis)

  # Second-order difference matrix on the n_basis spline columns
  D     <- diff(diag(n_basis), differences = 2)   # (n_basis - 2) x n_basis
  P_raw <- lambda * crossprod(D)                  # n_basis x n_basis

  # Enforce exact symmetry numerically
  P_raw <- (P_raw + t(P_raw)) / 2

  # Pad with a zero row/column for the intercept (column 1 of X)
  p      <- n_basis + 1L
  P_full <- matrix(0, p, p)
  P_full[2:p, 2:p] <- P_raw

  # Minimal ridge so the spline block stays positive-definite
  min_eig <- min(eigen(P_raw, only.values = TRUE, symmetric = TRUE)$values)
  if (min_eig < 1e-10) {
    ridge <- max(1e-8, -min_eig + 1e-8)
    P_full[2:p, 2:p] <- P_full[2:p, 2:p] + diag(ridge, n_basis)
  }

  P_full
}
