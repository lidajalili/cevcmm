#===============================================================================
# build_vcmm_design: construct the B-spline VCMM design and penalty
#
# Given covariates X (n x K), a time/index vector t (length n), and
# B-spline settings, builds:
#   X_design = cbind(1,  B * x_1,  B * x_2,  ...,  B * x_K)
#   penalty  = block-diag(0,  P_basis,  P_basis,  ...,  P_basis)
# where B is the n x n_basis B-spline basis matrix.
#
# This is the multi-covariate version of the simulation code's
# build_design_matrix(): the simulation handles the K = 1 case;
# here we support K >= 1 to match the generality of the paper's
# equation (1).
#===============================================================================

#' Build the design matrix and penalty for a VCMM
#'
#' Constructs the fixed-effects design matrix and matching penalty for
#' a VCMM with one or more varying coefficients. Given covariates
#' \code{X} (n by K) and a time/index vector \code{t} (length n), builds
#' the cubic-by-default B-spline basis \code{B} at \code{t} and
#' assembles
#' \itemize{
#'   \item \code{X_design = cbind(1, B * X[,1], B * X[,2], ..., B * X[,K])},
#'     dimension n by (1 + K * n_basis); the first column is the
#'     intercept, then each covariate gets its own n_basis-column block.
#'   \item \code{penalty}: a (1 + K * n_basis) by (1 + K * n_basis)
#'     block-diagonal second-order difference penalty matrix with the
#'     intercept unpenalised.
#' }
#' This is the natural multi-covariate generalisation of the
#' single-covariate design used in the simulation code, matching the
#' paper's general VCMM specification.
#'
#' @param X Numeric n by K matrix (or length-n vector if K = 1) of
#'   covariates that get varying coefficients in \code{t}.
#' @param t Numeric vector of length n. The variable in which the
#'   coefficients vary smoothly (time, index, location, etc.).
#' @param n_basis Integer (\eqn{\geq} \code{degree} + 1) or \code{NULL}.
#'   Number of B-spline basis columns per varying coefficient. If
#'   \code{NULL} (default), chosen as \code{max(floor(n^(1/3)) + 4, 10)}
#'   matching the simulation code's rule of thumb.
#' @param degree Integer. B-spline degree (default 3 = cubic).
#' @param lambda Non-negative numeric. Smoothing parameter for the
#'   penalty (default 1).
#' @param normalize_t Logical. If \code{TRUE} (default), \code{t} is
#'   linearly mapped to \code{[0, 1]} before building the basis. Set
#'   \code{FALSE} only if \code{t} is already in \code{[0, 1]}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{X_design}: numeric n by (1 + K * n_basis) design matrix.
#'   \item \code{penalty}: numeric (1 + K * n_basis) by (1 + K * n_basis)
#'     penalty matrix.
#'   \item \code{B_spline}: numeric n by n_basis basis matrix at \code{t}.
#'   \item \code{internal_knots}, \code{boundary_knots}: spline knot
#'     locations, recorded so predictions at new \code{t} can use the
#'     same basis.
#'   \item \code{degree}, \code{n_basis}, \code{K}, \code{lambda}:
#'     scalars recording how the design was built.
#'   \item \code{normalize_t}, \code{t_min}, \code{t_max}: how to remap
#'     new \code{t} values at prediction time.
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
#' n <- 200
#' t <- runif(n)
#' x1 <- runif(n); x2 <- runif(n)
#'
#' # Single varying coefficient
#' d1 <- build_vcmm_design(X = x1, t = t, n_basis = 8)
#' dim(d1$X_design)   # 200 x 9   (1 + 1*8)
#'
#' # Two varying coefficients
#' d2 <- build_vcmm_design(X = cbind(x1, x2), t = t, n_basis = 8)
#' dim(d2$X_design)   # 200 x 17  (1 + 2*8)
build_vcmm_design <- function(X,
                              t,
                              n_basis     = NULL,
                              degree      = 3L,
                              lambda      = 1,
                              normalize_t = TRUE) {

  # ----- Input validation -------------------------------------------------
  if (is.vector(X) && is.numeric(X)) X <- matrix(X, ncol = 1L)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix or vector.", call. = FALSE)
  }
  if (anyNA(X)) {
    stop("NA values are not allowed in `X`.", call. = FALSE)
  }
  n <- nrow(X)
  K <- ncol(X)

  if (!is.numeric(t) || length(t) != n) {
    stop(sprintf("`t` must be numeric of length %d (= nrow(X)).", n),
         call. = FALSE)
  }
  if (anyNA(t)) {
    stop("NA values are not allowed in `t`.", call. = FALSE)
  }
  if (!is.numeric(degree) || length(degree) != 1L ||
      degree < 1L || degree != as.integer(degree)) {
    stop("`degree` must be a single positive integer.", call. = FALSE)
  }
  degree <- as.integer(degree)

  # Auto-pick n_basis using the simulation code's rule of thumb
  if (is.null(n_basis)) {
    n_basis <- max(floor(n^(1 / 3)) + 4L, 10L)
  }
  if (!is.numeric(n_basis) || length(n_basis) != 1L ||
      n_basis < (degree + 1L) || n_basis != as.integer(n_basis)) {
    stop(sprintf("`n_basis` must be a single integer >= degree + 1 (= %d).",
                 degree + 1L),
         call. = FALSE)
  }
  n_basis <- as.integer(n_basis)

  # ----- Normalise t to [0, 1] if requested -------------------------------
  t_min <- min(t)
  t_max <- max(t)
  if (normalize_t) {
    if (t_max <= t_min) {
      stop("`t` has zero range; cannot normalise.", call. = FALSE)
    }
    t_use <- (t - t_min) / (t_max - t_min)
  } else {
    if (any(t < 0) || any(t > 1)) {
      warning("`t` not in [0, 1] but normalize_t = FALSE; bs() may extrapolate.",
              call. = FALSE)
    }
    t_use <- t
  }

  # ----- B-spline basis (n x n_basis) -------------------------------------
  # df = n_basis with intercept = FALSE yields exactly n_basis columns
  B <- splines::bs(t_use,
                   df             = n_basis,
                   degree         = degree,
                   intercept      = FALSE,
                   Boundary.knots = c(0, 1))
  internal_knots <- attr(B, "knots")
  boundary_knots <- attr(B, "Boundary.knots")
  B <- unclass(B)
  attributes(B) <- list(dim = dim(B))   # keep it a plain matrix

  # ----- Build X_design = cbind(1, B * x_1, ..., B * x_K) -----------------
  # B * x_k means: row i of B scaled by X[i, k]
  blocks <- vector("list", K)
  for (k in seq_len(K)) {
    blocks[[k]] <- B * X[, k]
  }
  X_design <- cbind(1, do.call(cbind, blocks))

  # ----- Matching block-diagonal penalty ---------------------------------
  penalty <- build_penalty_matrix(n_basis = n_basis,
                                  lambda  = lambda,
                                  n_blocks = K)

  list(
    X_design       = X_design,
    penalty        = penalty,
    B_spline       = B,
    internal_knots = internal_knots,
    boundary_knots = boundary_knots,
    degree         = degree,
    n_basis        = n_basis,
    K              = K,
    lambda         = lambda,
    normalize_t    = normalize_t,
    t_min          = t_min,
    t_max          = t_max
  )
}
