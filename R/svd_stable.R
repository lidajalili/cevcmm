#===============================================================================
# SVD-stabilised matrix inversion
#
# Implements Algorithm 2 of Lin & Jalili (2026): split-merge SVD for large
# matrices plus a master dispatcher that auto-selects between solve() and
# SVD pseudo-inverse based on dimension and condition number.
#===============================================================================

#-------------------------------------------------------------------------------
# Internal: split-merge SVD (Algorithm 2 helper)
#-------------------------------------------------------------------------------

#' Split-merge SVD via row partitioning
#'
#' Computes an SVD of a matrix \eqn{X} by partitioning its rows into
#' \eqn{s} blocks, taking the SVD of each block, and merging the results.
#' More numerically stable than direct \code{svd()} for ill-conditioned
#' large matrices.
#'
#' This is an internal helper for \code{invert_matrix()} and
#' \code{svd_pseudo_inverse()}; end users typically should not call it
#' directly.
#'
#' @param X Numeric matrix.
#' @param s Integer. Number of row partitions (default 10).
#' @param verbose Logical. Currently unused; retained for API compatibility.
#'
#' @return A list with elements \code{u}, \code{d}, \code{v}, matching the
#'   structure returned by base \code{\link[base]{svd}}.
#'
#' @keywords internal
split_merge_svd_row <- function(X, s = 10, verbose = FALSE) {
  m <- nrow(X)
  n <- ncol(X)

  if (s > m) s <- m

  # Equal-ish row partitions
  base_size  <- floor(m / s)
  remainder  <- m %% s
  part_sizes <- rep(base_size, s)
  if (remainder > 0) part_sizes[1:remainder] <- part_sizes[1:remainder] + 1

  row_indices <- vector("list", s)
  start <- 1
  for (i in seq_len(s)) {
    end              <- start + part_sizes[i] - 1
    row_indices[[i]] <- start:end
    start            <- end + 1
  }

  # SVD per partition
  U_list <- D_list <- V_list <- vector("list", s)
  for (i in seq_len(s)) {
    Xi <- X[row_indices[[i]], , drop = FALSE]
    if (nrow(Xi) == 1) {
      nrm <- sqrt(sum(Xi^2))
      if (nrm > 1e-14) {
        U_list[[i]] <- matrix(1, 1, 1)
        D_list[[i]] <- nrm
        V_list[[i]] <- t(Xi) / nrm
      } else {
        U_list[[i]] <- matrix(0, 1, 1)
        D_list[[i]] <- 0
        V_list[[i]] <- matrix(0, n, 1)
      }
    } else {
      sv          <- svd(Xi)
      U_list[[i]] <- sv$u
      D_list[[i]] <- sv$d
      V_list[[i]] <- sv$v
    }
  }

  # Assemble Y = stacked V_i D_i and take its SVD
  Y_parts <- vector("list", s)
  for (i in seq_len(s)) {
    nz <- D_list[[i]] > 1e-14
    if (any(nz)) {
      Y_parts[[i]] <- t(V_list[[i]][, nz, drop = FALSE] %*%
                          diag(D_list[[i]][nz], nrow = sum(nz)))
    }
  }
  Y_parts <- Y_parts[!vapply(Y_parts, is.null, logical(1))]
  if (length(Y_parts) == 0) {
    return(list(u = matrix(0, m, 1), d = 0, v = matrix(0, n, 1)))
  }
  Y    <- do.call(rbind, Y_parts)
  sv_Y <- svd(Y)

  # Block-diagonal U_tilde
  tc <- sum(vapply(U_list, ncol, integer(1)))
  tr <- sum(vapply(U_list, nrow, integer(1)))
  U_tilde <- matrix(0, tr, tc)
  ro <- co <- 0
  for (i in seq_len(s)) {
    ri <- nrow(U_list[[i]]); ci <- ncol(U_list[[i]])
    U_tilde[(ro + 1):(ro + ri), (co + 1):(co + ci)] <- U_list[[i]]
    ro <- ro + ri; co <- co + ci
  }

  # Final U
  min_d   <- min(ncol(U_tilde), nrow(sv_Y$u))
  U_final <- U_tilde[, 1:min_d, drop = FALSE] %*% sv_Y$u[1:min_d, , drop = FALSE]

  if (nrow(U_final) != m) {
    warning("split_merge_svd_row: dimension mismatch, falling back to base svd().",
            call. = FALSE)
    return(svd(X))
  }

  list(u = U_final, d = sv_Y$d, v = sv_Y$v)
}

#-------------------------------------------------------------------------------
# Internal: SVD pseudo-inverse
#-------------------------------------------------------------------------------

#' SVD-based Moore-Penrose pseudo-inverse
#'
#' Computes \eqn{A^{+}} via SVD, optionally using the split-merge variant
#' for large matrices. Singular values below
#' \code{.Machine$double.eps * d[1] * max(dim(A))} are treated as zero.
#'
#' This is an internal helper for \code{invert_matrix()}; end users should
#' call \code{invert_matrix()} instead.
#'
#' @param A Numeric square matrix.
#' @param use_split_merge Logical. If \code{TRUE} and \code{nrow(A) >= 100},
#'   use \code{split_merge_svd_row()}; otherwise use base \code{svd()}.
#' @param verbose Logical. If \code{TRUE}, print method and condition number.
#'
#' @return A list with elements \code{inverse} (the pseudo-inverse matrix),
#'   \code{condition_number}, \code{effective_rank}, and \code{method}
#'   (character: \code{"standard SVD"} or \code{"split-merge SVD"}).
#'
#' @keywords internal
svd_pseudo_inverse <- function(A, use_split_merge = TRUE, verbose = FALSE) {
  n <- nrow(A)

  use_sm <- use_split_merge && n >= 100
  if (use_sm) {
    s_adj  <- max(2, min(10, floor(n / 10)))
    sv_res <- split_merge_svd_row(A, s = s_adj)
    method <- "split-merge SVD"
  } else {
    sv_res <- svd(A)
    method <- "standard SVD"
  }

  U <- sv_res$u; d <- sv_res$d; V <- sv_res$v
  tol   <- .Machine$double.eps * d[1] * max(dim(A))
  d_inv <- ifelse(d > tol, 1 / d, 0)

  cond_num <- d[1] / max(d[length(d)], 1e-16)
  if (verbose) {
    cat(sprintf("  [SVD Inverse] method: %s | cond: %.2e | eff-rank: %d/%d\n",
                method, cond_num, sum(d > tol), length(d)))
  }

  list(
    inverse          = V %*% (d_inv * t(U)),
    condition_number = cond_num,
    effective_rank   = sum(d > tol),
    method           = method
  )
}

#-------------------------------------------------------------------------------
# Public: master matrix-inversion dispatcher
#-------------------------------------------------------------------------------

#' Numerically stable matrix inversion with automatic method selection
#'
#' Inverts a square matrix \eqn{A} using a dispatch rule that balances speed
#' and numerical stability:
#' \itemize{
#'   \item If \code{q < 100}: attempt base \code{\link[base]{solve}}; if the
#'     condition number exceeds \eqn{10^8} or \code{solve()} fails, fall
#'     back to SVD pseudo-inverse.
#'   \item If \code{q >= 100}: skip \code{solve()} and use the SVD
#'     pseudo-inverse directly, via the split-merge variant
#'     \code{split_merge_svd_row()}.
#' }
#' This implements the SVD-stabilised step of Algorithm 2 in
#' Lin & Jalili (2026), used inside the SS and CSL estimators when the
#' aggregated Gram matrices are dense or ill-conditioned.
#'
#' @param A Numeric square matrix to invert.
#' @param q Optional integer. Routing dimension used to pick the inversion
#'   strategy (defaults to \code{nrow(A)}). Pass an explicit value if you
#'   know \code{A} is a curvature block with a meaningful dimension that
#'   differs from its row count.
#' @param verbose Logical. If \code{TRUE}, print the chosen method and
#'   condition number.
#'
#' @return A numeric matrix with the same dimensions as \code{A}.
#'
#' @references
#' Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
#' Varying Coefficient Mixed-Effects Models.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' A <- crossprod(matrix(rnorm(50), 10, 5)) + diag(5)
#' A_inv <- invert_matrix(A)
#' max(abs(A %*% A_inv - diag(5)))  # ~ machine epsilon
invert_matrix <- function(A, q = NULL, verbose = FALSE) {
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop(sprintf("`A` must be square; got %d x %d.", nrow(A), ncol(A)),
         call. = FALSE)
  }
  if (is.null(q)) q <- nrow(A)

  if (q < 100) {
    cn <- tryCatch(kappa(A), error = function(e) Inf)
    if (cn < 1e8) {
      result <- tryCatch(solve(A), error = function(e) NULL)
      if (!is.null(result)) {
        if (verbose) cat(sprintf("  [Invert] solve() | cond: %.2e\n", cn))
        return(result)
      }
    }
    if (verbose) cat(sprintf("  [Invert] kappa=%.2e - falling back to SVD\n", cn))
  }

  # SVD path: q >= 100, or solve() failed
  res <- svd_pseudo_inverse(A, use_split_merge = (q >= 100), verbose = verbose)
  res$inverse
}
