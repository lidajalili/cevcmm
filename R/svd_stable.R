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
# Internal: RSpectra-based truncated SVD pseudo-inverse (Day 18)
#-------------------------------------------------------------------------------

#' Truncated SVD pseudo-inverse via RSpectra (internal)
#'
#' Computes the Moore-Penrose pseudo-inverse using RSpectra::svds to
#' extract only the top-k singular triplets above a truncation tolerance.
#' Useful for large q (>= 200) ill-conditioned matrices whose effective
#' rank is much smaller than the matrix dimension; in that regime the
#' iterative Lanczos approach inside RSpectra is much faster than a
#' dense LAPACK svd().
#'
#' The function iteratively grows k starting from 50, doubling each
#' round, until either the smallest singular value found falls below
#' the LAPACK-style truncation tolerance, or k reaches an internal cap
#' (default approximately half of the matrix dimension n), or k reaches
#' n minus 1, which is the largest value that RSpectra::svds supports.
#'
#' If the matrix turns out to have near-full rank (no truncation point
#' found within the cap), this function returns NULL so the caller can
#' fall back to dense LAPACK SVD.
#'
#' @param A Numeric square matrix.
#' @param rel_tol Optional relative truncation tolerance. If NULL
#'   (default), uses the LAPACK convention based on machine epsilon
#'   times the largest singular value times the matrix dimension.
#' @param k_max Maximum number of singular values to extract before
#'   giving up. Defaults to approximately half of the matrix size.
#' @param verbose Logical. If TRUE, print the iterations and final
#'   effective rank.
#'
#' @return A list with the same shape as svd_pseudo_inverse (inverse,
#'   condition_number, effective_rank, method) on success; NULL if the
#'   truncation point is not found within k_max, signalling the caller
#'   to fall back to LAPACK.
#'
#' @keywords internal
#' @noRd
svd_pseudo_inverse_rspectra <- function(A,
                                        rel_tol = NULL,
                                        k_max   = NULL,
                                        verbose = FALSE) {
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("svd_pseudo_inverse_rspectra() requires the RSpectra package.\n",
         "  install.packages('RSpectra')",
         call. = FALSE)
  }
  
  n <- nrow(A)
  if (is.null(k_max)) k_max <- min(n - 1L, n %/% 2L)
  k_max <- as.integer(max(50L, k_max))    # ensure cap is at least 50
  
  # Initial probe size: 50 (or as large as possible)
  k     <- min(50L, n - 1L)
  sv    <- NULL
  done  <- FALSE
  
  while (!done) {
    sv <- RSpectra::svds(A, k = k, nu = k, nv = k)
    
    tol_eff <- if (is.null(rel_tol)) {
      .Machine$double.eps * sv$d[1] * max(dim(A))
    } else {
      rel_tol * sv$d[1]
    }
    
    if (verbose) {
      cat(sprintf(
        "  [RSpectra] probe k=%d, d[1]=%.3e, d[k]=%.3e, tol=%.3e\n",
        k, sv$d[1], sv$d[k], tol_eff
      ))
    }
    
    if (sv$d[k] < tol_eff) {
      # Found the truncation point: smallest probed singular value is
      # already below tolerance, so any further singular values would be
      # discarded.
      done <- TRUE
    } else if (k >= n - 1L || k >= k_max) {
      # Hit a hard cap. Either RSpectra can't take k = n (only up to n-1)
      # or the user-specified k_max was reached. Caller falls back to
      # LAPACK because we never found the truncation point.
      if (verbose) {
        cat(sprintf(
          "  [RSpectra] no truncation point found within k_max=%d; falling back\n",
          k_max
        ))
      }
      return(NULL)
    } else {
      # Need more singular values; double k.
      k <- min(2L * k, k_max, n - 1L)
    }
  }
  
  # Build the truncated pseudo-inverse from singular values above tol.
  tol_eff <- if (is.null(rel_tol)) {
    .Machine$double.eps * sv$d[1] * max(dim(A))
  } else {
    rel_tol * sv$d[1]
  }
  keep   <- sv$d > tol_eff
  d_keep <- sv$d[keep]
  U_keep <- sv$u[, keep, drop = FALSE]
  V_keep <- sv$v[, keep, drop = FALSE]
  
  if (verbose) {
    cat(sprintf("  [RSpectra] final: kept %d / %d probed (rank-deficient by %d)\n",
                length(d_keep), k, n - length(d_keep)))
  }
  
  list(
    inverse          = V_keep %*% (1 / d_keep * t(U_keep)),
    condition_number = sv$d[1] / max(d_keep[length(d_keep)], 1e-16),
    effective_rank   = length(d_keep),
    method           = "RSpectra truncated SVD"
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
#'   \item If \code{q < 100}: try Cholesky (via
#'     \code{invert_spd_cpp()}), falling through to LU
#'     (\code{invert_general_cpp()}) and finally SVD pseudo-inverse if
#'     the matrix is not positive-definite.
#'   \item If \code{q >= 100}: route through the SVD pseudo-inverse path.
#'     By default this is the dense LAPACK split-merge variant
#'     (\code{svd_pseudo_inverse()}, paper Algorithm 2). When the user
#'     opts in (see the \code{method} argument below), the iterative
#'     truncated SVD from \code{RSpectra::svds()} is used instead --
#'     faster when the matrix has effective rank much smaller than
#'     \code{q}, slower otherwise.
#' }
#' The Cholesky fast path is roughly 2-3x faster than the original
#' \code{kappa(A) + solve(A)} R path on VCMM K matrices. See the Day 17
#' and Day 18 validation scripts in \code{inst/validation/} for the
#' bit-equivalence and timing details.
#'
#' If you need to reproduce the original R-only path (e.g. for a
#' bit-equivalence test), pass \code{use_cpp = FALSE}.
#'
#' @param A Numeric square matrix to invert.
#' @param q Optional integer. Routing dimension used to pick the inversion
#'   strategy (defaults to \code{nrow(A)}). Pass an explicit value if you
#'   know \code{A} is a curvature block with a meaningful dimension that
#'   differs from its row count.
#' @param verbose Logical. If \code{TRUE}, print the chosen method.
#' @param use_cpp Logical. If \code{TRUE} (default since Day 17), use the
#'   RcppArmadillo Cholesky/LU backend. If \code{FALSE}, use the original
#'   pure-R path with \code{kappa(A)} plus \code{solve(A)}.
#' @param method Character string, one of \code{"auto"} (default),
#'   \code{"lapack"}, or \code{"rspectra"}. Only affects the
#'   \code{q >= 100} branch. With \code{"auto"}, the routing checks
#'   \code{getOption("cevcmm.use_rspectra", FALSE)}: if \code{TRUE} and
#'   the RSpectra package is installed, the truncated-SVD path is used;
#'   otherwise the dense LAPACK split-merge SVD is used. With
#'   \code{"lapack"}, always use the dense LAPACK path (the original
#'   Algorithm 2 behaviour). With \code{"rspectra"}, force the truncated
#'   SVD via RSpectra; if the matrix turns out to be full-rank, the
#'   function silently falls back to LAPACK so results are always
#'   correct.
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
invert_matrix <- function(A,
                          q       = NULL,
                          verbose = FALSE,
                          use_cpp = TRUE,
                          method  = c("auto", "lapack", "rspectra")) {
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop(sprintf("`A` must be square; got %d x %d.", nrow(A), ncol(A)),
         call. = FALSE)
  }
  method <- match.arg(method)
  if (is.null(q)) q <- nrow(A)
  
  if (!isTRUE(use_cpp)) {
    return(.invert_matrix_R_legacy(A, q = q, verbose = verbose))
  }
  
  # ------------------------- C++ fast path --------------------------------
  if (q < 100) {
    result <- tryCatch(invert_spd_cpp(A), error = function(e) NULL)
    if (!is.null(result)) {
      if (verbose) cat(sprintf("  [Invert] Cholesky (q=%d)\n", q))
      return(result)
    }
    
    result <- tryCatch(invert_general_cpp(A), error = function(e) NULL)
    if (!is.null(result)) {
      if (verbose) {
        cat(sprintf("  [Invert] LU (q=%d, Cholesky declined)\n", q))
      }
      return(result)
    }
    
    if (verbose) {
      cat(sprintf("  [Invert] both Cholesky and LU failed at q=%d; SVD fallback\n",
                  q))
    }
  }
  
  # ------------------------- SVD path -------------------------------------
  # Day 18: optionally use RSpectra truncated SVD for q >= 100.
  use_rspectra <- switch(
    method,
    "rspectra" = TRUE,
    "lapack"   = FALSE,
    "auto"     = isTRUE(getOption("cevcmm.use_rspectra", FALSE)) &&
      requireNamespace("RSpectra", quietly = TRUE)
  )
  
  if (q >= 100 && use_rspectra) {
    res <- svd_pseudo_inverse_rspectra(A, verbose = verbose)
    if (!is.null(res)) return(res$inverse)
    # RSpectra returned NULL: matrix has near-full rank, fall through to
    # dense SVD below.
    if (verbose) {
      cat("  [Invert] RSpectra found no truncation; falling back to LAPACK\n")
    }
  }
  
  if (q >= 100) {
    res <- svd_pseudo_inverse(A, use_split_merge = TRUE, verbose = verbose)
    return(res$inverse)
  } else {
    return(pinv_cpp(A))
  }
}

# Internal: pure-R legacy implementation used when use_cpp = FALSE.
# Preserved verbatim so the Day-17 bit-equivalence test has a stable
# reference and so the package still works if the compiled .so is
# absent.
.invert_matrix_R_legacy <- function(A, q, verbose) {
  if (q < 100) {
    cn <- tryCatch(kappa(A), error = function(e) Inf)
    if (cn < 1e8) {
      result <- tryCatch(solve(A), error = function(e) NULL)
      if (!is.null(result)) {
        if (verbose) {
          cat(sprintf("  [Invert] solve() | cond: %.2e\n", cn))
        }
        return(result)
      }
    }
    if (verbose) {
      cat(sprintf("  [Invert] kappa=%.2e - falling back to SVD\n", cn))
    }
  }
  res <- svd_pseudo_inverse(A, use_split_merge = (q >= 100), verbose = verbose)
  res$inverse
}
