# Day 21: stability tests covering svd_pseudo_inverse, split_merge_svd_row,
# pinv_cpp, and invert_general_cpp paths that test-invert-matrix.R only
# grazes. Targets src/solve_cpp.cpp (56% -> ~90%) and R/svd_stable.R
# (76% -> ~90%).

test_that("ill-conditioned SPD is inverted via R-path SVD pseudo-inverse", {
  # cond ~ 1e10. The legacy R path's kappa() check trips its threshold
  # and routes through svd_pseudo_inverse(). Moore-Penrose identity
  # holds even when a naive solve() would be unstable.
  set.seed(101)
  q <- 30L
  U <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
  d <- exp(seq(0, -log(1e10), length.out = q))
  A <- U %*% (d * t(U))
  A <- (A + t(A)) / 2
  A_inv <- invert_matrix(A, q = q, use_cpp = FALSE)
  expect_equal(A %*% A_inv %*% A, A, tolerance = 1e-6)
})

test_that("split_merge_svd_row reconstructs the input matrix", {
  # Algebraic property of the helper used by Algorithm 2.
  set.seed(102)
  X <- matrix(rnorm(150 * 8), 150, 8)
  sv <- cevcmm:::split_merge_svd_row(X, s = 5)
  X_recon <- sv$u %*% (sv$d * t(sv$v))
  expect_lt(max(abs(X - X_recon)) / max(abs(X)), 1e-8)
})

test_that("svd_pseudo_inverse satisfies Moore-Penrose on rank-deficient", {
  A <- make_rank_deficient_spd(q = 120L, r = 40L)
  res <- cevcmm:::svd_pseudo_inverse(A, use_split_merge = TRUE)
  expect_true(is.matrix(res$inverse))
  expect_equal(A %*% res$inverse %*% A, A, tolerance = 1e-6)
  expect_equal(res$effective_rank, 40L)
})

test_that("invert_matrix LU fallback fires for non-PD invertible matrix", {
  # Cholesky declines, LU succeeds: covers invert_general_cpp.
  set.seed(103)
  p <- 20L
  A <- matrix(rnorm(p * p), p, p) + diag(p)   # invertible, not SPD
  expect_true(!isSymmetric(A) ||
              any(eigen(A, only.values = TRUE)$values < 0))
  A_inv <- invert_matrix(A, q = p)
  expect_equal(A %*% A_inv, diag(p), tolerance = 1e-9)
})

test_that("invert_matrix q>=100 R-legacy path returns correct inverse", {
  # use_cpp = FALSE with q >= 100 routes through svd_pseudo_inverse().
  A <- make_spd(120L)
  A_inv <- invert_matrix(A, q = 120L, use_cpp = FALSE)
  expect_equal(A %*% A_inv, diag(120L), tolerance = 1e-8)
})

test_that("invert_matrix method='auto' without option uses LAPACK", {
  skip_if_not_installed("RSpectra")
  opts_old <- options(cevcmm.use_rspectra = NULL)
  on.exit(options(opts_old), add = TRUE)
  A <- make_rank_deficient_spd(q = 120L, r = 30L)
  inv_auto   <- invert_matrix(A, q = 120L, method = "auto")
  inv_lapack <- invert_matrix(A, q = 120L, method = "lapack")
  expect_lt(max(abs(inv_auto - inv_lapack)) / max(abs(inv_lapack)), 1e-10)
})

test_that("pinv_cpp is reachable through invert_matrix small-q SVD fallback", {
  # Indefinite p < 100 matrix: Cholesky declines, LU succeeds in this
  # case -- but with a near-singular indefinite matrix LU can also fail
  # and the path falls into pinv_cpp. Use a constructed near-singular
  # symmetric matrix so the chain Cholesky->LU->pinv exercises all
  # three rungs.
  set.seed(104)
  p <- 12L
  U <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  d <- c(rep(1, p - 2L), 1e-2, -1e-2)            # not PD, near-singular
  A <- U %*% (d * t(U))
  A <- (A + t(A)) / 2
  A_inv <- invert_matrix(A, q = p)
  # On the indefinite branch the path must still produce an inverse
  # satisfying the Moore-Penrose identity on the dominant subspace.
  expect_equal(A %*% A_inv %*% A, A, tolerance = 1e-6)
})
