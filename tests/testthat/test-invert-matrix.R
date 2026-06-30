# Mirrors inst/validation/day17_rcpp_solve_validation.R.

test_that("invert_matrix produces a correct inverse on SPD input", {
  A <- make_spd(10L)
  A_inv <- invert_matrix(A)
  expect_equal(A %*% A_inv, diag(10L), tolerance = 1e-10)
})

test_that("invert_matrix R legacy and C++ Cholesky agree on SPD", {
  for (p in c(10L, 30L, 50L)) {
    A <- make_spd(p)
    inv_R   <- invert_matrix(A, q = p, use_cpp = FALSE)
    inv_cpp <- invert_matrix(A, q = p, use_cpp = TRUE)
    tol <- max(1e-12, p * 1e-13)
    expect_equal(inv_R, inv_cpp, tolerance = tol,
                 label = sprintf("p = %d", p))
  }
})

test_that("invert_matrix falls back to LU on non-PD input", {
  A_indef <- make_indef(p = 10L)
  # Pre-condition: matrix really is not PD
  expect_true(
    any(eigen(A_indef, symmetric = TRUE, only.values = TRUE)$values < 0)
  )
  # Should still produce a valid inverse via LU fallback
  A_inv <- invert_matrix(A_indef, q = 10L)
  expect_equal(A_indef %*% A_inv, diag(10L), tolerance = 1e-9)
})

test_that("invert_matrix rejects non-square or non-numeric input", {
  expect_error(invert_matrix(matrix(1:6, 2, 3)), "must be square")
  expect_error(invert_matrix(matrix("x", 2, 2)), "must be a numeric matrix")
})

test_that("invert_matrix q >= 100 path returns a valid inverse", {
  A <- make_spd(120L)
  A_inv <- invert_matrix(A, q = 120L)
  expect_equal(A %*% A_inv, diag(120L), tolerance = 1e-8)
})

# Day 18 RSpectra path (skipped if RSpectra not installed) ---------------------

test_that("RSpectra path agrees with LAPACK on rank-deficient input", {
  skip_if_not_installed("RSpectra")

  A <- make_rank_deficient_spd(q = 150L, r = 30L)
  inv_lapack   <- invert_matrix(A, q = 150L, method = "lapack")
  inv_rspectra <- invert_matrix(A, q = 150L, method = "rspectra")

  # Both compute the Moore-Penrose pseudo-inverse of the same retained
  # subspace; relative agreement should be tight.
  rel_diff <- max(abs(inv_lapack - inv_rspectra)) / max(abs(inv_lapack))
  expect_lt(rel_diff, 1e-5)

  # Both should satisfy the Moore-Penrose identity A %*% A^+ %*% A == A
  expect_equal(A %*% inv_lapack %*% A,   A, tolerance = 1e-8)
  expect_equal(A %*% inv_rspectra %*% A, A, tolerance = 1e-8)
})

test_that("RSpectra falls back to LAPACK on full-rank input", {
  skip_if_not_installed("RSpectra")

  A <- make_rank_deficient_spd(q = 120L, r = 120L, cond_kept = 1e3)
  inv_lapack <- invert_matrix(A, q = 120L, method = "lapack")
  inv_auto   <- invert_matrix(A, q = 120L, method = "rspectra")

  # After fallback the two paths produce identical output
  rel_diff <- max(abs(inv_lapack - inv_auto)) / max(abs(inv_lapack))
  expect_lt(rel_diff, 1e-8)
})

test_that("cevcmm.use_rspectra option toggles auto routing", {
  skip_if_not_installed("RSpectra")
  opts_old <- options(cevcmm.use_rspectra = TRUE)
  on.exit(options(opts_old), add = TRUE)

  A <- make_rank_deficient_spd(q = 150L, r = 30L)
  inv_auto <- invert_matrix(A, q = 150L, method = "auto")
  expect_equal(A %*% inv_auto %*% A, A, tolerance = 1e-8)
})
