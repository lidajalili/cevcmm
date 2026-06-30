# Mirrors inst/validation/day16_rcpp_suffstats_validation.R.

test_that("compute_sufficient_stats returns a well-formed vcmm_ss object", {
  d  <- gen_diag_data(N = 200L, q = 3L)
  X  <- matrix(d$x, ncol = 1L)
  ss <- compute_sufficient_stats(d$y, X, d$Z)

  expect_s3_class(ss, "vcmm_ss")
  expect_named(ss, c("a", "b", "C", "ZtZ", "Zty", "XtZ", "n_obs"))
  expect_equal(ss$n_obs, d$N)
  expect_type(ss$a, "double")
  expect_true(is.matrix(ss$b))
  expect_equal(dim(ss$b), c(1L, 1L))
  expect_equal(dim(ss$C), c(1L, 1L))
  expect_equal(dim(ss$ZtZ), c(d$q, d$q))
  expect_equal(dim(ss$Zty), c(d$q, 1L))
  expect_equal(dim(ss$XtZ), c(1L, d$q))
})

test_that("compute_sufficient_stats R and C++ paths agree", {
  d <- gen_diag_data(N = 500L, q = 5L)
  X <- matrix(d$x, ncol = 1L)

  ss_R   <- compute_sufficient_stats(d$y, X, d$Z, use_cpp = FALSE)
  ss_cpp <- compute_sufficient_stats(d$y, X, d$Z, use_cpp = TRUE)

  # Tolerance scales with N (BLAS summation order error).
  tol <- max(1e-12, d$N * 1e-13)
  expect_equal(ss_R$a,   ss_cpp$a,   tolerance = tol)
  expect_equal(ss_R$b,   ss_cpp$b,   tolerance = tol)
  expect_equal(ss_R$C,   ss_cpp$C,   tolerance = tol)
  expect_equal(ss_R$ZtZ, ss_cpp$ZtZ, tolerance = tol)
  expect_equal(ss_R$Zty, ss_cpp$Zty, tolerance = tol)
  expect_equal(ss_R$XtZ, ss_cpp$XtZ, tolerance = tol)
  expect_equal(ss_R$n_obs, ss_cpp$n_obs)
})

test_that("vcmm_ss summaries are additive across nodes", {
  d <- gen_diag_data(N = 300L, q = 4L)
  X <- matrix(d$x, ncol = 1L)

  ss_full <- compute_sufficient_stats(d$y, X, d$Z)

  splits <- list(1:100, 101:200, 201:300)
  parts  <- lapply(splits, function(i)
    compute_sufficient_stats(
      d$y[i],
      X[i, , drop = FALSE],
      d$Z[i, , drop = FALSE]
    )
  )
  ss_sum <- Reduce(`+`, parts)

  expect_equal(ss_full$a,     ss_sum$a,     tolerance = 1e-10)
  expect_equal(ss_full$b,     ss_sum$b,     tolerance = 1e-10)
  expect_equal(ss_full$C,     ss_sum$C,     tolerance = 1e-10)
  expect_equal(ss_full$ZtZ,   ss_sum$ZtZ,   tolerance = 1e-10)
  expect_equal(ss_full$Zty,   ss_sum$Zty,   tolerance = 1e-10)
  expect_equal(ss_full$XtZ,   ss_sum$XtZ,   tolerance = 1e-10)
  expect_equal(ss_full$n_obs, ss_sum$n_obs)
})

test_that("compute_sufficient_stats rejects malformed inputs", {
  expect_error(
    compute_sufficient_stats(numeric(0), matrix(0, 0, 1), matrix(0, 0, 1)),
    NA  # empty inputs are fine, just produce trivial output
  )
  expect_error(
    compute_sufficient_stats(1:5, matrix(0, 4, 1), matrix(0, 5, 1))
  )
})
