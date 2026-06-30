# Mirrors inst/validation/day10_distributed_validation.R.
# Focuses on the +.vcmm_ss S3 method: vcmm_ss objects sum across nodes
# producing an aggregate identical to the single-node summary. The full
# end-to-end distributed refit pipeline (which exercises internal design
# construction) lives in inst/validation/day10; the regression value of
# THIS test file is the additivity property, which is the actual math
# guarantee the distributed API depends on.

test_that("+.vcmm_ss adds two summaries to match the full-data summary", {
  d <- gen_diag_data(N = 200L, q = 3L)
  X <- matrix(d$x, ncol = 1L)

  ss1 <- compute_sufficient_stats(d$y[1:100], X[1:100, , drop = FALSE],
                                  d$Z[1:100, , drop = FALSE])
  ss2 <- compute_sufficient_stats(d$y[101:200], X[101:200, , drop = FALSE],
                                  d$Z[101:200, , drop = FALSE])
  ss_sum <- ss1 + ss2

  expect_s3_class(ss_sum, "vcmm_ss")

  ss_full <- compute_sufficient_stats(d$y, X, d$Z)
  expect_equal(ss_full$a,     ss_sum$a,     tolerance = 1e-10)
  expect_equal(ss_full$b,     ss_sum$b,     tolerance = 1e-10)
  expect_equal(ss_full$C,     ss_sum$C,     tolerance = 1e-10)
  expect_equal(ss_full$ZtZ,   ss_sum$ZtZ,   tolerance = 1e-10)
  expect_equal(ss_full$Zty,   ss_sum$Zty,   tolerance = 1e-10)
  expect_equal(ss_full$XtZ,   ss_sum$XtZ,   tolerance = 1e-10)
  expect_equal(ss_full$n_obs, ss_sum$n_obs)
})

test_that("+.vcmm_ss aggregates a 3-way split via Reduce", {
  d <- gen_diag_data(N = 300L, q = 4L)
  X <- matrix(d$x, ncol = 1L)
  splits <- list(1:100, 101:200, 201:300)

  parts <- lapply(splits, function(i)
    compute_sufficient_stats(
      d$y[i], X[i, , drop = FALSE], d$Z[i, , drop = FALSE]
    )
  )
  ss_sum <- Reduce(`+`, parts)

  expect_s3_class(ss_sum, "vcmm_ss")
  expect_equal(ss_sum$n_obs, 300L)

  ss_full <- compute_sufficient_stats(d$y, X, d$Z)
  expect_equal(ss_full$a, ss_sum$a, tolerance = 1e-10)
  expect_equal(ss_full$C, ss_sum$C, tolerance = 1e-10)
})

test_that("+.vcmm_ss preserves the C++ vs R agreement", {
  # If the R path and the C++ path of compute_sufficient_stats both
  # produce additive summaries, their aggregates must agree to the same
  # tolerance as the single-call ones. Catches a class of bugs where a
  # backend silently produces non-additive output.
  d <- gen_diag_data(N = 200L, q = 3L)
  X <- matrix(d$x, ncol = 1L)
  i <- 1:100; j <- 101:200

  sum_R <- compute_sufficient_stats(d$y[i], X[i, , drop = FALSE],
                                    d$Z[i, , drop = FALSE], use_cpp = FALSE) +
           compute_sufficient_stats(d$y[j], X[j, , drop = FALSE],
                                    d$Z[j, , drop = FALSE], use_cpp = FALSE)
  sum_C <- compute_sufficient_stats(d$y[i], X[i, , drop = FALSE],
                                    d$Z[i, , drop = FALSE], use_cpp = TRUE) +
           compute_sufficient_stats(d$y[j], X[j, , drop = FALSE],
                                    d$Z[j, , drop = FALSE], use_cpp = TRUE)

  expect_equal(sum_R$a,     sum_C$a,     tolerance = 1e-10)
  expect_equal(sum_R$b,     sum_C$b,     tolerance = 1e-10)
  expect_equal(sum_R$C,     sum_C$C,     tolerance = 1e-10)
  expect_equal(sum_R$ZtZ,   sum_C$ZtZ,   tolerance = 1e-10)
  expect_equal(sum_R$Zty,   sum_C$Zty,   tolerance = 1e-10)
  expect_equal(sum_R$XtZ,   sum_C$XtZ,   tolerance = 1e-10)
})
