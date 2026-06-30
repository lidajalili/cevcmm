# Mirrors inst/validation/day8_kronecker_validation.R.

test_that("kronecker covariance fit converges", {
  k <- gen_kron_data(N = 600L, G = 8L)
  fit <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
              method             = "csl",
              re_cov             = "kronecker",
              n_groups           = k$G,
              Sigma_spatial_init = k$Sigma_spatial,
              control            = vcmm_control(
                sigma_eps       = 0.5,
                sigma_alpha     = sqrt(0.5),
                update_variance = TRUE))

  expect_true(fit$converged)
  expect_true(is.finite(fit$sigma_eps))
  expect_lt(abs(fit$sigma_eps - 0.5), 0.15)
})

test_that("kronecker re_cov_state carries Sigma_2x2 with right shape", {
  k <- gen_kron_data(N = 600L, G = 8L)
  fit <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
              method             = "csl",
              re_cov             = "kronecker",
              n_groups           = k$G,
              Sigma_spatial_init = k$Sigma_spatial,
              control            = vcmm_control(
                sigma_eps       = 0.5,
                sigma_alpha     = sqrt(0.5),
                update_variance = TRUE))

  expect_equal(fit$re_cov_state$type, "kronecker")
  expect_true(is.matrix(fit$re_cov_state$Sigma_left))
  expect_equal(dim(fit$re_cov_state$Sigma_left), c(2L, 2L))
  # Sigma must be PD
  expect_true(all(eigen(fit$re_cov_state$Sigma_left,
                        symmetric = TRUE,
                        only.values = TRUE)$values > 0))
})
