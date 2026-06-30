# Mirrors inst/validation/day9_separable_validation_v3.R.

test_that("separable covariance fit converges with arbitrary q_left", {
  s <- gen_separable_data(N = 800L, G = 10L, q_left = 3L)
  fit <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t,
              method       = "csl",
              re_cov       = "separable",
              n_groups     = s$G,
              q_left       = s$q_left,
              Omega_G_init = s$Omega_G,
              control      = vcmm_control(
                sigma_eps       = 0.5,
                sigma_alpha     = sqrt(0.5),
                update_variance = TRUE))

  expect_true(fit$converged)
  expect_equal(fit$re_cov_state$type, "separable")
  expect_equal(dim(fit$re_cov_state$Sigma_left),
               c(s$q_left, s$q_left))
})

test_that("separable with q_left=2 is backward-compatible with kronecker", {
  # Day 9 backward-compatibility: when q_left=2, separable should
  # behave identically to a Day 8 kronecker fit on the same data.
  k <- gen_kron_data(N = 600L, G = 8L)

  fit_kron <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                   method             = "csl",
                   re_cov             = "kronecker",
                   n_groups           = k$G,
                   Sigma_spatial_init = k$Sigma_spatial,
                   control            = vcmm_control(
                     sigma_eps       = 0.5,
                     sigma_alpha     = sqrt(0.5),
                     update_variance = TRUE))

  fit_sep <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                  method       = "csl",
                  re_cov       = "separable",
                  n_groups     = k$G,
                  q_left       = 2L,
                  Omega_G_init = k$Sigma_spatial,
                  control      = vcmm_control(
                    sigma_eps       = 0.5,
                    sigma_alpha     = sqrt(0.5),
                    update_variance = TRUE))

  expect_true(fit_kron$converged)
  expect_true(fit_sep$converged)
  # Beta vectors should agree to machine precision
  expect_equal(fit_kron$beta, fit_sep$beta, tolerance = 1e-8)
})
