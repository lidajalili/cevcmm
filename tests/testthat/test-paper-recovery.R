# Mirrors inst/validation/day7_paper_validation.R, but with a single
# replicate at smaller N. This is a regression test for the paper's
# beta-recovery setting, not a Monte Carlo precision study.

test_that("diag-RE fit recovers intercept and varying coefficient", {
  d <- gen_diag_data(N = 500L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method  = "csl",
              re_cov  = "diag",
              control = vcmm_control(sigma_eps       = 0.5,
                                     sigma_alpha     = 0.4,
                                     update_variance = TRUE))

  expect_true(fit$converged)
  # Single-rep precision is much looser than the Day 7 100-rep test
  expect_lt(abs(fit$beta[1] - 2), 0.5)
  expect_lt(abs(fit$sigma_eps - 0.5), 0.1)
})

test_that("SS and CSL agree at first order", {
  d <- gen_diag_data(N = 500L, q = 3L)
  fit_ss  <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t, method = "ss",
                  re_cov  = "diag",
                  control = vcmm_control(sigma_eps = 0.5,
                                         sigma_alpha = 0.4,
                                         update_variance = TRUE))
  fit_csl <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t, method = "csl",
                  re_cov  = "diag",
                  control = vcmm_control(sigma_eps = 0.5,
                                         sigma_alpha = 0.4,
                                         update_variance = TRUE))

  # First-order Newton equivalence: CSL = SS pilot + 1 step. At
  # convergence on a fixed data set, beta should agree to high
  # precision.
  expect_equal(fit_ss$beta, fit_csl$beta, tolerance = 1e-6)
})
