# Mirrors inst/validation/day13_plot_validation.R.
# We just check that plot.vcmm_fit runs without erroring under each
# 'which' value -- visual correctness is verified by the validation
# script's PDF artifact, not here.

.with_pdf <- function(expr) {
  pdf(tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

test_that("plot works for which = 1 (varying coefficient curves)", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  expect_silent(.with_pdf(plot(fit, which = 1)))
})

test_that("plot works for which = 2 with data", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  expect_silent(.with_pdf(
    plot(fit, which = 2, data = list(y = d$y, X = d$x, Z = d$Z, t = d$t))
  ))
})

test_that("plot works for which = 3 (QQ / Sigma heatmaps)", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit_d <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
                method = "csl", re_cov = "diag",
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = 0.4,
                                       update_variance = TRUE))
  expect_silent(.with_pdf(plot(fit_d, which = 3)))

  k <- gen_kron_data(N = 400L, G = 6L)
  fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                method = "csl", re_cov = "kronecker",
                n_groups = k$G,
                Sigma_spatial_init = k$Sigma_spatial,
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = sqrt(0.5),
                                       update_variance = TRUE))
  expect_silent(.with_pdf(plot(fit_k, which = 3)))
})

test_that("plot rejects invalid which values", {
  d <- gen_diag_data(N = 150L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  expect_error(.with_pdf(plot(fit, which = 4)))
})
