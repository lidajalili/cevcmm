# Mirrors inst/validation/day11_s3_methods_validation.R.

.diag_fit_for_methods <- function() {
  d <- gen_diag_data(N = 300L, q = 3L)
  list(
    fit = vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
               method  = "csl",
               re_cov  = "diag",
               control = vcmm_control(sigma_eps = 0.5,
                                      sigma_alpha = 0.4,
                                      update_variance = TRUE)),
    d   = d
  )
}

test_that("vcmm_fit has expected class and structure", {
  out <- .diag_fit_for_methods()
  expect_s3_class(out$fit, "vcmm_fit")
  expect_true(out$fit$converged)
})

test_that("nobs returns the sample size", {
  out <- .diag_fit_for_methods()
  expect_equal(nobs(out$fit), out$d$N)
})

test_that("coef returns named numeric vector", {
  out <- .diag_fit_for_methods()
  cf <- coef(out$fit)
  expect_type(cf, "double")
  expect_true(length(cf) >= 1L)
  expect_equal(names(cf)[1], "(Intercept)")
})

test_that("fixef returns list with intercept and varying parts", {
  out <- .diag_fit_for_methods()
  fx <- fixef(out$fit)
  expect_true(is.list(fx))
  expect_named(fx, c("intercept", "varying"))
  expect_equal(fx$intercept, coef(out$fit)[1], tolerance = 1e-12)
})

test_that("ranef returns vector for diag, matrix for kron/sep", {
  out_d <- .diag_fit_for_methods()
  re_d  <- ranef(out_d$fit)
  expect_type(re_d, "double")

  k <- gen_kron_data(N = 600L, G = 8L)
  fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                method             = "csl",
                re_cov             = "kronecker",
                n_groups           = k$G,
                Sigma_spatial_init = k$Sigma_spatial,
                control            = vcmm_control(
                  sigma_eps       = 0.5,
                  sigma_alpha     = sqrt(0.5),
                  update_variance = TRUE))
  re_k <- ranef(fit_k)
  expect_true(is.matrix(re_k))
  expect_equal(dim(re_k), c(k$G, 2L))
})

test_that("vcov returns p x p positive-semidefinite matrix", {
  out <- .diag_fit_for_methods()
  V <- vcov(out$fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
  expect_true(all(eigen(V, symmetric = TRUE,
                        only.values = TRUE)$values >= -1e-10))
})

test_that("varying_coef returns curve and SE bands of correct shape", {
  out <- .diag_fit_for_methods()
  t_new <- seq(0, 1, length.out = 21L)
  vc <- varying_coef(out$fit, t_new = t_new, k = 1L, se.fit = TRUE)
  expect_named(vc, c("fit", "se.fit"))
  expect_equal(length(vc$fit), length(t_new))
  expect_equal(length(vc$se.fit), length(t_new))
  expect_true(all(vc$se.fit >= 0))
})

test_that("summary works for all re_cov modes", {
  d <- gen_diag_data(N = 300L, q = 3L)
  k <- gen_kron_data(N = 500L, G = 8L)
  s <- gen_separable_data(N = 600L, G = 8L, q_left = 3L)

  fit_d <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
                method = "csl", re_cov = "diag",
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = 0.4,
                                       update_variance = TRUE))
  fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                method = "csl", re_cov = "kronecker",
                n_groups = k$G,
                Sigma_spatial_init = k$Sigma_spatial,
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = sqrt(0.5),
                                       update_variance = TRUE))
  fit_s <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t,
                method = "csl", re_cov = "separable",
                n_groups = s$G, q_left = s$q_left,
                Omega_G_init = s$Omega_G,
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = sqrt(0.5),
                                       update_variance = TRUE))

  for (f in list(fit_d, fit_k, fit_s)) {
    sm <- summary(f)
    expect_s3_class(sm, "vcmm_summary")
    # print should not error
    expect_silent(invisible(capture.output(print(sm))))
  }
})
