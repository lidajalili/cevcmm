# Mirrors inst/validation/day12_predict_validation.R.

test_that("predict default (subject-specific) returns N-length vector", {
  d <- gen_diag_data(N = 300L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  yhat <- predict(fit)
  expect_type(yhat, "double")
  expect_equal(length(yhat), d$N)
  expect_true(all(is.finite(yhat)))
})

test_that("include_random=FALSE drops the RE contribution", {
  d <- gen_diag_data(N = 300L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  yhat_subj <- predict(fit, include_random = TRUE)
  yhat_marg <- predict(fit, include_random = FALSE)
  # The two should differ if random effects are non-trivial
  expect_false(isTRUE(all.equal(yhat_subj, yhat_marg)))
  # Marginal MSE should be larger (less precise) than subject-specific
  mse_subj <- mean((d$y - yhat_subj)^2)
  mse_marg <- mean((d$y - yhat_marg)^2)
  expect_gt(mse_marg, mse_subj)
})

test_that("predict with se.fit returns non-negative SEs", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  out <- predict(fit, se.fit = TRUE)
  expect_named(out, c("fit", "se.fit"))
  expect_equal(length(out$fit), d$N)
  expect_equal(length(out$se.fit), d$N)
  expect_true(all(out$se.fit >= 0))
})

test_that("logLik returns finite value with df and nobs attributes", {
  d <- gen_diag_data(N = 300L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_true(!is.null(attr(ll, "df")))
  expect_equal(attr(ll, "nobs"), d$N)

  # AIC and BIC should follow without error
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

test_that("logLik df differs across re_cov modes", {
  d <- gen_diag_data(N = 300L, q = 3L)
  k <- gen_kron_data(N = 400L, G = 6L)

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

  # kron has more variance-component parameters (3 in Sigma_2x2) than
  # diag (1 sigma_alpha), so df should be larger
  expect_gt(attr(logLik(fit_k), "df"), attr(logLik(fit_d), "df"))
})
