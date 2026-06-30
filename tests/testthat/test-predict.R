# Mirrors inst/validation/day12_predict_validation.R.
# All predict() calls supply newdata explicitly because the package
# convention is to require it (preventing accidental in-sample
# prediction that's easy to confuse with cross-validation results).

.predict_newdata <- function(d) {
  list(t = d$t, X = d$x, Z = d$Z)
}

test_that("predict default (subject-specific) returns N-length vector", {
  d <- gen_diag_data(N = 300L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  yhat <- predict(fit, newdata = .predict_newdata(d))
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

  nd <- .predict_newdata(d)
  yhat_subj <- predict(fit, newdata = nd, include_random = TRUE)
  yhat_marg <- predict(fit, newdata = nd, include_random = FALSE)

  expect_false(isTRUE(all.equal(yhat_subj, yhat_marg)))
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

  out <- predict(fit, newdata = .predict_newdata(d), se.fit = TRUE)
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

  expect_gt(attr(logLik(fit_k), "df"), attr(logLik(fit_d), "df"))
})
