# Day 21 patch C: defensive-branch coverage to push from 84.78% to ~88%.
# Hits the user-supplied-pilot path in fit_csl, the input-validation
# branches of predict.vcmm_fit / varying_coef, the vector-X coercion in
# predict, and the verbose-output paths in invert_matrix.

# ---- fit_csl: user-supplied pilot -----------------------------------------

test_that("fit_csl accepts a user-supplied pilot from fit_ss", {
  d <- gen_diag_data(N = 200L, q = 3L)
  design <- build_vcmm_design(X = d$x, t = d$t)
  ss <- compute_sufficient_stats(d$y, design$X_design, d$Z)
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                       update_variance = TRUE, max_iter = 5L)
  my_pilot <- fit_ss(ss, design$penalty, ctrl)
  fit <- fit_csl(ss, design$penalty, ctrl, pilot = my_pilot)
  expect_s3_class(fit, "vcmm_fit")
  expect_equal(fit$method, "CSL")
})

test_that("fit_csl rejects malformed user-supplied pilot", {
  d <- gen_diag_data(N = 200L, q = 3L)
  design <- build_vcmm_design(X = d$x, t = d$t)
  ss <- compute_sufficient_stats(d$y, design$X_design, d$Z)
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4)

  # not a vcmm_fit at all
  expect_error(
    fit_csl(ss, design$penalty, ctrl,
            pilot = list(beta = 1, alpha = 1)),
    "vcmm_fit")

  # wrong beta length
  bad_pilot <- fit_ss(ss, design$penalty, ctrl)
  bad_pilot$beta <- c(bad_pilot$beta, 0)
  expect_error(fit_csl(ss, design$penalty, ctrl, pilot = bad_pilot),
               "pilot\\$beta")

  # wrong alpha length
  bad_pilot2 <- fit_ss(ss, design$penalty, ctrl)
  bad_pilot2$alpha <- c(bad_pilot2$alpha, 0)
  expect_error(fit_csl(ss, design$penalty, ctrl, pilot = bad_pilot2),
               "pilot\\$alpha")
})

# ---- predict.vcmm_fit: defensive validation branches ---------------------

test_that("predict.vcmm_fit rejects malformed newdata", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  expect_error(predict(fit),                                "newdata")
  expect_error(predict(fit, newdata = 42),                  "list")
  # missing t
  expect_error(predict(fit, newdata = list(X = d$x)),       "numeric")
  # missing X
  expect_error(predict(fit, newdata = list(t = d$t)),       "varying-coefficient")
  # nrow(X) doesn't match length(t)
  expect_error(
    predict(fit, newdata = list(t = d$t, X = c(d$x, 1))),
    "does not match")
  # NA in X
  X_bad <- d$x; X_bad[1] <- NA
  expect_error(
    predict(fit, newdata = list(t = d$t, X = X_bad)),
    "NA")
  # ncol(Z) mismatch
  expect_error(
    predict(fit, newdata = list(t = d$t, X = d$x,
                                Z = matrix(0, d$N, 99))),
    "does not match")
  # NA in Z
  Z_bad <- d$Z; Z_bad[1, 1] <- NA
  expect_error(
    predict(fit, newdata = list(t = d$t, X = d$x, Z = Z_bad)),
    "NA")
})

test_that("predict accepts plain vector X for K = 1 fit", {
  # Exercises the is.vector(X_new) coercion in predict.vcmm_fit.
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))
  pred <- predict(fit, newdata = list(t = d$t, X = d$x, Z = d$Z))
  expect_equal(length(pred), d$N)
})

# ---- varying_coef: defensive validation branches -------------------------

test_that("varying_coef rejects bad inputs", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

  expect_error(varying_coef("not a fit", t_new = 0.5),   "vcmm_fit")
  expect_error(varying_coef(fit, t_new = "bad"),          "numeric")
  expect_error(varying_coef(fit, t_new = 0.5, k = 99L),   "lie in")
})

# ---- invert_matrix verbose paths -----------------------------------------

test_that("invert_matrix verbose=TRUE prints diagnostic info", {
  A <- make_spd(30L)
  # C++ path: Cholesky line
  expect_output(invert_matrix(A, q = 30L, verbose = TRUE, use_cpp = TRUE),
                "Cholesky")
  # R legacy path: solve() line
  expect_output(invert_matrix(A, q = 30L, verbose = TRUE, use_cpp = FALSE),
                "solve")
})
