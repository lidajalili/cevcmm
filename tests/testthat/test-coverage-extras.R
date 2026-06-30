# Day 21 patch: small extras to push coverage from 80% past 85%.
# Targets the gaps that didn't move on the first pass:
#   - fit_ss.R print kron/sep branch (uncovered: most tests use csl)
#   - predict.R logLik kron determinant + kron se.fit
#   - pinv_cpp + invert_general_cpp via direct ::: calls
#   - print methods for vcmm_ss / vcmm_accumulator / vcmm_control

test_that("SS print method renders for all three re_cov modes", {
  # print.vcmm_fit's kron/sep branch (lines ~275-308 in fit_ss.R) is
  # uncovered until method="ss" runs with those re_cov modes and we
  # print the result on the console.
  k <- gen_kron_data(N = 300L, G = 6L)
  fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                method = "ss", re_cov = "kronecker",
                n_groups = k$G,
                Sigma_spatial_init = k$Sigma_spatial,
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = sqrt(0.5),
                                       update_variance = TRUE))
  expect_output(print(fit_k), "Sigma_2x2")

  s <- gen_separable_data(N = 400L, G = 6L, q_left = 3L)
  fit_s <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t,
                method = "ss", re_cov = "separable",
                n_groups = s$G, q_left = s$q_left,
                Omega_G_init = s$Omega_G,
                control = vcmm_control(sigma_eps = 0.5,
                                       sigma_alpha = sqrt(0.5),
                                       update_variance = TRUE))
  expect_output(print(fit_s), "Sigma_q")
})

test_that("SS with update_variance = FALSE leaves variances at init", {
  d <- gen_diag_data(N = 200L, q = 3L)
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                       update_variance = FALSE)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "ss", re_cov = "diag", control = ctrl)
  expect_equal(fit$sigma_eps,   0.5, tolerance = 1e-12)
  expect_equal(fit$sigma_alpha, 0.4, tolerance = 1e-12)
})

test_that("predict se.fit and logLik exercise the kron branch", {
  k <- gen_kron_data(N = 400L, G = 6L)
  fit <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
              method = "csl", re_cov = "kronecker",
              n_groups = k$G,
              Sigma_spatial_init = k$Sigma_spatial,
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))
  # se.fit for kron exercises the joint W K_inv W' branch in predict.R
  pred <- predict(fit, newdata = list(t = k$t, X = k$x, Z = k$Z),
                  se.fit = TRUE)
  expect_named(pred, c("fit", "se.fit"))
  expect_true(all(pred$se.fit >= 0))
  # logLik with kron re_cov_state hits the kron determinant branch
  ll <- logLik(fit)
  expect_true(is.finite(as.numeric(ll)))
})

test_that("predict marginal (include_random = FALSE) with se.fit", {
  d <- gen_diag_data(N = 200L, q = 3L)
  fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method = "csl", re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))
  # Hits the V_beta-only se branch in predict.R (no Z contribution)
  pred <- predict(fit, newdata = list(t = d$t, X = d$x, Z = d$Z),
                  se.fit = TRUE, include_random = FALSE)
  expect_true(all(pred$se.fit >= 0))
})

test_that("pinv_cpp returns Moore-Penrose pseudo-inverse (direct call)", {
  # Note: routing a rank-1 SPSD matrix through invert_matrix() is
  # unreliable -- arma::inv() at the LU stage can return a numerically
  # garbage inverse for near-singular input without flagging
  # rank-deficiency. To guarantee coverage of pinv_cpp() we call it
  # directly. arma::pinv() uses SVD and handles rank-deficient input
  # correctly under all Armadillo versions.
  set.seed(801)
  v <- rnorm(10)
  A <- v %*% t(v)
  A <- (A + t(A)) / 2
  A_inv <- cevcmm:::pinv_cpp(A)
  expect_equal(A %*% A_inv %*% A, A, tolerance = 1e-8)
})

test_that("pinv_cpp accepts a custom positive tolerance", {
  # Hits the tol > 0 branch (vs the default tol = -1.0 branch).
  A <- diag(c(1, 0.5, 1e-15))
  A_inv <- cevcmm:::pinv_cpp(A, tol = 1e-10)
  expect_true(is.matrix(A_inv))
  # Singular value 1e-15 is below the cutoff 1e-10 and gets truncated,
  # so the third inverse entry should be 0 rather than 1e15.
  expect_lt(abs(A_inv[3, 3]), 1)
})

test_that("invert_general_cpp errors on a truly singular matrix", {
  # Hits the Rcpp::stop() branch in invert_general_cpp(). The all-ones
  # matrix is rank 1, so LU factorisation has zero pivots and
  # arma::inv() returns false.
  A_sing <- matrix(1, 5, 5)
  expect_error(cevcmm:::invert_general_cpp(A_sing), "singular")
})

test_that("print methods for vcmm_ss / accumulator / control render", {
  d <- gen_diag_data(N = 100L, q = 2L)
  X <- matrix(d$x, ncol = 1L)
  ss <- compute_sufficient_stats(d$y, X, d$Z)
  acc <- init_accumulator(p = 1L, q = 2L)
  acc <- accumulate_stats(acc, ss)

  expect_output(print(ss),  "vcmm_ss")
  expect_output(print(acc), "vcmm_accumulator")
  expect_output(print(vcmm_control()), "vcmm_control")
})
