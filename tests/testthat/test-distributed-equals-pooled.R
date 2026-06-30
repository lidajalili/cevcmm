# Day 21: end-to-end distributed-vs-pooled equivalence. This is the
# math guarantee of the distributed API (Theorem 1 of Lin & Jalili,
# 2026): the per-node summaries aggregated additively must yield a fit
# identical (to BLAS summation noise) to a pooled vcmm() call on the
# same data. Lifts R/distributed.R from 10% to ~80% in one file.

.split_indices <- function(N, K, seed = 7L) {
  set.seed(seed)
  split(seq_len(N), sample.int(K, N, replace = TRUE))
}

test_that("diag re_cov: 3-node fit_from_summaries == pooled vcmm()", {
  d <- gen_diag_data(N = 600L, q = 4L)
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                       update_variance = TRUE)

  fit_pooled <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
                     method = "ss", re_cov = "diag", control = ctrl)

  design <- build_vcmm_design(X = d$x, t = d$t)
  Xd <- design$X_design
  splits <- .split_indices(d$N, 3L)
  summaries <- lapply(splits, function(ii)
    node_summary(d$y[ii], Xd[ii, , drop = FALSE], d$Z[ii, , drop = FALSE])
  )
  fit_dist <- fit_from_summaries(summaries, penalty = design$penalty,
                                  control = ctrl, method = "ss",
                                  re_cov = "diag")

  expect_equal(fit_pooled$beta,      fit_dist$beta,      tolerance = 1e-8)
  expect_equal(fit_pooled$alpha,     fit_dist$alpha,     tolerance = 1e-8)
  expect_equal(fit_pooled$sigma_eps, fit_dist$sigma_eps, tolerance = 1e-8)
})

test_that("kronecker re_cov: distributed == pooled with rowsum_constant", {
  k <- gen_kron_data(N = 600L, G = 8L)
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = sqrt(0.5),
                       update_variance = TRUE)

  fit_pooled <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                     method = "ss", re_cov = "kronecker",
                     n_groups = k$G,
                     Sigma_spatial_init = k$Sigma_spatial,
                     control = ctrl)

  design <- build_vcmm_design(X = k$x, t = k$t)
  Xd <- design$X_design
  splits <- .split_indices(k$N, 3L)
  summaries <- lapply(splits, function(ii)
    node_summary(k$y[ii], Xd[ii, , drop = FALSE], k$Z[ii, , drop = FALSE])
  )
  # rowSums(k$Z) = 2 for OD indicator design; pass through so the
  # identifiability shift matches vcmm()'s automatic re-centering.
  fit_dist <- fit_from_summaries(
    summaries, penalty = design$penalty, control = ctrl,
    method = "ss", re_cov = "kronecker", n_groups = k$G,
    Sigma_spatial_init = k$Sigma_spatial,
    rowsum_constant = 2)

  expect_equal(fit_pooled$beta,  fit_dist$beta,  tolerance = 1e-7)
  expect_equal(fit_pooled$alpha, fit_dist$alpha, tolerance = 1e-7)
})

test_that("separable re_cov: distributed == pooled (csl)", {
  s <- gen_separable_data(N = 800L, G = 10L, q_left = 3L)
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = sqrt(0.5),
                       update_variance = TRUE)

  fit_pooled <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t,
                     method = "csl", re_cov = "separable",
                     n_groups = s$G, q_left = s$q_left,
                     Omega_G_init = s$Omega_G,
                     control = ctrl)

  design <- build_vcmm_design(X = s$x, t = s$t)
  Xd <- design$X_design
  splits <- .split_indices(s$N, 4L)
  summaries <- lapply(splits, function(ii)
    node_summary(s$y[ii], Xd[ii, , drop = FALSE], s$Z[ii, , drop = FALSE])
  )
  fit_dist <- fit_from_summaries(
    summaries, penalty = design$penalty, control = ctrl,
    method = "csl", re_cov = "separable",
    n_groups = s$G, q_left = s$q_left,
    Omega_G_init = s$Omega_G)

  expect_equal(fit_pooled$beta, fit_dist$beta, tolerance = 1e-6)
})

test_that("fit_from_summaries accepts vcmm_accumulator and pre-summed input", {
  d <- gen_diag_data(N = 400L, q = 3L)
  design <- build_vcmm_design(X = d$x, t = d$t)
  Xd <- design$X_design
  ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                       update_variance = TRUE)

  p <- ncol(Xd); q <- ncol(d$Z)
  acc <- init_accumulator(p, q)
  for (b in 1:4) {
    ii <- ((b - 1L) * 100L + 1L):(b * 100L)
    ss <- compute_sufficient_stats(d$y[ii], Xd[ii, , drop = FALSE],
                                    d$Z[ii, , drop = FALSE])
    acc <- accumulate_stats(acc, ss)
  }

  fit_acc <- fit_from_summaries(acc, penalty = design$penalty,
                                 control = ctrl, method = "ss",
                                 re_cov = "diag")
  ss_pooled <- compute_sufficient_stats(d$y, Xd, d$Z)
  fit_pre   <- fit_from_summaries(ss_pooled, penalty = design$penalty,
                                   control = ctrl, method = "ss",
                                   re_cov = "diag")

  expect_equal(fit_acc$beta, fit_pre$beta, tolerance = 1e-10)
})
