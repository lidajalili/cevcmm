# Mirrors inst/validation/day10_distributed_validation.R.
# Distributed bit-equivalence: split a dataset across N nodes, aggregate
# the summaries, and verify the refit matches the single-node fit.

test_that("3-node split aggregation matches single-node fit", {
  k <- gen_kron_data(N = 900L, G = 10L)

  # Single-node reference fit
  fit_full <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
                   method             = "csl",
                   re_cov             = "kronecker",
                   n_groups           = k$G,
                   Sigma_spatial_init = k$Sigma_spatial,
                   control            = vcmm_control(
                     sigma_eps       = 0.5,
                     sigma_alpha     = sqrt(0.5),
                     update_variance = TRUE))

  # Split into three random nodes
  set.seed(42L)
  node_id <- sample.int(3L, k$N, replace = TRUE)

  # Build the design (same as vcmm() would internally) and split per node
  design <- .build_design(X = k$x, t = k$t, control = vcmm_control())
  X_full <- design$X
  Z_full <- k$Z

  summaries <- lapply(seq_len(3L), function(g) {
    idx <- which(node_id == g)
    compute_sufficient_stats(
      k$y[idx], X_full[idx, , drop = FALSE], Z_full[idx, , drop = FALSE]
    )
  })

  # Aggregate via the + method
  ss_agg <- Reduce(`+`, summaries)

  # Refit from aggregate
  fit_dist <- fit_from_summaries(
    ss_list            = list(ss_agg),
    design             = design,
    re_cov             = "kronecker",
    n_groups           = k$G,
    Sigma_spatial_init = k$Sigma_spatial,
    method             = "csl",
    control            = vcmm_control(sigma_eps       = 0.5,
                                      sigma_alpha     = sqrt(0.5),
                                      update_variance = TRUE)
  )

  # Bit-equivalence: max abs diff should be well within numerical noise
  expect_lt(max(abs(fit_full$beta      - fit_dist$beta)),      1e-9)
  expect_lt(max(abs(fit_full$alpha     - fit_dist$alpha)),     1e-8)
  expect_lt(abs(fit_full$sigma_eps     - fit_dist$sigma_eps),  1e-9)
})

test_that("list-of-summaries path matches single-aggregate path", {
  d <- gen_diag_data(N = 400L, q = 3L)
  design <- .build_design(X = d$x, t = d$t, control = vcmm_control())
  X_full <- design$X

  splits <- list(1:150, 151:300, 301:400)
  summaries <- lapply(splits, function(i)
    compute_sufficient_stats(
      d$y[i], X_full[i, , drop = FALSE], d$Z[i, , drop = FALSE]
    )
  )

  fit_list <- fit_from_summaries(
    ss_list = summaries,
    design  = design,
    re_cov  = "diag",
    method  = "csl",
    control = vcmm_control(sigma_eps = 0.5,
                           sigma_alpha = 0.4,
                           update_variance = TRUE)
  )

  ss_agg <- Reduce(`+`, summaries)
  fit_agg <- fit_from_summaries(
    ss_list = list(ss_agg),
    design  = design,
    re_cov  = "diag",
    method  = "csl",
    control = vcmm_control(sigma_eps = 0.5,
                           sigma_alpha = 0.4,
                           update_variance = TRUE)
  )

  expect_equal(fit_list$beta,      fit_agg$beta,      tolerance = 1e-12)
  expect_equal(fit_list$sigma_eps, fit_agg$sigma_eps, tolerance = 1e-12)
})
