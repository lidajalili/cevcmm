# Day 21: multi-replicate Monte Carlo recovery of Sigma_2x2 and Sigma_q.
# Single-replicate tests in test-kron-covariance.R and test-separable.R
# only check that the fit runs; this file checks that the EM-corrected
# estimator points at the truth. Bumps coverage in covariance.R
# (.estimate_sigma_left_em, .apply_em_correction_if_kronecker) and gives
# the package its first numerical-accuracy assertion for the variance
# component.

test_that("kronecker Sigma_2x2 recovered within MC noise (5 reps)", {
  reps <- 5L
  Sigma_2x2_true <- matrix(c(0.5, 0.2, 0.2, 0.5), 2L, 2L)
  G <- 12L
  Sigma_spatial_true <- outer(seq_len(G), seq_len(G),
                              function(i, j) exp(-abs(i - j) / 4))

  Sigma_est <- array(0, c(2L, 2L, reps))
  for (r in seq_len(reps)) {
    set.seed(200L + r)
    N <- 1200L
    q <- 2L * G
    origin_id <- sample.int(G, N, replace = TRUE)
    dest_id   <- sample.int(G, N, replace = TRUE)
    Z <- matrix(0, N, q)
    Z[cbind(seq_len(N), origin_id)]   <- 1
    Z[cbind(seq_len(N), G + dest_id)] <- 1
    alpha_true <- as.vector(crossprod(
      chol(kronecker(Sigma_2x2_true, Sigma_spatial_true)),
      rnorm(q)))
    t <- runif(N); x <- runif(N)
    y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% alpha_true) +
         rnorm(N, sd = 0.5)

    fit <- vcmm(y, X = x, Z = Z, t = t,
                method             = "csl",
                re_cov             = "kronecker",
                n_groups           = G,
                Sigma_spatial_init = Sigma_spatial_true,
                control            = vcmm_control(
                  sigma_eps       = 0.5,
                  sigma_alpha     = sqrt(0.5),
                  update_variance = TRUE))
    Sigma_est[, , r] <- fit$re_cov_state$Sigma_left
  }

  mean_est <- apply(Sigma_est, c(1, 2), mean)
  expect_lt(abs(mean_est[1, 1] - 0.5), 0.15)
  expect_lt(abs(mean_est[2, 2] - 0.5), 0.15)
  expect_lt(abs(mean_est[1, 2] - 0.2), 0.15)
  # PD-ness preserved by .project_pd in every replicate
  for (r in seq_len(reps)) {
    eig <- eigen(Sigma_est[, , r], symmetric = TRUE,
                 only.values = TRUE)$values
    expect_true(all(eig > 0))
  }
})

test_that("separable Sigma_q diagonals recovered (3 reps)", {
  reps <- 3L
  G <- 10L
  q_left <- 3L
  Sigma_q_true <- outer(seq_len(q_left), seq_len(q_left),
                        function(i, j) 0.4 * 0.3^abs(i - j))
  Omega_G_true <- outer(seq_len(G), seq_len(G),
                        function(i, j) exp(-abs(i - j) / 5))

  diag_est <- matrix(0, q_left, reps)
  for (r in seq_len(reps)) {
    set.seed(300L + r)
    N <- 1500L
    q <- G * q_left
    group_id <- sample.int(G, N, replace = TRUE)
    Z <- matrix(0, N, q)
    for (k in seq_len(q_left)) {
      Z[cbind(seq_len(N), (k - 1L) * G + group_id)] <- rnorm(N)
    }
    alpha_true <- as.vector(crossprod(
      chol(kronecker(Sigma_q_true, Omega_G_true)),
      rnorm(q)))
    t <- runif(N); x <- runif(N)
    y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% alpha_true) +
         rnorm(N, sd = 0.5)

    fit <- vcmm(y, X = x, Z = Z, t = t,
                method       = "csl",
                re_cov       = "separable",
                n_groups     = G,
                q_left       = q_left,
                Omega_G_init = Omega_G_true,
                control      = vcmm_control(
                  sigma_eps       = 0.5,
                  sigma_alpha     = sqrt(0.4),
                  update_variance = TRUE))
    diag_est[, r] <- diag(fit$re_cov_state$Sigma_left)
  }

  expect_true(all(rowMeans(diag_est) > 0))
  expect_true(all(abs(rowMeans(diag_est) - 0.4) < 0.25))
})
