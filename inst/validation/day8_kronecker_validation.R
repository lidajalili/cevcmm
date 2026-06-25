#===============================================================================
# Day 8: Kronecker covariance validation (revised)
#
# Tests that vcmm() with re_cov = "kronecker" recovers Sigma_2x2 from
# realistic OD-style data when Sigma_spatial is supplied as known
# (parametric exponential kernel).
#
# Design rationale:
#   - True OD migration data has Z with TWO nonzeros per row: one for
#     the origin district, one for the destination district. We use that
#     here (the earlier version had only one nonzero per row, which is
#     not the OD design).
#   - Sigma_spatial is given as the true exponential kernel (i.e.,
#     supplied by the user as domain knowledge). The iteration estimates
#     only Sigma_2x2.
#   - We use M Monte Carlo replicates and AVERAGE the Sigma_2x2_hat
#     across them to absorb sampling variability from a small G.
#
# Setting:
#   G                   = 20 districts (n_groups)
#   q                   = 40 random effects (2 per district)
#   N                   = 3000 OD flows per replicate
#   M                   = 30 Monte Carlo replicates
#   beta_0              = 2,  beta_1(t) = sin(2 pi t)
#   sigma_eps_true      = 0.5
#   Sigma_2x2_true      = [[0.50, 0.20], [0.20, 0.50]]   (rho = 0.4)
#   Sigma_spatial_true  = exp(-|i - j| / 5)               (passed in as init)
#
# Passing criteria:
#   1. R CMD check still passes
#   2. Mean Sigma_2x2_hat (across replicates) is close to truth
#      ||mean(Sigma_2x2_hat) - Sigma_2x2_true||_F < 0.10
#   3. Mean bias of beta_0 is small: |mean(beta_0_hat) - 2| < 0.05
#   4. Convergence on every replicate
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

# ----- Truth ------------------------------------------------------------------
G              <- 20L
q              <- 2L * G
N              <- 3000L
M              <- 30L
sigma_eps_true <- 0.5
beta0_true     <- 2.0
beta1_true_fn  <- function(tt) sin(2 * pi * tt)

Sigma_2x2_true <- matrix(c(0.50, 0.20,
                           0.20, 0.50), 2L, 2L)
Sigma_spatial_true <- outer(seq_len(G), seq_len(G),
                            function(i, j) exp(-abs(i - j) / 5))
Sigma_alpha_true <- kronecker(Sigma_2x2_true, Sigma_spatial_true)
L_alpha          <- chol(Sigma_alpha_true)

# ----- Storage ----------------------------------------------------------------
Sigma_2x2_hats <- array(0, dim = c(M, 2L, 2L))
beta0_hats     <- numeric(M)
converged_vec  <- logical(M)

# ----- Monte Carlo loop -------------------------------------------------------
cat(sprintf("[Day 8] %d replicates at N=%d, G=%d, q=%d...\n", M, N, G, q))
set.seed(8L)
pb <- utils::txtProgressBar(min = 0, max = M, style = 3L)

for (rep in seq_len(M)) {
  # Draw alpha from N(0, Sigma_alpha_true)
  alpha_true <- as.vector(crossprod(L_alpha, rnorm(q)))

  # Generate OD flows: each obs has origin + destination
  origin_id <- sample.int(G, N, replace = TRUE)
  dest_id   <- sample.int(G, N, replace = TRUE)

  Z_mat <- matrix(0, N, q)
  Z_mat[cbind(seq_len(N), origin_id)]       <- 1   # origin block: cols 1..G
  Z_mat[cbind(seq_len(N), G + dest_id)]     <- 1   # dest block:   cols G+1..2G

  t_vec <- runif(N)
  x_vec <- runif(N)

  y_vec <- beta0_true +
           beta1_true_fn(t_vec) * x_vec +
           as.vector(Z_mat %*% alpha_true) +
           rnorm(N, sd = sigma_eps_true)

  ctrl <- vcmm_control(sigma_eps       = sigma_eps_true,
                       sigma_alpha     = sqrt(0.5),
                       update_variance = TRUE,
                       max_iter        = 50L)

  fit <- vcmm(y_vec,
              X                  = x_vec,
              Z                  = Z_mat,
              t                  = t_vec,
              method             = "csl",
              re_cov             = "kronecker",
              n_groups           = G,
              Sigma_spatial_init = Sigma_spatial_true,
              control            = ctrl)

  Sigma_2x2_hats[rep, , ] <- fit$re_cov_state$Sigma_2x2
  beta0_hats[rep]         <- fit$beta[1]
  converged_vec[rep]      <- isTRUE(fit$converged)
  utils::setTxtProgressBar(pb, rep)
}
close(pb)

# ----- Summary ---------------------------------------------------------------
cat("\n\n========== Recovery vs Truth ==========\n\n")

Sigma_2x2_mean <- apply(Sigma_2x2_hats, c(2L, 3L), mean)
cat("Sigma_2x2 (true):\n")
print(round(Sigma_2x2_true, 4))
cat("\nSigma_2x2 (mean across", M, "replicates):\n")
print(round(Sigma_2x2_mean, 4))

err_2x2 <- norm(Sigma_2x2_mean - Sigma_2x2_true, type = "F")
rel_2x2 <- err_2x2 / norm(Sigma_2x2_true, type = "F")
cat(sprintf("\n||mean(Sigma_2x2_hat) - Sigma_2x2_true||_F = %.4f  (relative = %.4f)\n",
            err_2x2, rel_2x2))

# OD correlation (true 0.4)
rho_hats <- vapply(seq_len(M), function(r) {
  S <- Sigma_2x2_hats[r, , ]
  S[1, 2] / sqrt(S[1, 1] * S[2, 2])
}, numeric(1))
rho_true <- 0.4
cat(sprintf("\nOD correlation:  true = %+.4f,  mean = %+.4f,  sd = %.4f\n",
            rho_true, mean(rho_hats), sd(rho_hats)))

# beta_0 recovery
cat(sprintf("\nbeta_0:   true = %.4f,  mean_hat = %.4f,  mean_bias = %+.4f\n",
            beta0_true, mean(beta0_hats), mean(beta0_hats) - beta0_true))

cat(sprintf("\nConvergence: %d / %d replicates\n", sum(converged_vec), M))

# ----- Pass / fail ------------------------------------------------------------
pass <- c(
  all_converged       = all(converged_vec),
  Sigma_2x2_close     = err_2x2 < 0.10,
  od_correlation_sign = sign(mean(rho_hats)) == sign(rho_true),
  beta_0_unbiased     = abs(mean(beta0_hats) - beta0_true) < 0.05
)

cat("\n========== Pass / fail ==========\n")
for (nm in names(pass)) {
  cat(sprintf("  %-25s %s\n", nm, if (pass[[nm]]) "PASS" else "FAIL"))
}
cat(sprintf("\nOverall: %s\n",
            if (all(pass)) "ALL CHECKS PASSED" else "SOME CHECKS FAILED"))

invisible(list(
  Sigma_2x2_hats = Sigma_2x2_hats,
  beta0_hats     = beta0_hats,
  converged      = converged_vec
))
