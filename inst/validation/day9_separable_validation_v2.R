#===============================================================================
# Day 9: Generalized Kronecker / Separable covariance validation
#
# Runs TWO checks in one script:
#
#   PART A -- Backward compatibility (Day 8 OD setting)
#     re_cov = "kronecker", q_left = 2 (default), G = 20, N = 3000,
#     Sigma_2x2_init / Sigma_spatial_init aliases must still work.
#     Pass criterion: relative Frobenius error of Sigma_2x2_hat < 0.08
#     (Day 8 used 0.05 at M=100; M=30 here implies ~80% more MC noise
#      on the mean, so the threshold relaxes proportionally.)
#
#   PART B -- New separable setting
#     re_cov = "separable", q_left = 5, G = 30, q_total = 150, N = 5000,
#     Z is block-dense (q_left nonzeros per row, one per "effect type"),
#     Sigma_q_true is AR(1)-style 5x5, Omega_G_true is exponential 30x30.
#     Pass criterion: relative Frobenius error of Sigma_q_hat < 0.15
#     (looser than Part A because q_left^2 = 25 entries vs 4, so the
#      per-entry MC SE accumulates over more entries in Frobenius norm).
#
# A small number of replicates (M = 30) is enough to detect a systematic
# bias; we are NOT trying to certify uncertainty bands on Sigma_q to 1%.
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

# ============================================================================
# PART A -- Backward compatibility: Day 8 OD setting
# ============================================================================
cat("\n========== PART A: backward-compat Day 8 OD setting ==========\n")

G_A              <- 20L
q_A              <- 2L * G_A
N_A              <- 3000L
M_A              <- 30L          # smaller M -- this is a regression check
sigma_eps_A      <- 0.5
beta0_A          <- 2.0
beta1_fn_A       <- function(tt) sin(2 * pi * tt)

Sigma_2x2_true_A <- matrix(c(0.50, 0.20,
                             0.20, 0.50), 2L, 2L)
Sigma_spatial_true_A <- outer(seq_len(G_A), seq_len(G_A),
                              function(i, j) exp(-abs(i - j) / 5))
Sigma_alpha_true_A   <- kronecker(Sigma_2x2_true_A, Sigma_spatial_true_A)
L_alpha_A            <- chol(Sigma_alpha_true_A)

S2_hats_A      <- array(0, dim = c(M_A, 2L, 2L))
beta0_hats_A   <- numeric(M_A)
converged_A    <- logical(M_A)

cat(sprintf("[A] %d replicates at N=%d, G=%d, q=%d  (q_left = 2, default)\n",
            M_A, N_A, G_A, q_A))
set.seed(8L)
pb_A <- utils::txtProgressBar(min = 0, max = M_A, style = 3L)

for (rep in seq_len(M_A)) {
  alpha_true <- as.vector(crossprod(L_alpha_A, rnorm(q_A)))

  origin_id <- sample.int(G_A, N_A, replace = TRUE)
  dest_id   <- sample.int(G_A, N_A, replace = TRUE)

  Z_mat <- matrix(0, N_A, q_A)
  Z_mat[cbind(seq_len(N_A), origin_id)]      <- 1
  Z_mat[cbind(seq_len(N_A), G_A + dest_id)]  <- 1

  t_vec <- runif(N_A)
  x_vec <- runif(N_A)

  y_vec <- beta0_A +
           beta1_fn_A(t_vec) * x_vec +
           as.vector(Z_mat %*% alpha_true) +
           rnorm(N_A, sd = sigma_eps_A)

  ctrl <- vcmm_control(sigma_eps       = sigma_eps_A,
                       sigma_alpha     = sqrt(0.5),
                       update_variance = TRUE,
                       max_iter        = 50L)

  ## NOTE: using the Day 8 alias names (Sigma_spatial_init) to verify
  ## backward compatibility. Day 9 also accepts Sigma_right_init.
  fit <- vcmm(y_vec,
              X                  = x_vec,
              Z                  = Z_mat,
              t                  = t_vec,
              method             = "csl",
              re_cov             = "kronecker",
              n_groups           = G_A,
              Sigma_spatial_init = Sigma_spatial_true_A,
              control            = ctrl)

  ## NOTE: still works because Sigma_2x2 is kept as a legacy alias
  ## inside re_cov_state when q_left = 2.
  S2_hats_A[rep, , ] <- fit$re_cov_state$Sigma_2x2
  beta0_hats_A[rep]  <- fit$beta[1L]
  converged_A[rep]   <- isTRUE(fit$converged)
  utils::setTxtProgressBar(pb_A, rep)
}
close(pb_A)

S2_mean_A  <- apply(S2_hats_A, c(2L, 3L), mean)
err_F_A    <- norm(S2_mean_A - Sigma_2x2_true_A, type = "F")
rel_F_A    <- err_F_A / norm(Sigma_2x2_true_A, type = "F")
rho_hat_A  <- vapply(seq_len(M_A), function(r) {
  S <- S2_hats_A[r, , ]; S[1, 2] / sqrt(S[1, 1] * S[2, 2])
}, numeric(1L))

cat("\n[A] Sigma_2x2 (true):\n");  print(round(Sigma_2x2_true_A, 4))
cat("[A] Sigma_2x2 (mean over", M_A, "reps):\n"); print(round(S2_mean_A, 4))
cat(sprintf("[A] ||err||_F = %.4f  (relative = %.4f)\n", err_F_A, rel_F_A))
cat(sprintf("[A] OD corr:  true = %+.4f,  mean = %+.4f\n",
            0.4, mean(rho_hat_A)))
cat(sprintf("[A] beta_0:   true = %.4f,  mean_hat = %.4f\n",
            beta0_A, mean(beta0_hats_A)))
cat(sprintf("[A] Convergence: %d / %d\n", sum(converged_A), M_A))

pass_A <- c(
  A_all_converged   = all(converged_A),
  A_Sigma_2x2_close = rel_F_A < 0.08,
  A_od_correlation  = abs(mean(rho_hat_A) - 0.4) < 0.05,
  A_beta_0_unbiased = {
    ## Adaptive 2-sigma test: the identifiability post-processing in vcmm.R
    ## absorbs 2 * mean(alpha) into beta_0 under indicator-Z designs, and
    ## mean(alpha) is noisy because alpha has spatially-correlated prior.
    ## At G = 20, this gives per-replicate SD of beta_0 ~ 0.80, so any
    ## fixed-tolerance test on the mean is under-powered at small M.
    mc_se_A <- sd(beta0_hats_A) / sqrt(M_A)
    cat(sprintf("[A] beta_0 per-rep SD = %.4f, MC SE on mean = %.4f, 2-sigma threshold = %.4f\n",
                sd(beta0_hats_A), mc_se_A, 2 * mc_se_A))
    abs(mean(beta0_hats_A) - beta0_A) < 2 * mc_se_A
  }
)

# ============================================================================
# PART B -- New: separable Sigma_q ⊗ Omega_G with q_left = 5
# ============================================================================
cat("\n\n========== PART B: separable Sigma_q ⊗ Omega_G (q_left = 5) ==========\n")

G_B          <- 30L
q_left_B     <- 5L
q_B          <- q_left_B * G_B    # 150
N_B          <- 5000L
M_B          <- 30L
sigma_eps_B  <- 0.5

# True Sigma_q: AR(1)-ish 5x5 with diag ~ 0.5
rho_q        <- 0.3
Sigma_q_true <- outer(seq_len(q_left_B), seq_len(q_left_B),
                      function(i, j) 0.5 * rho_q^abs(i - j))
# True Omega_G: exponential 30x30 (held FIXED in the fit)
Omega_G_true <- outer(seq_len(G_B), seq_len(G_B),
                      function(i, j) exp(-abs(i - j) / 6))
Sigma_alpha_true_B <- kronecker(Sigma_q_true, Omega_G_true)
L_alpha_B          <- chol(Sigma_alpha_true_B)

beta0_B    <- 2.0
beta1_fn_B <- function(tt) sin(2 * pi * tt)

Sq_hats_B    <- array(0, dim = c(M_B, q_left_B, q_left_B))
beta0_hats_B <- numeric(M_B)
converged_B  <- logical(M_B)

cat(sprintf("[B] %d replicates at N=%d, G=%d, q_left=%d, q_total=%d\n",
            M_B, N_B, G_B, q_left_B, q_B))
set.seed(9L)
pb_B <- utils::txtProgressBar(min = 0, max = M_B, style = 3L)

for (rep in seq_len(M_B)) {
  # Draw alpha from N(0, Sigma_q ⊗ Omega_G) under COLUMN-stacking:
  # alpha[(k-1)*G + g] = M[g, k]. Equivalent to sampling M = G x q_left.
  alpha_true <- as.vector(crossprod(L_alpha_B, rnorm(q_B)))

  # Block-dense Z: each row in some group g has q_left nonzero entries,
  # one per "effect type" k at column (k-1)*G + g. The nonzero values are
  # iid standard normal (the within-row covariate weights).
  group_id <- sample.int(G_B, N_B, replace = TRUE)
  W        <- matrix(rnorm(N_B * q_left_B), N_B, q_left_B)

  Z_mat <- matrix(0, N_B, q_B)
  for (k in seq_len(q_left_B)) {
    col_off <- (k - 1L) * G_B
    idx_lin <- cbind(seq_len(N_B), col_off + group_id)
    Z_mat[idx_lin] <- W[, k]
  }

  t_vec <- runif(N_B)
  x_vec <- runif(N_B)

  y_vec <- beta0_B +
           beta1_fn_B(t_vec) * x_vec +
           as.vector(Z_mat %*% alpha_true) +
           rnorm(N_B, sd = sigma_eps_B)

  ctrl <- vcmm_control(sigma_eps       = sigma_eps_B,
                       sigma_alpha     = sqrt(0.5),
                       update_variance = TRUE,
                       max_iter        = 50L)

  ## New API: re_cov = "separable" with q_left + Omega_G_init alias.
  fit <- vcmm(y_vec,
              X            = x_vec,
              Z            = Z_mat,
              t            = t_vec,
              method       = "csl",
              re_cov       = "separable",
              n_groups     = G_B,
              q_left       = q_left_B,
              Omega_G_init = Omega_G_true,
              control      = ctrl)

  ## Sigma_q is stored under that legacy alias for re_cov = "separable".
  Sq_hats_B[rep, , ] <- fit$re_cov_state$Sigma_q
  beta0_hats_B[rep]  <- fit$beta[1L]
  converged_B[rep]   <- isTRUE(fit$converged)
  utils::setTxtProgressBar(pb_B, rep)
}
close(pb_B)

Sq_mean_B <- apply(Sq_hats_B, c(2L, 3L), mean)
err_F_B   <- norm(Sq_mean_B - Sigma_q_true, type = "F")
rel_F_B   <- err_F_B / norm(Sigma_q_true, type = "F")

cat("\n[B] Sigma_q (true):\n");  print(round(Sigma_q_true, 4))
cat("[B] Sigma_q (mean over", M_B, "reps):\n"); print(round(Sq_mean_B, 4))
cat(sprintf("[B] ||err||_F = %.4f  (relative = %.4f)\n", err_F_B, rel_F_B))
cat(sprintf("[B] beta_0:   true = %.4f,  mean_hat = %.4f\n",
            beta0_B, mean(beta0_hats_B)))
cat(sprintf("[B] Convergence: %d / %d\n", sum(converged_B), M_B))

pass_B <- c(
  B_all_converged   = all(converged_B),
  B_Sigma_q_close   = rel_F_B < 0.15,
  B_beta_0_unbiased = {
    mc_se_B <- sd(beta0_hats_B) / sqrt(M_B)
    cat(sprintf("[B] beta_0 per-rep SD = %.4f, MC SE on mean = %.4f, 2-sigma threshold = %.4f\n",
                sd(beta0_hats_B), mc_se_B, 2 * mc_se_B))
    abs(mean(beta0_hats_B) - beta0_B) < 2 * mc_se_B
  }
)

# ============================================================================
# Combined pass / fail summary
# ============================================================================
cat("\n========== Day 9 pass / fail ==========\n")
all_pass <- c(pass_A, pass_B)
for (nm in names(all_pass)) {
  cat(sprintf("  %-22s %s\n", nm, if (all_pass[[nm]]) "PASS" else "FAIL"))
}
cat(sprintf("\nOverall: %s\n",
            if (all(all_pass)) "ALL CHECKS PASSED" else "SOME CHECKS FAILED"))

invisible(list(
  partA = list(Sigma_2x2_hats = S2_hats_A, beta0_hats = beta0_hats_A,
               converged = converged_A),
  partB = list(Sigma_q_hats = Sq_hats_B,   beta0_hats = beta0_hats_B,
               converged = converged_B)
))
