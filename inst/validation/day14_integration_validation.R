#===============================================================================
# Day 14: End-to-end integration test
#
# Single canonical script that exercises every public function in the
# package, in the order a user would call them:
#
#   simulate -> vcmm() -> S3 methods -> predict -> distributed -> plot
#
# Designed to be:
#   * Concise (under one screen of summary output).
#   * The kind of script a CRAN reviewer or a CV-viewer can run and
#     immediately see "yes, the package works end-to-end".
#   * Self-checking: prints PACKAGE IS HEALTHY iff every sanity gate
#     passes.
#
# Wall-clock budget: ~10 seconds.
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
# (1) Simulate -- OD migration setting, the package's headline use case
# ============================================================================
set.seed(14L)
N <- 1500L
G <- 12L
q <- 2L * G
sigma_eps_true   <- 0.5
sigma_alpha_true <- sqrt(0.5)

Sigma_2x2_true     <- matrix(c(0.50, 0.20, 0.20, 0.50), 2L, 2L)
Sigma_spatial_true <- outer(seq_len(G), seq_len(G),
                            function(i, j) exp(-abs(i - j) / 4))
Sigma_alpha_true   <- kronecker(Sigma_2x2_true, Sigma_spatial_true)
alpha_true <- as.vector(crossprod(chol(Sigma_alpha_true), rnorm(q)))

origin_id <- sample.int(G, N, replace = TRUE)
dest_id   <- sample.int(G, N, replace = TRUE)
Z <- matrix(0, N, q)
Z[cbind(seq_len(N), origin_id)]    <- 1
Z[cbind(seq_len(N), G + dest_id)]  <- 1

t   <- runif(N)
x   <- runif(N)
y   <- 2 + sin(2 * pi * t) * x +
       as.vector(Z %*% alpha_true) +
       rnorm(N, sd = sigma_eps_true)

# ============================================================================
# (2) Fit -- single call, CSL with kronecker covariance
# ============================================================================
ctrl <- vcmm_control(sigma_eps       = sigma_eps_true,
                     sigma_alpha     = sigma_alpha_true,
                     update_variance = TRUE,
                     max_iter        = 50L)

fit <- vcmm(y, X = x, Z = Z, t = t,
            method             = "csl",
            re_cov             = "kronecker",
            n_groups           = G,
            Sigma_spatial_init = Sigma_spatial_true,
            control            = ctrl)

# ============================================================================
# (3) Inspect -- every S3 method, one call each
# ============================================================================
n   <- nobs(fit)
b   <- coef(fit)
fe  <- fixef(fit)
re  <- ranef(fit)
V_b <- vcov(fit, which = "beta")
ll  <- logLik(fit)
aic <- AIC(fit)
bic <- BIC(fit)
sm  <- summary(fit)

# Varying-coefficient evaluation at a t-grid
t_grid <- seq(0, 1, length.out = 51L)
vc <- varying_coef(fit, t_new = t_grid, se.fit = TRUE)
b1_hat <- vc$fit[, 1L]
b1_se  <- vc$se.fit[, 1L]

# Recovery quality of the varying coefficient
b1_true_grid <- sin(2 * pi * t_grid)
mise_b1      <- mean((b1_hat - b1_true_grid)^2)

# ============================================================================
# (4) Predict -- both subject-specific and marginal
# ============================================================================
yhat_subject <- predict(fit, newdata = list(t = t, X = x, Z = Z))
yhat_marg    <- predict(fit, newdata = list(t = t, X = x, Z = Z),
                        include_random = FALSE)
mspe_subject <- mean((y - yhat_subject)^2)
mspe_marg    <- mean((y - yhat_marg)^2)

# Held-out test (paper's MSPE metric)
test_idx     <- sample.int(N, size = floor(0.2 * N))
fit_train    <- vcmm(y[-test_idx], X = x[-test_idx],
                     Z = Z[-test_idx, ], t = t[-test_idx],
                     method = "csl", re_cov = "kronecker",
                     n_groups = G,
                     Sigma_spatial_init = Sigma_spatial_true,
                     control = ctrl)
yhat_test    <- predict(fit_train,
                        newdata = list(t = t[test_idx],
                                       X = x[test_idx],
                                       Z = Z[test_idx, ]))
mspe_test    <- mean((y[test_idx] - yhat_test)^2)

# ============================================================================
# (5) Plot -- to PDF, so CI works headlessly
# ============================================================================
pdf_path <- "inst/validation/day14_integration_plots.pdf"
if (!dir.exists(dirname(pdf_path))) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
}
grDevices::pdf(pdf_path, width = 10, height = 5)
plot(fit, data = list(y = y, X = x, Z = Z, t = t), ask = FALSE)
grDevices::dev.off()

# ============================================================================
# (6) Distributed -- split into 3 nodes, aggregate, refit, verify bit-equiv
# ============================================================================
n_nodes  <- 3L
node_id  <- sample.int(n_nodes, N, replace = TRUE)
design   <- build_vcmm_design(X = x, t = t)
Xd       <- design$X_design

summaries <- lapply(seq_len(n_nodes), function(s) {
  idx <- which(node_id == s)
  node_summary(y[idx], Xd[idx, , drop = FALSE], Z[idx, , drop = FALSE])
})

fit_dist <- fit_from_summaries(
  summaries,
  penalty            = design$penalty,
  control            = ctrl,
  method             = "csl",
  re_cov             = "kronecker",
  n_groups           = G,
  Sigma_spatial_init = Sigma_spatial_true,
  rowsum_constant    = 2
)

dist_beta_diff  <- max(abs(fit$beta  - fit_dist$beta))
dist_alpha_diff <- max(abs(fit$alpha - fit_dist$alpha))

# ============================================================================
# (7) Sanity-check the other two re_cov modes also run end-to-end
# ============================================================================
# Diag mode
y_d   <- 2 + sin(2 * pi * t) * x + rnorm(N, sd = 0.5) +
         as.vector(matrix(rnorm(N * 5L), N, 5L) %*% rnorm(5L, sd = 0.4))
Z_d   <- matrix(rnorm(N * 5L), N, 5L)
fit_d <- vcmm(y_d, X = x, Z = Z_d, t = t, method = "csl",
              re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                                     update_variance = TRUE))
diag_ok <- isTRUE(fit_d$converged) && is.finite(logLik(fit_d))

# Separable mode
ql_s <- 3L; G_s <- 8L
g_id <- sample.int(G_s, N, replace = TRUE)
W_s  <- matrix(rnorm(N * ql_s), N, ql_s)
Z_s  <- matrix(0, N, ql_s * G_s)
for (kk in seq_len(ql_s)) {
  Z_s[cbind(seq_len(N), (kk - 1L) * G_s + g_id)] <- W_s[, kk]
}
Omega_s <- outer(seq_len(G_s), seq_len(G_s),
                 function(i, j) exp(-abs(i - j) / 3))
fit_s <- vcmm(y, X = x, Z = Z_s, t = t, method = "csl",
              re_cov       = "separable",
              n_groups     = G_s,
              q_left       = ql_s,
              Omega_G_init = Omega_s,
              control      = ctrl)
sep_ok <- isTRUE(fit_s$converged) && is.finite(logLik(fit_s))

# ============================================================================
# (8) Health summary
# ============================================================================
cat("\n=========================================================\n")
cat(" cevcmm END-TO-END INTEGRATION SUMMARY\n")
cat("=========================================================\n\n")

cat(sprintf("Data:           N = %d, G = %d, q = %d\n", n, G, q))
cat(sprintf("Truth:          sigma_eps = %.2f, beta_0 = 2.00, beta_1(t) = sin(2*pi*t)\n",
            sigma_eps_true))

cat("\n[Fit] vcmm(method = 'csl', re_cov = 'kronecker')\n")
cat(sprintf("  converged       : %s in %d iteration(s)\n",
            fit$converged, fit$iterations))
cat(sprintf("  sigma_eps_hat   : %.4f\n", fit$sigma_eps))
cat(sprintf("  beta_0_hat      : %.4f  (true 2.00)\n", fit$beta[1L]))
cat(sprintf("  MISE(beta_1(t)) : %.5f\n", mise_b1))

cat("\n[S3 methods] (every call succeeded)\n")
cat(sprintf("  nobs   = %d  | coef has %d entries  | ranef is %d x %d\n",
            n, length(b), nrow(re), ncol(re)))
cat(sprintf("  vcov(beta) : %d x %d  | logLik = %.2f (df=%d)\n",
            nrow(V_b), ncol(V_b), as.numeric(ll), attr(ll, "df")))
cat(sprintf("  AIC = %.2f  | BIC = %.2f\n", aic, bic))
cat(sprintf("  varying_coef returned curve of length %d with SE band\n",
            length(b1_hat)))

cat("\n[predict]\n")
cat(sprintf("  MSPE training (subject-specific): %.4f  (sigma_eps^2 = %.4f)\n",
            mspe_subject, fit$sigma_eps^2))
cat(sprintf("  MSPE training (marginal)        : %.4f  (ratio %.1fx)\n",
            mspe_marg, mspe_marg / mspe_subject))
cat(sprintf("  Held-out test MSPE (paper metric): %.4f\n", mspe_test))

cat("\n[Distributed]\n")
cat(sprintf("  Split into %d nodes, aggregated, refit.\n", n_nodes))
cat(sprintf("  max |beta_full - beta_dist|  : %.2e\n", dist_beta_diff))
cat(sprintf("  max |alpha_full - alpha_dist|: %.2e\n", dist_alpha_diff))

cat("\n[Cross-mode]\n")
cat(sprintf("  re_cov = 'diag'      : converged = %s, logLik finite = %s\n",
            fit_d$converged, is.finite(logLik(fit_d))))
cat(sprintf("  re_cov = 'separable' : converged = %s, logLik finite = %s\n",
            fit_s$converged, is.finite(logLik(fit_s))))

cat(sprintf("\n[plot] Wrote diagnostic PDF: %s (%.0f KB)\n",
            pdf_path, file.size(pdf_path) / 1024))

# ----- The gate ---------------------------------------------------------
gates <- c(
  kron_converged   = isTRUE(fit$converged),
  diag_converged   = diag_ok,
  sep_converged    = sep_ok,
  mspe_in_bounds   = mspe_subject < 1.3 * fit$sigma_eps^2,
  mise_in_bounds   = mise_b1 < 0.05,
  held_out_in_bounds = mspe_test < 3 * fit_train$sigma_eps^2,
  distributed_bit_equivalent = max(dist_beta_diff, dist_alpha_diff) < 1e-7,
  plot_artifact_exists = file.exists(pdf_path) && file.size(pdf_path) > 5000,
  loglik_finite    = is.finite(as.numeric(ll)),
  aic_bic_finite   = is.finite(aic) && is.finite(bic)
)

cat("\n[Gates]\n")
for (nm in names(gates)) {
  cat(sprintf("  %-30s %s\n", nm, if (gates[[nm]]) "PASS" else "FAIL"))
}

cat("\n=========================================================\n")
if (all(gates)) {
  cat(" PACKAGE IS HEALTHY\n")
} else {
  cat(" INTEGRATION ISSUES DETECTED\n")
}
cat("=========================================================\n")

invisible(list(fit = fit, fit_dist = fit_dist,
               fit_d = fit_d, fit_s = fit_s,
               gates = gates))
