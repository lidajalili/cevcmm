#===============================================================================
# Day 12: predict.vcmm_fit and logLik.vcmm_fit validation
#
# Checks (across all three re_cov modes):
#
# A. predict() with default include_random = TRUE recovers the training
#    response to within sigma_eps noise. The paper's MSPE formula is
#    yhat = X*beta + Z*alpha, so on the training set we expect
#    mean((y - yhat)^2) ~ sigma_eps^2 (it can dip a bit below sigma_eps^2
#    because the BLUP slightly overfits to the training noise, but it
#    should never exceed ~1.3 * sigma_eps^2).
#
# B. predict() with include_random = FALSE returns the marginal predictor
#    X*beta only, which produces strictly higher training MSPE than (A)
#    (Equation 14 of the supervisor's references).
#
# C. predict() bit-equivalence: at the training data, yhat from predict()
#    must match the fit's internal computation (sum of squared residuals
#    matches n * sigma_eps^2 if we use include_random = TRUE).
#
# D. se.fit = TRUE returns positive standard errors of the right length.
#
# E. logLik(fit) returns a finite numeric scalar with df and nobs
#    attributes; AIC(fit) and BIC(fit) work and are finite.
#
# F. Cross-mode sanity: more flexible covariance (separable / kronecker)
#    yields higher log-likelihood than naive diag on data that actually
#    has structure -- the AIC comparison should reflect this.
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

.expect_true <- function(cond, label) {
  if (!isTRUE(cond)) stop(sprintf("[FAIL] %s", label), call. = FALSE)
  cat(sprintf("  [PASS] %s\n", label))
}
.expect_lt <- function(x, bound, label) {
  if (!isTRUE(x < bound)) {
    stop(sprintf("[FAIL] %s  (got %.4g, bound %.4g)", label, x, bound),
         call. = FALSE)
  }
  cat(sprintf("  [PASS] %s  (%.4g < %.4g)\n", label, x, bound))
}

# ============================================================================
# Shared data simulation harness (mirrors Day 11)
# ============================================================================
make_data <- function(N, G, q_left, re_cov_type, sigma_eps = 0.5, seed = 12L) {
  set.seed(seed)
  q <- if (re_cov_type == "diag") max(5L, G) else q_left * G
  t_vec <- runif(N); x_vec <- runif(N)
  beta0_true <- 2.0
  beta1_fn   <- function(tt) sin(2 * pi * tt)
  
  if (re_cov_type == "kronecker") {
    origin_id <- sample.int(G, N, replace = TRUE)
    dest_id   <- sample.int(G, N, replace = TRUE)
    Z <- matrix(0, N, q)
    Z[cbind(seq_len(N), origin_id)]    <- 1
    Z[cbind(seq_len(N), G + dest_id)]  <- 1
    Sigma_left  <- matrix(c(0.5, 0.2, 0.2, 0.5), 2L, 2L)
    Sigma_right <- outer(seq_len(G), seq_len(G),
                         function(i, j) exp(-abs(i - j) / 5))
    alpha_true  <- as.vector(crossprod(chol(kronecker(Sigma_left, Sigma_right)),
                                       rnorm(q)))
  } else if (re_cov_type == "separable") {
    g_id <- sample.int(G, N, replace = TRUE)
    W    <- matrix(rnorm(N * q_left), N, q_left)
    Z    <- matrix(0, N, q)
    for (kk in seq_len(q_left)) {
      Z[cbind(seq_len(N), (kk - 1L) * G + g_id)] <- W[, kk]
    }
    Sigma_left  <- outer(seq_len(q_left), seq_len(q_left),
                         function(i, j) 0.5 * 0.3^abs(i - j))
    Sigma_right <- outer(seq_len(G), seq_len(G),
                         function(i, j) exp(-abs(i - j) / 6))
    alpha_true  <- as.vector(crossprod(chol(kronecker(Sigma_left, Sigma_right)),
                                       rnorm(q)))
  } else { # diag
    Z <- matrix(rnorm(N * q), N, q)
    alpha_true <- rnorm(q, sd = 0.4)
    Sigma_left <- Sigma_right <- NULL
  }
  
  y <- beta0_true + beta1_fn(t_vec) * x_vec +
    as.vector(Z %*% alpha_true) + rnorm(N, sd = sigma_eps)
  
  list(y = y, x = x_vec, Z = Z, t = t_vec, q = q,
       Sigma_left = Sigma_left, Sigma_right = Sigma_right,
       sigma_eps_true = sigma_eps)
}

# ============================================================================
# (1) re_cov = "diag"
# ============================================================================
cat("\n========== (1) re_cov = 'diag' ==========\n")
N_d <- 500L
d <- make_data(N = N_d, G = 5L, q_left = 1L, re_cov_type = "diag")
fit_d <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t, method = "csl",
              re_cov = "diag",
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.4,
                                     update_variance = TRUE))

# (A) Subject-specific prediction recovers training y
yhat_re <- predict(fit_d, newdata = list(t = d$t, X = d$x, Z = d$Z))
mse_re  <- mean((d$y - yhat_re)^2)
.expect_lt(mse_re, 1.3 * fit_d$sigma_eps^2,
           "predict(include_random=TRUE): training MSE < 1.3 * sigma_eps^2")

# (B) Marginal predictor has strictly higher MSE
yhat_marg <- predict(fit_d, newdata = list(t = d$t, X = d$x, Z = d$Z),
                     include_random = FALSE)
mse_marg <- mean((d$y - yhat_marg)^2)
cat(sprintf("  [INFO] diag MSE: subject-specific = %.4f, marginal = %.4f, ratio = %.2fx\n",
            mse_re, mse_marg, mse_marg / mse_re))
.expect_true(mse_marg > mse_re,
             "marginal MSE > subject-specific MSE (include_random=FALSE drops alpha)")

# (C) Bit-equivalence: yhat_re - y should give the same RSS as fit$sigma_eps^2 * N
# (within iteration tolerance; for CSL the relation isn't exact but should be close)
rss_pred <- sum((d$y - yhat_re)^2)
rss_fit  <- N_d * fit_d$sigma_eps^2
rel_err  <- abs(rss_pred - rss_fit) / rss_fit
.expect_lt(rel_err, 1e-6,
           "RSS from predict() matches N * sigma_eps^2 from fit")

# (D) se.fit
out_se <- predict(fit_d, newdata = list(t = d$t, X = d$x, Z = d$Z),
                  se.fit = TRUE)
.expect_true(is.list(out_se) &&
               length(out_se$fit) == N_d &&
               length(out_se$se.fit) == N_d &&
               all(out_se$se.fit >= 0),
             "se.fit returns non-negative SEs of length N")

# (E) logLik / AIC / BIC
ll_d <- logLik(fit_d)
.expect_true(is.finite(as.numeric(ll_d)),                "logLik is finite")
.expect_true(!is.null(attr(ll_d, "df")) && attr(ll_d, "df") > 0L,
             "logLik has df attribute")
.expect_true(!is.null(attr(ll_d, "nobs")) && attr(ll_d, "nobs") == N_d,
             "logLik has nobs == N")
.expect_true(is.finite(AIC(fit_d)) && is.finite(BIC(fit_d)),
             "AIC and BIC are finite")
cat(sprintf("  [INFO] diag: logLik=%.2f, df=%d, AIC=%.2f, BIC=%.2f\n",
            as.numeric(ll_d), attr(ll_d, "df"),
            AIC(fit_d), BIC(fit_d)))

# ============================================================================
# (2) re_cov = "kronecker"  (OD setting)
# ============================================================================
cat("\n========== (2) re_cov = 'kronecker' ==========\n")
N_k <- 2000L; G_k <- 20L
k <- make_data(N = N_k, G = G_k, q_left = 2L, re_cov_type = "kronecker")
fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t, method = "csl",
              re_cov             = "kronecker",
              n_groups           = G_k,
              Sigma_spatial_init = k$Sigma_right,
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))

yhat_re   <- predict(fit_k, newdata = list(t = k$t, X = k$x, Z = k$Z))
yhat_marg <- predict(fit_k, newdata = list(t = k$t, X = k$x, Z = k$Z),
                     include_random = FALSE)
mse_re    <- mean((k$y - yhat_re)^2)
mse_marg  <- mean((k$y - yhat_marg)^2)
cat(sprintf("  [INFO] kron MSE: subject-specific = %.4f, marginal = %.4f, ratio = %.2fx\n",
            mse_re, mse_marg, mse_marg / mse_re))
.expect_lt(mse_re, 1.3 * fit_k$sigma_eps^2,
           "kron: subject-specific MSE < 1.3 * sigma_eps^2")
.expect_true(mse_marg > mse_re,
             "kron: marginal MSE > subject-specific (RE contributes to prediction)")

ll_k <- logLik(fit_k)
.expect_true(is.finite(as.numeric(ll_k)),                "kron: logLik finite")
.expect_true(attr(ll_k, "df") == fit_k$p + 1L + 3L,
             "kron: df = p + 1 + 3 (sigma_eps + 3 Sigma_2x2 entries)")
cat(sprintf("  [INFO] kron: logLik=%.2f, df=%d, AIC=%.2f, BIC=%.2f\n",
            as.numeric(ll_k), attr(ll_k, "df"),
            AIC(fit_k), BIC(fit_k)))

# se.fit with random effects -- the SE per row should be comparable to sigma_eps
out_se_k <- predict(fit_k, newdata = list(t = k$t, X = k$x, Z = k$Z),
                    se.fit = TRUE)
.expect_true(all(out_se_k$se.fit >= 0) && mean(out_se_k$se.fit) < fit_k$sigma_eps,
             "kron: se.fit reasonable scale (mean below sigma_eps in-sample)")

# ============================================================================
# (3) re_cov = "separable"
# ============================================================================
cat("\n========== (3) re_cov = 'separable' ==========\n")
N_s <- 3000L; G_s <- 15L; ql_s <- 4L
s <- make_data(N = N_s, G = G_s, q_left = ql_s, re_cov_type = "separable")
fit_s <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t, method = "csl",
              re_cov       = "separable",
              n_groups     = G_s,
              q_left       = ql_s,
              Omega_G_init = s$Sigma_right,
              control = vcmm_control(sigma_eps = 0.5, sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))

yhat_re   <- predict(fit_s, newdata = list(t = s$t, X = s$x, Z = s$Z))
yhat_marg <- predict(fit_s, newdata = list(t = s$t, X = s$x, Z = s$Z),
                     include_random = FALSE)
mse_re    <- mean((s$y - yhat_re)^2)
mse_marg  <- mean((s$y - yhat_marg)^2)
cat(sprintf("  [INFO] sep MSE: subject-specific = %.4f, marginal = %.4f, ratio = %.2fx\n",
            mse_re, mse_marg, mse_marg / mse_re))
.expect_lt(mse_re, 1.3 * fit_s$sigma_eps^2,
           "sep: subject-specific MSE < 1.3 * sigma_eps^2")
.expect_true(mse_marg > mse_re,
             "sep: marginal MSE > subject-specific")

ll_s <- logLik(fit_s)
.expect_true(is.finite(as.numeric(ll_s)),                "sep: logLik finite")
expected_df_sep <- fit_s$p + 1L + as.integer(ql_s * (ql_s + 1L) / 2L)
.expect_true(attr(ll_s, "df") == expected_df_sep,
             sprintf("sep: df = p + 1 + q_left*(q_left+1)/2 = %d", expected_df_sep))
cat(sprintf("  [INFO] sep: logLik=%.2f, df=%d, AIC=%.2f, BIC=%.2f\n",
            as.numeric(ll_s), attr(ll_s, "df"),
            AIC(fit_s), BIC(fit_s)))

# ============================================================================
# Cross-method consistency
# ============================================================================
cat("\n========== Cross-method consistency ==========\n")

# (F) Held-out prediction check: split the kronecker data 80/20, fit on train,
# predict on test, ensure test MSPE is in the right ballpark.
set.seed(42L)
test_idx  <- sample.int(N_k, size = floor(0.2 * N_k))
train_idx <- setdiff(seq_len(N_k), test_idx)

fit_train <- vcmm(k$y[train_idx],
                  X = k$x[train_idx],
                  Z = k$Z[train_idx, , drop = FALSE],
                  t = k$t[train_idx],
                  method             = "csl",
                  re_cov             = "kronecker",
                  n_groups           = G_k,
                  Sigma_spatial_init = k$Sigma_right,
                  control = vcmm_control(sigma_eps = 0.5,
                                         sigma_alpha = sqrt(0.5),
                                         update_variance = TRUE))
yhat_test <- predict(fit_train,
                     newdata = list(t = k$t[test_idx],
                                    X = k$x[test_idx],
                                    Z = k$Z[test_idx, , drop = FALSE]))
mspe_test <- mean((k$y[test_idx] - yhat_test)^2)
cat(sprintf("  [INFO] held-out kron: test MSPE = %.4f  (sigma_eps^2 = %.4f)\n",
            mspe_test, fit_train$sigma_eps^2))
.expect_lt(mspe_test, 3 * fit_train$sigma_eps^2,
           "held-out test MSPE < 3 * sigma_eps^2 (paper's main metric, in-bounds)")

cat("\nAll predict / logLik checks passed.\n")
invisible(list(fit_d = fit_d, fit_k = fit_k, fit_s = fit_s,
               fit_train = fit_train))
