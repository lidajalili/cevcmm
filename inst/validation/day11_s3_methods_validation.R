#===============================================================================
# Day 11: S3 methods validation
#
# Smoke-tests every S3 method on the vcmm_fit object across all three
# re_cov modes. For each mode, checks:
#   - nobs(fit) returns N as an integer
#   - coef(fit) is a named numeric vector of length p
#   - fixef(fit) is a list with intercept (scalar) + varying (m x K matrix)
#   - ranef(fit) has the right shape:
#       diag      : named numeric vector of length q
#       kronecker : G x 2 matrix with origin/dest columns
#       separable : G x q_left matrix with k1..k_q_left columns
#   - varying_coef(fit, t_new) returns a (length(t_new) x K) matrix
#   - summary(fit) and print(summary(fit)) execute without error and
#     show the appropriate variance-components block per re_cov mode
#   - vcov(fit) returns the expected p x p (or other) matrix
#   - coef(fit)[1] equals fixef(fit)$intercept
#   - prediction consistency: at t_new = t (in-sample), varying_coef gives
#     the same beta_k(t) that's implicit in the fit
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

# ----- Reusable check helpers ----------------------------------------------
.expect_true <- function(cond, label) {
  if (!isTRUE(cond)) {
    stop(sprintf("[FAIL] %s", label), call. = FALSE)
  }
  cat(sprintf("  [PASS] %s\n", label))
}

.expect_equal <- function(actual, expected, tol, label) {
  diff <- max(abs(actual - expected))
  if (diff > tol) {
    stop(sprintf("[FAIL] %s  (diff = %.3e > %.3e)", label, diff, tol),
         call. = FALSE)
  }
  cat(sprintf("  [PASS] %s  (max diff = %.3e)\n", label, diff))
}

# ============================================================================
# Shared simulation harness
# ============================================================================
make_data <- function(N, G, q_left, re_cov_type, seed = 11L) {
  set.seed(seed)
  q <- if (re_cov_type == "diag") max(5L, G) else q_left * G
  t_vec <- runif(N)
  x_vec <- runif(N)
  beta0_true <- 2.0
  beta1_fn   <- function(tt) sin(2 * pi * tt)

  # Build Z based on the structure
  if (re_cov_type == "kronecker") {
    # Indicator-Z (OD)
    origin_id <- sample.int(G, N, replace = TRUE)
    dest_id   <- sample.int(G, N, replace = TRUE)
    Z <- matrix(0, N, q)
    Z[cbind(seq_len(N), origin_id)]     <- 1
    Z[cbind(seq_len(N), G + dest_id)]   <- 1
    Sigma_2x2     <- matrix(c(0.5, 0.2, 0.2, 0.5), 2L, 2L)
    Sigma_spatial <- outer(seq_len(G), seq_len(G),
                           function(i, j) exp(-abs(i - j) / 5))
    L <- chol(kronecker(Sigma_2x2, Sigma_spatial))
    alpha_true <- as.vector(crossprod(L, rnorm(q)))
  } else if (re_cov_type == "separable") {
    # Block-dense Z
    g_id <- sample.int(G, N, replace = TRUE)
    W    <- matrix(rnorm(N * q_left), N, q_left)
    Z    <- matrix(0, N, q)
    for (kk in seq_len(q_left)) {
      Z[cbind(seq_len(N), (kk - 1L) * G + g_id)] <- W[, kk]
    }
    Sigma_q <- outer(seq_len(q_left), seq_len(q_left),
                     function(i, j) 0.5 * 0.3^abs(i - j))
    Omega_G <- outer(seq_len(G), seq_len(G),
                     function(i, j) exp(-abs(i - j) / 6))
    L <- chol(kronecker(Sigma_q, Omega_G))
    alpha_true <- as.vector(crossprod(L, rnorm(q)))
  } else { # diag
    Z <- matrix(rnorm(N * q), N, q)
    alpha_true <- rnorm(q, sd = 0.4)
  }

  y <- beta0_true + beta1_fn(t_vec) * x_vec +
       as.vector(Z %*% alpha_true) + rnorm(N, sd = 0.5)

  list(y = y, x = x_vec, Z = Z, t = t_vec, q = q,
       alpha_true = alpha_true,
       Sigma_2x2 = if (re_cov_type == "kronecker") Sigma_2x2 else NULL,
       Sigma_spatial = if (re_cov_type == "kronecker") Sigma_spatial else NULL,
       Sigma_q = if (re_cov_type == "separable") Sigma_q else NULL,
       Omega_G = if (re_cov_type == "separable") Omega_G else NULL)
}

# ============================================================================
# (1) re_cov = "diag"
# ============================================================================
cat("\n========== (1) re_cov = 'diag' ==========\n")
N_d <- 500L
d   <- make_data(N = N_d, G = 5L, q_left = 1L, re_cov_type = "diag")
fit_d <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
              method  = "csl",
              re_cov  = "diag",
              control = vcmm_control(sigma_eps   = 0.5,
                                     sigma_alpha = 0.4,
                                     update_variance = TRUE))

.expect_true(inherits(fit_d, "vcmm_fit"),                "fit is vcmm_fit")
.expect_true(identical(nobs(fit_d), as.integer(N_d)),    "nobs(fit) == N")

b_d <- coef(fit_d)
.expect_true(is.numeric(b_d) && !is.null(names(b_d)),    "coef(fit) is named numeric")
.expect_true(length(b_d) == fit_d$p,                     "length(coef) == p")
.expect_true(identical(names(b_d)[1L], "(Intercept)"),   "first beta name = (Intercept)")

fe_d <- fixef(fit_d)
.expect_true(is.list(fe_d) && all(c("intercept", "varying") %in% names(fe_d)),
                                                         "fixef returns list(intercept, varying)")
.expect_equal(fe_d$intercept, unname(b_d[1L]), 1e-12,
                                                         "fixef$intercept == coef[1]")
.expect_true(is.matrix(fe_d$varying) &&
             nrow(fe_d$varying) == fit_d$design$n_basis &&
             ncol(fe_d$varying) == fit_d$design$K,
                                                         "fixef$varying shape m x K")

re_d <- ranef(fit_d)
.expect_true(is.numeric(re_d) && length(re_d) == d$q &&
             !is.matrix(re_d),
                                                         "ranef(diag) is flat named vector")

vc_d <- varying_coef(fit_d, t_new = seq(0, 1, length.out = 11L))
.expect_true(is.matrix(vc_d) && nrow(vc_d) == 11L && ncol(vc_d) == fit_d$design$K,
                                                         "varying_coef shape (length(t_new), K)")

vcov_d <- vcov(fit_d)
.expect_true(is.matrix(vcov_d) &&
             all(dim(vcov_d) == c(fit_d$p, fit_d$p)),    "vcov(fit) is p x p")

s_d <- summary(fit_d)
.expect_true(inherits(s_d, "vcmm_summary"),              "summary returns vcmm_summary")
cat("\n  [print(summary(fit)) output below]\n")
cat("  ----------------------------------------\n")
print(s_d)
cat("  ----------------------------------------\n")

# ============================================================================
# (2) re_cov = "kronecker"  (OD setting)
# ============================================================================
cat("\n========== (2) re_cov = 'kronecker' ==========\n")
N_k <- 2000L; G_k <- 20L
k <- make_data(N = N_k, G = G_k, q_left = 2L, re_cov_type = "kronecker")
fit_k <- vcmm(k$y, X = k$x, Z = k$Z, t = k$t,
              method             = "csl",
              re_cov             = "kronecker",
              n_groups           = G_k,
              Sigma_spatial_init = k$Sigma_spatial,
              control = vcmm_control(sigma_eps   = 0.5,
                                     sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))

.expect_true(identical(nobs(fit_k), as.integer(N_k)),    "nobs(fit) == N")

re_k <- ranef(fit_k)
.expect_true(is.matrix(re_k) &&
             all(dim(re_k) == c(G_k, 2L)),               "ranef(kronecker) is G x 2 matrix")
.expect_true(identical(colnames(re_k), c("origin", "dest")),
                                                         "ranef colnames = origin/dest")

s_k <- summary(fit_k)
.expect_true(!is.null(s_k$re_cov_state) &&
             identical(s_k$re_cov_state$type, "kronecker"),
                                                         "summary carries re_cov_state$type='kronecker'")
cat("\n  [print(summary(fit)) output below]\n")
cat("  ----------------------------------------\n")
print(s_k)
cat("  ----------------------------------------\n")

# ============================================================================
# (3) re_cov = "separable"
# ============================================================================
cat("\n========== (3) re_cov = 'separable' ==========\n")
N_s <- 3000L; G_s <- 15L; ql_s <- 4L
s <- make_data(N = N_s, G = G_s, q_left = ql_s, re_cov_type = "separable")
fit_s <- vcmm(s$y, X = s$x, Z = s$Z, t = s$t,
              method       = "csl",
              re_cov       = "separable",
              n_groups     = G_s,
              q_left       = ql_s,
              Omega_G_init = s$Omega_G,
              control = vcmm_control(sigma_eps   = 0.5,
                                     sigma_alpha = sqrt(0.5),
                                     update_variance = TRUE))

re_s <- ranef(fit_s)
.expect_true(is.matrix(re_s) &&
             all(dim(re_s) == c(G_s, ql_s)),             "ranef(separable) is G x q_left matrix")
.expect_true(identical(colnames(re_s), paste0("k", seq_len(ql_s))),
                                                         "ranef colnames = k1..k_qleft")

s_s <- summary(fit_s)
.expect_true(!is.null(s_s$re_cov_state) &&
             identical(s_s$re_cov_state$type, "separable"),
                                                         "summary carries re_cov_state$type='separable'")
cat("\n  [print(summary(fit)) output below]\n")
cat("  ----------------------------------------\n")
print(s_s)
cat("  ----------------------------------------\n")

# ============================================================================
# Cross-method consistency
# ============================================================================
cat("\n========== Cross-method consistency ==========\n")
# coef[1] == fixef$intercept across all three modes
for (lbl in list(c("diag", "fit_d"), c("kronecker", "fit_k"), c("separable", "fit_s"))) {
  fit <- get(lbl[2L])
  .expect_equal(coef(fit)[1L], fixef(fit)$intercept, 1e-12,
                sprintf("coef[1] == fixef$intercept (%s)", lbl[1L]))
}

# varying_coef at sample t's should equal X_design[, basis cols] %*% beta[basis]
# i.e. the same thing the fit "saw" internally.
fit <- fit_k
ds  <- fit$design
t_in <- k$t[1:5]
vc_in <- varying_coef(fit, t_new = t_in)
# Reproduce by hand: use stored basis to evaluate at t_in (normalized)
t_u  <- (t_in - ds$t_min) / (ds$t_max - ds$t_min)
B_in <- splines::bs(t_u, degree = ds$degree,
                    knots = ds$internal_knots,
                    Boundary.knots = ds$boundary_knots,
                    intercept = FALSE)
B_in <- unclass(B_in); attributes(B_in) <- list(dim = dim(B_in))
beta_basis <- coef(fit)[-1L]
vc_hand <- B_in %*% beta_basis  # K = 1 here, so single column
.expect_equal(as.numeric(vc_in), as.numeric(vc_hand), 1e-12,
              "varying_coef matches manual B %*% beta_basis")

cat("\nAll S3 method checks passed.\n")
invisible(list(fit_diag = fit_d, fit_kron = fit_k, fit_sep = fit_s))
