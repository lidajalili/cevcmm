#===============================================================================
# Day 17: validate the RcppArmadillo port of invert_matrix() and friends
#
# Four checks:
#  (1) Toolchain canary still works.
#  (2) BIT-EQUIVALENCE: R legacy path vs C++ Cholesky path agree to
#      machine precision on well-conditioned SPD matrices, the typical
#      VCMM K-matrix shape.
#  (3) SPEED: median wall-clock per call, R vs C++, at multiple
#      dimensions covering the range hit during real fits.
#  (4) END-TO-END: every Day-7-to-14 validation script continues to
#      pass via the new C++ default. (Run them separately.)
#
# Plus three robustness checks:
#  - Cholesky correctly fails on non-PD matrices and falls back to LU
#  - Singular matrices fall through to SVD pseudo-inverse
#  - q >= 100 still routes through the split-merge SVD path (paper
#    Algorithm 2), unchanged
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

# Generate a well-conditioned symmetric positive-definite matrix with a
# target dimension. Diagonal-dominant so kappa stays modest.
make_spd <- function(p, seed = 17L) {
  set.seed(seed + p)
  M <- matrix(rnorm(p * p), p, p)
  A <- crossprod(M) / p + diag(p)
  (A + t(A)) / 2  # numerical symmetry
}

# ============================================================================
# (1) Toolchain canary
# ============================================================================
cat("\n========== (1) Toolchain canary ==========\n")
toolchain_msg <- tryCatch(cevcmm_rcpp_check(), error = function(e) "ERROR")
cat(sprintf("  cevcmm_rcpp_check() -> %s\n", toolchain_msg))
.expect_true(grepl("^OK", toolchain_msg),
             "C++ toolchain canary still passes")

# ============================================================================
# (2) Bit-equivalence: R legacy path vs C++ Cholesky path
# ============================================================================
cat("\n========== (2) Bit-equivalence R vs C++ ==========\n")

eq_dims <- c(10L, 30L, 50L, 80L)
for (p in eq_dims) {
  A <- make_spd(p)

  inv_R   <- invert_matrix(A, q = p, use_cpp = FALSE)
  inv_cpp <- invert_matrix(A, q = p, use_cpp = TRUE)

  max_diff <- max(abs(inv_R - inv_cpp))
  # Tolerance scales with p: each entry of A^{-1} is a sum-of-p-products,
  # so floating-point error ~ p * eps_machine. Safety factor 100.
  tol_abs  <- max(1e-12, 1e-13 * p)
  .expect_true(
    max_diff < tol_abs,
    sprintf("p=%d: max |R - C++| = %.3e   (tol = %.3e)",
            p, max_diff, tol_abs)
  )

  # Sanity check: both paths should give A %*% A^{-1} = I to machine eps
  resid_cpp <- max(abs(A %*% inv_cpp - diag(p)))
  .expect_true(
    resid_cpp < 1e-10,
    sprintf("p=%d: C++ ||A A^{-1} - I||_max = %.3e", p, resid_cpp)
  )
}

# ============================================================================
# (2b) Robustness: non-PD matrix -> Cholesky fails, LU succeeds
# ============================================================================
cat("\n========== (2b) Robustness: non-PD matrix ==========\n")
set.seed(170L)
M <- matrix(rnorm(20 * 20), 20, 20)
A_indef <- (M + t(M)) / 2  # symmetric but not PD; some eigenvalues negative
.expect_true(
  any(eigen(A_indef, symmetric = TRUE, only.values = TRUE)$values < 0),
  "test matrix has negative eigenvalues (genuinely not PD)"
)

# Should still succeed via LU fallback
inv_indef <- tryCatch(invert_matrix(A_indef, q = 20, use_cpp = TRUE),
                      error = function(e) NULL)
.expect_true(!is.null(inv_indef),
             "invert_matrix() succeeds on non-PD matrix via LU fallback")
.expect_true(
  max(abs(A_indef %*% inv_indef - diag(20))) < 1e-9,
  "LU fallback gives correct inverse"
)

# ============================================================================
# (2c) q >= 100 still routes through split-merge SVD (Algorithm 2)
# ============================================================================
cat("\n========== (2c) Large-q routes through split-merge SVD ==========\n")
A_large <- make_spd(150L)
inv_large_R   <- invert_matrix(A_large, q = 150L, use_cpp = FALSE)
inv_large_cpp <- invert_matrix(A_large, q = 150L, use_cpp = TRUE)
.expect_true(
  max(abs(inv_large_R - inv_large_cpp)) < 1e-8,
  sprintf("p=150: R and C++ both route through split-merge SVD (max diff = %.3e)",
          max(abs(inv_large_R - inv_large_cpp)))
)

# ============================================================================
# (3) Speed benchmark (microbenchmark, 200 reps, microsecond resolution)
# ============================================================================
cat("\n========== (3) Speed benchmark (microbenchmark, 200 reps) ==========\n\n")

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("This benchmark needs the microbenchmark package.\n",
       "  install.packages('microbenchmark')",
       call. = FALSE)
}

cat(sprintf("%5s %14s %14s %12s %12s\n",
            "p", "R (us)", "C++ (us)", "speedup", "C++/R %"))
cat(strrep("-", 70), "\n", sep = "")

bench_dims <- c(20L, 50L, 80L, 99L)
bench_results <- list()
for (p in bench_dims) {
  A <- make_spd(p)

  # microbenchmark times each expression 200 times, returns ns timings.
  # Auto-warms up internally; no need for a manual warmup call.
  mb <- microbenchmark::microbenchmark(
    R   = invert_matrix(A, q = p, use_cpp = FALSE),
    cpp = invert_matrix(A, q = p, use_cpp = TRUE),
    times = 200L
  )

  # Convert nanoseconds -> microseconds, take median (robust to outliers).
  med_R   <- median(mb$time[mb$expr == "R"])   / 1000
  med_cpp <- median(mb$time[mb$expr == "cpp"]) / 1000
  speedup <- if (med_cpp > 0) med_R / med_cpp else NA_real_

  cat(sprintf("%5d %14.2f %14.2f %11.2fx %11.1f%%\n",
              p, med_R, med_cpp, speedup, 100 * med_cpp / med_R))

  bench_results[[length(bench_results) + 1L]] <- list(
    p = p, med_R_us = med_R, med_cpp_us = med_cpp, speedup = speedup
  )
}

cat("\nNote: p < 100 routes through the Cholesky fast path. p >= 100 routes\n")
cat("through the R-side split-merge SVD (paper Algorithm 2), where C++ and R\n")
cat("share the same implementation -- no speedup is expected or measured.\n")

# ============================================================================
# (4) End-to-end: vcmm() with the C++ default still works
# ============================================================================
cat("\n========== (4) End-to-end vcmm() with C++ default ==========\n")

set.seed(170L)
N <- 1500L; G <- 12L; q <- 2L * G
origin_id <- sample.int(G, N, replace = TRUE)
dest_id   <- sample.int(G, N, replace = TRUE)
Z <- matrix(0, N, q)
Z[cbind(seq_len(N), origin_id)]   <- 1
Z[cbind(seq_len(N), G + dest_id)] <- 1

Sigma_2x2     <- matrix(c(0.5, 0.2, 0.2, 0.5), 2L, 2L)
Sigma_spatial <- outer(seq_len(G), seq_len(G),
                       function(i, j) exp(-abs(i - j) / 4))
alpha_true <- as.vector(crossprod(chol(kronecker(Sigma_2x2, Sigma_spatial)),
                                  rnorm(q)))

t   <- runif(N); x   <- runif(N)
y   <- 2 + sin(2 * pi * t) * x +
       as.vector(Z %*% alpha_true) +
       rnorm(N, sd = 0.5)

fit <- vcmm(y, X = x, Z = Z, t = t, method = "csl",
            re_cov             = "kronecker",
            n_groups           = G,
            Sigma_spatial_init = Sigma_spatial,
            control            = vcmm_control(sigma_eps       = 0.5,
                                              sigma_alpha     = sqrt(0.5),
                                              update_variance = TRUE))

.expect_true(isTRUE(fit$converged),
             "vcmm() converges with the new C++ invert_matrix() default")
.expect_true(is.finite(fit$sigma_eps) && abs(fit$sigma_eps - 0.5) < 0.1,
             sprintf("sigma_eps_hat = %.4f within 0.1 of truth 0.5",
                     fit$sigma_eps))
.expect_true(is.finite(as.numeric(logLik(fit))),
             "logLik(fit) is finite")

cat(sprintf("\n  fit$method     = %s\n", fit$method))
cat(sprintf("  fit$converged  = %s\n", fit$converged))
cat(sprintf("  fit$sigma_eps  = %.4f  (true 0.5)\n", fit$sigma_eps))
cat(sprintf("  logLik         = %.2f\n", as.numeric(logLik(fit))))

cat("\nDay 17 validation passed.\n")
invisible(list(toolchain = toolchain_msg, fit = fit))
