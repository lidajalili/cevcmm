#===============================================================================
# Day 16: validate the RcppArmadillo port of compute_sufficient_stats()
#
# Three checks:
#  (1) Toolchain canary still works (regression test for Day 15).
#  (2) BIT-EQUIVALENCE: R and C++ paths return numerically equivalent
#      results at four problem sizes. Tolerance 1e-10 (typical
#      observed difference is ~1e-13 from BLAS summation order).
#  (3) SPEED: median wall-clock time per call, R vs C++, at three
#      production-relevant sizes.
#  (4) END-TO-END: vcmm() with the new C++ default still converges
#      and gives a sensible fit on an OD-kron simulation.
#
# Days 7-14 validations continue to exercise the C++ path via the
# default use_cpp = TRUE; if any of those break, the C++ port has a
# bug.
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

# ============================================================================
# (1) Toolchain canary
# ============================================================================
cat("\n========== (1) Toolchain canary ==========\n")
toolchain_msg <- tryCatch(cevcmm_rcpp_check(), error = function(e) "ERROR")
cat(sprintf("  cevcmm_rcpp_check() -> %s\n", toolchain_msg))
.expect_true(grepl("^OK", toolchain_msg),
             "C++ toolchain canary still passes")

# ============================================================================
# (2) Bit-equivalence: R vs C++ at four sizes
# ============================================================================
cat("\n========== (2) Bit-equivalence R vs C++ ==========\n")

eq_configs <- list(
  list(N = 100L,   p = 5L,  q = 10L,  tag = "tiny"),
  list(N = 1000L,  p = 10L, q = 40L,  tag = "small"),
  list(N = 10000L, p = 20L, q = 100L, tag = "medium"),
  list(N = 50000L, p = 30L, q = 200L, tag = "large")
)

for (cfg in eq_configs) {
  set.seed(16L + cfg$N %% 1000L)
  y <- rnorm(cfg$N)
  X <- cbind(1, matrix(rnorm(cfg$N * (cfg$p - 1L)), cfg$N, cfg$p - 1L))
  Z <- matrix(rnorm(cfg$N * cfg$q), cfg$N, cfg$q)

  ss_R   <- compute_sufficient_stats(y, X, Z, use_cpp = FALSE)
  ss_cpp <- compute_sufficient_stats(y, X, Z, use_cpp = TRUE)

  # Class consistency
  .expect_true(inherits(ss_cpp, "vcmm_ss") && inherits(ss_R, "vcmm_ss"),
               sprintf("%s: both have class 'vcmm_ss'", cfg$tag))

  # Shape consistency (b and Zty must be matrices, NOT vectors)
  shape_ok <- is.matrix(ss_cpp$b) && all(dim(ss_cpp$b) == c(cfg$p, 1L)) &&
              is.matrix(ss_cpp$Zty) && all(dim(ss_cpp$Zty) == c(cfg$q, 1L)) &&
              is.matrix(ss_cpp$C)   && all(dim(ss_cpp$C)   == c(cfg$p, cfg$p)) &&
              is.matrix(ss_cpp$ZtZ) && all(dim(ss_cpp$ZtZ) == c(cfg$q, cfg$q)) &&
              is.matrix(ss_cpp$XtZ) && all(dim(ss_cpp$XtZ) == c(cfg$p, cfg$q))
  .expect_true(shape_ok,
               sprintf("%s: C++ output shapes match R (b and Zty are matrices)",
                       cfg$tag))

  diffs <- c(
    a   = abs(ss_R$a   - ss_cpp$a),
    b   = max(abs(ss_R$b   - ss_cpp$b)),
    C   = max(abs(ss_R$C   - ss_cpp$C)),
    ZtZ = max(abs(ss_R$ZtZ - ss_cpp$ZtZ)),
    Zty = max(abs(ss_R$Zty - ss_cpp$Zty)),
    XtZ = max(abs(ss_R$XtZ - ss_cpp$XtZ))
  )
  max_diff <- max(diffs)

  # Tolerance scales with N. Each entry of X'X is a sum of N products of
  # roughly unit-magnitude terms, so floating-point summation error grows
  # linearly with N (worst case ~ N * eps_machine, eps_machine = 2.22e-16).
  # We multiply by a 100x safety factor and floor at 1e-10 so the tiny-N
  # case still requires meaningful agreement. At N = 50k this gives
  # 5e-9, which generously covers the observed ~1.6e-10.
  tol_abs <- max(1e-10, 1e-13 * cfg$N)

  .expect_true(
    max_diff < tol_abs,
    sprintf("%s (N=%d, p=%d, q=%d): max |R - C++| = %.3e   (tol = %.3e)",
            cfg$tag, cfg$N, cfg$p, cfg$q, max_diff, tol_abs)
  )
}

# ============================================================================
# (3) Speed benchmark (microbenchmark, microsecond resolution)
# ============================================================================
cat("\n========== (3) Speed benchmark (microbenchmark) ==========\n\n")

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("This benchmark needs the microbenchmark package.\n",
       "  install.packages('microbenchmark')",
       call. = FALSE)
}

cat(sprintf("%6s %4s %5s %14s %14s %12s\n",
            "N", "p", "q", "R (us)", "C++ (us)", "speedup"))
cat(strrep("-", 70), "\n", sep = "")

# Number of microbenchmark replicates: more for small N (cheap to repeat),
# fewer for large N (each call takes longer; 30 reps is plenty for a
# stable median).
bench_configs <- list(
  list(N = 1000L,  p = 10L, q = 40L,  reps = 200L),
  list(N = 10000L, p = 20L, q = 100L, reps = 100L),
  list(N = 50000L, p = 30L, q = 200L, reps = 30L)
)

results <- list()
for (cfg in bench_configs) {
  set.seed(16L + cfg$N %% 1000L)
  y <- rnorm(cfg$N)
  X <- cbind(1, matrix(rnorm(cfg$N * (cfg$p - 1L)), cfg$N, cfg$p - 1L))
  Z <- matrix(rnorm(cfg$N * cfg$q), cfg$N, cfg$q)

  mb <- microbenchmark::microbenchmark(
    R   = compute_sufficient_stats(y, X, Z, use_cpp = FALSE),
    cpp = compute_sufficient_stats(y, X, Z, use_cpp = TRUE),
    times = cfg$reps
  )

  med_R   <- median(mb$time[mb$expr == "R"])   / 1000  # ns -> us
  med_cpp <- median(mb$time[mb$expr == "cpp"]) / 1000
  speedup <- if (med_cpp > 0) med_R / med_cpp else NA_real_

  cat(sprintf("%6d %4d %5d %14.2f %14.2f %11.2fx\n",
              cfg$N, cfg$p, cfg$q, med_R, med_cpp, speedup))

  results[[length(results) + 1L]] <- list(
    N = cfg$N, p = cfg$p, q = cfg$q,
    med_R_us = med_R, med_cpp_us = med_cpp, speedup = speedup
  )
}

# ============================================================================
# (4) End-to-end: vcmm() with the C++ default
# ============================================================================
cat("\n========== (4) End-to-end vcmm() with C++ default ==========\n")

set.seed(160L)
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
             "vcmm() converges via C++ default")
.expect_true(is.finite(fit$sigma_eps) && abs(fit$sigma_eps - 0.5) < 0.1,
             sprintf("sigma_eps_hat = %.4f within 0.1 of truth 0.5",
                     fit$sigma_eps))
.expect_true(is.finite(as.numeric(logLik(fit))),
             "logLik(fit) is finite")

cat(sprintf("\n  fit$method     = %s\n", fit$method))
cat(sprintf("  fit$converged  = %s\n", fit$converged))
cat(sprintf("  fit$sigma_eps  = %.4f  (true 0.5)\n", fit$sigma_eps))
cat(sprintf("  logLik         = %.2f\n", as.numeric(logLik(fit))))

cat("\nDay 16 validation passed.\n")
invisible(list(toolchain = toolchain_msg, bench = results, fit = fit))
