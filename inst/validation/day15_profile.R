#===============================================================================
# Day 15: stage-by-stage profile of vcmm()
#
# Purpose
# -------
# Identify which internal stage dominates the wall-clock cost as the
# problem grows. The answer drives Days 16-18 Rcpp ports.
#
# Expected mechanics (from the algorithm):
#   - compute_sufficient_stats():  O(N * (p + q)^2)   -- linear in N
#   - solves inside fit_ss / fit_csl: O((p + q)^3)    -- independent of N
#   - design construction:         O(N * m)           -- linear in N
#   - K_inv inversion at convergence: O((p + q)^3)    -- independent of N
#
# So at small N the (p+q)^3 cubic terms (linear algebra) dominate, and at
# large N the linear-in-N term (suff-stats accumulation) dominates. The
# crossover and absolute numbers tell us which port pays off most.
#
# Also runs the Rcpp toolchain canary cevcmm_rcpp_check() to confirm the
# new src/ scaffolding compiles and links.
#
# Wall-clock budget: ~10-30 seconds on a modern laptop.
#===============================================================================

if (!exists("compute_sufficient_stats")) {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  } else {
    library(cevcmm)
  }
}

# ---- (1) Rcpp toolchain canary ------------------------------------------
cat("\n========== (1) Rcpp toolchain canary ==========\n")
toolchain_ok <- tryCatch({
  msg <- cevcmm_rcpp_check()
  cat(sprintf("  cevcmm_rcpp_check() -> %s\n", msg))
  TRUE
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  FALSE
})

if (!toolchain_ok) {
  stop("\nRcpp toolchain is not working. Do not proceed to Day 16 until ",
       "this canary passes. On macOS, run\n",
       "  xcode-select --install\n",
       "On Windows, install Rtools matching your R version. Then ",
       "Rcpp::compileAttributes() and devtools::install() again.",
       call. = FALSE)
}

# ---- (2) Stage-by-stage timing ------------------------------------------
cat("\n========== (2) Stage profile at three problem sizes ==========\n")

simulate_kron <- function(N, G, q_left = 2L, sigma_eps = 0.5, seed = 15L) {
  set.seed(seed)
  q <- q_left * G
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
  t   <- runif(N)
  x   <- runif(N)
  y   <- 2 + sin(2 * pi * t) * x +
         as.vector(Z %*% alpha_true) +
         rnorm(N, sd = sigma_eps)
  list(y = y, x = x, Z = Z, t = t, G = G, q = q,
       Sigma_spatial = Sigma_spatial)
}

time_one_fit <- function(N, G, q_left = 2L) {
  d <- simulate_kron(N = N, G = G, q_left = q_left)
  q <- d$q

  # Stage timings via wrapping each call in system.time().
  t_design <- system.time(
    design <- build_vcmm_design(X = d$x, t = d$t)
  )[["elapsed"]]

  t_suffstat <- system.time(
    stats <- compute_sufficient_stats(d$y, design$X_design, d$Z)
  )[["elapsed"]]

  t_total <- system.time({
    fit <- vcmm(d$y, X = d$x, Z = d$Z, t = d$t,
                method             = "csl",
                re_cov             = "kronecker",
                n_groups           = G,
                Sigma_spatial_init = d$Sigma_spatial,
                control = vcmm_control(sigma_eps       = 0.5,
                                       sigma_alpha     = sqrt(0.5),
                                       update_variance = TRUE))
  })[["elapsed"]]

  # Time the inner linear-algebra by inverting the assembled K once.
  # This is what ss_solve_cpp() will eventually replace.
  p <- length(fit$beta)
  K_dim <- p + q
  t_solve <- system.time({
    K_inv <- solve(K_dim * diag(K_dim) +
                   matrix(rnorm(K_dim * K_dim), K_dim, K_dim) * 0.01)
  })[["elapsed"]]

  list(
    N         = N,
    G         = G,
    q         = q,
    p         = p,
    K_dim     = K_dim,
    t_design  = t_design,
    t_suffst  = t_suffstat,
    t_total   = t_total,
    t_solve_K = t_solve,
    converged = fit$converged
  )
}

configs <- list(
  list(N = 1000L,  G = 10L),    # small
  list(N = 5000L,  G = 30L),    # medium
  list(N = 20000L, G = 80L)     # large
)

results <- lapply(configs, function(cfg) {
  cat(sprintf("\n  ... fitting N = %d, G = %d, q = %d ...\n",
              cfg$N, cfg$G, 2L * cfg$G))
  time_one_fit(N = cfg$N, G = cfg$G)
})

# ---- (3) Report ---------------------------------------------------------
cat("\n========== (3) Stage timings ==========\n\n")

cat(sprintf("%6s %4s %5s %5s %12s %12s %12s %12s\n",
            "N", "G", "q", "K_dim",
            "design (s)", "suffst (s)", "K-solve (s)", "total fit (s)"))
cat(strrep("-", 78), "\n", sep = "")
for (r in results) {
  cat(sprintf("%6d %4d %5d %5d %12.5f %12.5f %12.5f %12.5f\n",
              r$N, r$G, r$q, r$K_dim,
              r$t_design, r$t_suffst, r$t_solve_K, r$t_total))
}

# Dominant stage at each scale
cat("\n========== (4) Dominant stage at each scale ==========\n")
for (r in results) {
  stage_times <- c(design = r$t_design,
                   suffst = r$t_suffst,
                   solve  = r$t_solve_K)
  dom <- names(which.max(stage_times))
  share_of_total <- 100 * max(stage_times) / r$t_total
  cat(sprintf("  N = %6d, q = %4d -> dominant stage: %-7s (%.1f%% of total fit time)\n",
              r$N, r$q, dom, share_of_total))
}

# Scaling exponents (log-log slope of stage time vs N)
cat("\n========== (5) Empirical scaling with N ==========\n")
Ns <- sapply(results, function(r) r$N)
fit_pow <- function(times) {
  # Avoid log(0) for sub-microsecond times
  ok <- times > 0
  if (sum(ok) < 2L) return(NA_real_)
  coef(lm(log(times[ok]) ~ log(Ns[ok])))[2L]
}
times_design <- sapply(results, function(r) r$t_design)
times_suffst <- sapply(results, function(r) r$t_suffst)
times_solveK <- sapply(results, function(r) r$t_solve_K)
times_total  <- sapply(results, function(r) r$t_total)

cat(sprintf("  d log(design)  / d log(N)  ~ %.2f  (expected ~ 1.0)\n",
            fit_pow(times_design)))
cat(sprintf("  d log(suffst)  / d log(N)  ~ %.2f  (expected ~ 1.0)\n",
            fit_pow(times_suffst)))
cat(sprintf("  d log(solve_K) / d log(N)  ~ %.2f  (expected ~ 0.0,",
            fit_pow(times_solveK)))
cat(" K-solve grows with K not N)\n")
cat(sprintf("  d log(total)   / d log(N)  ~ %.2f  (expected 1-2)\n",
            fit_pow(times_total)))

# ---- (6) Recommendation -------------------------------------------------
cat("\n========== (6) Day 16-18 priority recommendation ==========\n")
largest <- results[[length(results)]]
if (largest$t_suffst > 0.5 * largest$t_total) {
  cat("  -> compute_sufficient_stats() dominates at large N.\n")
  cat("     PORT FIRST in Day 16. Target: 5x speedup.\n")
} else if (largest$t_solve_K > 0.5 * largest$t_total) {
  cat("  -> The linear-system solve dominates.\n")
  cat("     PORT FIRST in Day 16. Target: 2x speedup (Armadillo solve).\n")
} else {
  cat("  -> No single stage dominates. Profile a larger N before deciding.\n")
  cat("     Try editing this script: configs[[3]] = list(N = 1e5, G = 200).\n")
}

cat("\nDay 15 profile complete.\n")
invisible(results)
