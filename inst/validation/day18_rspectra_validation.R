#===============================================================================
# Day 18: validate the RSpectra truncated-SVD path in invert_matrix()
#
# Five checks:
#  (1) RSpectra availability and toolchain canary.
#  (2) ACCURACY: LAPACK SVD pseudo-inverse vs RSpectra truncated SVD on
#      a low-effective-rank matrix. Both should agree to ~1e-10 because
#      both are computing the same Moore-Penrose pseudo-inverse of the
#      same retained rank.
#  (3) SPEED: RSpectra should be faster than LAPACK when effective rank
#      is much less than q (Lanczos vs full SVD trade-off).
#  (4) FALL-BACK: full-rank matrix -> RSpectra correctly returns NULL and
#      invert_matrix() falls through to LAPACK silently.
#  (5) END-TO-END: vcmm() works with options(cevcmm.use_rspectra = TRUE).
#
# The RSpectra path is OFF by default. Users opt in via
#   options(cevcmm.use_rspectra = TRUE)
# or by passing method = "rspectra" to invert_matrix() directly.
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

# Build a q x q SPD matrix with chosen effective rank r and a controlled
# condition number for the retained part. The bottom (q - r) singular
# values are essentially zero (~1e-14 relative), giving a clean
# truncation point.
make_rank_deficient_spd <- function(q, r, seed = 18L, cond_kept = 1e6) {
  set.seed(seed + q + r)
  U <- qr.Q(qr(matrix(rnorm(q * q), q, q)))   # random orthogonal basis
  d <- numeric(q)
  d[seq_len(r)] <- exp(seq(0, -log(cond_kept), length.out = r))  # geom decay
  d[(r + 1):q] <- d[r] * 1e-14   # bottom block: numerically zero
  A <- U %*% (d * t(U))
  (A + t(A)) / 2                  # numerical symmetrization
}

# ============================================================================
# (1) Availability + toolchain canary
# ============================================================================
cat("\n========== (1) Availability ==========\n")
.expect_true(requireNamespace("RSpectra", quietly = TRUE),
             "RSpectra package is installed")
.expect_true(requireNamespace("microbenchmark", quietly = TRUE),
             "microbenchmark package is installed")
toolchain_msg <- tryCatch(cevcmm_rcpp_check(), error = function(e) "ERROR")
.expect_true(grepl("^OK", toolchain_msg),
             "C++ toolchain canary still passes")

# ============================================================================
# (2) Accuracy: LAPACK vs RSpectra on rank-deficient matrices
# ============================================================================
cat("\n========== (2) Accuracy: LAPACK vs RSpectra ==========\n")

eq_configs <- list(
  list(q = 200L, r = 40L,  tag = "q=200, rank=40"),
  list(q = 500L, r = 100L, tag = "q=500, rank=100"),
  list(q = 800L, r = 150L, tag = "q=800, rank=150")
)

for (cfg in eq_configs) {
  A <- make_rank_deficient_spd(cfg$q, cfg$r)

  inv_lapack   <- invert_matrix(A, q = cfg$q, method = "lapack")
  inv_rspectra <- invert_matrix(A, q = cfg$q, method = "rspectra")

  # Both should be valid pseudo-inverses of the same effective rank r.
  # Compare via the "Moore-Penrose identity" residual on the kept subspace:
  #   A * A^+ * A should equal A (to within numerical precision)
  resid_lapack   <- max(abs(A %*% inv_lapack   %*% A - A)) / max(abs(A))
  resid_rspectra <- max(abs(A %*% inv_rspectra %*% A - A)) / max(abs(A))

  .expect_true(
    resid_lapack < 1e-8,
    sprintf("%s: LAPACK   relative ||A A^+ A - A|| = %.3e",
            cfg$tag, resid_lapack)
  )
  .expect_true(
    resid_rspectra < 1e-8,
    sprintf("%s: RSpectra relative ||A A^+ A - A|| = %.3e",
            cfg$tag, resid_rspectra)
  )

  # Direct comparison: both Moore-Penrose pseudo-inverses on the same rank
  # should be element-wise close (the small singular vectors live in the
  # truncated subspace and cancel out).
  diff_norm <- max(abs(inv_lapack - inv_rspectra)) / max(abs(inv_lapack))
  .expect_true(
    diff_norm < 1e-6,
    sprintf("%s: max relative |LAPACK - RSpectra| = %.3e",
            cfg$tag, diff_norm)
  )
}

# ============================================================================
# (3) Speed: LAPACK vs RSpectra
# ============================================================================
cat("\n========== (3) Speed benchmark (microbenchmark) ==========\n\n")
cat(sprintf("%5s %6s %14s %16s %12s\n",
            "q", "rank", "LAPACK (ms)", "RSpectra (ms)", "speedup"))
cat(strrep("-", 70), "\n", sep = "")

speed_configs <- list(
  list(q = 200L, r = 40L,  reps = 30L),
  list(q = 500L, r = 100L, reps = 20L),
  list(q = 800L, r = 150L, reps = 10L)
)

speed_results <- list()
for (cfg in speed_configs) {
  A <- make_rank_deficient_spd(cfg$q, cfg$r)

  mb <- microbenchmark::microbenchmark(
    lapack   = invert_matrix(A, q = cfg$q, method = "lapack"),
    rspectra = invert_matrix(A, q = cfg$q, method = "rspectra"),
    times = cfg$reps
  )

  med_lap <- median(mb$time[mb$expr == "lapack"])   / 1e6   # ns -> ms
  med_rsp <- median(mb$time[mb$expr == "rspectra"]) / 1e6
  speedup <- if (med_rsp > 0) med_lap / med_rsp else NA_real_

  cat(sprintf("%5d %6d %14.2f %16.2f %11.2fx\n",
              cfg$q, cfg$r, med_lap, med_rsp, speedup))

  speed_results[[length(speed_results) + 1L]] <- list(
    q = cfg$q, r = cfg$r,
    lapack_ms = med_lap, rspectra_ms = med_rsp, speedup = speedup
  )
}

# ============================================================================
# (4) Fall-back: full-rank matrix -> RSpectra returns NULL gracefully
# ============================================================================
cat("\n========== (4) Full-rank fallback ==========\n")

# A 200x200 matrix that is genuinely full-rank (no truncation point
# below k_max). RSpectra should iterate up, fail to find a truncation
# point, return NULL, and invert_matrix() should silently fall back to
# LAPACK.
set.seed(184L)
A_full <- make_rank_deficient_spd(200L, 200L, cond_kept = 1e3)

inv_lapack <- invert_matrix(A_full, q = 200L, method = "lapack")
inv_auto   <- invert_matrix(A_full, q = 200L, method = "rspectra")

# After fallback, the two should agree to LAPACK precision.
.expect_true(
  max(abs(inv_lapack - inv_auto)) / max(abs(inv_lapack)) < 1e-8,
  sprintf("full-rank matrix: rspectra fallback matches lapack (rel diff = %.3e)",
          max(abs(inv_lapack - inv_auto)) / max(abs(inv_lapack)))
)

# ============================================================================
# (5) End-to-end: vcmm() with options(cevcmm.use_rspectra = TRUE)
# ============================================================================
cat("\n========== (5) End-to-end vcmm() with RSpectra option ==========\n")

# Build a separable (q_left = 4) problem with q_total = 100 -> SVD path
set.seed(180L)
N <- 1500L; G <- 25L; q_left <- 4L; q <- G * q_left
group_id <- sample.int(G, N, replace = TRUE)
Z <- matrix(0, N, q)
for (k in seq_len(q_left)) {
  Z[cbind(seq_len(N), (k - 1L) * G + group_id)] <- rnorm(N)
}
Omega_G <- outer(seq_len(G), seq_len(G),
                 function(i, j) exp(-abs(i - j) / 5))

t   <- runif(N); x   <- runif(N)
y   <- 2 + sin(2 * pi * t) * x + rnorm(N, sd = 0.5)

# Fit with the option ON
opts_old <- options(cevcmm.use_rspectra = TRUE)
on.exit(options(opts_old), add = TRUE)

fit_rspectra <- vcmm(y, X = x, Z = Z, t = t, method = "csl",
                     re_cov       = "separable",
                     n_groups     = G,
                     q_left       = q_left,
                     Omega_G_init = Omega_G,
                     control      = vcmm_control(sigma_eps       = 0.5,
                                                 sigma_alpha     = sqrt(0.5),
                                                 update_variance = TRUE))

# Fit with the option OFF (back to default LAPACK path)
options(cevcmm.use_rspectra = FALSE)
fit_lapack <- vcmm(y, X = x, Z = Z, t = t, method = "csl",
                   re_cov       = "separable",
                   n_groups     = G,
                   q_left       = q_left,
                   Omega_G_init = Omega_G,
                   control      = vcmm_control(sigma_eps       = 0.5,
                                               sigma_alpha     = sqrt(0.5),
                                               update_variance = TRUE))

.expect_true(isTRUE(fit_rspectra$converged),
             "vcmm() converges with cevcmm.use_rspectra = TRUE")
.expect_true(isTRUE(fit_lapack$converged),
             "vcmm() converges with cevcmm.use_rspectra = FALSE")

# beta should match between the two paths (both produce the same
# Moore-Penrose pseudo-inverse on the retained subspace)
beta_diff <- max(abs(fit_rspectra$beta - fit_lapack$beta)) /
             max(abs(fit_lapack$beta))
.expect_true(
  beta_diff < 1e-6,
  sprintf("end-to-end beta: max relative diff RSpectra vs LAPACK = %.3e",
          beta_diff)
)

cat(sprintf("\n  fit_rspectra$converged = %s, sigma_eps = %.4f\n",
            fit_rspectra$converged, fit_rspectra$sigma_eps))
cat(sprintf("  fit_lapack$converged   = %s, sigma_eps = %.4f\n",
            fit_lapack$converged, fit_lapack$sigma_eps))

cat("\nDay 18 validation passed.\n")
invisible(list(
  toolchain     = toolchain_msg,
  speed_results = speed_results,
  fit_rspectra  = fit_rspectra,
  fit_lapack    = fit_lapack
))
