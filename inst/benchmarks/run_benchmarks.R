#===============================================================================
# cevcmm Benchmark Suite (Day 19)
#
# Runs four benchmark families that document the Day 15-18 performance
# interventions and an end-to-end fit-time sweep:
#
#   1. compute_sufficient_stats() : R vs C++
#                                   (Day 16: RcppArmadillo BLAS-bound port)
#   2. invert_matrix()            : R legacy vs C++ Cholesky-first
#                                   (Day 17: Cholesky-first dispatch)
#   3. SVD pseudo-inverse         : LAPACK split-merge vs RSpectra
#                                   (Day 18: opt-in truncated SVD)
#   4. vcmm() end-to-end fit      : wall-clock across N, mode
#                                   (NEW: paper-supplement headline)
#
# Outputs:
#   * inst/benchmarks/bench_compute_sufficient_stats.csv
#   * inst/benchmarks/bench_invert_matrix.csv
#   * inst/benchmarks/bench_rspectra.csv  (skipped if RSpectra absent)
#   * inst/benchmarks/bench_full_fit.csv
#   * inst/benchmarks/sessionInfo.txt
#
# Usage from R:
#   source("inst/benchmarks/run_benchmarks.R")
# From shell:
#   Rscript inst/benchmarks/run_benchmarks.R
#
# Optional environment variables:
#   CEVCMM_BENCH_OUTDIR  = directory to write CSVs (default: inst/benchmarks)
#   CEVCMM_BENCH_QUICK   = "1" to skip the largest configs (under 2 min)
#
# Total wall-clock: ~5 min on Apple Silicon (~10 min on Intel; longer on
# OpenBLAS-only systems because Day 16 wins less without Accelerate).
#===============================================================================

# ---- (1) Setup ---------------------------------------------------------------

# Load the package. Try devtools::load_all first (development workflow),
# then fall back to library(cevcmm) (installed-package workflow).
if (file.exists("DESCRIPTION") &&
    requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all(quiet = TRUE))
} else {
  library(cevcmm)
}

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("The benchmark suite requires the microbenchmark package.\n",
       "  install.packages('microbenchmark')",
       call. = FALSE)
}

# Output directory: respect env var, else default to inst/benchmarks
.outdir <- Sys.getenv("CEVCMM_BENCH_OUTDIR", unset = "")
if (!nzchar(.outdir)) {
  .outdir <- if (dir.exists(file.path("inst", "benchmarks"))) {
    file.path("inst", "benchmarks")
  } else if (dir.exists("benchmarks")) {
    "benchmarks"
  } else {
    "."
  }
}
if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)

.quick <- identical(Sys.getenv("CEVCMM_BENCH_QUICK"), "1")

# Pretty printing helpers
.banner <- function(label) {
  bar <- strrep("-", 70L)
  cat("\n", bar, "\n", label, "\n", bar, "\n", sep = "")
}
.tick <- function(label, t0) {
  cat(sprintf("  done in %.1fs : %s\n",
              as.numeric(Sys.time() - t0, units = "secs"), label))
}

# ---- (2) Section 1 : compute_sufficient_stats --------------------------------

bench_compute_ss <- function() {
  configs <- list(
    list(N =  1000L,  p = 10L, q =  40L, reps = 200L),
    list(N = 10000L,  p = 20L, q = 100L, reps = 100L),
    list(N = 50000L,  p = 30L, q = 200L, reps =  30L)
  )
  if (.quick) configs <- configs[1:2]

  rows <- list()
  for (cfg in configs) {
    set.seed(16L + cfg$N %% 1000L)
    y <- rnorm(cfg$N)
    X <- cbind(1, matrix(rnorm(cfg$N * (cfg$p - 1L)), cfg$N, cfg$p - 1L))
    Z <- matrix(rnorm(cfg$N * cfg$q), cfg$N, cfg$q)

    mb <- microbenchmark::microbenchmark(
      R   = compute_sufficient_stats(y, X, Z, use_cpp = FALSE),
      cpp = compute_sufficient_stats(y, X, Z, use_cpp = TRUE),
      times = cfg$reps
    )

    for (lvl in c("R", "cpp")) {
      tns <- mb$time[mb$expr == lvl]
      rows[[length(rows) + 1L]] <- data.frame(
        N = cfg$N, p = cfg$p, q = cfg$q,
        path = ifelse(lvl == "cpp", "C++", "R"),
        median_us = median(tns) / 1e3,
        mad_us    = mad(tns)    / 1e3,
        min_us    = min(tns)    / 1e3,
        max_us    = max(tns)    / 1e3,
        reps      = cfg$reps,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---- (3) Section 2 : invert_matrix -------------------------------------------

bench_invert_matrix <- function() {
  make_spd <- function(p, seed = 17L) {
    set.seed(seed + p)
    M <- matrix(rnorm(p * p), p, p)
    A <- crossprod(M) / p + diag(p)
    (A + t(A)) / 2
  }

  dims <- if (.quick) c(20L, 50L) else c(20L, 50L, 80L, 99L)
  rows <- list()

  for (p in dims) {
    A <- make_spd(p)
    mb <- microbenchmark::microbenchmark(
      R   = invert_matrix(A, q = p, use_cpp = FALSE),
      cpp = invert_matrix(A, q = p, use_cpp = TRUE),
      times = 200L
    )

    for (lvl in c("R", "cpp")) {
      tns <- mb$time[mb$expr == lvl]
      rows[[length(rows) + 1L]] <- data.frame(
        p         = p,
        path      = ifelse(lvl == "cpp", "C++ Cholesky", "R legacy"),
        median_us = median(tns) / 1e3,
        mad_us    = mad(tns)    / 1e3,
        min_us    = min(tns)    / 1e3,
        max_us    = max(tns)    / 1e3,
        reps      = 200L,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---- (4) Section 3 : RSpectra (LAPACK vs truncated SVD) ----------------------

bench_rspectra <- function() {
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    return(NULL)
  }

  # Synthetic rank-deficient SPD matrix with controlled effective rank.
  make_rank_deficient_spd <- function(q, r, seed = 18L, cond_kept = 1e6) {
    set.seed(seed + q + r)
    U <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
    d <- numeric(q)
    d[seq_len(r)] <- exp(seq(0, -log(cond_kept), length.out = r))
    # Safe no-op when r == q (Day-18 footgun fix).
    d[seq_len(q - r) + r] <- d[r] * 1e-14
    A <- U %*% (d * t(U))
    (A + t(A)) / 2
  }

  configs <- list(
    list(q = 200L, r =  40L, reps = 30L),
    list(q = 500L, r = 100L, reps = 20L),
    list(q = 800L, r = 150L, reps = 10L)
  )
  if (.quick) configs <- configs[1:2]

  rows <- list()
  for (cfg in configs) {
    A <- make_rank_deficient_spd(cfg$q, cfg$r)

    mb <- microbenchmark::microbenchmark(
      lapack   = invert_matrix(A, q = cfg$q, method = "lapack"),
      rspectra = invert_matrix(A, q = cfg$q, method = "rspectra"),
      times = cfg$reps
    )

    for (lvl in c("lapack", "rspectra")) {
      tns <- mb$time[mb$expr == lvl]
      rows[[length(rows) + 1L]] <- data.frame(
        q              = cfg$q,
        effective_rank = cfg$r,
        path           = ifelse(lvl == "rspectra",
                                "RSpectra truncated",
                                "LAPACK split-merge"),
        median_ms = median(tns) / 1e6,
        mad_ms    = mad(tns)    / 1e6,
        min_ms    = min(tns)    / 1e6,
        max_ms    = max(tns)    / 1e6,
        reps      = cfg$reps,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---- (5) Section 4 : end-to-end vcmm() fit walltimes -------------------------

bench_full_fit <- function() {
  # OD-kron data generator (matches Day 14's recipe so numbers are
  # directly comparable to that script's "PACKAGE IS HEALTHY" run).
  gen_kron <- function(N, G, seed) {
    q <- 2L * G
    set.seed(seed)
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
    t <- runif(N); x <- runif(N)
    y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% alpha_true) +
         rnorm(N, sd = 0.5)
    list(y = y, X = x, Z = Z, t = t,
         Sigma_spatial_init = Sigma_spatial, G = G, q = q)
  }

  configs <- list(
    list(N =  1000L, G = 10L),
    list(N =  5000L, G = 20L),
    list(N = 20000L, G = 40L),
    list(N = 50000L, G = 80L)
  )
  if (.quick) configs <- configs[1:3]

  rows <- list()
  for (cfg in configs) {
    cat(sprintf("    fitting N=%d, G=%d ...\n", cfg$N, cfg$G))
    d <- gen_kron(cfg$N, cfg$G, seed = 19L + cfg$N %% 1000L)

    # Three replicates per config; system.time is fine here since each
    # fit is well above microbenchmark's nanosecond regime.
    reps <- 3L
    walltimes <- numeric(reps)
    fit <- NULL
    for (rep in seq_len(reps)) {
      t0 <- Sys.time()
      fit <- vcmm(d$y, X = d$X, Z = d$Z, t = d$t, method = "csl",
                  re_cov             = "kronecker",
                  n_groups           = d$G,
                  Sigma_spatial_init = d$Sigma_spatial_init,
                  control            = vcmm_control(
                    sigma_eps       = 0.5,
                    sigma_alpha     = sqrt(0.5),
                    update_variance = TRUE))
      walltimes[rep] <- as.numeric(Sys.time() - t0, units = "secs")
    }

    rows[[length(rows) + 1L]] <- data.frame(
      N             = cfg$N,
      G             = d$G,
      q             = d$q,
      re_cov        = "kronecker",
      method        = "csl",
      median_s      = median(walltimes),
      mad_s         = mad(walltimes),
      min_s         = min(walltimes),
      max_s         = max(walltimes),
      iterations    = fit$iterations,
      converged     = isTRUE(fit$converged),
      sigma_eps_hat = fit$sigma_eps,
      reps          = reps,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ---- (6) Main runner ---------------------------------------------------------

main <- function() {
  t_start <- Sys.time()
  cat("=== cevcmm benchmark suite ===\n")
  cat("Started at  : ", format(t_start), "\n", sep = "")
  cat("Output dir  : ", normalizePath(.outdir, mustWork = FALSE), "\n", sep = "")
  cat("Quick mode  : ", .quick, "\n", sep = "")

  .banner("[1/4] compute_sufficient_stats() : R vs C++")
  t0 <- Sys.time()
  res_ss <- bench_compute_ss()
  print(res_ss, row.names = FALSE)
  write.csv(res_ss,
            file.path(.outdir, "bench_compute_sufficient_stats.csv"),
            row.names = FALSE)
  .tick("wrote bench_compute_sufficient_stats.csv", t0)

  .banner("[2/4] invert_matrix() : R legacy vs C++ Cholesky")
  t0 <- Sys.time()
  res_inv <- bench_invert_matrix()
  print(res_inv, row.names = FALSE)
  write.csv(res_inv,
            file.path(.outdir, "bench_invert_matrix.csv"),
            row.names = FALSE)
  .tick("wrote bench_invert_matrix.csv", t0)

  .banner("[3/4] SVD pseudo-inverse : LAPACK split-merge vs RSpectra")
  t0 <- Sys.time()
  res_rs <- bench_rspectra()
  if (is.null(res_rs)) {
    cat("  RSpectra not installed; section skipped.\n")
  } else {
    print(res_rs, row.names = FALSE)
    write.csv(res_rs,
              file.path(.outdir, "bench_rspectra.csv"),
              row.names = FALSE)
    .tick("wrote bench_rspectra.csv", t0)
  }

  .banner("[4/4] vcmm() end-to-end wall-clock")
  t0 <- Sys.time()
  res_full <- bench_full_fit()
  print(res_full, row.names = FALSE)
  write.csv(res_full,
            file.path(.outdir, "bench_full_fit.csv"),
            row.names = FALSE)
  .tick("wrote bench_full_fit.csv", t0)

  # Capture session info so the shipped CSVs are reproducible.
  writeLines(c(
    paste0("# cevcmm benchmark sessionInfo (",
           format(Sys.time(), tz = "UTC"), " UTC)"),
    "",
    capture.output(sessionInfo()),
    "",
    "BLAS:",
    capture.output(extSoftVersion()["BLAS"])
  ), file.path(.outdir, "sessionInfo.txt"))

  total <- as.numeric(Sys.time() - t_start, units = "secs")
  cat(sprintf("\nDone in %.1f s. Wrote CSVs to %s\n",
              total, normalizePath(.outdir, mustWork = FALSE)))

  invisible(list(
    compute_ss = res_ss,
    invert     = res_inv,
    rspectra   = res_rs,
    full_fit   = res_full
  ))
}

# If sourced or Rscript'd, run main(). When dot-sourced from interactive R,
# also run main() so the user just gets the results.
main()
