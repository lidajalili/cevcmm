# cevcmm Benchmark Suite

Reproducible micro- and end-to-end benchmarks for the cevcmm package.

## What's here

| File | What it measures |
|---|---|
| `run_benchmarks.R` | Single entry point. Runs all four benchmark families and writes the CSVs. |
| `bench_compute_sufficient_stats.csv` | `compute_sufficient_stats()`: R vs C++ at three problem sizes. Documents the Day 16 RcppArmadillo port. |
| `bench_invert_matrix.csv` | `invert_matrix()`: R legacy (`kappa(A) + solve(A)`) vs C++ Cholesky-first dispatch. Documents the Day 17 port. |
| `bench_rspectra.csv` | SVD pseudo-inverse: LAPACK split-merge (paper Algorithm 2) vs `RSpectra::svds()` truncated SVD on rank-deficient matrices. Documents the Day 18 opt-in path. |
| `bench_full_fit.csv` | End-to-end `vcmm()` wall-clock for OD-kron data at N ∈ {1k, 5k, 20k, 50k}. |
| `sessionInfo.txt` | Hardware and software context for the shipped numbers. |

## How to reproduce

From the package root, in an R session:

```r
source("inst/benchmarks/run_benchmarks.R")
```

Or from the shell:

```sh
Rscript inst/benchmarks/run_benchmarks.R
```

Requires `microbenchmark` (always) and `RSpectra` (optional — the third
section is skipped if absent).

```r
install.packages(c("microbenchmark", "RSpectra"))
```

Total wall-clock: roughly 5 minutes on Apple Silicon with Accelerate
BLAS. Longer on Intel + OpenBLAS or reference BLAS.

For a faster sanity check (~2 minutes, drops the largest config in each
section):

```sh
CEVCMM_BENCH_QUICK=1 Rscript inst/benchmarks/run_benchmarks.R
```

## Reading the results

### `bench_compute_sufficient_stats.csv`

Columns: `N, p, q, path, median_us, mad_us, min_us, max_us, reps`.

The function computes the six per-node cross-products (`a, b, C, ZtZ,
Zty, XtZ`) that distributed nodes transmit to the server. The C++ port
(Day 16) replaces six R `crossprod()` calls with six RcppArmadillo
matrix products. On Apple Silicon with Accelerate BLAS, both paths
dispatch to the same LAPACK kernels, so the speedup comes purely from
reduced R-side dispatch overhead — typically **~1.4×**.

This is a BLAS-bound ceiling. The speedup is real but modest. On
systems with slower LAPACK (e.g. reference BLAS) the gap may be larger.

### `bench_invert_matrix.csv`

Columns: `p, path, median_us, mad_us, min_us, max_us, reps`.

The R legacy path runs `kappa(A) + solve(A)` — and `kappa(A)` is
**itself an SVD**, used defensively to decide whether to call `solve()`.
The C++ path (Day 17) eliminates this redundant probe and uses Cholesky
directly (with LU and SVD fallbacks). Speedups are typically **3-4× at
small p** and **~2× at p ≈ 100**, because the eliminated `kappa()`
dominates at small p but amortises into the solve at large p.

### `bench_rspectra.csv`

Columns: `q, effective_rank, path, median_ms, mad_ms, min_ms, max_ms,
reps`.

The Day 18 opt-in `RSpectra` path uses Lanczos iteration to compute
only the top-k singular triplets. It wins when effective rank is much
less than q (its scaling is O(k² · q) vs LAPACK's O(q³)) and loses or
breaks even when the matrix is close to full rank.

Pattern from the reference machine: **wins big at q = 200** (rank ≈
q/5), **loses ~15% at q = 500** (Lanczos workspace bumps memory
bandwidth limits), **wins again at q = 800** (LAPACK's O(q³) regrows
faster). The path is opt-in via `options(cevcmm.use_rspectra = TRUE)`
because of this non-monotonic profile — users with known low-rank
problems benefit; users with full-rank problems should leave it off.

### `bench_full_fit.csv`

Columns: `N, G, q, re_cov, method, median_s, mad_s, min_s, max_s,
iterations, converged, sigma_eps_hat, reps`.

End-to-end `vcmm()` wall-clock for kronecker covariance (OD migration
layout). Three replicates per configuration; the inner loop is
deterministic so spread is small. `iterations` and `sigma_eps_hat` are
reported as a smoke test that the fit converged to the truth — they
should be 1 and ≈ 0.5 respectively.

## Reference numbers (Apple Silicon M-series, R 4.5, Accelerate BLAS)

Approximate medians observed during Day 16–18 validation. Your numbers
will differ; what matters is the speedup ratio.

**compute_sufficient_stats (microseconds):**

| N | R (μs) | C++ (μs) | speedup |
|---|---|---|---|
| 1,000 | 85 | 62 | 1.37× |
| 10,000 | 1,650 | 1,116 | 1.48× |
| 50,000 | 20,010 | 14,794 | 1.35× |

**invert_matrix (microseconds):**

| p | R (μs) | C++ (μs) | speedup |
|---|---|---|---|
| 20 | 37 | 10 | 3.85× |
| 50 | 96 | 48 | 2.01× |
| 80 | 183 | 104 | 1.76× |
| 99 | 270 | 169 | 1.60× |

**RSpectra (milliseconds):**

| q | rank | LAPACK (ms) | RSpectra (ms) | speedup |
|---|---|---|---|---|
| 200 | 40 | 3.6 | 0.95 | 3.80× |
| 500 | 100 | 24 | 28 | 0.85× |
| 800 | 150 | 73 | 55 | 1.35× |

## Caveats

* All numbers are sensitive to the BLAS library. Apple Silicon's
  Accelerate is heavily optimised for small matrices; numbers on
  OpenBLAS or reference BLAS look different.
* `microbenchmark` reports nanosecond timings but the underlying clock
  on macOS has microsecond resolution. Treat sub-microsecond
  differences as noise.
* The `bench_full_fit.csv` walltimes include data generation overhead
  outside the `vcmm()` call but exclude the data-prep that's already
  hoisted out of the loop.
* `sessionInfo.txt` documents what produced the shipped CSVs. If a
  reviewer asks "what hardware?", this is the answer.

## When the numbers go stale

Re-run the suite and overwrite the CSVs whenever:

* A new perf intervention lands (Day 18+ ports).
* The CRAN reference machine changes.
* The user reports significantly different numbers on their setup.

The runner is deterministic given the same random seed and BLAS, so
year-over-year comparisons make sense if those don't change.
