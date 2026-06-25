# cevcmm — Development Roadmap

Target: **v0.1.0 release on GitHub + CRAN submission in ~4 weeks** (28 working days @ 4–6 hrs/day).

Mark each day with ✅ when complete, 🚧 when in progress, ⬜ when not started.
Commit this file at the end of every working day.

---

## Week 1 — Foundation & core refactor

| Day | Status | Deliverable | Done when |
|----:|:------:|-------------|-----------|
| 1   |  ✅    | Package skeleton: `DESCRIPTION`, `NAMESPACE`, `.Rbuildignore`, `.gitignore`, `R/cevcmm-package.R`, git init, GitHub repo | `R CMD check` runs with 0 errors |
| 2   | ✅     | `R/sufficient_stats.R`, `R/penalty.R`, `R/svd_stable.R` refactored from current modules with roxygen docs | All three load via `devtools::load_all()` |
| 3   | ✅     | `R/ss_estimator.R` refactored; `vcmm_control()` options object | SS gives a numeric result on tiny data |
| 4   | ✅     | `R/csl_estimator.R` refactored; pilot logic clean | CSL gives same first-order behavior as SS |
| 5   | 🚧     | `R/design.R` — explicit-spec interface `vcmm(y, X, Z, t, group)` (formula sugar deferred to v0.2) | Three test specs run end-to-end |
| 6   | ⬜     | `R/vcmm.R` — main wrapper with `auto` method selection (rule: `csl` if N×q > 1e7 or q > 100; else `ss`) | `vcmm(...)` runs end-to-end |
| 7   | ⬜     | Smoke test — replicate the SS validation setting (N=1000, q=1) | Estimates within MC noise of paper Table A1 |

## Week 2 — Variants, S3 methods, distributed API

| Day | Status | Deliverable | Done when |
|----:|:------:|-------------|-----------|
| 8   | ⬜     | `R/covariance.R` — Kronecker Σ₂ₓ₂ ⊗ Σ_spatial | OD simulation runs through package |
| 9   | ⬜     | Group-shared dense Σ_q ⊗ Ω_G via `re_cov="separable"` | Dense simulation runs through package |
| 10  | ⬜     | `R/distributed.R` — `node_summary()` + `fit_from_summaries()` | 3-node mock fit = single-node fit |
| 11  | ⬜     | S3: `print`, `summary`, `coef`, `fixef`, `ranef` | All five give sensible output |
| 12  | ⬜     | S3: `predict.vcmm` (with newdata), `vcov.vcmm`, `logLik.vcmm` | Predictions match held-out MSPE |
| 13  | ⬜     | `plot.vcmm` — varying-coef CI bands, residuals, RE QQ | Three plot types render |
| 14  | ⬜     | Integration day — full end-to-end on OD and dense | Both match current code numerically |

## Week 3 — Rcpp + tests

| Day | Status | Deliverable | Done when |
|----:|:------:|-------------|-----------|
| 15  | ⬜     | `src/Makevars`, profile to confirm hot paths | Profile shows accumulation + solve > 70% |
| 16  | ⬜     | Rcpp `accumulate_ss_cpp()` with RcppArmadillo | ≥5× faster on N=500k; output identical to 1e-12 |
| 17  | ⬜     | Rcpp `ss_solve_cpp()` iteration kernel | ≥2× faster; identical iterates |
| 18  | ⬜     | RSpectra truncated-SVD path; auto-route by q | SVD path runs on q=300 ill-cond test |
| 19  | ⬜     | Benchmark Rcpp vs R; document in `inst/benchmarks/` | Speedup table written |
| 20  | ⬜     | `tests/testthat/test-ss-recovers-MLE.R`, `test-csl-onestep-equivalence.R` | Tests pass; coverage ≥70% on core |
| 21  | ⬜     | `test-svd-stability.R`, `test-kronecker-recovery.R`, `test-distributed-equals-pooled.R`, `test-spec-parser.R` | Coverage ≥85% overall |

## Week 4 — Docs, vignettes, CRAN prep

| Day | Status | Deliverable | Done when |
|----:|:------:|-------------|-----------|
| 22  | ⬜     | Vignette 1: "Getting started" | Knits under 30s |
| 23  | ⬜     | Vignette 2: "Distributed fitting" | Knits; recovery vs pooled shown |
| 24  | ⬜     | Vignette 3: "OD migration" (small `inst/extdata/`) | Knits; reproduces one paper figure |
| 25  | ⬜     | All man pages, `README.Rmd`, `NEWS.md`, `pkgdown` config | `pkgdown::build_site()` builds |
| 26  | ⬜     | `R CMD check --as-cran` clean | 0 errors, 0 warnings, ≤1 NOTE |
| 27  | ⬜     | GitHub Actions CI, codecov, **release v0.1.0** (Zenodo DOI for CV) | Green badges; v0.1.0 tagged with binaries |
| 28  | ⬜     | CRAN submission via `devtools::release()` | Submission confirmation email |

---

## Known risks

1. **Day 5 formula parser** — deferred. Using explicit-spec API in v0.1; formula sugar in v0.2.
2. **Week 3 toolchain** — needs Rtools (Windows) / Xcode CLT (Mac). Verify before Day 15.
3. **CRAN review queue** — 1–3 weeks after Day 28. First submissions often bounce once for trivial fixes. Budget +1 week.

## CV-ready milestones

- **Day 27:** GitHub release v0.1.0 with Zenodo DOI — citeable, can list on CV/website immediately.
- **Day 28 + ~2 weeks:** CRAN acceptance — "available on CRAN" line on CV.
