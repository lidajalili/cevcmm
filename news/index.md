# Changelog

## cevcmm 0.1.3 (2026-07-16)

### CRAN reviewer feedback (round 2)

Addresses a follow-up NOTE from win-builder R-devel on cevcmm 0.1.2
asking that arXiv preprints be cited via their arXiv DOI form rather
than the plain URL form. Metadata-only change; no code changes.

- DESCRIPTION: change the arXiv reference from
  `<https://arxiv.org/abs/2511.12732>` to
  `<doi:10.48550/arXiv.2511.12732>`, as required by CRAN policy for
  arXiv preprints.

## cevcmm 0.1.2 (2026-07-07)

### CRAN reviewer feedback (round 1)

Addresses comments from CRAN reviewer Konstanze Lauseker on the initial
submission of cevcmm 0.1.1.

- DESCRIPTION: expand the “SVD” acronym on first use to “Singular Value
  Decomposition”.
- DESCRIPTION: reformat the reference to Jalili and Lin (2025) in the
  CRAN-preferred autolinking form (later refined to DOI form in 0.1.3).
- R/methods.R: add `\value` sections to the `fixef` and `ranef` S3
  generic function documentation, describing the return-value contract
  dispatched to method-specific documentation.
- R/plot.R: change `\dontrun{}` to `\donttest{}` in the `plot.vcmm_fit`
  example.
- .Rbuildignore: exclude `inst/validation` and `inst/benchmarks` from
  the CRAN tarball. These are development-time scripts (the Day 7-18
  validation harness and micro-benchmarks used during package
  development) that write files and change
  [`par()`](https://rdrr.io/r/graphics/par.html) /
  [`options()`](https://rdrr.io/r/base/options.html) without resetting.
  They remain available in the GitHub repository for reference but are
  no longer shipped to CRAN users.

## cevcmm 0.1.1 (2026-07-07)

### Bug fixes

- Add explicit `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` to
  `src/Makevars` for BLAS / LAPACK / Fortran linking on Linux. The
  equivalent flags were already present in `src/Makevars.win` for
  Windows. Fixes install failure with `undefined symbol: dpotrf_` on
  CRAN’s Debian pretest.

## cevcmm 0.1.0 (2026-07-01)

First public release.

### Features

#### Main fit function

- [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) is
  the main entry point: a single call that builds the spline design,
  computes sufficient statistics, and fits a Varying Coefficient
  Mixed-Effects Model. The default `method = "auto"` routes between the
  iterative sufficient-statistics estimator (SS) and the one-step
  communication-efficient surrogate-likelihood estimator (CSL) of Jalili
  and Lin (2025), based on the problem size.

#### Random-effect covariance

Three structures, all selectable via `re_cov`:

- `re_cov = "diag"`: independent random effects with shared variance
  $`\sigma_\alpha^2`$.
- `re_cov = "kronecker"`: $`\Sigma_\alpha = \Sigma_{\text{left}} \otimes
  \Sigma_{\text{right}}`$, with $`\Sigma_{\text{left}}`$ (`q_left` x
  `q_left`) estimated and $`\Sigma_{\text{right}}`$ (`G` x `G`) held at
  its initial value. Designed for origin-destination flow data.
- `re_cov = "separable"`: same internal Kronecker machinery as the
  kronecker mode, but with the `q_left` argument required (no default)
  and the right-side $`\Omega_G`$ held fixed at its initial value.

#### Distributed and streaming API

- [`node_summary()`](https://lidajalili.github.io/cevcmm/reference/node_summary.md)
  computes a sufficient-statistics summary on a slice of data; the
  resulting `vcmm_ss` object is independent of the slice size (a few KB
  regardless of $`n_k`$) and can be shipped between machines.
- `+.vcmm_ss` adds summaries together; `Reduce("+", list_of_summaries)`
  aggregates an arbitrary number.
- [`init_accumulator()`](https://lidajalili.github.io/cevcmm/reference/init_accumulator.md)
  plus repeated calls to
  [`accumulate_stats()`](https://lidajalili.github.io/cevcmm/reference/accumulate_stats.md)
  provides a streaming alternative for chunks that arrive over time.
- [`fit_from_summaries()`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md)
  accepts a single summary, a list of summaries, or an accumulator, and
  returns a fit bit-equivalent to a pooled
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) call
  on the same data (Theorem 1 of Jalili and Lin, 2025).
- For designs with constant `rowSums(Z)` (e.g. origin-destination
  indicators),
  [`fit_from_summaries()`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md)
  takes a `rowsum_constant` argument so the same identifiability shift
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md)
  applies automatically can be reproduced in the distributed setting.

#### S3 methods

Standard generics work on the `vcmm_fit` class:
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fixef()`](https://lidajalili.github.io/cevcmm/reference/fixef.md),
[`ranef()`](https://lidajalili.github.io/cevcmm/reference/ranef.md),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), plus the
package-specific helper
[`varying_coef()`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md)
for evaluating $`\hat\beta_k(t)`$ at new $`t`$ values.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the asymptotic
  Wald variance from the prior-augmented Hessian inverse.
  `which = "beta"` (default), `"alpha"`, or `"both"`.
- [`predict()`](https://rdrr.io/r/stats/predict.html) supports both
  subject-specific (with `Z`) and marginal predictions, with optional
  pointwise standard errors via `se.fit`.
- [`varying_coef()`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md)
  evaluates $`\hat\beta_k(t)`$ on an arbitrary grid with optional
  pointwise standard errors.
- [`plot.vcmm_fit()`](https://lidajalili.github.io/cevcmm/reference/plot.vcmm_fit.md)
  shows three diagnostic panels (`which = 1, 2, 3`): varying coefficient
  with CI band; residuals vs fitted and vs $`t`$; random-effect
  diagnostics that adapt to the chosen `re_cov` (Q-Q plot for diag,
  heatmap of the estimated $`\Sigma_{\text{left}}`$ for kronecker and
  separable).

#### Performance backends

- RcppArmadillo backends for the hot paths: sufficient-statistics
  cross-products, Cholesky / LU / Moore-Penrose inversion. Falls back
  cleanly to R-only implementations if the C++ build is unavailable.
- Optional truncated-SVD path via RSpectra (in `Suggests:`), enabled
  with `options(cevcmm.use_rspectra = TRUE)` for the dense
  high-dimensional case.

### Bundled data

- `inst/extdata/od_migration.csv`: a 3000-row simulated
  origin-destination migration dataset (10 regions over 30 years, all
  100 OD pairs observed annually). Used by the `od-migration` vignette.

### Vignettes

- `getting-started` — simulate, fit, and inspect on a simple
  diagonal-random-effects example.
- `distributed-fitting` — split data across multiple nodes and recover
  the same fit; covers both list-of-summaries and streaming-accumulator
  patterns.
- `od-migration` — Kronecker covariance for origin-destination flow data
  using the bundled CSV example.

### Reference

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed Effect Models: Methodology, Theory, and
Applications. arXiv:2511.12732; under review at *Journal of the American
Statistical Association*.
