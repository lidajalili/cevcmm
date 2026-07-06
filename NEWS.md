# cevcmm 0.1.0 (in development)

First public release.

## Features

### Main fit function

* `vcmm()` is the main entry point: a single call that builds the spline
  design, computes sufficient statistics, and fits a Varying Coefficient
  Mixed-Effects Model. The default `method = "auto"` routes between the
  iterative sufficient-statistics estimator (SS) and the one-step
  communication-efficient surrogate-likelihood estimator (CSL) of
  Jalili and Lin (2025), based on the problem size.

### Random-effect covariance

Three structures, all selectable via `re_cov`:

* `re_cov = "diag"`: independent random effects with shared variance
  $\sigma_\alpha^2$.
* `re_cov = "kronecker"`: $\Sigma_\alpha = \Sigma_{\text{left}} \otimes
  \Sigma_{\text{right}}$, with $\Sigma_{\text{left}}$ (`q_left` x
  `q_left`) estimated and $\Sigma_{\text{right}}$ (`G` x `G`) held at
  its initial value. Designed for origin-destination flow data.
* `re_cov = "separable"`: same internal Kronecker machinery as the
  kronecker mode, but with the `q_left` argument required (no default)
  and the right-side $\Omega_G$ held fixed at its initial value.

### Distributed and streaming API

* `node_summary()` computes a sufficient-statistics summary on a slice
  of data; the resulting `vcmm_ss` object is independent of the slice
  size (a few KB regardless of $n_k$) and can be shipped between
  machines.
* `+.vcmm_ss` adds summaries together; `Reduce("+", list_of_summaries)`
  aggregates an arbitrary number.
* `init_accumulator()` plus repeated calls to `accumulate_stats()`
  provides a streaming alternative for chunks that arrive over time.
* `fit_from_summaries()` accepts a single summary, a list of summaries,
  or an accumulator, and returns a fit bit-equivalent to a pooled
  `vcmm()` call on the same data (Theorem 1 of Jalili and Lin, 2025).
* For designs with constant `rowSums(Z)` (e.g. origin-destination
  indicators), `fit_from_summaries()` takes a `rowsum_constant`
  argument so the same identifiability shift `vcmm()` applies
  automatically can be reproduced in the distributed setting.

### S3 methods

Standard generics work on the `vcmm_fit` class: `print()`, `summary()`,
`coef()`, `fixef()`, `ranef()`, `vcov()`, `nobs()`, `logLik()`, `AIC()`,
`BIC()`, `predict()`, `plot()`, plus the package-specific helper
`varying_coef()` for evaluating $\hat\beta_k(t)$ at new $t$ values.

* `vcov()` returns the asymptotic Wald variance from the prior-augmented
  Hessian inverse. `which = "beta"` (default), `"alpha"`, or `"both"`.
* `predict()` supports both subject-specific (with `Z`) and marginal
  predictions, with optional pointwise standard errors via `se.fit`.
* `varying_coef()` evaluates $\hat\beta_k(t)$ on an arbitrary grid with
  optional pointwise standard errors.
* `plot.vcmm_fit()` shows three diagnostic panels (`which = 1, 2, 3`):
  varying coefficient with CI band; residuals vs fitted and vs $t$;
  random-effect diagnostics that adapt to the chosen `re_cov` (Q-Q
  plot for diag, heatmap of the estimated $\Sigma_{\text{left}}$ for
  kronecker and separable).

### Performance backends

* RcppArmadillo backends for the hot paths: sufficient-statistics
  cross-products, Cholesky / LU / Moore-Penrose inversion. Falls back
  cleanly to R-only implementations if the C++ build is unavailable.
* Optional truncated-SVD path via RSpectra (in `Suggests:`), enabled
  with `options(cevcmm.use_rspectra = TRUE)` for the dense
  high-dimensional case.

## Bundled data

* `inst/extdata/od_migration.csv`: a 3000-row simulated origin-destination
  migration dataset (10 regions over 30 years, all 100 OD pairs observed
  annually). Used by the `od-migration` vignette.

## Vignettes

* `getting-started` — simulate, fit, and inspect on a simple
  diagonal-random-effects example.
* `distributed-fitting` — split data across multiple nodes and recover
  the same fit; covers both list-of-summaries and streaming-accumulator
  patterns.
* `od-migration` — Kronecker covariance for origin-destination flow data
  using the bundled CSV example.

## Reference

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed Effect Models: Methodology, Theory, and
Applications. arXiv:2511.12732; under review at *Journal of the American
Statistical Association*.
