# Package index

## Main fit function

Single-call entry point for all three estimators and covariance modes.

- [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) :
  Fit a Varying Coefficient Mixed-Effects Model
- [`vcmm_control()`](https://lidajalili.github.io/cevcmm/reference/vcmm_control.md)
  [`print(`*`<vcmm_control>`*`)`](https://lidajalili.github.io/cevcmm/reference/vcmm_control.md)
  : Control parameters for VCMM fitting

## Distributed and streaming API

Per-node summaries that aggregate additively, plus the central fit
function that consumes them.

- [`node_summary()`](https://lidajalili.github.io/cevcmm/reference/node_summary.md)
  : Compute one node's sufficient-statistics summary
- [`compute_sufficient_stats()`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md)
  : Compute one batch or node sufficient statistics for a normal linear
  VCMM
- [`init_accumulator()`](https://lidajalili.github.io/cevcmm/reference/init_accumulator.md)
  : Initialize an empty sufficient-statistics accumulator
- [`accumulate_stats()`](https://lidajalili.github.io/cevcmm/reference/accumulate_stats.md)
  : Add one batch or node statistics to a running accumulator
- [`` `+`( ``*`<vcmm_ss>`*`)`](https://lidajalili.github.io/cevcmm/reference/plus-.vcmm_ss.md)
  : Additive aggregation of sufficient-statistics summaries
- [`fit_from_summaries()`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md)
  : Fit a VCMM from aggregated sufficient-statistics summaries

## Lower-level estimators

Direct access to the SS and CSL fit functions, called internally by
vcmm(). Advanced users only.

- [`fit_ss()`](https://lidajalili.github.io/cevcmm/reference/fit_ss.md)
  [`print(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/fit_ss.md)
  : Fit a VCMM via the sufficient-statistics estimator
- [`fit_csl()`](https://lidajalili.github.io/cevcmm/reference/fit_csl.md)
  : Fit a VCMM via the one-step CSL estimator

## Design and penalty

B-spline design matrix and difference-penalty construction.

- [`build_vcmm_design()`](https://lidajalili.github.io/cevcmm/reference/build_vcmm_design.md)
  : Build the design matrix and penalty for a VCMM
- [`build_penalty_matrix()`](https://lidajalili.github.io/cevcmm/reference/build_penalty_matrix.md)
  : Build a B-spline second-order difference penalty matrix

## Covariance utilities

Helpers for the Kronecker / separable covariance machinery.

- [`build_kronecker_precision()`](https://lidajalili.github.io/cevcmm/reference/build_kronecker_precision.md)
  : Build the Kronecker-structured random-effect precision matrix
- [`estimate_kronecker_components()`](https://lidajalili.github.io/cevcmm/reference/estimate_kronecker_components.md)
  : Moment-based estimator for the Kronecker left component

## Linear-algebra backends

Stabilised matrix inversion used inside the estimators.

- [`invert_matrix()`](https://lidajalili.github.io/cevcmm/reference/invert_matrix.md)
  : Numerically stable matrix inversion with automatic method selection
- [`svd_pseudo_inverse()`](https://lidajalili.github.io/cevcmm/reference/svd_pseudo_inverse.md)
  : SVD-based Moore-Penrose pseudo-inverse
- [`split_merge_svd_row()`](https://lidajalili.github.io/cevcmm/reference/split_merge_svd_row.md)
  : Split-merge SVD via row partitioning
- [`cevcmm_rcpp_check()`](https://lidajalili.github.io/cevcmm/reference/cevcmm_rcpp_check.md)
  : Verify the Rcpp / RcppArmadillo toolchain is wired up

## S3 methods on fitted models

Standard generics for inspecting and using a vcmm_fit.

- [`coef(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/coef.vcmm_fit.md)
  : Fixed-effects coefficient vector from a vcmm fit
- [`fixef()`](https://lidajalili.github.io/cevcmm/reference/fixef.md) :
  Extract fixed-effects from a fitted model object
- [`fixef(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/fixef.vcmm_fit.md)
  : Fixed effects of a VCMM, reshaped by varying-coefficient block
- [`ranef()`](https://lidajalili.github.io/cevcmm/reference/ranef.md) :
  Extract random effects from a fitted model object
- [`ranef(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/ranef.vcmm_fit.md)
  : Random effects of a VCMM, reshaped by re_cov structure
- [`vcov(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/vcov.vcmm_fit.md)
  : Variance-covariance matrix of the fixed-effects from a vcmm fit
- [`nobs(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/nobs.vcmm_fit.md)
  : Number of observations from a vcmm fit
- [`logLik(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/logLik.vcmm_fit.md)
  : Log-likelihood of a vcmm fit
- [`summary(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/summary.vcmm_fit.md)
  [`print(`*`<vcmm_summary>`*`)`](https://lidajalili.github.io/cevcmm/reference/summary.vcmm_fit.md)
  : Summarise a vcmm fit
- [`predict(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/predict.vcmm_fit.md)
  : Predictions from a fitted VCMM
- [`plot(`*`<vcmm_fit>`*`)`](https://lidajalili.github.io/cevcmm/reference/plot.vcmm_fit.md)
  : Diagnostic plots for a vcmm fit
- [`varying_coef()`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md)
  : Evaluate the varying coefficient(s) at new t values

## Package

- [`cevcmm`](https://lidajalili.github.io/cevcmm/reference/cevcmm-package.md)
  [`cevcmm-package`](https://lidajalili.github.io/cevcmm/reference/cevcmm-package.md)
  : cevcmm: Communication-Efficient Varying Coefficient Mixed-Effects
  Models
