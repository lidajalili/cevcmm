# Fit a VCMM via the one-step CSL estimator

Performs a single Newton refinement of a pilot estimator using the
aggregated sufficient statistics, implementing the one-step
communication-efficient surrogate-likelihood (CSL) estimator of Jalili
and Lin (2025, Section 3). For the normal linear VCMM the update is \$\$
\widehat{\theta}\_{CSL} = \widehat{\theta}\_0 - K^{-1} \\
g(\widehat{\theta}\_0), \$\$ where \\\theta = (\beta, \alpha)\\, the
gradient \\g\\ uses the full aggregated stats, and the Hessian \\K\\
also uses the full aggregated stats with prior augmentation. Theorem 3.1
of the paper shows that whenever the pilot estimator is
\\\sqrt{N}\\-consistent, the one-step CSL estimator is first-order
equivalent to the full SS estimator.

## Usage

``` r
fit_csl(
  stats,
  penalty,
  control = vcmm_control(),
  pilot = NULL,
  pilot_max_iter = 5L,
  pilot_tol_beta = 0.001,
  pilot_tol_alpha = 0.001,
  re_cov_state = NULL
)
```

## Arguments

- stats:

  A `vcmm_ss` or `vcmm_accumulator` object containing the aggregated
  sufficient statistics.

- penalty:

  A symmetric \\p \times p\\ penalty matrix from
  [`build_penalty_matrix()`](https://lidajalili.github.io/cevcmm/reference/build_penalty_matrix.md).

- control:

  A `vcmm_control` object with fitting options. The `sigma_eps` and
  `sigma_alpha` entries provide the initial variance values used by the
  internal pilot (when `pilot = NULL`).

- pilot:

  Optional `vcmm_fit` object to use as the pilot estimator. If `NULL`
  (default), an internal SS pilot is run.

- pilot_max_iter:

  Integer. Maximum iterations for the internal SS pilot (default 5).
  Ignored if `pilot` is supplied.

- pilot_tol_beta:

  Positive numeric. Loose tolerance for the internal SS pilot (default
  1e-3). Ignored if `pilot` is supplied.

- pilot_tol_alpha:

  Positive numeric. Loose tolerance for the internal SS pilot (default
  1e-3). Ignored if `pilot` is supplied.

- re_cov_state:

  Optional. An internal random-effects covariance state object (NULL for
  diagonal, or constructed for kronecker via
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) with
  `re_cov = "kronecker"`). Passed through to the internal SS pilot and
  reused for the Newton-step Hessian. Advanced users typically reach
  `re_cov_state` via
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md)
  rather than calling `fit_csl()` directly.

## Value

A list of class `"vcmm_fit"` with the same fields as
[`fit_ss()`](https://lidajalili.github.io/cevcmm/reference/fit_ss.md),
plus:

- `pilot`: the `vcmm_fit` pilot used.

- `pilot_elapsed_sec`: pilot fitting time.

- `step_elapsed_sec`: Newton step time.

The `method` field is `"CSL"`.

## Details

**Pilot estimator.** If `pilot = NULL` (default), an internal SS pilot
is run with loose tolerances and a small number of iterations
(`pilot_max_iter = 5`, `pilot_tol_beta = 1e-3`,
`pilot_tol_alpha = 1e-3`); these defaults match the dense simulation
study of the paper. Setting `pilot_max_iter = 1L` gives the cheapest
possible pilot (the OLS-like single-step pilot used in the
origin-destination simulation), still \\\sqrt{N}\\-consistent in the
normal linear case. Advanced users may pass any `vcmm_fit` object via
the `pilot` argument – e.g. a previously fitted
[`fit_ss()`](https://lidajalili.github.io/cevcmm/reference/fit_ss.md)
result.

**Hessian.** This implementation uses the full aggregated Hessian, not
the reference-node curvature approximation \\\tilde K\\ from the paper.
The two are first-order equivalent under the conditions of Theorem 3.1;
the full-aggregated form is the most accurate and is the natural default
for a single-server fit, which is the typical use case of this package.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
n <- 500; p <- 4; q <- 3
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
Z <- matrix(rnorm(n * q), n, q)
alpha_true <- rnorm(q, sd = 0.5)
y <- as.vector(
  X %*% c(2, 0.5, -0.3, 0.8) + Z %*% alpha_true + rnorm(n, sd = 0.5)
)

stats <- compute_sufficient_stats(y, X, Z)
P     <- build_penalty_matrix(n_basis = p - 1, lambda = 0.1)

# Default: internal SS pilot (5 loose iterations) then one Newton step
fit <- fit_csl(stats, P,
               vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
fit
#> <vcmm_fit>  Varying Coefficient Mixed-Effects Model fit
#>   method      : CSL
#>   n_obs       : 500
#>   p (fixed)   : 4
#>   q (random)  : 3
#>   RE cov      : diag
#>   pilot iter  : 3 (converged)
#>   newton step : 1
#>   sigma_eps   : 0.5000
#>   sigma_alpha : 0.5000
#>   elapsed     : 0.0010 sec (pilot <0.001s + newton <0.001s)

# Cheapest pilot: single SS step
fit_one <- fit_csl(stats, P,
                   vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5),
                   pilot_max_iter = 1L)

# Advanced: user-supplied pilot
my_pilot <- fit_ss(stats, P,
                   vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5,
                                max_iter = 3))
fit_user <- fit_csl(stats, P,
                    vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5),
                    pilot = my_pilot)
```
