# Fit a VCMM via the sufficient-statistics estimator

Iteratively solves for the fixed-effects coefficient vector `beta` and
the random-effects vector `alpha` using only the aggregated sufficient
summary, implementing Algorithm 1 of Lin and Jalili (2026) for the
normal linear VCMM. The raw response and design matrices are not needed:
everything reads off `stats`, which may come from a single call to
[`compute_sufficient_stats()`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md)
or from a streaming accumulator
([`init_accumulator()`](https://lidajalili.github.io/cevcmm/reference/init_accumulator.md)
plus repeated
[`accumulate_stats()`](https://lidajalili.github.io/cevcmm/reference/accumulate_stats.md)).

## Usage

``` r
fit_ss(stats, penalty, control = vcmm_control(), re_cov_state = NULL)

# S3 method for class 'vcmm_fit'
print(x, ...)
```

## Arguments

- stats:

  A `vcmm_ss` or `vcmm_accumulator` object containing the aggregated
  sufficient statistics.

- penalty:

  A symmetric \\p \times p\\ penalty matrix from
  [`build_penalty_matrix()`](https://lidajalili.github.io/cevcmm/reference/build_penalty_matrix.md).

- control:

  A `vcmm_control` object with fitting options. Pass
  [`vcmm_control()`](https://lidajalili.github.io/cevcmm/reference/vcmm_control.md)
  to use defaults.

- re_cov_state:

  Optional. An internal random-effects covariance state object (NULL for
  diagonal, or constructed for kronecker via
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) with
  `re_cov = "kronecker"`). When NULL, the prior precision is
  `(sigma_eps^2 / sigma_alpha^2) * I_q` matching the diagonal case.
  Advanced users typically reach `re_cov_state` via
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md)
  rather than calling `fit_ss()` directly.

- x:

  A `vcmm_fit` object.

- ...:

  Unused; present for S3 method consistency.

## Value

A list of class `"vcmm_fit"` with elements:

- `beta`: fitted fixed-effects vector, length p.

- `alpha`: fitted random-effects vector, length q.

- `sigma_eps`: final residual standard deviation.

- `sigma_alpha`: final random-effect standard deviation.

- `iterations`: number of iterations performed.

- `converged`: `TRUE` if convergence tolerances were met.

- `elapsed_sec`: wall-clock fitting time in seconds.

- `n_obs`, `p`, `q`: data and design dimensions.

- `method`: character, `"SS"`.

- `control`: the `vcmm_control` object used.

- `call`: the matched call.

## Details

At each iteration the algorithm solves:

- beta-update: `(C + P) beta = b - XtZ %*% alpha`

- alpha-update:
  `(ZtZ + (sigma_eps^2 / sigma_alpha^2) * I) alpha = Zty - t(XtZ) %*% beta`

Both linear systems are solved via
[`invert_matrix()`](https://lidajalili.github.io/cevcmm/reference/invert_matrix.md),
which automatically dispatches between
[`solve()`](https://rdrr.io/r/base/solve.html) and SVD pseudo-inverse
depending on dimension and condition number. Convergence is declared
when the relative change in both `beta` and `alpha` falls below their
respective tolerances.

If `control$update_variance` is `TRUE`, the residual variance is
re-estimated each iteration from the SS residual sum of squares, and the
random-effect variance is re-estimated as `mean(alpha^2)`.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
n <- 200; p <- 4; q <- 3
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
Z <- matrix(rnorm(n * q), n, q)
alpha_true <- rnorm(q, sd = 0.5)
y <- as.vector(
  X %*% c(2, 0.5, -0.3, 0.8) + Z %*% alpha_true + rnorm(n, sd = 0.5)
)

# Build sufficient statistics and penalty
stats <- compute_sufficient_stats(y, X, Z)
P     <- build_penalty_matrix(n_basis = p - 1, lambda = 0.1)

# Fit with fixed variances
fit <- fit_ss(stats, P,
              vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
fit
#> <vcmm_fit>  Varying Coefficient Mixed-Effects Model fit
#>   method      : SS
#>   n_obs       : 200
#>   p (fixed)   : 4
#>   q (random)  : 3
#>   RE cov      : diag
#>   iterations  : 5 (converged)
#>   sigma_eps   : 0.5000
#>   sigma_alpha : 0.5000
#>   elapsed     : 0.0010 sec
coef(fit)
#>     beta_1     beta_2     beta_3     beta_4 
#>  2.0268439  0.5047267 -0.2624521  0.8393673 
```
