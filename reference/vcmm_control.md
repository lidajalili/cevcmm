# Control parameters for VCMM fitting

Builds a validated options list controlling the iterative SS estimator
(and, later, the CSL and SVD-stabilised estimators). Pass the returned
object as the `control` argument of
[`fit_ss()`](https://lidajalili.github.io/cevcmm/reference/fit_ss.md).

## Usage

``` r
vcmm_control(
  max_iter = 200L,
  tol_beta = 1e-06,
  tol_alpha = 1e-06,
  sigma_eps = 1,
  sigma_alpha = 1,
  update_variance = FALSE,
  verbose = FALSE
)

# S3 method for class 'vcmm_control'
print(x, ...)
```

## Arguments

- max_iter:

  Integer. Maximum number of iterations (default 200).

- tol_beta:

  Positive numeric. Relative-change convergence tolerance for the
  fixed-effects coefficient vector `beta` (default 1e-6).

- tol_alpha:

  Positive numeric. Relative-change convergence tolerance for the
  random-effects vector `alpha` (default 1e-6).

- sigma_eps:

  Positive numeric. Initial residual standard deviation. If
  `update_variance = FALSE`, this value is held fixed throughout fitting
  (default 1).

- sigma_alpha:

  Positive numeric. Initial random-effect standard deviation. If
  `update_variance = FALSE`, this value is held fixed throughout fitting
  (default 1).

- update_variance:

  Logical. If `FALSE` (default), `sigma_eps` and `sigma_alpha` are held
  fixed at the supplied values – matching Algorithm 1 of Jalili and
  Lin (2025) as written. If `TRUE`, both are re-estimated at every
  iteration using the residual sum of squares formula (for `sigma_eps`)
  and a method-of-moments update (for `sigma_alpha`).

- verbose:

  Logical. If `TRUE`, print progress every 20 iterations (default
  `FALSE`).

- x:

  A `vcmm_control` object.

- ...:

  Unused.

## Value

A list of class `"vcmm_control"` containing the validated options.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
# Defaults: fix variances at 1, iterate up to 200 times.
ctrl <- vcmm_control()
ctrl
#> <vcmm_control>  fitting options
#>   max_iter        : 200
#>   tol_beta        : 1.00e-06
#>   tol_alpha       : 1.00e-06
#>   sigma_eps       : 1.0000
#>   sigma_alpha     : 1.0000
#>   update_variance : FALSE
#>   verbose         : FALSE

# Fix variances at user-supplied values.
ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5)

# Re-estimate variances each iteration.
ctrl <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5,
                     update_variance = TRUE)
```
