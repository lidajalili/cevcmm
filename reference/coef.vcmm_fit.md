# Fixed-effects coefficient vector from a vcmm fit

Returns \\\hat\beta\\ as a named numeric vector. Under the package's
design (`X_design = cbind(1, B*x_1, ..., B*x_K)`, see
[`build_vcmm_design`](https://lidajalili.github.io/cevcmm/reference/build_vcmm_design.md)),
entry 1 is the constant intercept and entries
`(1 + (k-1)*m + 1):(1 + k*m)` are the spline basis coefficients for the
\\k\\-th varying coefficient \\\beta_k(t)\\.

## Usage

``` r
# S3 method for class 'vcmm_fit'
coef(object, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- ...:

  Unused.

## Value

Named numeric vector of length \\p = 1 + K \cdot m\\.

## Details

To evaluate \\\beta_k(t)\\ at user-supplied `t` values, use
[`varying_coef`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md).
To get the same vector reshaped into a (basis x covariate) matrix with
the intercept split out, use
[`fixef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/fixef.vcmm_fit.md).
