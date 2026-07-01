# Fixed effects of a VCMM, reshaped by varying-coefficient block

Splits the coefficient vector returned by
[`coef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/coef.vcmm_fit.md)
into:

- `intercept`: the constant scalar \\\hat\beta_0\\.

- `varying`: an \\m \times K\\ matrix of B-spline basis coefficients,
  with row names `"basis1", ..., "basisM"` and column names
  `"X1", ..., "XK"`.

For \\K = 0\\ (no varying covariate; intercept-only model), the
`varying` slot is `NULL`.

## Usage

``` r
# S3 method for class 'vcmm_fit'
fixef(object, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- ...:

  Unused.

## Value

A two-element list.

## See also

[`coef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/coef.vcmm_fit.md),
[`varying_coef`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md),
[`ranef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/ranef.vcmm_fit.md).
