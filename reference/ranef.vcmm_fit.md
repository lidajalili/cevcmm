# Random effects of a VCMM, reshaped by re_cov structure

For `re_cov = "diag"`, returns a named numeric vector of length \\q\\.
For `"kronecker"` and `"separable"`, reshapes \\\hat\alpha\\ into a \\G
\times q\_{\mathrm{left}}\\ matrix (column-stacking convention used
throughout the package). Row names are `"g1", ..., "gG"`; column names
are `c("origin", "dest")` when `re_cov = "kronecker"` and `q_left = 2`,
and `"k1", ..., "kK"` otherwise.

## Usage

``` r
# S3 method for class 'vcmm_fit'
ranef(object, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- ...:

  Unused.

## Value

Numeric vector (`"diag"`) or numeric matrix (`"kronecker"` /
`"separable"`).

## See also

[`coef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/coef.vcmm_fit.md),
[`fixef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/fixef.vcmm_fit.md).
