# Summarise a vcmm fit

Produces a `vcmm_summary` object with a coefficient table for the
fixed-effects, the variance components, and basic fit diagnostics.
[`print()`](https://rdrr.io/r/base/print.html) renders it in the style
of `lm` / `lmer`.

## Usage

``` r
# S3 method for class 'vcmm_fit'
summary(object, ...)

# S3 method for class 'vcmm_summary'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object:

  A `vcmm_fit` object.

- ...:

  Unused.

- x:

  A `vcmm_summary` object.

- digits:

  Integer. Number of significant digits to display.

- signif.stars:

  Logical. If `TRUE`, append significance stars to the coefficient
  table.

## Value

A list of class `"vcmm_summary"`.

## Details

The standard errors are the square roots of the diagonal of
`vcov(object, which = "beta")`, treating \\\sigma\_\varepsilon\\ and
\\\sigma\_\alpha\\ (or \\\Sigma\_\alpha\\) as fixed at their fitted
values. They are first-order valid (Theorem 3.1) but do not adjust for
variance-component uncertainty.

For `re_cov %in% c("kronecker", "separable")`, the variance- components
block of the printed summary displays the estimated
\\\Sigma\_{\mathrm{left}}\\ (\\\Sigma\_{2\times 2}\\ for OD or
\\\Sigma_q\\ for group-shared dense) instead of a scalar
\\\sigma\_\alpha\\.
