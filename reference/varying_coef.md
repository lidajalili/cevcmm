# Evaluate the varying coefficient(s) at new t values

For a VCMM with varying coefficients \\\beta_k(t), k = 1, ..., K\\, this
function evaluates \\\hat\beta_k(t)\\ at any vector of `t_new` values,
using the same B-spline basis the fit was built on. Useful for
diagnostic plots, predictions, and reporting; the Day-13 plot method
uses it internally.

## Usage

``` r
varying_coef(object, t_new, k = NULL, se.fit = FALSE, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- t_new:

  Numeric vector at which to evaluate. Same scale as the original `t`;
  the package's normalisation is applied internally.

- k:

  Integer vector of which varying coefficients to evaluate (1-based).
  Default `NULL` = all `K`.

- se.fit:

  Logical. If `TRUE`, also returns pointwise standard errors.

- ...:

  Unused.

## Value

Either a numeric matrix (`length(t_new)` by `length(k)`, default), or a
list with components `fit` and `se.fit` when `se.fit = TRUE`.

## Details

The constant intercept \\\hat\beta_0\\ is *not* returned (it does not
vary in `t`); use
[`fixef`](https://lidajalili.github.io/cevcmm/reference/fixef.md)`(fit)\$intercept`
for that.

Pointwise standard errors are available via `se.fit = TRUE`, computed as
\$\$ \mathrm{SE}(\hat\beta_k(t)) = \sqrt{B(t)^\top
\widehat{\mathrm{Var}}(\beta\_{(k)})\\ B(t)} \$\$ where \\\beta\_{(k)}\\
is the basis-coefficient sub-vector for coefficient \\k\\ and the
covariance comes from `vcov(object, which = "beta")`.

## See also

[`fixef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/fixef.vcmm_fit.md),
[`coef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/coef.vcmm_fit.md).
