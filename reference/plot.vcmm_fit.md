# Diagnostic plots for a vcmm fit

Three panels, each requestable via `which`:

- 1:

  Varying-coefficient curves \\\hat\beta_k(t)\\ with pointwise Wald
  confidence bands.

- 2:

  Residual diagnostics: residuals vs fitted, residuals vs \\t\\.
  **Requires `data`** since a vcmm fit does not carry the raw training
  data.

- 3:

  Random-effects diagnostics: Normal Q-Q for `re_cov = "diag"`, or a
  heatmap of the \\G \times q\_{\mathrm{left}}\\ random-effects matrix
  together with a heatmap of the estimated \\\Sigma\_{\mathrm{left}}\\
  for `"kronecker"` / `"separable"`.

## Usage

``` r
# S3 method for class 'vcmm_fit'
plot(
  x,
  which = 1:3,
  data = NULL,
  t_grid = NULL,
  n_grid = 200L,
  conf_level = 0.95,
  ask = (length(which) > 1L) && interactive(),
  ...
)
```

## Arguments

- x:

  A `vcmm_fit` object.

- which:

  Integer vector subset of `1:3`. Default `1:3`.

- data:

  Optional list with components `y`, `X`, `Z`, `t` (the training data,
  or any data on which to compute residuals). Required when panel 2 is
  requested.

- t_grid:

  Numeric. Grid of \\t\\ values for panel 1. Default `NULL` means an
  evenly-spaced grid over the stored training range `[t_min, t_max]`.

- n_grid:

  Integer. Number of grid points if `t_grid` is `NULL`. Default 200.

- conf_level:

  Numeric in (0, 1). Confidence level for panel-1 bands. Default 0.95.

- ask:

  Logical. Passed to
  [`devAskNewPage()`](https://rdrr.io/r/grDevices/devAskNewPage.html)
  when multiple panels are requested in an interactive session.

- ...:

  Further arguments passed to base graphics calls.

## Value

Invisibly `NULL`. Called for side effects (plots).

## Details

Uses base R graphics, no ggplot2 dependency. Each requested panel is
drawn on its own figure (or set of subfigures via `par(mfrow)`); call
`par(mfrow = c(2, 2))` or similar before
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to combine
panels in one figure.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## See also

[`varying_coef`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md),
[`ranef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/ranef.vcmm_fit.md),
[`predict.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/predict.vcmm_fit.md).

## Examples

``` r
set.seed(1)
n <- 500
t <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 3), n, 3)
a <- rnorm(3, sd = 0.5)
y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% a) + rnorm(n, sd = 0.5)
fit <- vcmm(y, X = x, Z = Z, t = t,
            control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
if (FALSE) { # \dontrun{
plot(fit)                                         # all three panels
plot(fit, which = 1)                              # only varying coefs
plot(fit, which = 2, data = list(y = y, X = x, Z = Z, t = t))
plot(fit, which = 3)                              # ranef diagnostics
} # }
```
