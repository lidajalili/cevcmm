# Predictions from a fitted VCMM

Produces predicted responses at new `(X, t)` (and optionally `Z`)
values, using the paper's subject-specific predictor \$\$\hat y_i =
\tilde{\mathbf x}\_i^{\top} \hat{\tilde\beta} + \mathbf z_i^{\top}
\hat\alpha\$\$ by default (Jalili and Lin, 2025, Section 5). For the
alternative "new groups not seen in training" scenario, pass
`include_random = FALSE` to use the marginal predictor \\\hat y_i =
\tilde{\mathbf x}\_i^{\top} \hat{\tilde\beta}\\.

## Usage

``` r
# S3 method for class 'vcmm_fit'
predict(object, newdata, include_random = TRUE, se.fit = FALSE, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- newdata:

  A named list or data frame. See Details.

- include_random:

  Logical. If `TRUE` (default) and `newdata\$Z` is supplied, adds
  \\Z\hat\alpha\\ to the prediction. Set `FALSE` for the marginal
  predictor (new groups scenario).

- se.fit:

  Logical. If `TRUE`, also returns per-prediction standard errors.

- ...:

  Unused.

## Value

Either a numeric vector of length \\N\_{\mathrm{new}}\\, or (when
`se.fit = TRUE`) a list with components `fit` and `se.fit`.

## Details

**newdata format.** A named list (or data frame whose columns match
these names) containing:

- `t`: numeric vector of length \\N\_{\mathrm{new}}\\.

- `X`: numeric matrix \\N\_{\mathrm{new}} \times K\\ (or
  length-\\N\_{\mathrm{new}}\\ vector when \\K = 1\\) of varying-
  coefficient covariates, in the same column order used at fit time.

- `Z`: optional \\N\_{\mathrm{new}} \times q\\ random- effects design
  matrix. Must follow the same column-stacking convention as the
  training `Z` so that `Z %*% alpha` references the appropriate
  random-effect slots.

**Standard errors.** With `se.fit = TRUE` the per-prediction standard
error is the square root of \\\[W, Z\]\\\hat\sigma\_\varepsilon^2
K^{-1}\\\[W, Z\]^{\top}\\ when `include_random = TRUE` (joint
uncertainty of fixed and random effects), where \\W\\ is the
spline-expanded fixed-effects row; the random-effect block is omitted
when `include_random = FALSE`. These are confidence-interval SEs on the
mean; for prediction intervals add \\\hat\sigma\_\varepsilon^2\\ to the
variance.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## See also

[`varying_coef`](https://lidajalili.github.io/cevcmm/reference/varying_coef.md),
[`vcov.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/vcov.vcmm_fit.md).

## Examples

``` r
set.seed(1)
n <- 400
t <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 3), n, 3)
y <- 2 + sin(2 * pi * t) * x +
     as.vector(Z %*% rnorm(3, sd = 0.5)) + rnorm(n, sd = 0.5)

fit <- vcmm(y, X = x, Z = Z, t = t,
            control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))

# Default: subject-specific predictor (training-group prediction)
yhat_train <- predict(fit, newdata = list(t = t, X = x, Z = Z))
mean((y - yhat_train)^2)  # ~ sigma_eps^2 = 0.25
#> [1] 0.2853015

# Marginal predictor (new groups scenario)
yhat_marg <- predict(fit, newdata = list(t = t, X = x, Z = Z),
                      include_random = FALSE)
```
