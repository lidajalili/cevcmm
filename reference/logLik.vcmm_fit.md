# Log-likelihood of a vcmm fit

Returns the marginal log-likelihood \$\$\ell(\hat\beta,
\hat\sigma\_\varepsilon, \hat\Sigma\_\alpha) = -\tfrac{n}{2}\log(2\pi) -
\tfrac{1}{2}\log\|\Sigma_y\| - \tfrac{1}{2}(y - X\hat\beta)^{\top}
\Sigma_y^{-1}(y - X\hat\beta),\$\$ evaluated at the fitted parameter
values, with \\\Sigma_y = \sigma\_\varepsilon^2 I +
Z\\\Sigma\_\alpha\\Z^{\top}\\. The value is computed once at convergence
and cached on the fit object as `object$marginal_loglik`; this method
simply retrieves it and attaches `df` and `nobs` attributes so that
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work out of the box.

## Usage

``` r
# S3 method for class 'vcmm_fit'
logLik(object, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- ...:

  Unused.

## Value

An object of class `"logLik"`; numeric scalar with `df` and `nobs`
attributes.

## Details

Degrees of freedom counted are \\p\\ (fixed-effects, including all
spline basis coefficients) plus the number of free variance-component
parameters:

- `re_cov = "diag"`: 2 (\\\sigma\_\varepsilon\\, \\\sigma\_\alpha\\).

- `re_cov = "kronecker"` / `"separable"`: \\1 +
  q\_{\mathrm{left}}(q\_{\mathrm{left}} + 1)/2\\
  (\\\sigma\_\varepsilon\\ plus the free entries of
  \\\Sigma\_{\mathrm{left}}\\). \\\Sigma\_{\mathrm{right}}\\ is held
  fixed at its user-supplied value and contributes `0` df.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
n <- 400
t <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 3), n, 3)
y <- 2 + sin(2 * pi * t) * x +
     as.vector(Z %*% rnorm(3, sd = 0.5)) + rnorm(n, sd = 0.5)
fit <- vcmm(y, X = x, Z = Z, t = t,
            control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
logLik(fit)
#> 'log Lik.' -328.0214 (df=14)
AIC(fit)
#> [1] 684.0428
BIC(fit)
#> [1] 739.9233
```
