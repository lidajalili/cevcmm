# Variance-covariance matrix of the fixed-effects from a vcmm fit

Returns the asymptotic variance-covariance matrix of the fixed-effects
coefficient vector \\\hat\beta\\ from a fitted `vcmm_fit`. The matrix is
computed as \$\$ \widehat{\mathrm{Var}}(\hat\beta) =
\hat\sigma\_\varepsilon^2 \cdot \[K^{-1}\]\_{1:p,\\ 1:p}, \$\$ where
\\K\\ is the prior-augmented Hessian assembled at convergence and cached
in `object$K_inv`. This is the standard plug-in asymptotic-normal
variance estimator for the linear normal VCMM with fixed variance
components.

## Usage

``` r
# S3 method for class 'vcmm_fit'
vcov(object, which = c("beta", "alpha", "both"), ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- which:

  Character: `"beta"` (default), `"alpha"`, or `"both"`.

- ...:

  Unused.

## Value

A numeric matrix:

- `"beta"`: p by p.

- `"alpha"`: q by q.

- `"both"`: (p+q) by (p+q), joint.

## Details

Pass `which = "alpha"` for the random-effect block, `which = "both"` for
the full \\(p+q) \times (p+q)\\ joint matrix.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
n <- 300
t <- runif(n); x <- runif(n)
Z <- matrix(rnorm(n * 3), n, 3)
y <- 2 + sin(2 * pi * t) * x +
     as.vector(Z %*% rnorm(3, sd = 0.5)) + rnorm(n, sd = 0.5)

fit <- vcmm(y, X = x, Z = Z, t = t,
            control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))

V_beta <- vcov(fit)
dim(V_beta)
#> [1] 11 11
```
