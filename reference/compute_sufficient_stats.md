# Compute one batch or node sufficient statistics for a normal linear VCMM

Computes the six-component summary that fully encodes one node's
contribution to the normal linear VCMM likelihood. These statistics are
additive across nodes and batches, so a central server can recover the
full-data estimator by summing per-node summaries – without ever seeing
the raw response or design matrices.

## Usage

``` r
compute_sufficient_stats(y, X, Z, use_cpp = TRUE)
```

## Arguments

- y:

  Numeric response vector of length n.

- X:

  Numeric n by p fixed-effects design matrix (intercept plus spline
  basis columns).

- Z:

  Numeric n by q random-effects design matrix.

- use_cpp:

  Logical. If `TRUE` (default), dispatch to the RcppArmadillo backend
  `compute_sufficient_stats_cpp()`. If `FALSE`, use the pure-R
  [`crossprod()`](https://rdrr.io/r/base/crossprod.html) reference path.
  The two paths agree to within floating-point summation order
  (typically below 1e-13 relative).

## Value

A list of class `"vcmm_ss"` with elements `a`, `b`, `C`, `ZtZ`, `Zty`,
`XtZ`, and `n_obs` (the number of observations summarized).

## Details

The returned components are:

- `a`: `sum(y^2)`, a scalar.

- `b`: `crossprod(X, y)`, dimension p by 1.

- `C`: `crossprod(X)`, dimension p by p.

- `ZtZ`: `crossprod(Z)`, dimension q by q.

- `Zty`: `crossprod(Z, y)`, dimension q by 1.

- `XtZ`: `crossprod(X, Z)`, dimension p by q.

Since Day 16 the default backend is a RcppArmadillo implementation
(`use_cpp = TRUE`). Pass `use_cpp = FALSE` to fall back to the pure-R
reference path; this is used by the Day-16 bit-equivalence validation
and is also useful for debugging.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## See also

Other sufficient statistics:
[`accumulate_stats()`](https://lidajalili.github.io/cevcmm/reference/accumulate_stats.md),
[`init_accumulator()`](https://lidajalili.github.io/cevcmm/reference/init_accumulator.md)

## Examples

``` r
set.seed(1)
n <- 100; p <- 3; q <- 2
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
Z <- matrix(rnorm(n * q), n, q)
y <- rnorm(n)

ss <- compute_sufficient_stats(y, X, Z)
str(ss)
#> List of 7
#>  $ a    : num 136
#>  $ b    : num [1:3, 1] -3.91 14.87 15.14
#>  $ C    : num [1:3, 1:3] 100 10.89 -3.78 10.89 81.06 ...
#>  $ ZtZ  : num [1:2, 1:2] 106 11.4 11.4 97.6
#>  $ Zty  : num [1:2, 1] 3.44 2.54
#>  $ XtZ  : num [1:3, 1:2] 2.97 2.01 -4.97 5.16 -3.89 ...
#>  $ n_obs: int 100
#>  - attr(*, "class")= chr [1:2] "vcmm_ss" "list"
```
