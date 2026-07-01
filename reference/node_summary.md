# Compute one node's sufficient-statistics summary

Convenience alias for
[`compute_sufficient_stats`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md),
intended for distributed-computing workflows. Each compute node calls
this on its local data and transmits the small returned summary; the
central server aggregates summaries with `+` (or
`Reduce("+", summaries)`) and fits via
[`fit_from_summaries`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md).

## Usage

``` r
node_summary(y, X, Z)
```

## Arguments

- y:

  Numeric response vector of length n.

- X:

  Numeric n by p fixed-effects design matrix (intercept plus spline
  basis columns).

- Z:

  Numeric n by q random-effects design matrix.

## Value

A `vcmm_ss` object (additive via `+.vcmm_ss`).

## Details

**What's transmitted.** The returned object contains six fixed-size
arrays whose dimensions depend only on \\p = \mathrm{ncol}(X)\\ and \\q
= \mathrm{ncol}(Z)\\, never on the node-local sample size. For a typical
VCMM with \\p \approx 15\\ and \\q \approx 50\\, one summary is a few
kilobytes regardless of \\N_s\\.

**Basis alignment.** All nodes must build their local `X` (the
spline-expanded fixed-effects design) using the *same* basis
specification. The recommended workflow is to call
[`build_vcmm_design`](https://lidajalili.github.io/cevcmm/reference/build_vcmm_design.md)
once with the full `t`-range (or pre-agreed knots), broadcast the
resulting basis, and have each node evaluate it on local `t`. The
package itself does not enforce this; mis-aligned bases will silently
give incorrect fits.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## See also

[`fit_from_summaries`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md),
[`compute_sufficient_stats`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md),
[`build_vcmm_design`](https://lidajalili.github.io/cevcmm/reference/build_vcmm_design.md)

## Examples

``` r
set.seed(1)
n_per_node <- 100; p <- 5; q <- 3
X <- cbind(1, matrix(rnorm(n_per_node * (p - 1)), n_per_node, p - 1))
Z <- matrix(rnorm(n_per_node * q), n_per_node, q)
y <- rnorm(n_per_node)

gamma_1 <- node_summary(y, X, Z)
gamma_2 <- node_summary(y, X, Z)
gamma_pooled <- gamma_1 + gamma_2
gamma_pooled
#> <vcmm_ss>  one-batch sufficient statistics
#>   n_obs : 200
#>   p     : 5   (fixed-effects columns)
#>   q     : 3   (random-effects columns)
#>   a     : 238.4
```
