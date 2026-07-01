# Add one batch or node statistics to a running accumulator

Performs the additive aggregation `acc <- acc + stats` component by
component. After processing all batches, the accumulator holds the
full-data sufficient summary.

## Usage

``` r
accumulate_stats(acc, stats)
```

## Arguments

- acc:

  A `vcmm_accumulator` object from
  [`init_accumulator()`](https://lidajalili.github.io/cevcmm/reference/init_accumulator.md).

- stats:

  A `vcmm_ss` object from
  [`compute_sufficient_stats()`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md).

## Value

The updated accumulator (class `"vcmm_accumulator"`).

## See also

Other sufficient statistics:
[`compute_sufficient_stats()`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md),
[`init_accumulator()`](https://lidajalili.github.io/cevcmm/reference/init_accumulator.md)

## Examples

``` r
set.seed(1)
n_batch <- 50; p <- 3; q <- 2
acc <- init_accumulator(p, q)

for (b in 1:3) {
  X <- cbind(1, matrix(rnorm(n_batch * (p - 1)), n_batch, p - 1))
  Z <- matrix(rnorm(n_batch * q), n_batch, q)
  y <- rnorm(n_batch)
  ss <- compute_sufficient_stats(y, X, Z)
  acc <- accumulate_stats(acc, ss)
}

acc$n_obs  # 150
#> [1] 150
```
