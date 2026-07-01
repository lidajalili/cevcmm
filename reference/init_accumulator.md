# Initialize an empty sufficient-statistics accumulator

Allocates a zero-filled accumulator with the correct dimensions to
receive batched calls to
[`accumulate_stats()`](https://lidajalili.github.io/cevcmm/reference/accumulate_stats.md).
Use this once before a streaming loop over data batches or nodes.

## Usage

``` r
init_accumulator(p, q)
```

## Arguments

- p:

  Integer. Number of fixed-effects columns (intercept plus spline
  basis).

- q:

  Integer. Number of random-effects columns (length of alpha).

## Value

A list of class `"vcmm_accumulator"` with the same six matrix slots as
[`compute_sufficient_stats()`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md),
all initialised to zero, plus `n_obs = 0L`.

## See also

Other sufficient statistics:
[`accumulate_stats()`](https://lidajalili.github.io/cevcmm/reference/accumulate_stats.md),
[`compute_sufficient_stats()`](https://lidajalili.github.io/cevcmm/reference/compute_sufficient_stats.md)

## Examples

``` r
acc <- init_accumulator(p = 5, q = 3)
dim(acc$C)   # 5 x 5
#> [1] 5 5
dim(acc$ZtZ) # 3 x 3
#> [1] 3 3
acc$n_obs    # 0
#> [1] 0
```
