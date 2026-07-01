# Split-merge SVD via row partitioning

Computes an SVD of a matrix \\X\\ by partitioning its rows into \\s\\
blocks, taking the SVD of each block, and merging the results. More
numerically stable than direct
[`svd()`](https://rdrr.io/r/base/svd.html) for ill-conditioned large
matrices.

## Usage

``` r
split_merge_svd_row(X, s = 10, verbose = FALSE)
```

## Arguments

- X:

  Numeric matrix.

- s:

  Integer. Number of row partitions (default 10).

- verbose:

  Logical. Currently unused; retained for API compatibility.

## Value

A list with elements `u`, `d`, `v`, matching the structure returned by
base [`svd`](https://rdrr.io/r/base/svd.html).

## Details

This is an internal helper for
[`invert_matrix()`](https://lidajalili.github.io/cevcmm/reference/invert_matrix.md)
and
[`svd_pseudo_inverse()`](https://lidajalili.github.io/cevcmm/reference/svd_pseudo_inverse.md);
end users typically should not call it directly.
