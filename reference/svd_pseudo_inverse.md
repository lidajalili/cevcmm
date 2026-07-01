# SVD-based Moore-Penrose pseudo-inverse

Computes \\A^{+}\\ via SVD, optionally using the split-merge variant for
large matrices. Singular values below
`.Machine$double.eps * d[1] * max(dim(A))` are treated as zero.

## Usage

``` r
svd_pseudo_inverse(A, use_split_merge = TRUE, verbose = FALSE)
```

## Arguments

- A:

  Numeric square matrix.

- use_split_merge:

  Logical. If `TRUE` and `nrow(A) >= 100`, use
  [`split_merge_svd_row()`](https://lidajalili.github.io/cevcmm/reference/split_merge_svd_row.md);
  otherwise use base [`svd()`](https://rdrr.io/r/base/svd.html).

- verbose:

  Logical. If `TRUE`, print method and condition number.

## Value

A list with elements `inverse` (the pseudo-inverse matrix),
`condition_number`, `effective_rank`, and `method` (character:
`"standard SVD"` or `"split-merge SVD"`).

## Details

This is an internal helper for
[`invert_matrix()`](https://lidajalili.github.io/cevcmm/reference/invert_matrix.md);
end users should call
[`invert_matrix()`](https://lidajalili.github.io/cevcmm/reference/invert_matrix.md)
instead.
