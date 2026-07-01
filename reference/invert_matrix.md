# Numerically stable matrix inversion with automatic method selection

Inverts a square matrix \\A\\ using a dispatch rule that balances speed
and numerical stability:

- If `q < 100`: try Cholesky (via `invert_spd_cpp()`), falling through
  to LU (`invert_general_cpp()`) and finally SVD pseudo-inverse if the
  matrix is not positive-definite.

- If `q >= 100`: route through the SVD pseudo-inverse path. By default
  this is the dense LAPACK split-merge variant
  ([`svd_pseudo_inverse()`](https://lidajalili.github.io/cevcmm/reference/svd_pseudo_inverse.md),
  paper Algorithm 2). When the user opts in (see the `method` argument
  below), the iterative truncated SVD from
  [`RSpectra::svds()`](https://rdrr.io/pkg/RSpectra/man/svds.html) is
  used instead – faster when the matrix has effective rank much smaller
  than `q`, slower otherwise.

The Cholesky fast path is roughly 2-3x faster than the original
`kappa(A) + solve(A)` R path on VCMM K matrices. See the Day 17 and Day
18 validation scripts in `inst/validation/` for the bit-equivalence and
timing details.

## Usage

``` r
invert_matrix(
  A,
  q = NULL,
  verbose = FALSE,
  use_cpp = TRUE,
  method = c("auto", "lapack", "rspectra")
)
```

## Arguments

- A:

  Numeric square matrix to invert.

- q:

  Optional integer. Routing dimension used to pick the inversion
  strategy (defaults to `nrow(A)`). Pass an explicit value if you know
  `A` is a curvature block with a meaningful dimension that differs from
  its row count.

- verbose:

  Logical. If `TRUE`, print the chosen method.

- use_cpp:

  Logical. If `TRUE` (default since Day 17), use the RcppArmadillo
  Cholesky/LU backend. If `FALSE`, use the original pure-R path with
  `kappa(A)` plus `solve(A)`.

- method:

  Character string, one of `"auto"` (default), `"lapack"`, or
  `"rspectra"`. Only affects the `q >= 100` branch. With `"auto"`, the
  routing checks `getOption("cevcmm.use_rspectra", FALSE)`: if `TRUE`
  and the RSpectra package is installed, the truncated-SVD path is used;
  otherwise the dense LAPACK split-merge SVD is used. With `"lapack"`,
  always use the dense LAPACK path (the original Algorithm 2 behaviour).
  With `"rspectra"`, force the truncated SVD via RSpectra; if the matrix
  turns out to be full-rank, the function silently falls back to LAPACK
  so results are always correct.

## Value

A numeric matrix with the same dimensions as `A`.

## Details

If you need to reproduce the original R-only path (e.g. for a
bit-equivalence test), pass `use_cpp = FALSE`.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
A <- crossprod(matrix(rnorm(50), 10, 5)) + diag(5)
A_inv <- invert_matrix(A)
max(abs(A %*% A_inv - diag(5)))  # ~ machine epsilon
#> [1] 2.220446e-16
```
