# Build a B-spline second-order difference penalty matrix

Constructs a symmetric \\p \times p\\ penalty matrix
\\\mathbf{P}\_\lambda\\ for use with penalised B-spline varying
coefficients. For a single varying coefficient (`n_blocks = 1`,
default), \\p = n\\basis + 1\\; the intercept (row/column 1) is
unpenalised, and the remaining block is \\\lambda \mathbf{D}\_2^\top
\mathbf{D}\_2\\ where \\\mathbf{D}\_2\\ is the second-order difference
operator on `n_basis` coefficients.

## Usage

``` r
build_penalty_matrix(n_basis, lambda, n_blocks = 1L)
```

## Arguments

- n_basis:

  Integer (\\\geq 3\\). Number of B-spline basis columns per varying
  coefficient.

- lambda:

  Non-negative numeric. Smoothing parameter \\\lambda\\.

- n_blocks:

  Integer (\\\geq 1\\). Number of varying coefficients (covariates)
  sharing the same penalty structure. Default 1.

## Value

A symmetric `(1 + n_blocks * n_basis) x (1 + n_blocks * n_basis)`
numeric matrix. Row/column 1 is zero (intercept is unpenalised).

## Details

For `n_blocks > 1` (multiple covariates each with their own varying
coefficient), the result is block-diagonal: one zero entry for the
intercept, followed by `n_blocks` copies of the same `n_basis x n_basis`
penalty block. The total dimension is `1 + n_blocks * n_basis`.

A small ridge is added if any block is not positive-definite, to
guarantee numerical stability in downstream solves.

## References

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89–121.

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
# Single varying coefficient
P1 <- build_penalty_matrix(n_basis = 10, lambda = 1)
dim(P1)            # 11 x 11
#> [1] 11 11
P1[1, 1]           # 0 -- intercept is not penalised
#> [1] 0

# Three varying coefficients sharing the same penalty
P3 <- build_penalty_matrix(n_basis = 10, lambda = 1, n_blocks = 3)
dim(P3)            # 31 x 31  (1 + 3*10)
#> [1] 31 31
isSymmetric(P3)
#> [1] TRUE
```
