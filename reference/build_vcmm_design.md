# Build the design matrix and penalty for a VCMM

Constructs the fixed-effects design matrix and matching penalty for a
VCMM with one or more varying coefficients. Given covariates `X` (n by
K) and a time/index vector `t` (length n), builds the cubic-by-default
B-spline basis `B` at `t` and assembles

- `X_design = cbind(1, B * X[,1], B * X[,2], ..., B * X[,K])`, dimension
  n by (1 + K \* n_basis); the first column is the intercept, then each
  covariate gets its own n_basis-column block.

- `penalty`: a (1 + K \* n_basis) by (1 + K \* n_basis) block-diagonal
  second-order difference penalty matrix with the intercept unpenalised.

This is the natural multi-covariate generalisation of the
single-covariate design used in the simulation code, matching the
paper's general VCMM specification.

## Usage

``` r
build_vcmm_design(
  X,
  t,
  n_basis = NULL,
  degree = 3L,
  lambda = 1,
  normalize_t = TRUE
)
```

## Arguments

- X:

  Numeric n by K matrix (or length-n vector if K = 1) of covariates that
  get varying coefficients in `t`.

- t:

  Numeric vector of length n. The variable in which the coefficients
  vary smoothly (time, index, location, etc.).

- n_basis:

  Integer (\\\geq\\ `degree` + 1) or `NULL`. Number of B-spline basis
  columns per varying coefficient. If `NULL` (default), chosen as
  `max(floor(n^(1/3)) + 4, 10)` matching the simulation code's rule of
  thumb.

- degree:

  Integer. B-spline degree (default 3 = cubic).

- lambda:

  Non-negative numeric. Smoothing parameter for the penalty (default 1).

- normalize_t:

  Logical. If `TRUE` (default), `t` is linearly mapped to `[0, 1]`
  before building the basis. Set `FALSE` only if `t` is already in
  `[0, 1]`.

## Value

A list with elements:

- `X_design`: numeric n by (1 + K \* n_basis) design matrix.

- `penalty`: numeric (1 + K \* n_basis) by (1 + K \* n_basis) penalty
  matrix.

- `B_spline`: numeric n by n_basis basis matrix at `t`.

- `internal_knots`, `boundary_knots`: spline knot locations, recorded so
  predictions at new `t` can use the same basis.

- `degree`, `n_basis`, `K`, `lambda`: scalars recording how the design
  was built.

- `normalize_t`, `t_min`, `t_max`: how to remap new `t` values at
  prediction time.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
n <- 200
t <- runif(n)
x1 <- runif(n); x2 <- runif(n)

# Single varying coefficient
d1 <- build_vcmm_design(X = x1, t = t, n_basis = 8)
dim(d1$X_design)   # 200 x 9   (1 + 1*8)
#> [1] 200   9

# Two varying coefficients
d2 <- build_vcmm_design(X = cbind(x1, x2), t = t, n_basis = 8)
dim(d2$X_design)   # 200 x 17  (1 + 2*8)
#> [1] 200  17
```
