# Fit a VCMM from aggregated sufficient-statistics summaries

Server-side fit using only the small per-node summaries; no raw data
required. The mathematical guarantee (Theorem 1 of Lin and Jalili, 2026)
is that for a fixed partition of the data into nodes, the fit obtained
here is identical to the one a centralised
[`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) call
would produce on the pooled data, up to floating-point summation noise.

## Usage

``` r
fit_from_summaries(
  summaries,
  penalty,
  control = vcmm_control(),
  method = c("csl", "ss"),
  re_cov = c("diag", "kronecker", "separable"),
  n_groups = NULL,
  q_left = NULL,
  Sigma_left_init = NULL,
  Sigma_right_init = NULL,
  Sigma_2x2_init = NULL,
  Sigma_spatial_init = NULL,
  Sigma_q_init = NULL,
  Omega_G_init = NULL,
  rowsum_constant = NULL,
  ...
)
```

## Arguments

- summaries:

  Either a list of `vcmm_ss` objects (one per node, the typical case) OR
  a single `vcmm_ss` that has already been aggregated (e.g., via
  `Reduce("+", ...)`) OR a `vcmm_accumulator`.

- penalty:

  The \\p \times p\\ smoothing-penalty matrix used in the spline basis
  the nodes shared. Build it via the package's
  [`build_vcmm_design()`](https://lidajalili.github.io/cevcmm/reference/build_vcmm_design.md)
  (returned as `design$penalty`) or directly via
  [`build_penalty_matrix()`](https://lidajalili.github.io/cevcmm/reference/build_penalty_matrix.md).

- control:

  A
  [`vcmm_control`](https://lidajalili.github.io/cevcmm/reference/vcmm_control.md)
  object.

- method:

  Either `"csl"` (default) or `"ss"`.

- re_cov:

  Either `"diag"`, `"kronecker"`, or `"separable"`. See
  [`vcmm`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) for
  full semantics.

- n_groups:

  Required if `re_cov` is `"kronecker"` or `"separable"`.

- q_left:

  Required for `"separable"`; defaults to 2 for `"kronecker"`.

- Sigma_left_init, Sigma_right_init, Sigma_2x2_init, Sigma_spatial_init,
  Sigma_q_init, Omega_G_init:

  Same aliases accepted as in
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md).

- rowsum_constant:

  Optional numeric. If supplied and non-zero, apply the same
  identifiability re-centering that
  [`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md)
  applies post-fit. See Details.

- ...:

  Passed to
  [`fit_csl`](https://lidajalili.github.io/cevcmm/reference/fit_csl.md).

## Value

A `vcmm_fit` object.

## Details

**Identifiability re-centering.**
[`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md)
applies a post-fit re-centering when `rowSums(Z)` is constant
(indicator-Z OD designs and similar), absorbing
`rowsum_constant * mean(alpha_hat)` into `beta_0`. The server has no
access to `Z`, so the caller must signal that the original `Z` had
constant row sums by passing `rowsum_constant`. Default `NULL` means "no
re-centering" (appropriate for dense-Z and block-Z designs).

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## See also

[`node_summary`](https://lidajalili.github.io/cevcmm/reference/node_summary.md),
[`vcmm`](https://lidajalili.github.io/cevcmm/reference/vcmm.md).

## Examples

``` r
set.seed(1)
n <- 600
t  <- runif(n); x <- runif(n); Z <- matrix(rnorm(n * 4), n, 4)
a_true <- rnorm(4, sd = 0.3)
y <- 2 + sin(2 * pi * t) * x + as.vector(Z %*% a_true) + rnorm(n, sd = 0.5)

design <- build_vcmm_design(X = x, t = t)
Xd     <- design$X_design
ctrl   <- vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.3)

# Split into 3 nodes and aggregate
idx_node <- sample.int(3, n, replace = TRUE)
summaries <- lapply(seq_len(3), function(s) {
  ii <- which(idx_node == s)
  node_summary(y[ii], Xd[ii, , drop = FALSE], Z[ii, , drop = FALSE])
})
fit <- fit_from_summaries(summaries,
                           penalty = design$penalty,
                           control = ctrl)
```
