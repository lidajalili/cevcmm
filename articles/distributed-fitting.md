# Distributed fitting with cevcmm

## When you need this

Three scenarios call for the distributed API:

1.  **Data lives on $`K`$ different machines** — privacy or governance
    constraints prevent pooling, but each site can publish a small
    summary.
2.  **Data is on one machine but doesn’t fit in memory** — process chunk
    by chunk and accumulate.
3.  **Streaming data** — new batches arrive over time and refitting from
    scratch is expensive.

All three are solved with the same primitive: each chunk produces a
`vcmm_ss` summary; summaries are additive; the fit is computed from the
aggregate.

## The math guarantee

Theorem 1 of Jalili and Lin (2025): partition the data into $`K`$
chunks, compute a sufficient-statistics summary on each chunk, sum the
summaries, and fit. The result is **bit-equivalent** (up to BLAS
summation noise) to fitting on the pooled data. The summaries lose no
information — they *are* the sufficient statistics for the VCMM
likelihood.

The summary on chunk $`k`$ stores
``` math
\bigl(a_k,\; b_k,\; C_k,\; (Z^\top Z)_k,\; (Z^\top y)_k,\; (X^\top Z)_k,\; n_k\bigr)
```
where $`a = y^\top y`$, $`b = X^\top y`$, $`C = X^\top X`$. Adding
summaries across chunks is element-wise addition; the aggregate is the
same as the single summary you would have computed on the pooled data.

## Setup

``` r

library(cevcmm)
```

## Simulate

A modest example: $`N = 1200`$ observations, $`K = 1`$ varying
covariate, $`q = 4`$ independent random effects. Then we’ll split it
into 3 nodes.

``` r

set.seed(23)
N <- 1200L
q <- 4L

t <- runif(N)
x <- runif(N)
Z <- matrix(rnorm(N * q), N, q)

beta_1_fun <- function(u) sin(2 * pi * u)
alpha_true <- rnorm(q, sd = 0.5)
y <- 2 + beta_1_fun(t) * x + as.vector(Z %*% alpha_true) +
     rnorm(N, sd = 0.5)
```

## Pattern 1 — list of node summaries

This is the “data on $`K`$ machines” case. Each node:

1.  computes the same spline design matrix on its own slice of $`t`$ and
    $`X`$,
2.  calls
    [`node_summary()`](https://lidajalili.github.io/cevcmm/reference/node_summary.md)
    on its slice,
3.  ships the resulting `vcmm_ss` object (small — a few KB regardless of
    $`n_k`$) to the central node.

The central node sums the summaries and calls
[`fit_from_summaries()`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md).

``` r

# Assign each observation to one of 3 nodes
node_id <- sample.int(3L, N, replace = TRUE)
splits  <- split(seq_len(N), node_id)
lengths(splits)
#>   1   2   3 
#> 405 401 394
```

Build the spline design once. (In a real deployment, each node would
construct its piece of the design locally using the same `n_basis`,
`degree`, and knot placement — the package handles this automatically as
long as the inputs agree.)

``` r

design <- build_vcmm_design(X = x, t = t)
X_d    <- design$X_design
# N rows x p columns; p = intercept + spline-basis coefficients,
# with n_basis auto-picked from the sample size.
dim(X_d)
#> [1] 1200   15
```

Each node produces a summary:

``` r

summaries <- lapply(splits, function(idx) {
  node_summary(y[idx], X_d[idx, , drop = FALSE], Z[idx, , drop = FALSE])
})
summaries[[1]]
#> <vcmm_ss>  one-batch sufficient statistics
#>   n_obs : 405
#>   p     : 15   (fixed-effects columns)
#>   q     : 4   (random-effects columns)
#>   a     : 2461
```

The size of each summary is independent of $`n_k`$:

``` r

object.size(summaries[[1]])
#> 4776 bytes
object.size(y[splits[[1]]]) + object.size(X_d[splits[[1]], ])
#> 52104 bytes
```

Now fit from the summaries — pass the list directly:

``` r

ctrl <- vcmm_control(sigma_eps       = 0.5,
                     sigma_alpha     = 0.5,
                     update_variance = TRUE)

fit_dist <- fit_from_summaries(summaries,
                               penalty = design$penalty,
                               control = ctrl,
                               method  = "ss",
                               re_cov  = "diag")
fit_dist
#> <vcmm_fit>  Varying Coefficient Mixed-Effects Model fit
#>   method      : SS
#>   n_obs       : 1200
#>   p (fixed)   : 15
#>   q (random)  : 4
#>   RE cov      : diag
#>   iterations  : 5 (converged)
#>   sigma_eps   : 0.4845
#>   sigma_alpha : 0.7251
#>   elapsed     : 0.0010 sec
```

## Verify bit-equivalence against pooled fitting

``` r

fit_pooled <- vcmm(y, X = x, Z = Z, t = t,
                   method = "ss", re_cov = "diag",
                   control = ctrl)

# beta estimates
max(abs(fit_pooled$beta - fit_dist$beta))
#> [1] 1.065814e-14

# alpha estimates
max(abs(fit_pooled$alpha - fit_dist$alpha))
#> [1] 4.440892e-16

# variance components
all.equal(fit_pooled$sigma_eps,   fit_dist$sigma_eps)
#> [1] TRUE
all.equal(fit_pooled$sigma_alpha, fit_dist$sigma_alpha)
#> [1] TRUE
```

The maximum coefficient difference sits at the BLAS summation-noise
floor — typically near machine epsilon ($`\sim 10^{-15}`$). The
distributed fit *is* the pooled fit; only the path that produced it
differs.

## Pattern 2 — streaming accumulator

When you have many small chunks and don’t want to keep them all in
memory, use the accumulator instead of a list. The accumulator stores a
running sum of the sufficient statistics; chunks can be incorporated and
discarded as they arrive.

``` r

p <- ncol(X_d)
acc <- init_accumulator(p = p, q = q)

# Pretend chunks arrive one at a time
for (k in seq_along(splits)) {
  idx <- splits[[k]]
  ss_k <- compute_sufficient_stats(y[idx],
                                   X_d[idx, , drop = FALSE],
                                   Z[idx, , drop = FALSE])
  acc <- accumulate_stats(acc, ss_k)
}
acc
#> <vcmm_accumulator>  accumulated sufficient statistics
#>   n_obs : 1200
#>   p     : 15   (fixed-effects columns)
#>   q     : 4   (random-effects columns)
#>   a     : 7587
```

Fit from the accumulator — same call signature:

``` r

fit_acc <- fit_from_summaries(acc, penalty = design$penalty,
                              control = ctrl, method = "ss", re_cov = "diag")

# Same answer as Pattern 1
max(abs(fit_acc$beta - fit_dist$beta))
#> [1] 0
```

## Identifiability for OD designs

When `rowSums(Z)` is constant — for example, a row of $`Z`$ has exactly
one “origin” indicator and one “destination” indicator, so each row sums
to $`s = 2`$ — the model has a one-dimensional identifiability
redundancy. For any scalar $`c`$, replacing
``` math
\alpha \;\mapsto\; \alpha + c\mathbf{1}, \qquad
\beta_0 \;\mapsto\; \beta_0 - c\,s
```
leaves $`Z\alpha + \beta_0\mathbf{1}`$ — and hence the likelihood —
unchanged.
[`vcmm()`](https://lidajalili.github.io/cevcmm/reference/vcmm.md)
detects this automatically by inspecting the raw $`Z`$ matrix and
applies a centering shift so that $`\sum_j \hat\alpha_j = 0`$.

[`fit_from_summaries()`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md)
doesn’t auto-detect the row-sum from the sufficient statistics, so you
supply it explicitly via `rowsum_constant`:

``` r

# Skeleton — see vignette("od-migration") for a complete runnable example.
# G and Sigma_spatial come from the OD setup (number of regions and an
# initial spatial covariance matrix).
fit_dist_od <- fit_from_summaries(
  summaries, penalty = design$penalty, control = ctrl,
  method = "ss", re_cov = "kronecker",
  n_groups           = G,              # number of OD regions
  Sigma_spatial_init = Sigma_spatial,  # G x G initial spatial covariance
  rowsum_constant    = 2               # matches what vcmm() applies internally
)
```

Without `rowsum_constant`, the distributed fit’s $`\hat\beta_0`$ and
$`\hat\alpha`$ differ from the pooled fit’s by exactly the centering
shift. The fits are the same up to the redundancy; only the choice of
representative within the equivalence class differs.

See
[`vignette("od-migration", package = "cevcmm")`](https://lidajalili.github.io/cevcmm/articles/od-migration.md)
for a full OD example.

## Where to go next

- **OD migration with Kronecker covariance**:
  [`vignette("od-migration", package = "cevcmm")`](https://lidajalili.github.io/cevcmm/articles/od-migration.md).
- **Basic single-machine usage**:
  [`vignette("getting-started", package = "cevcmm")`](https://lidajalili.github.io/cevcmm/articles/getting-started.md).

## Reference

Jalili, L. and Lin, L.-H. (2025). *Scalable and Communication-Efficient
Varying Coefficient Mixed Effect Models: Methodology, Theory, and
Applications.* arXiv:2511.12732; under review at *Journal of the
American Statistical Association*.
