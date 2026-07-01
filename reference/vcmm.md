# Fit a Varying Coefficient Mixed-Effects Model

The main user-facing fit function. Builds the B-spline design and
penalty for one or more varying coefficients in `t`, computes the
aggregated sufficient statistics, and fits the model using either the
one-step CSL estimator (default) or the iterative SS estimator.

## Usage

``` r
vcmm(
  y,
  X,
  Z,
  t,
  method = c("auto", "csl", "ss"),
  re_cov = c("diag", "kronecker", "separable"),
  n_groups = NULL,
  q_left = NULL,
  Sigma_left_init = NULL,
  Sigma_right_init = NULL,
  Sigma_2x2_init = NULL,
  Sigma_spatial_init = NULL,
  Sigma_q_init = NULL,
  Omega_G_init = NULL,
  n_basis = NULL,
  degree = 3L,
  lambda = 1,
  control = vcmm_control(),
  normalize_t = TRUE,
  ...
)
```

## Arguments

- y:

  Numeric response vector of length \\n\\.

- X:

  Numeric \\n \times K\\ matrix (or length-\\n\\ vector if \\K = 1\\) of
  covariates that get varying coefficients in `t`.

- Z:

  Numeric \\n \times q\\ random-effects design matrix. For
  `re_cov = "kronecker"` or `"separable"`, \\q\\ must equal
  `q_left * n_groups`.

- t:

  Numeric vector of length \\n\\ in which the coefficients vary
  smoothly.

- method:

  Character: `"auto"` (default), `"csl"`, or `"ss"`.

- re_cov:

  Character: `"diag"` (default), `"kronecker"`, or `"separable"`. See
  Details.

- n_groups:

  Integer \\G\\. Required for `"kronecker"` and `"separable"`.

- q_left:

  Integer. The left (within) dimension of the Kronecker factor. For
  `re_cov = "kronecker"` defaults to 2 (OD setting); for
  `re_cov = "separable"` this is the per-group random-effect dimension
  and **must be supplied**.

- Sigma_left_init:

  Optional \\k \times k\\ initial left covariance (\\k =
  q\_{\mathrm{left}}\\). Aliases accepted: `Sigma_2x2_init` (for
  `"kronecker"` with `q_left = 2`) and `Sigma_q_init` (for
  `"separable"`). Default is `sigma_alpha^2 * I_k`.

- Sigma_right_init:

  Optional \\G \times G\\ initial right covariance. Aliases accepted:
  `Sigma_spatial_init` (for `"kronecker"`) and `Omega_G_init` (for
  `"separable"`). Default is \\I_G\\.

- Sigma_2x2_init:

  Legacy alias for `Sigma_left_init` when `re_cov = "kronecker"` and
  `q_left = 2`.

- Sigma_spatial_init:

  Legacy alias for `Sigma_right_init` when `re_cov = "kronecker"`.

- Sigma_q_init:

  Alias for `Sigma_left_init` when `re_cov = "separable"`.

- Omega_G_init:

  Alias for `Sigma_right_init` when `re_cov = "separable"`.

- n_basis:

  Integer or `NULL`. Number of B-spline basis functions per varying
  coefficient. `NULL` (default) auto-picks
  `max(floor(n^(1/3)) + 4, 10)`.

- degree:

  Integer. B-spline degree (default 3 = cubic).

- lambda:

  Non-negative numeric. Smoothing parameter (default 1).

- control:

  A `vcmm_control` object.

- normalize_t:

  Logical. If `TRUE` (default), `t` is linearly mapped to `[0, 1]`
  before building the basis.

- ...:

  Further arguments passed to
  [`fit_csl()`](https://lidajalili.github.io/cevcmm/reference/fit_csl.md)
  when CSL is used.

## Value

A `vcmm_fit` object as returned by
[`fit_ss()`](https://lidajalili.github.io/cevcmm/reference/fit_ss.md) or
[`fit_csl()`](https://lidajalili.github.io/cevcmm/reference/fit_csl.md),
augmented with `design` (basis metadata) and `re_cov`. When `re_cov` is
`"kronecker"` or `"separable"`, `re_cov_state` contains the canonical
fields `Sigma_left`, `Sigma_right`, plus legacy aliases
`Sigma_2x2`/`Sigma_spatial` (when `q_left = 2`) or `Sigma_q`/`Omega_G`
(for separable).

## Details

**Model.** For observations \\i = 1, \ldots, n\\ the fitted model is
\$\$ y_i = \beta_0(t_i) + \sum\_{k=1}^{K} x\_{ik}\\ \beta_k(t_i) +
z_i^\top \alpha + \varepsilon_i, \$\$ where each \\\beta_k(t)\\ is a
cubic B-spline with `n_basis` basis functions and a second-order
difference penalty, and \\\alpha \sim N(0, \Sigma\_\alpha)\\ with
structure chosen by `re_cov`.

**Method selection.** The default is `method = "auto"`, which picks
`"csl"` when \\N \cdot q \> 10^5\\ or \\q \> 50\\ and `"ss"` otherwise.

**Random-effects covariance.** Three structures:

- `re_cov = "diag"`: \\\alpha \sim N(0, \sigma\_\alpha^2 I_q)\\.

- `re_cov = "kronecker"`: \\\alpha \sim N(0, \Sigma\_{\mathrm{left}}
  \otimes \Sigma\_{\mathrm{right}})\\ with \\\Sigma\_{\mathrm{left}}\\
  of size `q_left x q_left` (default `q_left = 2`; OD-style with origin
  / destination blocks). User-facing names `Sigma_2x2_init`,
  `Sigma_spatial_init` are accepted as aliases.

- `re_cov = "separable"`: \\\alpha \sim N(0, \Sigma_q \otimes
  \Omega_G)\\ with `Sigma_q` of size `q_left x q_left` (required, no
  default) and `Omega_G` of size \\G \times G\\. User-facing names
  `Sigma_q_init`, `Omega_G_init` are accepted as aliases for
  `Sigma_left_init`, `Sigma_right_init`.

**Column-stacking convention.** For `re_cov = "kronecker"` or
`"separable"`, the random-effects vector is \\\alpha =
\mathrm{vec}\_{\mathrm{col}}(M)\\ where \\M \in \mathbb R^{G \times
q\_{\mathrm{left}}}\\, i.e.\\ `alpha[(k - 1) * G + g] = M[g, k]`. The
`Z` matrix must be constructed so that `Z %*% alpha` gives the correct
random-effect contribution. `ncol(Z)` must equal `q_left * n_groups`.

**Identifiability of the right component.** The right-side covariance
(`Sigma_spatial` / `Omega_G`) is not separately identifiable from a
single \\\hat\alpha\\, so it is held fixed at the user-supplied initial
value (default \\I_G\\). Supply a parametric kernel via
`Sigma_right_init` (e.g., `exp(-D / phi)` for OD; AR(1) for separable).
The left-side covariance (`Sigma_2x2` / `Sigma_q`) is updated every
iteration via the EM-style estimator (Theorem 1 \\M\_\eta\\ rule) when
`control$update_variance = TRUE`.

**Choosing `re_cov`.** The three modes specify the assumed covariance
structure of the random-effects vector \\\alpha\\; they do **not**
constrain what `Z`'s entries look like. The distribution of `Z` (binary
indicators, continuous values, mixed) does not enter this choice — only
the structure of \\\alpha\\ does. Pick by data shape:

- `"diag"`: independent random effects, one per group or subject, no
  assumed cross-group dependence. Simplest case. Use when each group
  contributes a single offset and groups are exchangeable.

- `"kronecker"`: origin-destination flow data. Every row of `Z`
  activates one origin indicator and one destination indicator, so
  `ncol(Z) = 2 * G` where \\G\\ is the number of regions. The 2x2 left
  block captures origin/destination dependence; the \\G \times G\\ right
  block captures spatial dependence between regions.

- `"separable"`: each group carries multiple correlated random effects
  (`q_left > 2`) — for example a random intercept plus a random slope
  per region. Use when the per-group random-effect dimension exceeds 2.

## References

Lin, L.-H. and Jalili, L. (2026). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.

## Examples

``` r
set.seed(1)
n <- 500
t  <- runif(n)
x  <- runif(n)
Z  <- matrix(rnorm(n * 3), n, 3)
alpha_true <- rnorm(3, sd = 0.5)
y  <- 2 + sin(2 * pi * t) * x +
      as.vector(Z %*% alpha_true) + rnorm(n, sd = 0.5)

fit <- vcmm(y, X = x, Z = Z, t = t,
            control = vcmm_control(sigma_eps = 0.5, sigma_alpha = 0.5))
fit
#> <vcmm_fit>  Varying Coefficient Mixed-Effects Model fit
#>   method      : SS
#>   n_obs       : 500
#>   p (fixed)   : 12
#>   q (random)  : 3
#>   RE cov      : diag
#>   iterations  : 6 (converged)
#>   sigma_eps   : 0.5000
#>   sigma_alpha : 0.5000
#>   elapsed     : <0.001 sec
```
