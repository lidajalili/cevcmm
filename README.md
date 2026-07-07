
<!-- README.md is generated from README.Rmd. Please edit that file, then run devtools::build_readme(). -->

# cevcmm

<!-- badges: start -->

[![R-CMD-check](https://github.com/lidajalili/cevcmm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lidajalili/cevcmm/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lidajalili/cevcmm/branch/main/graph/badge.svg)](https://app.codecov.io/gh/lidajalili/cevcmm?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/cevcmm)](https://CRAN.R-project.org/package=cevcmm)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**cevcmm** fits **Varying Coefficient Mixed-Effects Models (VCMMs)** at
scale, including settings where data is split across multiple machines
or doesn’t fit in memory. It implements the sufficient-statistics (SS)
and one-step communication-efficient surrogate-likelihood (CSL)
estimators of Jalili and Lin (2025), with diagonal, Kronecker, and
separable random-effect covariance structures.

The model is

$$y_i \;=\; \beta_0(t_i) \;+\; \sum_{k=1}^{K} x_{ik}\,\beta_k(t_i) \;+\; z_i^\top \alpha \;+\; \varepsilon_i,$$

where each $\beta_k(t)$ is a smooth function of $t$ (cubic B-splines
with a second-order-difference penalty),
$\alpha \sim N(0, \Sigma_\alpha)$, and
$\varepsilon_i \sim N(0, \sigma_\varepsilon^2)$.

## Installation

You can install the released version of cevcmm from CRAN with:

``` r
install.packages("cevcmm")
```

Or the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("lidajalili/cevcmm")
```

## Quick example

A 500-observation fit with one varying coefficient and three diagonal
random effects:

``` r
library(cevcmm)

set.seed(1)
N <- 500L
q <- 3L

t <- runif(N)
x <- runif(N)
Z <- matrix(rnorm(N * q), N, q)
y <- 2 + sin(2 * pi * t) * x +
     as.vector(Z %*% rnorm(q, sd = 0.4)) +
     rnorm(N, sd = 0.5)

fit <- vcmm(y, X = x, Z = Z, t = t,
            control = vcmm_control(sigma_eps       = 0.5,
                                   sigma_alpha     = 0.4,
                                   update_variance = TRUE))
fit
#> <vcmm_fit>  Varying Coefficient Mixed-Effects Model fit
#>   method      : SS
#>   n_obs       : 500
#>   p (fixed)   : 12
#>   q (random)  : 3
#>   RE cov      : diag
#>   iterations  : 6 (converged)
#>   sigma_eps   : 0.5076
#>   sigma_alpha : 0.6327
#>   elapsed     : 0.0010 sec
```

All standard S3 methods work on the result: `coef()`, `fixef()`,
`ranef()`, `vcov()`, `nobs()`, `logLik()`, `AIC()`, `BIC()`,
`predict()`, `summary()`, and `plot()`.

## Choosing the random-effects structure

The package supports three `re_cov` modes. They specify how the
random-effects vector $\alpha$ is assumed to covary — **not** what $Z$’s
entries look like. The distribution of $Z$ (binary indicators,
continuous values, or a mix) doesn’t enter the choice; only the
structure of $\alpha$ does.

| If your data looks like… | Use | Example |
|----|----|----|
| Independent random effects (one offset per group/subject, no cross-group dependence) | `re_cov = "diag"` | Patients in different clinics, classrooms in different schools |
| **Origin-destination flows** — each row has a “from” group and a “to” group | `re_cov = "kronecker"`, `q_left = 2` | Migration between regions, commuting between zones, trade between countries |
| Multiple correlated random effects per group (e.g., random intercept + random slope per region) | `re_cov = "separable"`, `q_left = number of effects per group` | Longitudinal data where each subject has a baseline and a trend, group-shared dense designs |

Each mode is demonstrated end-to-end in a vignette (links below).

## Three vignettes

Three vignettes ship with the package, covering the main use cases
end-to-end:

- **`vignette("getting-started", package = "cevcmm")`** — single-machine
  fit on a simple diagonal-random-effects example. Walks through every
  public S3 method.
- **`vignette("distributed-fitting", package = "cevcmm")`** — split data
  across $K$ nodes, compute per-node summaries with `node_summary()`,
  and recover a fit bit-equivalent to pooled `vcmm()` via
  `fit_from_summaries()`. Also covers the streaming-accumulator pattern
  for memory-constrained settings.
- **`vignette("od-migration", package = "cevcmm")`** — Kronecker
  covariance for origin-destination flow data, using a bundled simulated
  migration dataset.

## Citation

If you use cevcmm in published work, please cite:

> Jalili, L. and Lin, L.-H. (2025). *Scalable and
> Communication-Efficient Varying Coefficient Mixed Effect Models:
> Methodology, Theory, and Applications.* arXiv:2511.12732; under review
> at *Journal of the American Statistical Association*.

BibTeX:

``` bibtex
@article{jalili2025cevcmm,
  title   = {Scalable and Communication-Efficient Varying Coefficient
             Mixed Effect Models: Methodology, Theory, and Applications},
  author  = {Jalili, Lida and Lin, Li-Hsiang},
  journal = {arXiv preprint arXiv:2511.12732},
  year    = {2025},
  note    = {Under review at Journal of the American Statistical Association}
}
```

## License

Released under the GPL (\>= 3) license.
