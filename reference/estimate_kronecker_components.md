# Moment-based estimator for the Kronecker left component

Given a length-\\kG\\ random-effects estimate \\\hat\alpha =
\mathrm{vec}\_{\mathrm{col}}(M)\\ (column-stacked, \\M \in \mathbb R^{G
\times k}\\) and the right-side covariance `Sigma_right`, returns the
weighted moment estimate \$\$\hat\Sigma\_{\mathrm{left}} =
\tfrac{1}{G}\\ M^{\top}\\ \Sigma\_{\mathrm{right}}^{-1}\\ M.\$\$

## Usage

``` r
estimate_kronecker_components(alpha, n_groups, q_left = 2L, Sigma_right = NULL)
```

## Arguments

- alpha:

  Numeric vector of length \\kG\\.

- n_groups:

  Integer \\G\\.

- q_left:

  Integer \\k\\, the left (within) dimension. Defaults to 2 for backward
  compatibility with the OD setting.

- Sigma_right:

  Optional \\G \times G\\ positive-definite covariance. If `NULL`, the
  unweighted sample covariance `cov(alpha_mat)` is returned (unbiased
  only when rows are uncorrelated).

## Value

A \\k \times k\\ symmetric positive-definite matrix.

## Details

This is unbiased under the Kronecker model \\\alpha \sim N(0,
\Sigma\_{\mathrm{left}} \otimes \Sigma\_{\mathrm{right}})\\ when
`Sigma_right` is correct. When \\\hat\alpha\\ is the BLUP rather than
the true \\\alpha\\, apply the EM-style correction in
[`vcmm`](https://lidajalili.github.io/cevcmm/reference/vcmm.md) (handled
automatically by `fit_ss` / `fit_csl`).

Backwards-compatible alias: if `q_left = 2` (the default), this is the
same estimator as the previous `estimate_kronecker_components` for the
OD setting.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.
