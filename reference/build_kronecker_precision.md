# Build the Kronecker-structured random-effect precision matrix

For \\\Sigma\_\alpha = \Sigma\_{\mathrm{left}} \otimes
\Sigma\_{\mathrm{right}}\\ (column-stacking convention), returns \$\$
\sigma\_\varepsilon^2 \cdot (\Sigma\_{\mathrm{left}}^{-1} \otimes
\Sigma\_{\mathrm{right}}^{-1}) \$\$ which is added to `crossprod(Z)`
inside the random-effect block of the VCMM Hessian. This is the
structured analogue of `(sigma_eps^2 / sigma_alpha^2) * I_q` used under
`re_cov = "diag"`.

## Usage

``` r
build_kronecker_precision(Sigma_left, Sigma_right, sigma_eps)
```

## Arguments

- Sigma_left:

  A \\k \times k\\ positive-definite numeric matrix (\\k = 2\\ for
  OD-style; arbitrary \\k\\ for general separable).

- Sigma_right:

  A \\G \times G\\ positive-definite numeric matrix.

- sigma_eps:

  Positive numeric. Residual standard deviation.

## Value

A \\kG \times kG\\ numeric matrix.

## References

Jalili, L. and Lin, L.-H. (2025). Scalable and Communication-Efficient
Varying Coefficient Mixed-Effects Models.
