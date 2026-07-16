# Extract fixed-effects from a fitted model object

Generic in the style of
[`nlme::fixef`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) /
`lme4::fixef`, redefined here so cevcmm avoids a hard dependency on
either. If you also have nlme or lme4 loaded, call `cevcmm::fixef(fit)`
explicitly.

## Usage

``` r
fixef(object, ...)
```

## Arguments

- object:

  A model object.

- ...:

  Method-specific arguments.

## Value

The return value depends on the class of `object`; see the appropriate
method (e.g.
[`fixef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/fixef.vcmm_fit.md)
for `vcmm_fit` objects, which returns a two-element list with scalar
`intercept` and a matrix `varying` of B-spline basis coefficients).
