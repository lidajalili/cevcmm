# Extract random effects from a fitted model object

Generic in the style of
[`nlme::ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html) /
`lme4::ranef`, redefined here so cevcmm avoids a hard dependency on
either. If you also have nlme or lme4 loaded, call `cevcmm::ranef(fit)`
explicitly.

## Usage

``` r
ranef(object, ...)
```

## Arguments

- object:

  A model object.

- ...:

  Method-specific arguments.

## Value

The return value depends on the class of `object`; see the appropriate
method (e.g.
[`ranef.vcmm_fit`](https://lidajalili.github.io/cevcmm/reference/ranef.vcmm_fit.md)
for `vcmm_fit` objects, which returns a numeric vector, matrix, or list
reshaped to match the fitted `re_cov` structure).
