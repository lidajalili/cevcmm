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
