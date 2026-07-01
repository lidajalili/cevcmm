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
