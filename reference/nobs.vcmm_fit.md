# Number of observations from a vcmm fit

Standard [`stats::nobs`](https://rdrr.io/r/stats/nobs.html) generic.
Returns \\N\\, the total number of observations used to compute the fit
(or the sum of node sample sizes for fits produced via
[`fit_from_summaries`](https://lidajalili.github.io/cevcmm/reference/fit_from_summaries.md)).

## Usage

``` r
# S3 method for class 'vcmm_fit'
nobs(object, ...)
```

## Arguments

- object:

  A `vcmm_fit` object.

- ...:

  Unused.

## Value

Integer.
