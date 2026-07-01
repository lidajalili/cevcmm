# Additive aggregation of sufficient-statistics summaries

Defines the `+` method for `vcmm_ss` objects, so that
`Reduce("+", summaries)` aggregates a list of node summaries
element-wise across the six sufficient-statistics components.

## Usage

``` r
# S3 method for class 'vcmm_ss'
e1 + e2
```

## Arguments

- e1, e2:

  `vcmm_ss` objects.

## Value

A `vcmm_ss` object holding the component-wise sums.

## Details

Dimensions (\\p\\, \\q\\) must match between operands; an informative
error is thrown otherwise.
