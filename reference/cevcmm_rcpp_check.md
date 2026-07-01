# Verify the Rcpp / RcppArmadillo toolchain is wired up

Internal stub used by Day 15's profiling script to confirm the C++
compiler, Rcpp dispatch, and Armadillo linking all work before the real
ports on Days 16-18. Performs a trivial 3x3 identity-matrix trace to
exercise the linker.

## Usage

``` r
cevcmm_rcpp_check()
```

## Value

A character string of the form "OK (Armadillo X.Y.Z; trace(I_3) = 3)".
