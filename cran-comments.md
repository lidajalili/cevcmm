## Resubmission

This is a maintenance release of cevcmm 0.1.4 addressing the
gcc-UBSAN "additional issue" reported for cevcmm 0.1.3 on the CRAN
Fedora tests-gcc-SAN farm (2026-07-24):

  compute_sufficient_stats_cpp.cpp:54:14: runtime error:
    null pointer passed as argument 1, which is declared to
    never be null

## Root cause

`src/compute_sufficient_stats_cpp.cpp` copied the input response
vector `y` into an `n x 1` `arma::mat` via
`std::memcpy(y_mat.memptr(), y.memptr(), n * sizeof(double))`.

When `n == 0` (the empty-input edge case exercised by
`tests/testthat/test-defensive-branches.R`), an empty `arma::vec`
returns `NULL` from `.memptr()`. Passing `NULL` to `std::memcpy`
is undefined behavior per the C standard even when the byte count
is zero, because `std::memcpy` carries the `nonnull` attribute on
both pointers. gcc-UBSAN correctly flags this. Most C libraries
short-circuit the zero-byte case before touching the pointers, so
the issue is latent on most platforms (all six main CRAN check
flavors returned OK on 0.1.3) but observable under UBSAN
instrumentation.

## Fix

Replace the `std::memcpy` call with Armadillo's `arma::mat(arma::vec)`
copy constructor, which handles the `n == 0` case internally without
dereferencing any pointer and produces the same `n x 1` matrix for
all `n > 0`. No user-visible behavior change and no numerical
change; the fix is purely a memory-safety hardening.

## Test environments

* Local: macOS aarch64 (Apple Silicon), R 4.5.2
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (release, devel, oldrel-1)
* win-builder: R-devel (Status: 1 NOTE, unchanged from 0.1.3)

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE flags "Jalili" (maintainer's surname) and "VCMMs"
  (paper's method acronym) as possibly misspelled. Unchanged
  from 0.1.3.

## Downstream dependencies

None.
