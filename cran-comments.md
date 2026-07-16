## Resubmission

This is a resubmission of cevcmm 0.1.2 addressing the five items
raised by CRAN reviewer Konstanze Lauseker on cevcmm 0.1.1:

1. Expanded the "SVD" acronym in the DESCRIPTION as
   "Singular Value Decomposition".
2. Reformatted the Jalili and Lin (2025) reference in the
   CRAN-preferred autolinking form
   <https://arxiv.org/abs/2511.12732>.
3. Added `\value` sections to `fixef.Rd` and `ranef.Rd` (the S3
   generic function documentation), describing the return-value
   contract dispatched to method-specific documentation.
4. Changed `\dontrun{}` to `\donttest{}` in the `plot.vcmm_fit`
   example.
5. Excluded `inst/validation` and `inst/benchmarks` from the CRAN
   tarball via `.Rbuildignore`. These development-time scripts
   were the source of the `par()` / `options()` reset warnings
   and the `on.exit()` outside of a function. They remain in the
   GitHub repository for reference but are no longer shipped to
   CRAN users.

## Test environments

* Local: macOS aarch64 (Apple Silicon), R 4.5.2
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (release, devel, oldrel-1)
* win-builder: R-devel (Status: 1 NOTE, "New submission" only)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* NOTE flags "Jalili" (maintainer's surname) and "VCMMs"
  (paper's method acronym) as possibly misspelled.

## Downstream dependencies

None (new package).

## Additional notes for the reviewer

* The package implements the methodology of Jalili and Lin (2025),
  currently available as an arXiv preprint
  <https://arxiv.org/abs/2511.12732> and under review at the
  Journal of the American Statistical Association. The
  DESCRIPTION and citation entries reflect this status; we will
  update on acceptance.
* The bundled simulated dataset under
  `inst/extdata/od_migration.csv` is approximately 88 KB; we
  considered RDS but kept CSV for human-readable inspection and
  to keep the load path explicit in the OD-migration vignette.
* Vignettes are pre-built; total knit time on the local machine
  is well under 30 seconds.
