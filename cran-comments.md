## Resubmission

This is a resubmission of cevcmm 0.1.1 following the CRAN automated
pretest of cevcmm 0.1.0, which failed to install on Debian with
`undefined symbol: dpotrf_`. The Windows pretest passed with just
the harmless "New submission" NOTE.

The fix adds explicit `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS)
$(FLIBS)` to `src/Makevars` for BLAS/LAPACK/Fortran linking on
Linux. The equivalent flags were already present in
`src/Makevars.win` for Windows.

Verified locally on macOS with an updated gfortran install; all
236 tests pass.

## Test environments

* Local: macOS aarch64 (Apple Silicon), R 4.5.2
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (release, devel, oldrel-1)
* win-builder: R-devel and R-release (Status: 1 NOTE, "New
  submission" only)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* NOTE flags "Jalili" (maintainer's surname), "Scalable" (common
  word), and "VCMMs" (paper's method acronym) as possibly
  misspelled.

## Downstream dependencies

None (new package).

## Additional notes for the reviewer

* The package implements the methodology of Jalili and Lin (2025),
  currently available as an arXiv preprint (arXiv:2511.12732) and
  under review at the Journal of the American Statistical
  Association. The DESCRIPTION and citation entries reflect this
  status; we will update on acceptance.
* The bundled simulated dataset under `inst/extdata/od_migration.csv`
  is approximately 88 KB; we considered RDS but kept CSV for
  human-readable inspection and to keep the load path explicit in
  the OD-migration vignette.
* Vignettes are pre-built; total knit time on the local machine is
  well under 30 seconds.