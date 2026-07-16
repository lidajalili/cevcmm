## Resubmission

This is a minor metadata-only resubmission of cevcmm 0.1.3 following
a NOTE from win-builder R-devel on cevcmm 0.1.2 asking that arXiv
preprints be cited via their arXiv DOI form
(`<doi:10.48550/arXiv.YYMM.NNNNN>`) rather than the plain URL form
(`<https://arxiv.org/abs/YYMM.NNNNN>`). Only the DESCRIPTION
reference was changed; no code or documentation changes.

All five items from CRAN reviewer Konstanze Lauseker's earlier
feedback on cevcmm 0.1.1 remain addressed as in 0.1.2:

1. "SVD" expanded to "Singular Value Decomposition" in DESCRIPTION.
2. Paper reference reformatted for autolinking (now in DOI form).
3. `\value` sections added to `fixef.Rd` and `ranef.Rd`.
4. `\dontrun{}` changed to `\donttest{}` in `plot.vcmm_fit`.
5. `inst/validation` and `inst/benchmarks` excluded from the CRAN
   tarball via `.Rbuildignore`.

## Test environments

* Local: macOS aarch64 (Apple Silicon), R 4.5.2
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (release, devel, oldrel-1)
* win-builder: R-devel (Status: 1 NOTE, "New submission" only,
  after the DOI change)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* NOTE flags "Jalili" (maintainer's surname) and "VCMMs"
  (paper's method acronym) as possibly misspelled.

## Downstream dependencies

None (new package).

## Additional notes for the reviewer

* The paper by Jalili and Lin (2025) is currently available as an
  arXiv preprint (arXiv:2511.12732, DOI 10.48550/arXiv.2511.12732)
  and under review at the Journal of the American Statistical
  Association. The DESCRIPTION reference and package citation
  will be updated on journal acceptance.
* The bundled simulated dataset under
  `inst/extdata/od_migration.csv` is approximately 88 KB; we
  considered RDS but kept CSV for human-readable inspection and
  to keep the load path explicit in the OD-migration vignette.
* Vignettes are pre-built; total knit time on the local machine
  is well under 30 seconds.
