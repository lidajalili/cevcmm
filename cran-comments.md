## Test environments

* Local: macOS aarch64 (Apple Silicon), R 4.5.2
* GitHub Actions:
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (release, devel, oldrel-1)
* win-builder (release, devel)
* R-hub (linux, macos, windows)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Downstream dependencies

There are currently no downstream reverse dependencies for this
package.

## Additional notes for the reviewer

* The package implements the methodology of Lin and Jalili (2026),
  forthcoming in JASA. The reference is marked `(under review)` in
  the DESCRIPTION because the paper is still in peer review at
  submission time; we will update the DESCRIPTION upon acceptance.
* The bundled simulated dataset under `inst/extdata/od_migration.csv`
  is approximately 88 KB; we considered RDS but kept CSV for
  human-readable inspection and to keep the load path explicit in
  the OD-migration vignette.
* Vignettes are pre-built; total knit time on the local machine is
  well under 30 seconds.
