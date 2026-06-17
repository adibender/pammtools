# pammtools 0.8.0

This is a feature release. It supersedes the GitHub-only 0.7.5 and 0.7.6
development versions; the CRAN-published version is 0.7.4.

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked locally with `R CMD check --as-cran` (including `--run-donttest` and the
PDF manual). The only local NOTE is "HTML Tidy not found", which is specific to
the local machine (no `tidy` binary installed) and does not occur on CRAN's
check infrastructure.

## Test environments

* local: Linux Mint 22.1 (Ubuntu 24.04 base), R 4.6.0
* (recommended before upload: win-builder R-devel/R-release and R-hub
  linux / macos / windows / nosuggests)

## Major changes since 0.7.4

* New `scam` estimation engine for shape-constrained PAMMs (`pamm(engine =
  "scam")`).
* Interval-censored data support via a multiple-imputation workflow
  (`pamm_ic()`, `pamm_ic_cr()`, `add_inspections()`).
* The post-processing / confidence-interval machinery is now backend-pluggable
  (internal `get_hazard()`/`sim_hazard()` S3 seam), demonstrated by a new
  gradient-boosting (xgboost) vignette.
* Breaking change: `make_newdata()` output no longer contains internal PED
  columns; `add_cif()` now uses an exact closed-form CIF integral; simulation
  CIs now use type-6 quantiles. See NEWS.md for the full, categorized list.

## Reverse dependencies

We checked all 7 reverse dependencies (adjustedCurves, autoReg, contsurvplot,
KMunicate, relsurv, robber, weightedsurv; pammtools is a Suggested dependency
for all), comparing R CMD check results across the CRAN (0.7.4) and dev (0.8.0)
versions with `revdepcheck::revdep_check()`.

* 0 new problems
* 0 packages failed to check
