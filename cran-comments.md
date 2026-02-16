## Test environments
* local: Ubuntu Linux (x86_64), R 4.5.2
* GitHub Actions: macOS (release), Windows (release), Ubuntu (release, devel)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Andreas Bender <andreas.bender@stat.uni-muenchen.de>'
  This is a maintenance release addressing CRAN check NOTEs.

## Downstream dependencies
No known reverse dependencies with issues.

## Changes in this version
This release fixes the "no visible binding for global variable 'id'" NOTE
flagged in CRAN checks. It also fixes several bugs in competing risks data
transformation (factor status variables) and transition probability calculation,
and improves documentation.
