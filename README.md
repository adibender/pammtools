
<!-- badges: start -->

[![R-CMD-check](https://github.com/adibender/pammtools/workflows/R-CMD-check/badge.svg)](https://github.com/adibender/pammtools/actions)
[![cran
checks](https://badges.cranchecks.info/worst/pammtools.svg)](https://cran.r-project.org/web/checks/check_results_pammtools.html)
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![codecov.io](https://codecov.io/github/adibender/pammtools/coverage.svg?branch=master)](https://app.codecov.io/github/adibender/pammtools/branch/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-ago/pammtools)](https://cran.r-project.org/package=pammtools)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/pammtools)](https://cran.r-project.org/package=pammtools)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/license/mit)
<!-- badges: end -->

# **`pammtools`**: Piece-Wise Exponential Additive Mixed Modeling Tools

### Installation

Install from CRAN or GitHub using:

``` r
# CRAN
install.packages("pammtools")
# Development version
remotes::install_github("adibender/pammtools")
```

### Overview

**`pammtools`** facilitates the estimation of Piece-wise exponential
Additive Mixed Models (PAMMs) for time-to-event data. PAMMs can be
represented as generalized additive models and can therefore be
estimated using GAM software (e.g. **`mgcv`**), which, compared to other
packages for survival analysis, often offers more flexibility w.r.t. to
the specification of covariate effects (e.g. non-linear, time-varying
effects, cumulative effects, etc.). The package supports single-event
analysis, left-truncation, recurrent events, competing risks and
multi-state models.

To get started, see the
[Articles](https://adibender.github.io/pammtools/articles/) section.

<!-- Zaehlmarke VGWort -->

<img src="https://vg09.met.vgwort.de/na/f993b9c06b5249adb509a4df30d807a6" width="1" height="1" alt="">
