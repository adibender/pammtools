[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://travis-ci.org/adibender/pammtools.svg?branch=master)](https://travis-ci.org/adibender/pammtools)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/adibender/pammtools?branch=master&svg=true)](https://ci.appveyor.com/project/adibender/pammtools/branch/master)
[![codecov.io](https://codecov.io/github/adibender/pammtools/coverage.svg?branch=master)](https://codecov.io/github/adibender/pammtools/branch/master)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version-ago/pammtools)](https://cran.r-project.org/package=pammtools)
[![CRAN\_Download\_Badge](https://cranlogs.r-pkg.org/badges/pammtools)](https://cran.r-project.org/package=pammtools)
[![MIT license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)

## `pammtools`: Piece-Wise Exponential Additive Mixed Modeling Tools <img src="man/figures/logo.png" align="right" width="150px" />


### Installation

Install from CRAN or GitHub using:

```r
# CRAN
install.packages("pammtools")

# GitHub
devtools::install_github("adibender/pammtools")
```



### Overview

**`pammtools`** facilitates the estimation of Piece-wise exponential Additive Mixed Models (PAMMs) for time-to-event data. PAMMs can be represented as generalized additive models and can therefore be estimated using GAM software (e.g. **`mgcv`**), which, compared to other packages for survival analysis, often offers more flexibility w.r.t. to the specification of covariate effects (e.g. non-linear, time-varying effects, cumulative effects, etc.).

To get started, see the [Articles](https://adibender.github.io/pammtools/articles/) section.


<!-- An overview over the packages functionality is given in

- Andreas Bender and Fabian Scheipl, 2018: "pammtools: Piece-wise exponential
Additive Mixed Modeling tools", arXiv eprint, 2018, https://arxiv.org/abs/1806.01042

For a tutorial-like introduction to PAMMs see:

  - Andreas Bender, Andreas Groll, and Fabian Scheipl. 2018. “A Generalized Additive Model Approach to Time-to-Event Analysis.” Statistical Modelling. https://doi.org/10.1177/1471082X17748083.


A general framework for the representation and estimation of cumulative effects
(or exposure-lag-response associations) is described in:

- Andreas Bender, Fabian Scheipl, Wolfgang Hartl, Andrew G Day, Helmut Küchenhoff; "Penalized estimation of complex, non-linear exposure-lag-response associations", Biostatistics, , kxy003, https://doi.org/10.1093/biostatistics/kxy003
 -->
