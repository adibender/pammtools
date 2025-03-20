# pamtools 0.5.93
+ Maintnance (some tidyverse deprecations, link fixes, etc., smaller bugs)

# pammtools 0.5.92
+ Fixed competing risks data trafo in case of more than 2 causes

# pammtools 0.5.9
+ Fixes issue 154: direction argument to `geom_stepribbon`

# pammtools 0.5.8
+ removed argument `methods` from `pamm`. Can be specified via `...`. Fixes #200
+ adapted `warn_about_new_time_points` when original data not stored in model object. Fixes #203
+ Fixed issue where not all ped attributes were retained when applying dplyr functions #202

# pammtools 0.5.7
+ added staph data with recurrent events

# pammtools 0.5.6
+ maintenance fix
+ fixes to URLs and DOIs


# pammtools 0.5.4
+ updates to the `split_data` function that now accepts `Surv(start, stop, event)` type inputs, e.g., to construct left-truncated data.
+ Support and [vignette for left truncated data](https://adibender.github.io/pammtools/articles/left-truncation.html)
+ Support and [vignette for competing risks data](https://adibender.github.io/pammtools/articles/competing-risks.html)
+ Support and [vignette for recurrent events data](https://adibender.github.io/pammtools/articles/recurrent-events.html)


# pammtools 0.2.4
* CRAN fix. Discrepancy between man page and code.

# pammtools 0.2.3
* CRAN fix. Compliance with new dplyr version (1.0.0)

# pammtools 0.2.2
* CRAN fix, removed plyr dependency (see issue #141)
* `as_ped.ped` now also works for transformations with time-dependent covariates

# pammtools 0.2.1
* Adds a new interface for model estimation called `pamm`, which is a thin wrapper
around `mgcv::gam` with some arguments pre-set.
* Adds S3 method `predictSurvProb.pamm`
* Adds support and vignette for model evaluation using package **`pec`**
* Fixed bug when CIs were calculated simulation based and model contained factor variables
* Removed unnecessary dependencies in Imports/Suggests

# pammtools 0.1.15
* Interface for specification of data transformation in `as_ped` changed. The vertical bar `|` is no longer necessary to indicate concurrent or cumulative effects

# pammtools 0.1.14

* Support for new interface to tidyr

# pammtools 0.1.13

* Functions `get_hazard` and `add_hazard` also gain `reference` argument.
Allows to calculate (log-)hazard ratios.

* Introduces breaking changes to `add_term` function. Argument `relative` is replaced by `reference`, makes calculation of relative (log-)hazards, i.e. hazard ratios, more flexible. Argument `se.fit` is replaced by `ci`.



# pammtools 0.1.11

## bugs
* fixes bug in **`dplyr`** reverse dependency (see #101)
* fixes bug in tidiers for Aalen models (see #99)

## documentation
* Better documentation and functionality for `make_newdata`
* Added new vignette linking to tutorial paper (online only)

# pammtools 0.1.9
* maintenance update: fixes CRAN issues due to new RNG

# pammtools 0.1.8

## documentation
* Updates to cumulative effect vignette
* Updates to time-dependent covariate vignette (+ data transformation)
* Update citation information

## Features
* `concurrent` now has a `lag = 0` argument, can be set to positive integer values
* `as_ped` accepts multiple `concurrent` specials with different `lag` specifications

## Bug/Issue fixes
* Fixed bug caused by changes in **`checkmate`** [#73](https://github.com/adibender/pammtools/issues/73)
* Bug Fixes [#42](https://github.com/adibender/pammtools/issues/42), [#76](https://github.com/adibender/pammtools/issues/76), [#63](https://github.com/adibender/pammtools/issues/63), [#77](https://github.com/adibender/pammtools/issues/77)


# pammtools 0.1.7

* Further improved support for cumulative effects
* Added vignette on estimation and visualization of cumulative effect
* Updated vignette on convenience functions (now "Workflow and convenience functions")
* Other (minor) upgrades/updates to documentation/vignettes
* Updates homepage (via pkgdown)


# pammtools 0.1.3

## Minor changes

* Update documentation
* More tests/improved coverage
* Lag-lead column is adjusted in `make-newdata.fped`

## Bug fixes
- visualization functions `gg_laglead` and `gg_partial_ll` did not
calculate the lag-lead-window correctly when applied to `ped` data

# pammtools 0.1.0

## Features
* Better support for cumulative effects
* Lag-Lead matrix now contains quadrature weights
* Better support for visualization of cumulative effects


# pammtools 0.0.9

## Breaking changes

*  `make_newdata` loses arguments `expand` and `n` and
gains `...` where arbitrary covariate specifications can be placed, i.e.
e.g. `age=seq_range(age, n=20)`. Multiple such expression can be provided and
a data frame with one row for each combination of the evaluated expressions
will be returned. All variables not specified in \code{...} will be set to
respective mean or modus values. For data of class `ped` or `fped` `make_newdata` will try to specify time-dependent variables intelligently.


* `te_var` argument in `concurrent` and `cumulative` was renamed to
`tz_var`

* `te` arguments have been replaced by `tz` (time points at which `z` was observed) in all functions to avoid confusion with `mgcv::te`
(e.g., `gg_laglead`)


## Updates and new features

* Overall better support for cumulative effects

* Added convenience functions for work with cumulative effects, namely
    - `gg_partial` and
    - `gg_slice`

* Added helper functions to calculate and visualize Lag-lead windows
    - `get_laglead`
    - `gg_laglead`

* Added convenience `geom`s for piece-wise constant hazards (see examples in
`?geom_hazard`, cumulative hazards and survival probabilities (usually
`aes(x=time, y = surv_prob)`, but data set doesn't contain extra row for
`time = 0`), thus
    - `geom_stephazard` adds row (x=0, y = y[1]) to the data before plotting
    - `geom_hazard` adds row (x = 0, y = 0) before plotting (can also be used
    for cumulative hazard)
    - `geom_surv` add row (x = 0, y = 1) before plotting


# pammtools 0.0.8

* All data transformation is now handled using `as_ped` (see
[data transformation vignette](https://adibender.github.io/pammtools/articles/data-transformation.html))

* Data transformation now handles
    - standard time-to-event data
    - time-to-event data with concurrent effects of time-dependent covariates
    - time-to-event data with cumulative effects of time-dependent covariates

* Added functionality to flexibly simulate data from PEXP including cumulative effects, see `?sim_pexp`

* Added functionality to calculate Aalen-model style cumulative coefficients,
see `?cumulative_coefficient`


* Breaking change in `split_data` (`as_ped` now main data trafo function):
    - removed `max.end` argument
    - added `max_time` argument to introduce administrative censoring at
    `max_time` when no custom interval split points are provided


# pammtools 0.0.3

## pammtools 0.0.3.2
* More `tidyeval` adaptations
* consistent handling of "no visible global binding" NOTEs
* Release used in <br>
A. Bender, Groll A., Scheipl F., "A generalized additive model approach to
time-to-event analysis" (2017). Statistical Modelling (*to appear*)

## pammtools 0.0.3.1
* some adaptations to `tidyeval`
* Minor bug fixes


# pammtools 0.0.2

* Ported `pamm` package to `pammtools` due to naming conflicts with `PAMM`
package on CRAN
