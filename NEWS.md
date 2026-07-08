# pammtools 0.8.1

## Enhancements
* The cumulative post-processing functions (`add_cumu_hazard()`,
  `add_surv_prob()`, `add_cif()` -- for both their default and `pamm_ic`
  methods -- and `add_trans_prob()`) now **error** when `newdata` is not grouped
  so that the time variable is unique within each group (typically a forgotten
  `group_by()`), instead of silently accumulating cumulative quantities across
  distinct covariate profiles / causes / transitions and returning wrong
  curves. This generalises the guard already used for the interval-censoring
  pooling path (and mirrors the `stopifnot()` guard in the RMST example). All
  of these functions accept `check_grouping = FALSE` to opt out (e.g. for
  advanced workflows that intentionally accumulate over rows with varying
  covariates). The guard detects mis-grouping via repeated time values within
  a group: grids built with `make_newdata()` are always caught, but hand-built
  grids stacking profiles with *disjoint* time grids are indistinguishable
  from a single profile with time-varying covariates and pass undetected.

# pammtools 0.8.0

This release collects all changes since the last CRAN version (0.7.4),
including the previously GitHub-only 0.7.5 and 0.7.6 development versions.

## Breaking changes
* `make_newdata()` output no longer contains internal PED columns (`tstart`,
  `intlen`, `interval`, `offset`, `ped_status`). Output now contains `tend` +
  `id` + user covariates (plus `cause`/`transition` for competing risks /
  multi-state models). `ped_info()` output is unchanged. `intlen` is
  reconstructed on demand by downstream functions (`add_cumu_hazard()`,
  `add_surv_prob()`, `add_cif()`, `add_trans_prob()`) and dropped from
  user-facing output.
* `add_cif()` now uses the exact closed-form integral of the cumulative
  incidence function under piecewise-exponential hazards instead of the previous
  left-Riemann approximation. CIF estimates from existing user code change
  numerically; results are now invariant to the time grid passed to
  `make_newdata()`.
* Simulation-based confidence intervals (`ci_type = "sim"`) now use type-6
  empirical quantiles instead of the `stats::quantile()` default (type 7).
  Type-7 quantiles made these intervals systematically too narrow for small
  `nsim` (at the default `nsim = 100` the enclosed central mass is ~93% rather
  than the nominal 95%); type-6 removes this inward bias (#288). As a result,
  all simulation-based CI bounds change slightly (intervals widen at both ends)
  relative to versions <= 0.7.4.

## New features
* Full support for shape-constrained additive models fit with `scam::scam()`
  (#286): the post-processing workflow (`add_hazard()`, `add_cumu_hazard()`,
  `add_surv_prob()`, `add_term()`, `add_cif()`, `add_trans_prob()`,
  `get_cumu_coef()`, `get_cumu_eff()`, `tidy_fixed()`, `tidy_smooth()`,
  `gg_smooth()`, ...) now works for `scam` fits exactly as for `gam` fits,
  including delta-method and simulation-based confidence intervals. The
  calculations correctly use the re-parametrized coefficients
  (`$coefficients.t`) and their covariance (`$Vp.t`). `pamm()` gained
  `engine = "scam"`. See the new "Shape-constrained effects (scam)" article.
* Interval-censored time-to-event data are now supported via a multiple-
  imputation (MI) workflow. Data specified with `Surv(L, R, type = "interval2")`
  are detected automatically by `as_ped()`. The new `pamm_ic()` (single event)
  and `pamm_ic_cr()` (competing risks) fit a PAMM by repeatedly drawing exact
  event times from the model-based conditional hazard distribution on `(L, R]`
  and re-fitting the standard right-censored pipeline. Inference pools the
  imputations: `add_hazard()`, `add_cumu_hazard()`, `add_surv_prob()` and
  `add_cif()` gain `pamm_ic` methods that combine per-imputation posterior draws
  (within- plus between-imputation variance). The `iter` argument enables
  chained (refit-and-reimpute) imputation, recommended for sparsely inspected
  data. `add_inspections()` turns exact simulated times (e.g. from `sim_pexp()`)
  into interval-censored panel data for testing and coverage studies.
  `print()`/`summary()` of a `pamm_ic` report the *pooled* (Rubin-combined) fit.
  See the new "Interval-Censored Data" vignette.
* The post-processing / confidence-interval machinery is now backend-pluggable:
  all `add_*()` quantities and their delta-method and simulation-based CIs are
  derived from two internal S3 primitives, `get_hazard()` and `sim_hazard()`, so
  an alternative estimation backend only needs to provide methods for those two.
  The new "Defining a new backend: gradient boosting with xgboost" vignette
  demonstrates this end-to-end (and the Bayesian vignette was reworked to use
  the same unified interface).
* `gg_state_occupation()` is now exported.

## Enhancements
* `gg_smooth()` is now fully general across univariate smooth terms: a bare
  variable name selects every 1d smooth over that variable (main effect plus any
  `by`-variable or factor-smooth interaction term), `terms` is optional
  (defaulting to all univariate smooths), and 1d `ti()` as well as factor-smooth
  interactions (`bs = "fs"`, `bs = "sz"`) are supported. Factor-indexed smooths
  are drawn in a single facet with one curve per factor level, identified by a
  new `level` column in the `get_terms()` output. Random-effect smooths
  (`bs = "re"`, `bs = "mrf"`) and multivariate/tensor smooths are excluded
  (use `gg_re()` / `gg_tensor()`).
* `add_cif()` now supports arbitrary time points in `make_newdata()` (parity
  with `add_cumu_hazard()`); missing breakpoints are inserted internally so CIF
  estimates are independent of the chosen prediction grid.
* `add_surv_prob()`, `add_cif()`, `add_trans_prob()` and `add_cumu_hazard()` now
  include plotting boundary rows at `tend = 0` (or the selected `time_var`).
  Boundary values are set to their known limits, `S(0) = 1`, `CIF(0) = 0`,
  off-diagonal transition probabilities `P_rs(0) = 0`, and `cumu_hazard = 0`,
  with collapsed confidence-interval bounds when requested. Boundary rows are
  added only for continuous-time models (`gam`/`scam`/`pamm`).
* Simulation-based CIs now draw a single shared posterior coefficient sample
  across groups (consistent with `add_cif()` / `add_trans_ci()`); single-group
  results are unchanged.
* `get_trans_prob()` now supports non-integer (categorical) state labels
  (e.g. `"healthy->ill"`) in addition to integer-coded transitions.
* Transition probability calculation is faster due to a base R refactor.
* `expand_df()` preserves the `cause` column when `make_newdata()` is called with
  only `tend` and `cause`, fixing a competing-risks edge case.
* `predictSurvProb.pamm()` now respects non-default `id` column names and works
  when `trafo_args` are not attached to the fitted object.

## Bug fixes
* `add_trans_prob()` and `add_trans_ci()` no longer require the input data to be
  pre-sorted (#255, related to #227).
* Fixed transition probability matrix dimensions when transitions start from
  state 0 (off-by-one in state indexing).
* `add_trans_prob()` / `add_trans_ci()` / `get_trans_prob()` now consistently
  thread `time_var` and `interval_length`, fixing argument forwarding for
  nonstandard column names.
* `add_counterfactual_transitions()` now fully honors `from_col`, `to_col`,
  and `transition_col`.
* `gg_smooth()` / `get_terms()` now select smooth terms via the model's `mgcv`
  smooth metadata instead of unanchored `grep()`, fixing two errors reported in
  #283 (variable names matched by several smooths, and factor terms). Names that
  match no smooth are skipped with a warning rather than erroring.
* Fixed CIF cause mislabeling under non-alphabetical factor levels.
* Fixed interval-censored pooling regressions (backend-aware prediction,
  boundary rows, and a survival-probability length mismatch).
* Registered the `.glm`/`.pamm` methods for the internal
  `warn_about_new_time_points()` generic (previously "no applicable method").
* Pooled `pamm_ic` adders now warn when given under-grouped `newdata`.

## Deprecations
* The `trafo_args` argument of `pamm()` is deprecated; convert data with
  `as_ped()` before calling `pamm()`.

## Documentation
* Added derivation of the piecewise-exponential CIF integral to the
  competing-risks vignette.

# pammtools 0.7.4

## Bug fixes
* Fixed competing risks data transformation when status variable is a factor (#220, #216, #233)
* Fixed CIF calculation to use factor levels from newdata instead of model attribute (#245)
* Fixed cut point extraction for factor/character status variables
* Fixed transition probability matrix initialization for states starting at 0 or 1
* Fixed CRAN NOTE: added `id` to global variables for dplyr compatibility (#260)

## Enhancements
* Improved `add_trans_prob`: better documentation, proper examples, attribute attachment, and base R speedup
* Added warning in `pamm()` when data does not contain an offset column
* Added `broom` to Suggests

## Documentation
* Updated `add_trans_prob` help page with proper parameter descriptions and working example
* Added simulations vignette

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
