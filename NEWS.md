## pammtools 0.0.9

* Overall better support for cumulative effects

* **Breaking change** `make_newdata` loses arguments `expand` and `n` and
gains `...` where arbitrary covariate specifications can be placed, i.e.
e.g. `age=seq_range(age, n=20)`. Multiple such expression can be provided and
a data frame with one row for each combination of the evaluated expressions
will be returned. All variables not specified in \code{...} will be set to
respective mean or modus values. For data of class `ped` or `fped` `make_newdata` will try to specify time-dependent variables intelligently.

* Added convenience functions for work with cumulative effects, namely
    - `gg_partial` and
    - `gg_slice`

* Added helper functions to calculate and visualize Lag-lead windows
    - `get_laglead`
    - `gg_laglead`


# pammtools 0.0.8

* All data transformation is now handled using `as_ped` (see
[data transformation vignette](../articles/data-transformation.html))

* Data transformation now handles
    - standard time-to-event data
    - time-to-event data with concurrent effects of time-dependent covariates
    - time-to-event data with cumulative effects of time-dependent covariates

* Added functionality to flexibly simulate data from PEXP including cumulative effects, see `?sim_pexp`

* Added functionality to calculate Aaalen-model style cumulative coefficients,
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
