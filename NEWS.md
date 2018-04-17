# pammtools 0.0.8

* All data transformation is now hadled using `as_ped` (see
[data transformation vignette](../articles/data-transformation.html))

* Data transformation now handles
    - standard time-to-event data
    - time-to-event data with concurrent effects of time-dependent covariates
    - time-to-event data with cumulative effects of time-dependent covariates

* `make_newdata` handles creation of new data sets for all of the above cases

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
