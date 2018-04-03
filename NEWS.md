# pammtools 0.0.7

## pammtools 0.0.7.3

+ Added support for data transformation of TDCs for cumulative effects
(see `as_ped` and the [data transformation vignette](articles/data-transformation.html))

* Added functionality to flexibly simulate data from PEXP including cumulative effects, see `?sim_pexp`

* Added functionality to calculate Aaalen-model style cumulative coefficients,
see `?cumulative_coefficient`


* Breaking change in `split_data`:
    - removed `max.end` argument
    - added `max.time` argument to introduce administrative censoring at
    `max_time`  when no custom cut-points provided


* Readded details tags in vignettes (after pandoc update)


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
