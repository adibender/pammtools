#pammtools 0.0.6
<hr>

## pammtools 0.0.6.0
* Added functionality to calculate Aaalen-model style cumulative coefficients,
see `?cumulative_coefficient`
* Minor fixes in `sample_info` for ped objects

# pammtools 0.0.5
<hr>

## pammtools 0.0.5.0

* Breaking change in `split_data`:
    - rename `max.end` to `include_last` and only used when no custom cut-points provided
    - added `max.time` argument to introduce administrative censoring when no
    custom cut-points provided

# pammtools 0.0.4
<hr>

## pammtools 0.0.4.0
* Readded details tags in vignettes (after pandoc update)
* First prototype of functions that create necessary data to estimate cumulative effects


# pammtools 0.0.3
<hr>

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
<hr>

* Ported `pamm` package to `pammtools` due to naming conflicts with `PAMM`
package on CRAN
