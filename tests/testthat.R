Sys.setenv("R_TESTS" = "") # see https://github.com/hadley/testthat/issues/86
library(testthat)
library(checkmate)
library(dplyr)
library(purrr)
library(tidyr)
library(data.table)
library(mvna)
# library(pammtools)

test_check("pammtools")
