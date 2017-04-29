Sys.setenv("R_TESTS" = "") # see https://github.com/hadley/testthat/issues/86
library(testthat)
library(pam)

test_check("pam")
