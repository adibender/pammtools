context("Convenience functions for calculation of hazard and similar")

library(mgcv)
data("veteran", package="survival")
ped <- split_data(Surv(time, status)~ trt + age, data=veteran, 
  cut=c(0, 50, 100, 200, 300, 400), id="id")
pam <- gam(ped_status ~ s(tend, k=5) + trt, data=ped,
  family = poisson(), offset = offset)
pem <- glm(ped_status ~ 0 + interval + trt, data=ped,
  family = poisson(), offset = offset)

test_that("hazard functions work for PAM", {
	expect_is(haz <- add_hazard(ped_info(ped), pam), "data.frame")
	expect_identical(dim(haz), c(5L, 11L))
	expect_error(add_hazard(haz, pam))
	expect_identical(dim(add_hazard(haz, pam, overwrite = TRUE)), c(5L, 11L))
})

test_that("hazard functions work for PEM", {
  expect_is(haz <- add_hazard(ped_info(ped), pem), "data.frame")
  expect_identical(dim(haz), c(5L, 11L))
  expect_error(add_hazard(haz, pem))
  expect_identical(dim(add_hazard(haz, pem, overwrite = TRUE)), c(5L, 11L))
})


test_that("cumulative hazard functions work for PAM", {
	expect_is(haz <- add_cumhazard(ped_info(ped), pam), "data.frame")
	expect_identical(dim(haz), c(5L, 10L))
	expect_error(add_cumhazard(haz, pam))
	expect_identical(dim(add_cumhazard(haz, pam, overwrite = TRUE)), c(5L, 10L))
})

test_that("cumulative hazard functions work for PEM", {
  expect_is(haz <- add_cumhazard(ped_info(ped), pem), "data.frame")
  expect_identical(dim(haz), c(5L, 10L))
  expect_error(add_cumhazard(haz, pem))
  expect_identical(dim(add_cumhazard(haz, pem, overwrite = TRUE)), c(5L, 10L))
})

test_that("adding terms works for PAM", {
	expect_is(add_term(ped_info(ped), pam, term="trt"), "data.frame")
	expect_is(add_term(ped_info(ped), pam, term="trt", relative=TRUE),
		"data.frame")
	pam$model <- NULL
	expect_error(add_term(ped_info(ped), pam, term="trt", relative=TRUE))
})

test_that("adding terms works for PEM", {
  expect_is(add_term(ped_info(ped), pem, term="trt"), "data.frame")
  expect_is(add_term(ped_info(ped), pem, term="trt", relative=TRUE),
    "data.frame")
  pem$model <- NULL
  expect_error(add_term(ped_info(ped), pem, term="trt", relative=TRUE))
})

test_that("warns about / aborts for unknown intervals", {
  weird <- make_newdata(ped_info(ped), tend = 150, interval = "(1.4, 4]")
  expect_warning(add_hazard(weird, pam), "not used in original fit")
  expect_error(add_hazard(weird, pem), "not used in original fit")
})

test_that("works for nonstandard baseline arguments", {
  pseudonymous <- ped %>% dplyr::rename(stop = tend, int = interval)
  pseudonymous <-  pseudonymous %>% dplyr::mutate(length = stop - tstart)
  ped <- ped %>% dplyr::mutate(intlen = tend - tstart)

  p_pam <- gam(ped_status ~ s(stop, k=5) + trt, data=pseudonymous,
    family = poisson(), offset = offset)
  p_pem <- glm(ped_status ~ 0 + int + trt, data=pseudonymous,
    family = poisson(), offset = offset)
  expect_equal(
    add_hazard(pseudonymous[1:5,], p_pam, time_variable = "stop")$hazard,
    add_hazard(ped[1:5,], pam)$hazard)
  expect_equal(
    add_hazard(pseudonymous[1:5,], p_pem, time_variable = "int")$hazard,
    add_hazard(ped[1:5,], pem)$hazard)

  expect_equal(
    add_cumhazard(pseudonymous[1:5,], p_pam, time_variable = "stop",
      interval_length = dplyr::quo(length))$cumhazard,
    add_cumhazard(ped[1:5,], pam)$cumhazard)
  expect_equal(
    add_cumhazard(pseudonymous[1:5,], p_pem, time_variable = "int",
      interval_length = dplyr::quo(length))$cumhazard,
    add_cumhazard(ped[1:5,], pem)$cumhazard)
})


## test sensibility 

test_that("hazards and CI positiv for type response", {
  
  ped <- split_data(Surv(time, status)~ trt + age, data=veteran, id="id")
  pam <- gam(ped_status ~ s(tend, k=5) + trt, 
    data=ped, family = poisson(), offset = offset)
  haz_test <- add_hazard(ped_info(ped), pam) %>% 
    summarize_at(c("hazard", "lower", "upper"), funs(any(. < 0)))
  expect_equal(any(unlist(haz_test)), FALSE)

})