context("Convenience functions for calculation of hazard and similar")

data("veteran", package = "survival")
ped <- veteran %>% as_ped(Surv(time, status)~ trt + age,
  cut = c(0, 50, 100, 200, 300, 400), id = "id")
pam <- mgcv::gam(ped_status ~ s(tend, k = 5) + trt, data = ped,
  family = poisson(), offset = offset)
bam <- mgcv::bam(ped_status ~ s(tend, k = 5) + trt, data = ped,
  family = poisson(), offset = offset, method = "fREML", discrete = TRUE)
pem <- glm(ped_status ~ 0 + interval + trt, data = ped,
  family = poisson(), offset = offset)

test_that("hazard functions work for PAM", {

  expect_data_frame(haz <- add_hazard(ped_info(ped), bam), nrows = 5L,
    ncols = 11L)
  expect_data_frame(haz <- add_hazard(ped_info(ped), pam), nrows = 5L,
    ncols = 11L)
  expect_equal(all(haz$ci_lower < haz$hazard), TRUE)
  expect_equal(all(haz$ci_upper > haz$hazard), TRUE)
  expect_equal(round(haz$hazard, 3), c(0.009, 0.009, 0.007, 0.006, 0.005))
  expect_equal(round(haz$ci_lower, 3), c(0.007, 0.007, 0.006, 0.004, 0.003))
  expect_error(add_hazard(haz, pam))
  expect_data_frame(add_hazard(haz, pam, overwrite = TRUE),
    nrows = 5L, ncols = 11L)

  haz2 <- add_hazard(ped_info(ped), pam, type = "link")
  expect_equal(all(haz2$ci_lower < haz2$hazard), TRUE)
  expect_equal(all(haz2$ci_upper > haz2$hazard), TRUE)
  expect_equal(round(haz2$hazard, 2), c(-4.67, -4.75, -4.91, -5.07, -5.23))
  expect_equal(round(haz2$ci_lower, 2), c(-4.91, -4.94, -5.12, -5.42, -5.74))

  ## delta rule
  expect_data_frame(add_hazard(ped_info(ped), bam, ci_type = "delta"),
    nrows = 5L, ncols = 11L)
  haz3 <- add_hazard(ped_info(ped), pam, ci_type = "delta")
  expect_data_frame(haz3, nrows = 5L, ncols = 11L)
  expect_equal(round(haz3$hazard * 100, 2), c(.94, .87, .74, .63, .54))
  expect_equal(round(haz3$se * 100, 2), c(.11, .08, .08, .11, .14))
  expect_equal(round(haz3$ci_lower * 100, 2), c(.72, .70, .58, .41, .26))
  expect_equal(round(haz3$ci_upper * 100, 2), c(1.16, 1.03, .9, .85, .82))

  ## simulation based ci (0.95)
  haz4 <- add_hazard(ped_info(ped), pam, ci_type = "sim")

})

test_that("hazard functions work for PEM", {

  expect_data_frame(haz <- add_hazard(ped_info(ped), pem),
    nrows = 5L, ncols = 11L)
  expect_error(add_hazard(haz, pem))
  expect_data_frame(add_hazard(haz, pem, overwrite = TRUE),
    nrows = 5L, ncols = 11L)

})


test_that("cumulative hazard functions work for PAM", {

  expect_data_frame(add_cumu_hazard(ped_info(ped), bam, ci = FALSE),
   nrows = 5L, ncols = 8L)
 expect_data_frame(haz <- add_cumu_hazard(ped_info(ped), pam, ci = FALSE),
   nrows = 5L, ncols = 8L)
  expect_data_frame(haz <- add_cumu_hazard(ped_info(ped), pam),
    nrows = 5L, ncols = 10L)
  expect_equal(round(haz$cumu_hazard, 2), c(0.47, 0.90, 1.64, 2.27, 2.81))
  expect_equal(round(haz$cumu_lower, 2), c(0.37, 0.73, 1.32, 1.77, 2.09))
  expect_equal(all(diff(haz$cumu_hazard) >= 0), TRUE)
  # overwrite works
  expect_data_frame(add_cumu_hazard(haz, pam, overwrite = TRUE),
    nrows = 5L, ncols = 10L)

  # error on wrong input
  expect_error(add_cumu_hazard(haz, pam))

  ## test that cumu_hazard works for grouped data
  grouped_haz <- ped %>% group_by(trt) %>%
    ped_info() %>%
    add_cumu_hazard(pam)
  expect_data_frame(grouped_haz, nrows = 10L, ncols = 10L)
  expect_equal(round(grouped_haz$cumu_hazard, 2),
    c(0.46, 0.88, 1.60, 2.21, 2.74, 0.48, 0.93, 1.69, 2.33, 2.89))

  ## delta method
  haz2 <- ped_info(ped) %>% add_cumu_hazard(pam, ci_type = "delta")
  expect_equal(round(haz2$cumu_upper, 2), c(.58, 1.09, 1.94, 2.72, 3.48))
  expect_equal(round(haz2$cumu_lower, 2), c(.36, .71, 1.34, 1.82, 2.14))

  suppressWarnings(RNGversion("3.5.0"))
  ## sim CI (0.95)
  set.seed(123)
  haz3 <- ped_info(ped) %>% add_cumu_hazard(pam, ci_type = "sim")
  expect_equal(round(haz3$cumu_upper, 2), c(.58, 1.11, 1.95, 2.75, 3.59))
  expect_equal(round(haz3$cumu_lower, 2), c(.38, .76, 1.42, 1.92, 2.28))

  ## check that hazard columns are not deleted
  newdata <- ped_info(ped) %>% add_hazard(pam) %>%
    add_cumu_hazard(pam)
  expect_data_frame(newdata, nrows = 5L, ncols = 14L)
  newdata <- ped_info(ped) %>% add_hazard(pam, ci = FALSE) %>%
    add_cumu_hazard(pam)
  expect_data_frame(newdata, nrows = 5L, ncols = 11L)

})

test_that("cumulative hazard functions work for PEM", {

  expect_data_frame(haz <- add_cumu_hazard(ped_info(ped), pem),
    nrows = 5L, ncols = 10L)
  expect_error(add_cumu_hazard(haz, pem))
  expect_data_frame(add_cumu_hazard(haz, pem, overwrite = TRUE),
    nrows = 5L, ncols = 10L)

})

test_that("adding terms works for PAM", {
  expect_data_frame(term <- add_term(ped_info(ped), bam, term = "trt"),
    nrows = 5L, ncols = 10L)
  expect_data_frame(term <- add_term(ped_info(ped), pam, term = "trt"),
    nrows = 5L, ncols = 10L)
  expect_data_frame(add_term(ped_info(ped), pam, term = "trt", relative = TRUE),
    nrows = 5L, ncols = 10L)
  pam$model <- NULL
  expect_error(add_term(ped_info(ped), pam, term = "trt", relative = TRUE))
})

test_that("adding terms works for PEM", {
expect_data_frame(term <- add_term(ped_info(ped), pem, term = "trt"),
    nrows = 5L, ncols = 10L)
  expect_data_frame(add_term(ped_info(ped), pem, term = "trt", relative = TRUE),
    nrows = 5L, ncols = 10L)
  pem$model <- NULL
  expect_error(add_term(ped_info(ped), pem, term = "trt", relative = TRUE))
})

test_that("warns about / aborts for unknown intervals", {
  weird <- make_newdata(ped_info(ped), tend = c(150), interval = c("(1.4, 4]"))
  expect_warning(add_hazard(weird, pam), "not used in original fit")
  expect_error(add_hazard(weird, pem), "not used in original fit")
})

test_that("works for nonstandard baseline arguments", {
  pseudonymous <- ped %>% dplyr::rename(stop = tend, int = interval)
  pseudonymous <-  pseudonymous %>% dplyr::mutate(length = stop - tstart)
  ped <- ped %>% dplyr::mutate(intlen = tend - tstart)

  p_pam <- mgcv::gam(ped_status ~ s(stop, k = 5) + trt, data = pseudonymous,
    family = poisson(), offset = offset)
  p_pem <- glm(ped_status ~ 0 + int + trt, data = pseudonymous,
    family = poisson(), offset = offset)
  expect_equal(
    add_hazard(pseudonymous[1:5, ], p_pam, time_var = "stop")$hazard,
    add_hazard(ped[1:5, ], pam)$hazard)
  expect_equal(
    add_hazard(pseudonymous[1:5, ], p_pem, time_var = "int")$hazard,
    add_hazard(ped[1:5, ], pem)$hazard)

  expect_equal(
    add_cumu_hazard(pseudonymous[1:5, ], p_pam, time_var = "stop",
      interval_length = length)$cumu_hazard,
    add_cumu_hazard(ped[1:5, ], pam)$cumu_hazard)
  expect_equal(
    add_cumu_hazard(pseudonymous[1:5, ], p_pem, time_var = "int",
      interval_length = length)$cumu_hazard,
    add_cumu_hazard(ped[1:5, ], pem)$cumu_hazard)
  expect_equal(
    add_cumu_hazard(pseudonymous[1:5, ], p_pem, time_var = "int",
      interval_length = "length")$cumu_hazard,
    add_cumu_hazard(ped[1:5, ], pem)$cumu_hazard)
})


## test surv_prob
test_that("survival probabilities functions work for PAM", {

  suppressWarnings(RNGversion("3.5.0"))

  expect_data_frame(add_surv_prob(ped_info(ped), bam, ci = FALSE),
    nrows = 5L, ncols = 8L)
  expect_data_frame(surv <- add_surv_prob(ped_info(ped), pam, ci = FALSE),
    nrows = 5L, ncols = 8L)
  expect_data_frame(
    surv <- add_surv_prob(ped_info(ped), pam), nrows = 5L, ncols = 10L)
  stest <- sapply(surv[, c("surv_prob", "surv_lower", "surv_upper")],
    function(z) {
      all(z >= 0 & z <= 1)
    })
  expect_identical(all(stest), TRUE)
  expect_identical(round(surv$surv_prob, 2), c(0.63, 0.41, 0.19, 0.10, 0.06))
  expect_identical(round(surv$surv_lower, 2), c(0.55, 0.33, 0.13, 0.05, 0.02))
  expect_identical(round(surv$surv_upper, 2), c(0.69, 0.48, 0.27, 0.17, 0.12))
  # check that overwrite works
  expect_data_frame(add_surv_prob(surv, pam, overwrite = TRUE),
    nrows = 5L, ncols = 10L)
  # error on wrong input
  expect_error(add_surv_prob(surv, pam))

  ## test that cumu_hazard works for grouped data
  grouped_surv <- ped %>% group_by(trt) %>%
    ped_info() %>%
    add_surv_prob(pam)
  expect_data_frame(grouped_surv, nrows = 10L, ncols = 10L)
  expect_equal(round(grouped_surv$surv_prob, 2),
    c(0.63, 0.42, 0.20, .11, .06, .62, 0.40, .19, .10, .06))

  ## delta CI
  surv2 <- add_surv_prob(ped_info(ped), pam, ci_type = "delta")
  expect_equal(round(surv2$surv_lower * 10, 2), c(5.56, 3.28, 1.36, .57, .2))
  expect_equal(round(surv2$surv_upper * 10, 2), c(6.95, 4.83, 2.51, 1.49, 1.01))

  # sim CI
  set.seed(123)
  surv3 <- add_surv_prob(ped_info(ped), pam, ci_type = "sim")
  expect_equal(round(surv3$surv_lower * 10, 2), c(5.59, 3.28, 1.42, .64, .28))
  expect_equal(round(surv3$surv_upper * 10, 2), c(6.86, 4.70, 2.42, 1.47, 1.02))

})


## test sensibility

test_that("hazards and CI positive for type response", {

  ped <- veteran %>% as_ped(Surv(time, status)~ trt + age, id = "id")
  pam <- mgcv::gam(ped_status ~ s(tend, k = 5) + trt,
    data = ped, family = poisson(), offset = offset)
  haz_test <- add_hazard(ped_info(ped), pam) %>%
    summarize_at(c("hazard", "ci_lower", "ci_upper"), funs(any(. < 0)))
  expect_equal(any(unlist(haz_test)), FALSE)

})
