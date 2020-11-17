context("Test as_ped functions")

test_that("Trafo works and attributes are appended", {
  # preparations
  data("veteran", package = "survival")
  veteran <- veteran[c(1:3, 135:137), ]
  ped <- as_ped(
    data    = veteran,
    formula = Surv(time, status)~ trt + age,
    cut     = c(0, 100, 400))
  # retransform to ped
  expect_data_frame(ped, nrow = 10L, ncols = 8L)
  expect_is(ped, "ped")
  expect_subset(c("ped_status", "tstart", "tend", "interval", "offset"),
    names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(is.ped(ped), TRUE)

  ped <- as_ped(
    data = veteran,
    formula = Surv(time, status)~ trt + age)
  expect_data_frame(ped, nrows = 21L, ncols = 8L)


})

test_that("Trafo works for list objects (with TDCs)", {
  data("patient")
  event_df  <- filter(patient, CombinedID %in% c(1110, 1116))
  ped <- as_ped(data = list(event_df), formula = Surv(survhosp, PatientDied)~ .,
    cut = 0:30, id = "CombinedID")
  expect_data_frame(ped, nrows = 40, ncols = 15)
  tdc_df    <- filter(daily, CombinedID  %in% c(1110, 1116))
  ## check nesting
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(survhosp, PatientDied) ~ . +
      cumulative(survhosp, Study_Day, caloriesPercentage, tz_var = "Study_Day") +
        cumulative(proteinGproKG, tz_var = "Study_Day"),
    cut     = 0:30,
    id  = "CombinedID")
  expect_subset("survhosp_Study_Day_mat", colnames(ped))
  expect_data_frame(ped, nrows = 40L, ncols = 20L)
  expect_identical(any(is.na(ped$caloriesPercentage_Study_Day)), FALSE)
  expect_identical(colnames(ped$Study_Day), paste0("Study_Day", 1:12))
  ped <- as_ped(
      data    = list(event_df, tdc_df),
      formula = Surv(survhosp, PatientDied) ~ . +
        cumulative(Study_Day, caloriesPercentage, tz_var = "Study_Day") +
          cumulative(proteinGproKG, tz_var = "Study_Day"),
      id  = "CombinedID")
  expect_data_frame(ped, nrows = 2L, ncols = 19L)

})


test_that("Trafo works for left truncated data", {

  mort_ped <- as_ped(Surv(tstart, exit, event) ~ ses, data = mort2)
  expect_data_frame(mort_ped, nrows = 3L, ncols = 7L)
  expect_identical(round(mort_ped$tstart, 2), c(0.00, 3.48, 3.48))
  expect_identical(round(mort_ped$tend, 2), c(3.48, 17.56, 17.56))
  expect_identical(round(mort_ped$offset, 2), c(1.25, 2.65, 2.65))
  expect_identical(mort_ped$ped_status, c(0, 0, 1))
  expect_identical(mort_ped$ses, factor(c("upper", "upper", "lower")))

})
