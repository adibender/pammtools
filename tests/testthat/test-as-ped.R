context("Test as_ped functions")

test_that("Trafo works and attributes are appended", {
  # preparations
  data("tumor")
  tumor <- tumor[c(1:3, 135:137), ]
  ped <- as_ped(
    data    = tumor,
    formula = Surv(days, status)~ complications + age,
    cut     = c(0, 100, 400))
  # retransform to ped
  expect_data_frame(ped, nrow = 12L, ncols = 8L)
  expect_is(ped, "ped")
  expect_subset(c("ped_status", "tstart", "tend", "interval", "offset"),
    names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(is.ped(ped), TRUE)

  ped <- as_ped(
    data = tumor,
    formula = Surv(days, status)~ complications + age)
  expect_data_frame(ped, nrows = 11L, ncols = 8L)


})

test_that("Trafo works for list objects (with TDCs)", {
  data("patient")
  event_df  <- filter(patient, CombinedID %in% c(1110, 1116))
  ped <- as_ped(data = list(event_df), formula = Surv(survhosp, PatientDied)~ .,
    cut = 0:30, id = "CombinedID")
  expect_data_frame(ped, nrows = 40, ncols = 15)
  tdc_df    <- filter(daily, CombinedID  %in% c(1110, 1116))
  ## check nesting
  expect_error(as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(survhosp, PatientDied) ~ .,
    cut     = 0:30,
    id  = "CombinedID"))
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

  mort2 <- mort %>% group_by(id) %>% slice(1) %>% filter(id %in% c(1:3))
  mort_ped <- as_ped(Surv(tstart, exit, event) ~ ses, data = mort2)
  expect_data_frame(mort_ped, nrows = 8L, ncols = 7L)
  expect_identical(round(mort_ped$tstart, 2), c(0.00, 3.48, 13.46, 17.56, 3.48, 13.46, 0.00, 3.48))
  expect_identical(round(mort_ped$tend, 2), c(3.48, 13.46, 17.56, 20.00, 13.46, 17.56, 3.48, 13.46))
  expect_identical(round(mort_ped$offset, 2), c(1.25, 2.30, 1.41, 0.89, 2.30, 1.41, 1.25, 2.30))
  expect_identical(mort_ped$ped_status, c(rep(0, 5), 1, 0, 0))
  expect_identical(mort_ped$ses, factor(rep(c("upper", "lower", "upper"), times = c(4,2,2))))

})


test_that("Trafo works for recurrent events data", {

  test_df <- data.frame(
    id     = c(1,1, 2,2,2),
    tstart = c(0, .5, 0, .8, 1.2),
    tstop  = c(.5, 3, .8, 1.2, 3),
    status = c(1, 0, 1, 1, 0),
    enum   = c(1, 2, 1, 2, 3),
    age    = c(50, 50, 24, 24, 24))
  # GAP timescale
  gap_df <- as_ped(
    data       = test_df,
    formula    = Surv(tstart, tstop, status)~ enum + age,
    transition = "enum",
    id         = "id",
    timescale  = "gap")

  expect_data_frame(gap_df, nrows = 9L, ncols = 8L)
  expect_identical(
    round(gap_df$tstart, 1),
    c(0.0, 0.4, 0.0, 0.4, 0.5, 0.0, 0.4, 0.5, 0.0))
  expect_identical(
    round(gap_df$tend, 1),
    c(0.4, 0.5, 0.4, 0.5, 0.8, 0.4, 0.5, 0.8, 0.4))
  expect_identical(
    gap_df$ped_status,
    c(0, 1, 0, 0, 1, 0, 0, 0, 1)
  )
  expect_identical(
    gap_df$enum,
    rep(c(1, 2), times = c(5, 4))
  )

  ## CALENDAR timescale
  cal_df <- as_ped(
    data       = test_df,
    formula    = Surv(tstart, tstop, status)~ age,
    id         = "id",
    transition = "enum",
    timescale  = "calendar")

  expect_data_frame(cal_df, nrows = 6L, ncols = 8L)
  expect_identical(
    round(cal_df$tstart, 1),
    c(0.0, 0.0, 0.5, 0.5, 0.8, 0.8))
  expect_identical(
    round(cal_df$tend, 1),
    c(0.5, 0.5, 0.8, 0.8, 1.2, 1.2))
  expect_identical(
    cal_df$ped_status,
    c(1, 0, 1, 0, 0, 1)
  )
  expect_identical(
    cal_df$enum,
    rep(c(1, 2), each = 3)
  )

})


test_that("Trafo works for multi-state data without recurrent events", {

  test_df <- data.frame(
    id     = c(1, 2,2),
    tstart = c(0, 0,1.2),
    tstop  = c(3, 1.2,3),
    status = c(1, 1,1),
    from = c(1, 1,2),
    to = c(3, 2,3),
    transition   = c("1->3", "1->2","2->3"),
    age    = c(24, 36,36))

  test_df <- test_df %>% add_counterfactual_transitions()

   # CALENDAR timescale
  cal_df <- as_ped(
    data       = test_df,
    formula    = Surv(tstart, tstop, status)~ .,
    transition = "transition",
    id         = "id",
    timescale  = "calendar")

   # according to code: order by transition -> id -> tstart
  expect_data_frame(cal_df, nrows = 7L, ncols = 10L)
  expect_identical(cal_df$transition,
    as.factor(c("1->2", "1->2","1->2","1->3", "1->3","1->3","2->3")))
  expect_identical(
    cal_df$id,
    c(1, 1, 2, 1, 1, 2, 2))
  expect_identical(
    round(cal_df$tstart, 1),
    c(0.0, 1.2, 0.0, 0.0, 1.2, 0.0, 1.2))
  expect_identical(
    round(cal_df$tend, 1),
    c(1.2, 3.0, 1.2, 1.2, 3.0, 1.2, 3.0))
  expect_identical(
    cal_df$ped_status,
    c(0, 0, 1, 0, 1, 0, 1)
  )
  expect_identical(
    cal_df$from,
    c(1,1,1, 1,1,1, 2)
  )
  expect_identical(
    cal_df$to,
    c(2,2,2, 3,3,3, 3)
  )

  }
)
