context("Test as_ped_cr functions")


test_that("Trafo works and attributes are appended.", {
  # preparations
  ped <- as_ped(
    data         = sir_adm,
    formula      = Surv(time, status) ~ age + pneu,
    cut          = c(0, 10, 100)
  )
  expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
  expect_is(ped, "ped_cr_union")
  expect_subset(c("ped_status", "tstart", "tend",
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(sum(as.numeric(ped$cause)), 36)

  # check that trafo can be recovered
  ped2 <- as_ped(ped, newdata = sir_adm)
  expect_equal(ped, ped2, check.attributes = FALSE)

  # check that list output identical for given cut points
  ped_list <- as_ped(
    data    = sir_adm,
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100),
    combine = FALSE)
  ped2 <- do.call(rbind, ped_list)
  expect_true(all.equal(do.call(rbind, ped_list), ped, check.attributes = FALSE))
  expect_identical(length(ped_list), 2L)
  expect_identical(class(ped_list), c("ped_cr_list", "ped_cr", "ped", "list"))
  expect_identical(names(attributes(ped_list)), c("class", "names", "trafo_args", "risks"))
  expect_identical(length(attr(ped_list, "trafo_args")$cut), 2L)

  # check that trafo can be recovered for ped list objects
  ped_list2 <- as_ped(ped_list, newdata = sir_adm)
  expect_equal(ped_list, ped_list2, check.attributes = FALSE)

  # test when split points not specified
  ped <- as_ped(data = sir_adm, formula = Surv(time, status) ~ .)
  expect_data_frame(ped, nrows = 56L, ncols = 10L)
  expect_equal(sum(as.numeric(ped$cause)), 84L)

  ped_list <- as_ped_cr(sir_adm, Surv(time, status) ~ ., combine = FALSE)
  expect_identical(attr(ped_list[[1]], "breaks"), c(4L, 10L, 24L, 37L, 101L))
  expect_identical(attr(ped_list[[2]], "breaks"), c(22L, 25L))
  ped_list2 <- as_ped(ped_list, newdata = sir_adm)
  expect_equal(ped_list, ped_list2, check.attributes = FALSE)

})

test_that("Trafo works for more than two risks.", {
  # preparations
  sir_adm$status[2] <- 3
  ped <- as_ped(
    data         = sir_adm,
    formula      = Surv(time, status) ~ age + pneu,
    cut          = c(0, 10, 100)
  )
  expect_data_frame(ped, nrow = 12L * 3L, ncols = 9L)
  expect_is(ped, "ped_cr_union")
  expect_subset(c("ped_status", "tstart", "tend",
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(sum(as.numeric(ped$cause)), 72)
  expect_equal(sum(ped$ped_status[ped$cause == 3L]), 1)

})

test_that("Trafo works for status as factor.", {
  # create data set with status as factor
  dat <- data.frame(
    status_fct = factor(c(0, 0, 0, 0, 1, 1, 1, 2, 2, 2)),
    time = seq(from = 10, to = 100, by = 10)
  )
  
  # create as ped with status as factor
  ped <- as_ped(
    data = dat,
    formula = Surv(time, status_fct) ~ .
  )
  
  expect_data_frame(ped, nrow = 50L, ncols = 7L)
  expect_is(ped, "ped_cr_union")
  expect_subset(c("ped_status", "tstart", "tend",
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "risks"), "factor")
  expect_equal(levels(ped$cause), c("1", "2"))
})

test_that("TDC can be transformed for cause-specific data sets (combine = FALSE)", {

  set.seed(1)
  event_df <- data.frame(
    id     = factor(1:10),
    status = c(0L, 1L, 2L, 0L, 1L, 2L, 1L, 0L, 0L, 2L),
    time   = seq(from = 10, to = 100, by = 10))
  tdc_df <- data.frame(
    id  = factor(c(1:10, 1:10)),
    var = rbinom(20, 1, 0.3),
    tt  = c(rep(0, 10), seq(5, 23, 2)))

  ped_list <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~ concurrent(var, tz_var = "tt"),
    id      = "id",
    combine = FALSE)

  expect_class(ped_list, "ped_cr_list")
  expect_identical(names(ped_list), c("cause = 1", "cause = 2"))
  expect_identical(attr(ped_list, "trafo_args")$combine, FALSE)

  for (ped_i in ped_list) {
    expect_true("var" %in% names(ped_i))
    # the covariate value in an interval is the last one observed strictly
    # before the interval end (lag-lead function of concurrent(), lag = 0)
    expected <- vapply(
      seq_len(nrow(ped_i)),
      function(i) {
        tdc_i <- tdc_df[tdc_df$id == ped_i$id[i] & tdc_df$tt < ped_i$tend[i], ]
        if (nrow(tdc_i) == 0) NA_integer_ else tdc_i$var[which.max(tdc_i$tt)]
      },
      integer(1))
    expect_equal(ped_i$var, expected)
  }

})

test_that("Interval breaks cover all causes when TDCs are present", {

  set.seed(1)
  event_df <- data.frame(
    id     = factor(1:10),
    status = c(0L, 1L, 2L, 0L, 1L, 2L, 1L, 0L, 0L, 2L),
    time   = seq(from = 10, to = 100, by = 10))
  tdc_df <- data.frame(
    id  = factor(c(1:10, 1:10)),
    var = rbinom(20, 1, 0.3),
    tt  = c(rep(0, 10), seq(5, 23, 2)))

  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~ concurrent(var, tz_var = "tt"),
    id      = "id")

  # cause 1 has events at 20, 50, 70, cause 2 at 30, 60, 100
  event_times <- sort(event_df$time[event_df$status != 0])
  expect_true(all(event_times %in% attr(ped, "breaks")))
  expect_equal(max(attr(ped, "breaks")), max(event_times))

  # every event of every cause is retained, at its own time
  events <- ped[ped$ped_status == 1, ]
  expect_equal(sort(events$tend), event_times)
  expect_equal(
    as.vector(table(events$cause)),
    as.vector(table(event_df$status[event_df$status != 0])))

})
