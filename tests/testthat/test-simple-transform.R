context("Simple transformation to PED data")

test_that("Transformation of regular survival data works", {
  ## preparations
  data("veteran", package = "survival")
  veteran <- veteran[c(1:3, 135:137), ]
  veteran$ID <- sample(1:100, nrow(veteran))
  ped_vet <- split_data(data = veteran, Surv(time, status) ~ trt + age,
    cut = c(0, 100, 400), id = "ID")
  expect_identical(unique(ped_vet$ID), veteran$ID)
  veteran$ID <- NULL

  ## tests
  expect_data_frame(ped <- veteran %>% as_ped(Surv(time, status)~ trt + age,
    cut = c(0, 100, 400)), nrows = 10L, ncols = 8L)
  expect_is(ped, "ped")
  expect_subset(c("ped_status", "tstart", "tend", "interval", "offset"),
    names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_data_frame(veteran %>% as_ped(Surv(time, status)~ trt + age,
    cut = c(0, 100, 400), id = "id"), nrows = 10L, ncols = 8L)
  ## no error, when id in data and additionally specified
  veteran$id <- seq_len(nrow(veteran))
  expect_data_frame(ped <- veteran %>% as_ped(Surv(time, status)~trt,
    cut = c(0, 100, 400), id = "id"), nrows = 10L, ncols = 7L)
  ## no error when id already in data but not specified
  expect_data_frame(veteran %>% as_ped(Surv(time, status)~trt,
    cut = c(0, 100, 400)), nrows = 10L, ncols = 7L)
  ## no error when id has different name and is specified accordingly
  veteran$id2 <- veteran$id
  expect_data_frame(ped <- veteran %>% as_ped(Surv(time, status)~.,
    cut = c(0, 100, 400)), nrows = 10L, ncols = 13L)
  expect_identical(attr(ped, "id_var"), "id")
  ## no additional id when different id specified
  veteran$id <- NULL
  expect_data_frame(ped <- veteran %>% as_ped(Surv(time, status)~.,
    cut = c(0, 100, 400), id = "id2"), nrows = 10L, ncols = 12L)
  expect_identical(attr(ped, "id_var"), "id2")
  veteran$id2 <- NULL
  # max_time
  ped <- veteran %>% as_ped(Surv(time, status)~., max_time = 400)
  expect_data_frame(ped, nrows = 21L, ncols = 12L)
  expect_identical(max(ped$tend), 400)
  expect_identical(nlevels(ped$interval), 6L)

  # inlcude_last
  veteran[2, "status"] <- 0
  ped <- veteran %>% as_ped(Surv(time, status)~.)
  expect_data_frame(ped, nrows = 20L, ncols = 12L)
  expect_identical(max(ped$tend), 378)

})


test_that("Error on wrong input", {
  ## preparations
  data("veteran", package = "survival")
  veteran <- veteran[c(1:3, 135:137), ]

  ## tests
  expect_error(as_ped(veteran, x ~ y, cut = c(0:5, 10, 40)))
  expect_error(as_ped(veteran, Surv(time2, status) ~., cut = c(0:5, 10, 40)))
  expect_error(as_ped(
    data = rename(veteran, ped_time = time),
    formula = Surv(ped_time, status) ~.))
  # already in data set ped_time

  ## error when specified id variable not unique
  veteran$id <- rep(1:2, 3)
  expect_error(
    as_ped(veteran, Surv(time, status) ~ trt, cut = c(0, 100, 400), id = "id"),
    regexp = "Specified ID variable.*")
})
