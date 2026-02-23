context("Simple transformation to PED data")

test_that("Transformation of regular survival data works", {
  ## preparations
  data("tumor")
  tumor <- tumor[c(1:3, 135:137), ]
  tumor$ID <- sample(1:100, nrow(tumor))
  ped_vet <- split_data(data = tumor, Surv(days, status) ~ complications + age,
    cut = c(0, 100, 400), id = "ID")
  expect_identical(unique(ped_vet$ID), tumor$ID)
  tumor$ID <- NULL

  ## tests
  expect_data_frame(ped <- tumor %>% as_ped(Surv(days, status)~ complications + age,
    cut = c(0, 100, 400)), nrows = 12L, ncols = 8L)
  expect_is(ped, "ped")
  expect_subset(c("ped_status", "tstart", "tend", "interval", "offset"),
    names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_data_frame(tumor %>% as_ped(Surv(days, status)~ complications + age,
    cut = c(0, 100, 400), id = "id"), nrows = 12L, ncols = 8L)
  ## no error, when id in data and additionally specified
  tumor$id <- seq_len(nrow(tumor))
  expect_data_frame(ped <- tumor %>% as_ped(Surv(days, status)~complications,
    cut = c(0, 100, 400), id = "id"), nrows = 12L, ncols = 7L)
  ## no error when id already in data but not specified
  expect_data_frame(tumor %>% as_ped(Surv(days, status)~complications,
    cut = c(0, 100, 400)), nrows = 12L, ncols = 7L)
  ## no error when id has different name and is specified accordingly
  tumor$id2 <- tumor$id
  expect_data_frame(ped <- tumor %>% as_ped(Surv(days, status)~.,
    cut = c(0, 100, 400)), nrows = 12L, ncols = 14L)
  expect_identical(attr(ped, "id_var"), "id")
  ## no additional id when different id specified
  tumor$id <- NULL
  expect_data_frame(ped <- tumor %>% as_ped(Surv(days, status)~.,
    cut = c(0, 100, 400), id = "id2"), nrows = 12L, ncols = 13L)
  expect_identical(attr(ped, "id_var"), "id2")
  tumor$id2 <- NULL
  # max_time
  ped <- tumor %>% as_ped(Surv(days, status)~., max_time = 400)
  expect_data_frame(ped, nrows = 11L, ncols = 13L)
  expect_identical(max(ped$tend), 400)
  expect_identical(nlevels(ped$interval), 2L)

  # inlcude_last
  tumor[6, "days"] <- 358
  ped <- tumor %>% as_ped(Surv(days, status)~.)
  expect_data_frame(ped, nrows = 11L, ncols = 13L)
  expect_identical(max(ped$tend), 358)

})


test_that("Error on wrong input", {
  ## preparations
  data("tumor")
  tumor <- tumor[c(1:3, 135:137), ]

  ## tests
  expect_error(as_ped(tumor, x ~ y, cut = c(0:5, 10, 40)))
  expect_error(as_ped(tumor, Surv(days2, status) ~., cut = c(0:5, 10, 40)))
  expect_error(as_ped(
    data = rename(tumor, ped_time = time),
    formula = Surv(ped_time, status) ~.))
  # already in data set ped_time

  ## error when specified id variable not unique
  tumor$id <- rep(1:2, 3)
  expect_error(
    as_ped(tumor, Surv(days, status) ~ complications, cut = c(0, 100, 400), id = "id"),
    regexp = "Specified ID variable.*")
})

test_that("as_ped conserves total follow-up and event counts", {
  df <- data.frame(
    id = 1:4,
    time = c(2, 4, 3, 1),
    status = c(1, 0, 1, 1),
    x = c(-0.3, 0.8, 1.1, 0.2)
  )

  ped <- as_ped(
    df,
    survival::Surv(time, status) ~ x,
    id = "id",
    cut = c(0, 0.5, 1, 2, 4, 5)
  )

  observed <- as.data.frame(ped) |>
    dplyr::group_by(.data$id) |>
    dplyr::summarise(
      total_interval_length = sum(exp(.data$offset)),
      total_events = sum(.data$ped_status),
      .groups = "drop"
    ) |>
    dplyr::left_join(df, by = "id")

  expect_equal(observed$total_interval_length, observed$time, tolerance = 1e-12)
  expect_equal(observed$total_events, observed$status)
})

test_that("as_ped with left truncation conserves risk time per subject", {
  df <- data.frame(
    id = 1:3,
    tstart = c(0, 1, 0.5),
    tstop = c(2, 3, 1.5),
    status = c(1, 0, 1),
    x = c(0.1, 0.2, 0.3)
  )

  ped <- as_ped(
    df,
    survival::Surv(tstart, tstop, status) ~ x,
    id = "id",
    cut = seq(0, 3, by = 0.5)
  )

  observed <- as.data.frame(ped) |>
    dplyr::group_by(.data$id) |>
    dplyr::summarise(
      total_interval_length = sum(exp(.data$offset)),
      total_events = sum(.data$ped_status),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      dplyr::mutate(df, duration = .data$tstop - .data$tstart),
      by = "id"
    )

  expect_equal(
    observed$total_interval_length,
    observed$duration,
    tolerance = 1e-12
  )
  expect_equal(observed$total_events, observed$status)
})

test_that("split_data is invariant to ordering and duplicates in cut points", {
  df <- data.frame(
    id = 1:3,
    time = c(1, 2, 3),
    status = c(1, 0, 1)
  )

  out_sorted <- split_data(
    survival::Surv(time, status) ~ 1,
    data = df,
    cut = c(0, 1, 2, 3)
  )
  out_unsorted <- split_data(
    survival::Surv(time, status) ~ 1,
    data = df,
    cut = c(3, 0, 2, 1, 1)
  )

  expect_equal(
    as.data.frame(out_sorted),
    as.data.frame(out_unsorted)
  )
})
