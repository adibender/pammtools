context("Transformation invariants")

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
