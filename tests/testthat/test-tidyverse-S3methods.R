context("Tidyverse methods for specific classes")

test_that("ped class is preserved after dplyr operations", {
  data("tumor")
  tumor <- dplyr::slice(tumor, 2:3)
  ped <- as_ped(
    data = tumor,
    formula = Surv(days, status) ~ complications + age,
    cut = c(0, 100, 400),
    id = "id"
  )

  expect_is(filter(ped, id == 1), "ped")
  expect_is(slice(ped, 1), "ped")
  expect_is(arrange(ped, desc(id)), "ped")
  expect_is(select(ped, id), "ped")
  expect_is(rename(ped, ID = id), "ped")
  expect_is(mutate(ped, id = id + 1), "ped")
  expect_is(transmute(ped, id = id + 1), "ped")
  expect_is(sample_n(ped, 1), "ped")
  expect_is(sample_frac(ped, 0.5), "ped")
  expect_is(right_join(distinct(ped, id, interval), ped), "ped")
})

test_that("right_join.ped keeps unmatched right-hand rows", {

  data("tumor")
  ped <- as_ped(
    data    = dplyr::slice(tumor, 2:3),
    formula = Surv(days, status) ~ age,
    cut     = c(0, 100, 400),
    id      = "id")

  lookup <- data.frame(id = c(1, 999), marker = c("known", "new"))
  joined <- right_join(ped, lookup, by = "id")

  expect_true(any(joined$id == 999))
  expect_true(any(is.na(joined$ped_status[joined$id == 999])))

})

test_that("sample_frac.ped samples a fraction of rows", {

  data("tumor")
  ped <- as_ped(
    data    = dplyr::slice(tumor, 2:3),
    formula = Surv(days, status) ~ complications + age,
    cut     = c(0, 100, 400),
    id      = "id")

  set.seed(123)
  sampled <- sample_frac(ped, size = 0.5)

  expect_gt(nrow(sampled), 0)
  expect_lt(nrow(sampled), nrow(ped))

})


test_that("attributes are preserved", {
  # recurrent events data
  test_df <- data.frame(
    id = c(1, 1, 2, 2, 2),
    tstart = c(0, .5, 0, .8, 1.2),
    tstop = c(.5, 3, .8, 1.2, 3),
    status = c(1, 0, 1, 1, 0),
    enum = c(1, 2, 1, 2, 3),
    age = c(50, 50, 24, 24, 24)
  )
  # GAP timescale
  gap_df <- as_ped(
    data = test_df,
    formula = Surv(tstart, tstop, status) ~ enum + age,
    transition = "enum",
    id = "id",
    timescale = "gap"
  )

  expect_subset(
    names(attributes(gap_df)),
    c(
      "names",
      "row.names",
      "class",
      "breaks",
      "id_var",
      "intvars",
      "trafo_args",
      "time_var"
    )
  )
  expect_subset(
    names(attributes(mutate(gap_df, age = 10))),
    c(
      "names",
      "row.names",
      "class",
      "breaks",
      "id_var",
      "intvars",
      "trafo_args",
      "time_var"
    )
  )
  expect_subset(
    names(attributes(group_by(gap_df, id))),
    c(
      "names",
      "row.names",
      "class",
      "breaks",
      "id_var",
      "intvars",
      "trafo_args",
      "time_var",
      "groups"
    )
  )
  expect_subset(
    names(attributes(ungroup(group_by(gap_df, id)))),
    c(
      "names",
      "row.names",
      "class",
      "breaks",
      "id_var",
      "intvars",
      "trafo_args",
      "time_var"
    )
  )
})

test_that("inner_join and full_join preserve ped attributes", {
  data("tumor")
  ped <- as_ped(
    data = dplyr::slice(tumor, 2:4),
    formula = Surv(days, status) ~ age,
    cut = c(0, 100, 400),
    id = "id"
  )

  ped_keys <- distinct(ped, id, interval)
  joined_inner <- inner_join(ped, ped_keys, by = c("id", "interval"))
  joined_full <- full_join(ped, ped_keys, by = c("id", "interval"))

  expect_is(joined_inner, "ped")
  expect_is(joined_full, "ped")
  expect_identical(attr(joined_inner, "id_var"), attr(ped, "id_var"))
  expect_identical(attr(joined_full, "id_var"), attr(ped, "id_var"))
})

test_that("inner and full joins match dplyr row semantics", {
  data("tumor")
  ped <- as_ped(
    data = dplyr::slice(tumor, 2:4),
    formula = survival::Surv(days, status) ~ age,
    cut = c(0, 100, 400),
    id = "id"
  )

  keys <- dplyr::distinct(ped, id, interval)
  out_inner <- dplyr::inner_join(ped, keys, by = c("id", "interval"))
  out_full <- dplyr::full_join(ped, keys, by = c("id", "interval"))

  ref_inner <- dplyr::inner_join(
    as.data.frame(ped),
    keys,
    by = c("id", "interval")
  )
  ref_full <- dplyr::full_join(
    as.data.frame(ped),
    keys,
    by = c("id", "interval")
  )

  expect_equal(nrow(out_inner), nrow(ref_inner))
  expect_equal(nrow(out_full), nrow(ref_full))
  expect_identical(attr(out_inner, "id_var"), attr(ped, "id_var"))
  expect_identical(attr(out_full, "id_var"), attr(ped, "id_var"))
})

test_that("right_join.ped keeps unmatched right-hand keys", {
  data("tumor")
  ped <- as_ped(
    data = dplyr::slice(tumor, 2:4),
    formula = survival::Surv(days, status) ~ age,
    cut = c(0, 100, 400),
    id = "id"
  )

  keys <- dplyr::distinct(ped, id, interval) |>
    dplyr::slice(1:2)
  keys <- dplyr::bind_rows(
    keys,
    data.frame(
      id = 999L,
      interval = factor(
        levels(ped$interval)[1],
        levels = levels(ped$interval)
      )
    )
  )

  out_right <- dplyr::right_join(ped, keys, by = c("id", "interval"))
  ref_right <- dplyr::right_join(
    as.data.frame(ped),
    keys,
    by = c("id", "interval")
  )

  expect_equal(nrow(out_right), nrow(ref_right))
  expect_true(any(out_right$id == 999L))
  expect_true(any(is.na(out_right$ped_status[out_right$id == 999L])))
})

test_that("ped attributes are sanitized after structural column drops", {

  data("tumor")
  ped <- as_ped(
    data = dplyr::slice(tumor, 2:3),
    formula = Surv(days, status) ~ age,
    cut = c(0, 100, 400),
    id = "id"
  )
  ped_no_id <- dplyr::select(ped, -id)

  expect_null(attr(ped_no_id, "id_var"))
  expect_false("id" %in% attr(ped_no_id, "intvars"))
  expect_data_frame(sample_info(ped_no_id), nrows = 1L)

})
