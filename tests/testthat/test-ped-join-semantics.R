context("ped join semantics")

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
