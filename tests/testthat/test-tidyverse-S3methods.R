context("Tidyverse methods for specific classes")

test_that("ped class is preserved after dplyr operations", {
  data("veteran", package = "survival")
  veteran <- dplyr::slice(veteran, 1:2)
  ped <- as_ped(
    data    = veteran,
    formula = Surv(time, status) ~ trt + age,
    cut     = c(0, 100, 400),
    id      = "id")

  expect_is(group_by_(ped, "interval"), "ped")
  expect_is(distinct_(ped, "interval"), "ped")
  expect_is(filter(ped, id == 1), "ped")
  expect_is(slice(ped, 1), "ped")
  expect_is(arrange(ped, desc(id)), "ped")
  expect_is(select(ped, id), "ped")
  expect_is(select_(ped, "id"), "ped")
  expect_is(rename(ped, ID = id), "ped")
  expect_is(rename_(ped, "ID" = "id"), "ped")
  expect_is(mutate(ped, id = id + 1), "ped")
  expect_is(transmute(ped, id = id + 1), "ped")
  expect_is(sample_n(ped, 1), "ped")
  expect_is(sample_frac(ped, 0.5), "ped")
  expect_is(right_join(distinct(ped, id, interval), ped), "ped")
})



test_that("ped class is preserved after tidyr operations", {

  data("pbc", package = "survival")
  tdc_df <- pbcseq %>%
    filter(id %in% 1:5) %>%
    select(id, day, chol)
  class(tdc_df) <- c("ped", class(tdc_df))

  ## standard version
  temp <- tidyr::fill(tdc_df, chol)
  expect_is(temp, "ped")
  expect_identical(sum(is.na(temp)), 0L)

  ## non standard evaluation
  temp <- tidyr::fill(tdc_df, c("chol"))
  expect_is(temp, "ped")
  expect_identical(sum(is.na(temp)), 0L)

})

test_that("nested_fdf class is preserved after tidyr operations", {

  expect_is(sim_df <- filter(simdf_elra, id %in% c(1:2)), "nested_fdf")
  expect_is(arrange(sim_df, id), "nested_fdf")
  expect_is(group_by(sim_df, id), "nested_fdf")
  expect_is(group_by_(sim_df, "id"), "nested_fdf")
  expect_is(distinct(sim_df, id), "nested_fdf")
  expect_is(distinct_(sim_df, "id"), "nested_fdf")
  expect_is(sample_n(sim_df, 2), "nested_fdf")
  expect_is(sample_frac(sim_df, .1), "nested_fdf")
  expect_is(select(sim_df, id), "nested_fdf")
  expect_is(select_(sim_df, "id"), "nested_fdf")
  expect_is(mutate(sim_df, id = id + 1), "nested_fdf")
  expect_is(rename(sim_df, ID = id), "nested_fdf")
  expect_is(rename_(sim_df, "ID" = "id"), "nested_fdf")
  expect_is(summarise(sim_df, id = mean(id)), "nested_fdf")
  expect_is(left_join(sim_df, distinct(sim_df, id)), "nested_fdf")
  expect_is(right_join(distinct(sim_df, id), sim_df), "nested_fdf")
  expect_is(full_join(distinct(sim_df, id), sim_df), "nested_fdf")
  expect_is(inner_join(distinct(sim_df, id), sim_df), "nested_fdf")
  sim_df$id[1] <- NA
  expect_is(fill(sim_df, id, .direction = "up"), "nested_fdf")

})
