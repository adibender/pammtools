context("Tidyverse methods for specific classes")

test_that("ped class is preserved after dplyr operations", {

  data("tumor")
  tumor <- dplyr::slice(tumor, 2:3)
  ped <- as_ped(
    data    = tumor,
    formula = Surv(days, status) ~ complications + age,
    cut     = c(0, 100, 400),
    id      = "id")

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
  expect_is(distinct(sim_df, id), "nested_fdf")
  expect_is(sample_n(sim_df, 2), "nested_fdf")
  expect_is(sample_frac(sim_df, .1), "nested_fdf")
  expect_is(select(sim_df, id), "nested_fdf")
  expect_is(mutate(sim_df, id = id + 1), "nested_fdf")
  expect_is(rename(sim_df, ID = id), "nested_fdf")
  expect_is(summarise(sim_df, id = mean(id)), "nested_fdf")
  expect_is(left_join(sim_df, distinct(sim_df, id)), "nested_fdf")
  expect_is(right_join(distinct(sim_df, id), sim_df), "nested_fdf")
  expect_is(full_join(distinct(sim_df, id), sim_df), "nested_fdf")
  expect_is(inner_join(distinct(sim_df, id), sim_df), "nested_fdf")
  sim_df$id[1] <- NA
  expect_is(fill(sim_df, id, .direction = "up"), "nested_fdf")

})


test_that("attributes are preserved", {
  # recurrent events data
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

  expect_subset(names(attributes(gap_df)), c("names", "row.names", "class",
    "breaks", "id_var", "intvars", "trafo_args", "time_var"))
  expect_subset(
    names(attributes(mutate(gap_df, age = 10))),
    c("names", "row.names", "class", "breaks", "id_var", "intvars", "trafo_args",
    "time_var"))
  expect_subset(
    names(attributes(group_by(gap_df, id))),
    c("names", "row.names", "class", "breaks", "id_var", "intvars", "trafo_args",
    "time_var", "groups"))
  expect_subset(
    names(attributes(ungroup(group_by(gap_df, id)))),
    c("names", "row.names", "class", "breaks", "id_var", "intvars", "trafo_args",
    "time_var"))

})
