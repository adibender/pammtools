context("Transformation with TDC")


test_that("Concurrent TDC are transformed correctly", {
  data("pbc", package = "survival")
  # default case with lag = 0
  event_df <- filter(pbc, id %in% 1:3) %>% mutate(status = 1L*(status == 1))
  tdc_df <- filter(pbcseq, id %in% 1:3) %>% select(id, day, bili, protime)
  time <- sort(unique(event_df$time))[1:2]
  tz <- sort(unique(tdc_df$day))
  tz <- tz[tz <= max(time)][-1]
  expect_error(as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~. + concurrent(bili, protime, tz_var = "day"),
    id      = "id"), "No events in data")
  event_df <- filter(pbc, id %in% 1:3) %>% mutate(status = status == 2) %>%
    select(id, time, status, trt, age, bili, spiders)
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~. +
      concurrent(bili, protime, tz_var = "day"), id = "id")
  expect_equal(unique(ped$tend), c(176, 182, 192, 364, 365, 400, 743, 768, 1012))
  expect_equal(ped$bili,
    c(rep(14.5, 3), rep(21.3, 3), rep(1.1, 2), rep(0.8, 3), rep(1, 3),
      1.9, 1.4, rep(1.1, 3), rep(1.5, 3), rep(1.8, 2)))
  # lag != 0
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~. +
      concurrent(bili, protime, tz_var = "day", lag = 10),
    id = "id")
  expect_equal(
    unique(ped$tend),
    sort(c(time, tz + 10)))
  expect_equal(ped$bili,
    c(rep(14.5, 3), rep(21.3, 3), rep(1.1, 2), rep(0.8, 3), rep(1, 3),
      1.9, 1.4, rep(1.1, 3), rep(1.5, 3), rep(1.8, 2)))
  # unequal lags
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~. +
      concurrent(bili, tz_var = "day", lag = 10) +
      concurrent(protime, tz_var = "day", lag = 0),
    id = "id")
  expect_data_frame(ped, nrows = 40, ncols = 11)
  expect_equal(sum(ped$ped_status), 2)
  expect_equal(sort(unique(ped$tend)), sort(unique(c(time, tz, tz+10))))
  expect_equal(ped$bili,
               c(rep(14.5, 5), rep(21.3, 5), rep(1.1, 4), rep(0.8, 5), rep(1, 5),
                 1.9, rep(1.4, 3), rep(1.1, 5), rep(1.5, 4), rep(1.8, 3)))
  expect_equal(ped$protime,
               c(rep(12.2, 4), rep(11.2, 6), rep(10.6, 2), rep(11, 5),
                 rep(11.6, 6), rep(10.6, 2), rep(12, 11), rep(13.3, 4)))
# when maxtime is set
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status)~. + concurrent(bili, protime, tz_var = "day"),
    id      = "id",
    max_time = 1400)
  expect_equal(unique(ped$tend), sort(c(time, tz, 1400)))
  expect_equal(ped$bili,
    c(rep(14.5, 3), rep(21.3, 3), rep(1.1, 2), rep(0.8, 3), rep(1.0, 3), rep(1.9, 2),
      1.4, rep(1.1, 3), rep(1.5, 3), rep(1.8, 2)))
})

test_that("Covariate matrices are created correctly", {
  data <- simdf_elra %>% filter(id %in% c(1:2))
  time <- 0:2
  tz <- data %>% dplyr::pull("tz2") %>% unlist() %>% unique() %>% sort()
  nz <- length(tz)
  attr(data, "id_tseq") <- rep(1:3, 2)
  attr(data, "id_tz_seq") <- rep(1:2, times = c(3, 3))
  my_ll_fun <- function(t, tz) ( (t - tz) >= 0 & (t - tz) <= 5)
  expect_class(my_ll_fun, "function")
  Tmat <- make_time_mat(data, nz)
  TEmat <- make_z_mat(data, "tz2", nz)
  Ltmat <- make_latency_mat(data, tz)
  LLmat <- make_lag_lead_mat(data, tz, ll_fun = my_ll_fun)
  expect_equal(dim(Tmat), c(6, 11))
  expect_equal(dim(TEmat), c(6, 11))
  expect_equal(dim(Ltmat), c(6, 11))
  expect_equal(dim(LLmat), c(6, 11))
  expect_equal(all(Tmat[1, ] == 0), TRUE)
  expect_equal(all(Tmat[2, ] == 1), TRUE)
  expect_equal(all(Tmat[3, ] == 2), TRUE)
  expect_equal(all(TEmat[, 1] == -5), TRUE)
  expect_equal(all(TEmat[, 11] == 5), TRUE)
  expect_equal(Ltmat[1, ], c(5:0, rep(0, 5)))
  expect_equal(Ltmat[3, ], c(7:0, rep(0, 3)))
  expect_equal(LLmat[1, ], c(rep(1, 6), rep(0, 5)))
  expect_equal(LLmat[3, ], c(rep(0, 2), rep(1, 6), rep(0, 3)))
  expect_equal(max(Ltmat * LLmat), 5)
  ped <- as_ped(data,
      Surv(time, status) ~ . +
      cumulative(z.tz2, latency(tz2), tz_var = "tz2",
        ll_fun = function(t, tz) (t - tz) >= 0 & (t - tz) <= 5),
      cut = 0:2)
  expect_equal(max(ped$tz2_latency * ped$LL), 5)

})
