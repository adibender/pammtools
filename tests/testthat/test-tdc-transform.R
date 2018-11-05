context("Transformation with TDC")


test_that("Concurrent TDC are transformed correctly", {
  data("pbc", package = "survival")
  event_df <- filter(pbc, id %in% 1:3) %>% mutate(status = status == 1)
  tdc_df <- filter(pbcseq, id %in% 1:3) %>% select(id, day, bili, protime)
  expect_error(as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~. | concurrent(bili, protime, tz_var = "day"),
    id      = "id"), "No events in data")
  event_df <- filter(pbc, id %in% 1:3) %>% mutate(status = status == 2)
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status) ~. | concurrent(bili, protime, tz_var = "day"),
    id      = "id")
  expect_equal(unique(ped$tend)[1:7], c(176, 182, 192, 364, 365, 400, 743))
  ped <- as_ped(
    data    = list(event_df, tdc_df),
    formula = Surv(time, status)~. | concurrent(bili, protime, tz_var = "day"),
    id      = "id",
    cut = c(0, 3000))
  expect_equal(unique(ped$tend)[1:7], c(176, 182, 192, 364, 365, 743, 768))


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
      Surv(time, status) ~ . |
      cumulative(z.tz2, latency(tz2), tz_var = "tz2",
        ll_fun = function(t, tz) (t - tz) >= 0 & (t - tz) <= 5),
      cut = 0:2)
  expect_equal(max(ped$tz2_latency * ped$LL), 5)

})

test_that("split_tdc works correctly", {
  data("pbc", package = "survival")

  event_df <- pbc %>%
    filter(id <= 5) %>%
    mutate(event = 1L * (status == 2)) %>%
    select(id, time, event, sex, bili)
  tdc_df <- pbcseq %>%
    filter(id <= 5) %>%
    select(id, day, bili)

  expect_warning(pbc_ped <- split_tdc(Surv(time, event)~., event_df, tdc_df,
    tz_var = "day", status_var = "event"), "is deprecated")


  expect_data_frame(pbc_ped, nrows = 93L, ncols = 8L)
  expect_is(pbc_ped, "ped")
  expect_subset(c("breaks", "id_var", "intvars"), names(attributes(pbc_ped)))

})
