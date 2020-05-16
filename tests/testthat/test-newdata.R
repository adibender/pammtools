context("create newdata")

test_that("creating newdata works on ungrouped data", {

  iris2 <- iris %>%  group_by(Species) %>% slice(1:2) %>% ungroup()
  expect_data_frame(
    make_newdata(iris2),
    any.missing = FALSE, nrows = 1L, ncols = 5L)
  expect_equal(colnames(make_newdata(iris2)), colnames(iris2))
  expect_data_frame(
    make_newdata(iris2, Sepal.Length = c(5)),
    any.missing = FALSE, nrows = 1L, ncols = 5L)
  expect_equal(make_newdata(iris2, Sepal.Length = c(5))$Sepal.Length, 5)
  expect_data_frame(
    make_newdata(iris2, Sepal.Length = c(5, 6)),
    any.missing = FALSE, nrows = 2L, ncols = 5L)
  expect_data_frame(
    make_newdata(iris2, Sepal.Length = seq_range(Sepal.Length, 2)),
    any.missing = FALSE, nrows = 2L, ncols = 5L)
  expect_equal(
    make_newdata(iris2, Sepal.Length = seq_range(Sepal.Length, 2))$Sepal.Length,
    c(4.9, 7.0))
})


test_that("creating newdata fails on ungrouped data", {

  iris2 <- iris %>% group_by(Species) %>% slice(2) %>% ungroup()

  expect_warning(make_newdata(iris2, Sepal.length = c(5)))
  expect_error(make_newdata(iris2, Sepal.Length = 5))
  expect_error(make_newdata(iris2, Sepal.Length = seq_range(Sepal.length, 2)))
  expect_warning(make_newdata(iris2, Sepal.length = seq_range(Sepal.Length, 2)))

})


test_that("make_newdata works for PED data", {

  ped <- simdf_elra %>%
    slice(1:6) %>%
    as_ped(Surv(time, status)~x1 + x2, cut = seq(0, 10, by = 5))
  mdf <- ped %>% make_newdata(x1 = seq_range(x1, 2))
  expect_data_frame(mdf, nrows = 2L, ncols = 9L)
  expect_equal(mdf$tend, c(5, 5))
  expect_equal(mdf$x1, c(-2.43, 2.54), tolerance = 1e-2)
  expect_message(make_newdata(ped, tend = c(2.5)))
  mdf <- ped %>% make_newdata(tend = c(10), x1 = seq_range(x1, 2))
  expect_data_frame(mdf, nrows = 2L, ncols = 9L)
  mdf <- ped %>% make_newdata(x1 = seq_range(x1, 2), x2 = seq_range(x2, 2))
  expect_data_frame(mdf, nrows = 4L, ncols = 9L)
  mdf <- ped %>% make_newdata(tend = unique(tend), x2 = seq_range(x2, 2))
  expect_data_frame(mdf, nrows = 4L, ncols = 9L)

})


test_that("make_newdata works for PED with matrix columns", {

  ped_simdf <- simdf_elra %>% as_ped(
    Surv(time, status) ~ x1 + x2 +
      cumulative(time, latency(tz1), z.tz1, tz_var = "tz1",
        ll_fun = function(t, tz) t >= tz + 2) +
      cumulative(latency(tz2), z.tz2, tz_var = "tz2"),
    cut = 0:10)

  ## sample info
  expect_data_frame(sdf <- sample_info(ped_simdf), nrows = 1, ncols = 2)
  expect_equal(sdf$x1, 0.0718, tolerance = 1e-3)
  expect_equal(sdf$x2, 3.043, tolerance = 1e-3)

  ## ped info
  pinf <- ped_info(ped_simdf)
  expect_data_frame(pinf, nrows = 10L, ncols = 7L)
  expect_equal(pinf$x1[1], 0.0718, tolerance = 1e-3)
  expect_equal(pinf$x2[2], 3.043, tolerance = 1e-3)

  # make newdata
  nd1 <- ped_simdf %>% make_newdata(x1 = c(0.05))
  expect_data_frame(nd1, nrows = 1L, ncols = 16L)
  expect_equal(nd1$tstart, 0)
  expect_equal(nd1$tend, 1)
  expect_equal(nd1$x1, 0.05)
  expect_equal(nd1$x2, 2.65, tolerance = 1e-3)
  expect_equal(nd1$z.tz1_tz1, -0.370, 1e-3)

  nd2 <- ped_simdf %>% make_newdata(x1 = seq_range(x1, 2))
  expect_data_frame(nd2, nrows = 2L, ncols = 16L)
  expect_equal(nd2$x1[1], min(unlist(simdf_elra$x1)))
  expect_equal(nd2$x1[2], max(unlist(simdf_elra$x1)))

  nd3 <- ped_simdf %>% make_newdata(tend = unique(tend))
  expect_data_frame(nd3, nrows = 10L, ncols = 16L)
  expect_equal(nd3$tend, 1:10)

  nd4 <- ped_simdf %>% make_newdata(tz1_latency = c(0:5))
  expect_data_frame(nd4, nrows = 6L, ncols = 16L)
  expect_equal(nd4$tz1_latency, 0:5)

  nd5 <- ped_simdf %>%
    make_newdata(
      tend = c(1:10),
      tz1_latency = seq(1:5))
  expect_data_frame(nd5, nrows = 50L, ncols = 16L)
  expect_equal(nd5$tend, rep(1:10, 5L))
  expect_equal(nd5$tz1_latency, rep(1:5, each = 10L))
  expect_equal(nd5$LL_tz1, c(rep(0, 10), rep(1, nrow(nd5) - 10)))

})

test_that("Errors are thrown", {

  expect_error(combine_df(data.frame(x = 1), x = 2))

})
