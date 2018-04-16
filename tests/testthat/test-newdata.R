context("create newdata")


test_that("creating newdata works on ungrouped data", {
	set.seed(123)
	iris2 <- iris %>%  group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_data_frame(make_newdata(iris2), any.missing=FALSE, nrows=1, ncols=5)
	expect_equal(colnames(make_newdata(iris2)), colnames(iris2))
	expect_data_frame(make_newdata(iris2, Sepal.Length=c(5)), any.missing=FALSE,
		nrows=1, ncols=5)
	expect_equal(make_newdata(iris2, Sepal.Length=c(5))$Sepal.Length, 5)
	expect_data_frame(make_newdata(iris2, Sepal.Length=c(5, 6)),
		any.missing=FALSE, nrows=2, ncols=5)
	expect_data_frame(make_newdata(iris2, Sepal.Length=seq_range(Sepal.Length, 2)),
		any.missing=FALSE, nrows=2, ncols=5)
	expect_equal(make_newdata(iris2, Sepal.Length=seq_range(Sepal.Length, 2))$Sepal.Length,
		c(4.4, 7.1))
})


test_that("creating newdata fails on ungrouped data", {
	set.seed(123)
	iris2 <- iris %>% group_by(Species) %>% sample_n(2) %>% ungroup()

  expect_warning(make_newdata(iris2, Sepal.length=c(5)))
	expect_error(make_newdata(iris2, Sepal.Length=5))
  expect_error(make_newdata(iris2, Sepal.Length=seq_range(Sepal.length, 2)))
	expect_warning(make_newdata(iris2, Sepal.length=seq_range(Sepal.Length, 2)))

})


test_that("make_newdata works for PED data", {
  ped <- simdf_elra %>% slice(1:6) %>% as_ped(Surv(time, status)~x1 + x2,
    cut = seq(0, 10, by=5))
  mdf <- ped %>% make_newdata(x1 = seq_range(x1, 2))
  expect_data_frame(mdf, nrows = 2L, ncols=8)
  expect_equal(mdf$tend, c(5, 5))
  expect_equal(mdf$x1, c(-2.43, 2.54), tolerance = 1e-2)
  mdf <- ped %>% make_newdata(tend = c(10), x1 = seq_range(x1, 2))
  expect_data_frame(mdf, nrows = 2L, ncols = 8L)
  mdf <- ped %>% make_newdata(x1 = seq_range(x1, 2), x2 = seq_range(x2, 2))
  expect_data_frame(mdf, nrows=4L, ncols = 8L)
  mdf <- ped %>% make_newdata(tend=unique(tend), x2 = seq_range(x2, 2))
  expect_data_frame(mdf, nrows = 4L, ncols = 8L)

})


test_that("make_newdata works for PED with matrix columns", {
  # library(mgcv)
  ped_simdf <- simdf_elra %>% as_ped(
    Surv(time, status)~ x1 + x2|
      cumulative(time, latency(te1), z.te1, te_var="te1") +
      cumulative(latency(te2), z.te2, te_var="te2"),
    cut = 0:10)

  expect_data_frame(sdf <- sample_info(ped_simdf), nrows=1, ncols=2)
  expect_equal(sdf$x1, 0.0718, tolerance = 1e-3)
  expect_equal(sdf$x2, 3.043, tolerance = 1e-3)

})
