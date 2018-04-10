context("create newdata")


test_that("creating newdata works on ungrouped data", {
	set.seed(123)
	iris2 <- iris %>%  group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_data_frame(make_newdata(iris2), any.missing=FALSE, nrows=1, ncols=5)
	expect_equal(colnames(make_newdata(iris2)), colnames(iris2))
	expect_data_frame(make_newdata(iris2, Sepal.Length=5), any.missing=FALSE,
		nrows=1, ncols=5)
	expect_equal(make_newdata(iris2, Sepal.Length=5)$Sepal.Length, 5)
	expect_data_frame(make_newdata(iris2, Sepal.Length=c(5, 6)),
		any.missing=FALSE, nrows=2, ncols=5)
	expect_data_frame(make_newdata(iris2, expand="Sepal.Length", n = 2),
		any.missing=FALSE, nrows=2, ncols=5)
	expect_equal(make_newdata(iris2, expand="Sepal.Length", n = 2)$Sepal.Length,
		c(4.4, 7.1))
})


test_that("creating newdata fails on ungrouped data", {
	set.seed(123)
	iris2 <- iris %>% group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_error(make_newdata(iris2, Sepal.length=5))
	expect_error(make_newdata(iris2, expand="Sepal.length"))
	expect_error(make_newdata(iris2, expand="Sepal.Length", n = -3))

})

test_that("make_newdata.ped warns about intervals", {
  int_df <-  data.frame(time = 1:2, status = c(1,1)) %>%
    as_ped(Surv(time, status)~., id = "id")
  expect_warning(
    make_newdata(int_df, tstart = 1.5),
    "Setting interval borders")
  expect_warning(
    make_newdata(int_df, expand = "tend", n = 3),
    "Setting interval borders")
  expect_warning(
    make_newdata(int_df, tstart = 1.5, expand = "tend", n = 3),
    "Setting interval borders")
})


test_that("make_newdata works for PED with matrix columns", {
  # library(mgcv)
  ped_simdf <- simdf_elra %>% as_ped(
    Surv(time, status)~ x1 + x2|
      cumulative(time, latency(te1), z.te1, te_var="te1") +
      cumulative(latency(te2), z.te2, te_var="te2"),
    cut = 0:10)
  # mod_simdf <- gam(
  #   formula = ped_status ~ ti(tend, mc=1) + s(x1) + s(x2) +
  #     ti(time_te1, te1_latency, z.te1_te1, by = LL_te1, mc=c(1,1,1)) +
  #     ti(te2_latency, z.te2_te2, by = LL_te2, mc=c(1,1)),
  #   data = ped_simdf,
  #   family = poisson(),
  #   offset = offset)

  expect_data_frame(sdf <- sample_info(ped_simdf), nrows=1, ncols=2)
  expect_equal(sdf$x1, 0.0718, tolerance = 1e-3)
  expect_equal(sdf$x2, 3.043, tolerance = 1e-3)

  ndf <- mk_ndf(ped_simdf, x1=seq_range(x1, n=2), z.te1_te1 = c(0, 3, 5))
  # lp_mat <- predict(mod_simdf, newdata = ndf, type = "lpmatrix")


})
