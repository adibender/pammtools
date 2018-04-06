context("Interal info and median and modus information")

data("veteran", package="survival")
ped <- veteran %>% as_ped(Surv(time, status)~ trt + age,
	cut=seq(0,400, by=100), id="id")

ped <- filter(ped, id %in% c(1:3, 135:137))


test_that("Interval infos correct", {
	expect_data_frame(int_info(1:2), nrows=2L, ncols=5L)
	expect_equal(names(int_info(1:2)), c("tstart", "tend", "intlen", "intmid", "interval"))
	expect_equal(levels(int_info(1:2)$interval), c("(0,1]", "(1,2]"))
})

test_that("Interval info returned for ped objects", {
	expect_data_frame(int_info(ped), nrows=4L, ncols=5L, types=c("numeric", "factor"))
})

test_that("Sample info returned for data frame", {
	expect_data_frame(si <- sample_info(veteran), nrows=1L, ncols=8L)
	expect_equal(colnames(si), colnames(veteran))
	expect_data_frame(si <- veteran %>% group_by(trt, status) %>% sample_info(),
		nrows=4L, ncols=8L)
	expect_equal(colnames(si), colnames(veteran))
})

test_that("Sample info returned for ped objects", {
	expect_data_frame(sample_info(ped), nrows=1, ncols=2)
})

test_that("Sample info returned for grouped ped objects", {
	expect_data_frame(group_by(ped, trt) %>% sample_info(),
		nrows=2, ncols=2)
})


test_that("ped info returned for (grouped) ped objects", {
	# normal
	expect_data_frame(ped_info(ped), nrows=4L, ncols=7L)
	#grouped
	expect_data_frame(group_by(ped, trt) %>% ped_info(), nrows=8L, ncols=7L)
	# without covariates
	expect_data_frame(ped_info(select(ped, -trt, -age)), nrows=4L, ncols=5L)
})

test_that("riskset info returned for (grouped) ped objects", {
  expect_data_frame(riskset_info(ped), nrows=4L, ncols=4L)
  expect_data_frame(group_by(ped, trt) %>% riskset_info(), nrow=8L, ncols=5L)
})
