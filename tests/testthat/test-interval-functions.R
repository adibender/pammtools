context("Interval info and median and modus information")

data("tumor")
ped <- tumor[1:200, ] %>% as_ped(Surv(days, status)~ complications + age,
	cut=seq(0,400, by=100))

ped <- filter(ped, id %in% c(1:3, 135:137))


test_that("Interval infos correct", {
	expect_data_frame(int_info(1:2), nrows=2L, ncols=5L)
	expect_data_frame(int_info(2:1), nrows=2L, ncols=5L)
	expect_data_frame(int_info(data.frame(x1 = c(1,0), x2=c(2, 1))), nrows = 2L, ncols = 5L)
	expect_equal(names(int_info(1:2)), c("tstart", "tend", "intlen", "intmid", "interval"))
	expect_equal(levels(int_info(1:2)$interval), c("(0,1]", "(1,2]"))
})

test_that("Interval info returned for ped objects", {
	expect_data_frame(int_info(ped), nrows=4L, ncols=5L, types=c("numeric", "factor"))
})

test_that("Sample info returned for data frame", {
	expect_data_frame(si <- sample_info(tumor), nrows=1L, ncols=9L)
	expect_equal(colnames(si), colnames(tumor))
	expect_data_frame(si <- tumor %>% group_by(complications, status) %>% sample_info(),
		nrows=4L, ncols=9L)
	expect_equal(colnames(si), colnames(tumor))
})

test_that("Sample info returned for ped objects", {
	expect_data_frame(sample_info(ped), nrows=1, ncols=2)
})

test_that("Sample info returned for grouped ped objects", {
	expect_data_frame(group_by(ped, complications) %>% sample_info(),
		nrows=2, ncols=2)
})


test_that("ped info returned for (grouped) ped objects", {
	# normal
	expect_data_frame(ped_info(ped), nrows=4L, ncols=7L)
	#grouped
	expect_data_frame(group_by(ped, complications) %>% ped_info(), nrows=8L, ncols=7L)
	# without covariates
	expect_data_frame(ped_info(select(ped, -complications, -age)), nrows=4L, ncols=5L)
})
