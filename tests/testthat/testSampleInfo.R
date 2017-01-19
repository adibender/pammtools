context("Interal info and median and modus information")

data("leuk2", package="bpcp")
leuk.ped <- split_data(Surv(time, status)~., data=leuk2, cut=c(0:5, 10, 40), id="id")

test_that("Sample info returned for data frame", {
	expect_is(si <- sample_info(leuk2), "data.frame")
	expect_equal(nrow(si), 1L)
	expect_equal(ncol(si), 4L)
})

test_that("Sample info returned for ped objects", {
	expect_is(si2 <- sample_info(leuk.ped), "data.frame")
	expect_equal(nrow(si2), 1L)
	expect_equal(ncol(si2), 2L)
})


test_that("Interval infos correct", {
	expect_equal(nrow(int_info(1:2)), 2)
	expect_equal(ncol(int_info(1:2)), 5)
	expect_equal(names(int_info(1:2)), c("tstart", "tend", "intlen", "intmid", "interval"))
	expect_equal(levels(int_info(1:2)$interval), c("(0,1]", "(1,2]"))
})

test_that("Interval info returned for ped objects", {
	expect_is(int_info(leuk.ped), "data.frame")
	expect_equal(nrow(int_info(leuk.ped)), 7)
	expect_equal(ncol(int_info(leuk.ped)), 5)
})