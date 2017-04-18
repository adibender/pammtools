context("Interal info and median and modus information")

data("leuk2", package="bpcp")
leuk2$x <- rnorm(nrow(leuk2))
leuk.ped <- split_data(Surv(time, status)~., data=leuk2, cut=c(0:5, 10, 40), id="id")


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

test_that("Sample info returned for data frame", {
	expect_is(si <- sample_info(leuk2), "data.frame")
	expect_equal(dim(si), c(1L, 5L))
})

test_that("Sample info returned for ped objects", {
	expect_is(si2 <- sample_info(leuk.ped), "data.frame")
	expect_equal(dim(si2), c(1L, 3L))
})

test_that("Sample info returned for grouped ped objects", {
	expect_is(si2.ped <- group_by(leuk.ped, treatment) %>% sample_info(), "tbl_df")
	expect_equal(dim(si2.ped), c(2, 3))
})


test_that("ped info returned for (grouped) ped objects", {
	expect_is(si.ped <- ped_info(leuk.ped), "data.frame")
	expect_equal(dim(si.ped), c(7, 8))
	expect_is(si2.ped <- group_by(leuk.ped, treatment) %>% ped_info(), "data.frame")
	expect_equal(dim(si2.ped), c(14, 8))
})