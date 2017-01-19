context("Extract median and modus information")

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