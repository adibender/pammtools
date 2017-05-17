context("Simple transformation to PED data")

data("veteran", package="survival")
veteran <- veteran[c(1:3, 135:137), ]

test_that("Output as expected without id", {
	expect_is(ped <- split_data(Surv(time, status)~ trt + age, data=veteran, 
		cut=c(0, 100, 400)), "ped")
	expect_equal(dim(ped), c(10L, 7L))
	expect_is(ped, "data.frame")
	expect_true(all(c("ped_status", "tstart", "tend", "interval", "offset") %in% 
		names(ped)))
	expect_is(attr(ped, "cut"), "numeric")
	expect_is(attr(ped, "intvars"), "character")
})



test_that("Output as expected with id", {
	expect_is(ped <- split_data(Surv(time, status)~ trt +age, data=veteran, 
		cut=c(0, 100, 400), id="id"), "ped")
	expect_equal(dim(ped),c(10L, 8L))
})

test_that("ID kept when id variable excluded in formula", {
	expect_is(ped <- split_data(Surv(time, status)~trt, data=veteran, 
		cut=c(0, 100, 400), id="id"), "ped")
	expect_equal(ncol(ped), 7L)
})


test_that("Error on wrong input", {
	expect_error(split_data())
	expect_error(split_data(x~y, data=veteran, cut=c(0:5, 10, 40)))
	expect_error(split_data(Surv(time2, status)~., data=veteran, cut=c(0:5, 10, 40)))
	expect_error(split_data(Surv(ped_time, status)~., 
		data=rename(veteran, ped_time=time))) # already in data set ped_time
})