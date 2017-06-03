context("Simple transformation to PED data")

test_that("Output as expected without and with id", {
	## preparations 
	data("veteran", package="survival")
	veteran <- veteran[c(1:3, 135:137), ]

	## tests 
	expect_data_frame(ped <- split_data(Surv(time, status)~ trt + age, 
		data=veteran, cut=c(0, 100, 400)), nrows=10L, ncols=8L)
	expect_is(ped, "ped")
	expect_subset(c("ped_status", "tstart", "tend", "interval", "offset"), names(ped))
	expect_is(attr(ped, "cut"), "numeric")
	expect_is(attr(ped, "intvars"), "character")
	expect_is(attr(ped, "id_var"), "character")
	expect_equal(attr(ped, "id_var"), "id")
	expect_data_frame(split_data(Surv(time, status)~ trt +age, 
		data=veteran, cut=c(0, 100, 400), id="id"), nrows=10L, ncols=8L)
	## no error, when id in data and additionally specified
	veteran$id <- seq_len(nrow(veteran))
	expect_data_frame(ped <- split_data(Surv(time, status)~trt, 
		data=veteran, cut=c(0, 100, 400), id="id"), nrows=10L, ncols=7L)
	## no error when id already in data but not specified
	expect_data_frame(split_data(Surv(time, status)~trt, 
		data=veteran, cut=c(0, 100, 400)), nrows=10L, ncols=7L)
	## no error when id has different name and is specified accordingly 
	veteran$id2 = veteran$id
	expect_data_frame(ped <- split_data(Surv(time, status)~., 
		data=veteran, cut=c(0, 100, 400)), nrows=10L, ncols=13L)
	expect_identical(attr(ped, "id_var"), "id")
	## no additional id when different id specified
	veteran$id <- NULL
	expect_data_frame(ped <- split_data(Surv(time, status)~., 
		data=veteran, cut=c(0, 100, 400), id="id2"), nrows=10L, ncols=12L)
	expect_identical(attr(ped, "id_var"), "id2")
})


test_that("Error on wrong input", {
	## preparations 
	data("veteran", package="survival")
	veteran <- veteran[c(1:3, 135:137), ]

	## tests
	expect_error(split_data(x~y, data=veteran, cut=c(0:5, 10, 40)))
	expect_error(split_data(Surv(time2, status)~., data=veteran, cut=c(0:5, 10, 40)))
	expect_error(split_data(Surv(ped_time, status)~., 
		data=rename(veteran, ped_time=time))) # already in data set ped_time

	## error when specified id variable not unique
	veteran$id <- rep(1:2, 3)
	expect_error(
		split_data(Surv(time, status)~trt, data=veteran, cut=c(0, 100, 400), id="id"), 
		regexp="Specified ID variable.*")
})