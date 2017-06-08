context("Test pamm wrapper function")


test_that("pamm function works correctly", {
	
	data("veteran", package="survival")
	ped <- split_data(Surv(time, status)~ trt + karno, veteran[1:20,])
	pam <- pamm(ped_status ~ s(tend, k=3) + trt + karno, data=ped)
	expect_is(pam, "pamm")
	expect_is(summary(pam), "summary.gam")
	expect_data_frame(int_info(pam), nrows=18L, ncols=5L)
	expect_identical(is.pamm(pam), TRUE)

})