context("Convenience functiions for calculation of hazard and similar")

data("leuk2", package="bpcp")
leuk2$x <- rnorm(nrow(leuk2))
leuk.ped <- split_data(Surv(time, status)~., data=leuk2, cut=c(0:5, 10, 40), id="id")
fit <- mgcv::gam(status ~ s(tend, k=5) + treatment, data=leuk.ped)

test_that("hazard functions work", {
	expect_is(haz <- add_hazard(ped_info(leuk.ped), fit), "data.frame")
	expect_identical(dim(haz), c(7L, 12L))
	expect_error(add_hazard(haz, fit))
	expect_identical(dim(add_hazard(haz, fit, overwrite = TRUE)), c(7L, 12L))
})


test_that("cumulative hazard functions work", {
	expect_is(haz <- add_cumhazard(ped_info(leuk.ped), fit), "data.frame")
	expect_identical(dim(haz), c(7L, 11L))
	expect_error(add_cumhazard(haz, fit))
	expect_identical(dim(add_cumhazard(haz, fit, overwrite = TRUE)), c(7L, 11L))
})