context("mgcv convenience functions")

test_that("mgcv convenience works", {
	
	library(mgcv)
	g <- gam(Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data=iris)
	expect_data_frame(s1d <- tidy_smooth(g), nrows=200, ncols=7)

})