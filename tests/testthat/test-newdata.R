context("create newdata")


test_that("creating newdata works on ungrouped data", {
	library(checkmate)
	library(magrittr)
	library(dplyr)
	iris %<>% group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_data_frame(make_newdata(iris), any.missing=FALSE, nrows=1, ncols=5)
	expect_equal(colnames(make_newdata(iris)), colnames(iris))
	expect_data_frame(make_newdata(iris, Sepal.Length=5), any.missing=FALSE, 
		nrows=1, ncols=5)
	expect_equal(make_newdata(iris, Sepal.Length=5)$Sepal.Length, 5)
	expect_data_frame(make_newdata(iris, Sepal.Length=c(5, 6)), 
		any.missing=FALSE, nrows=2, ncols=5)
	expect_data_frame(make_newdata(iris, expand="Sepal.Length", length.out=2), 
		any.missing=FALSE, nrows=2, ncols=5)
	expect_equal(make_newdata(iris, expand="Sepal.Length", length.out=2)$Sepal.Length, 
		c(5.1, 6.4))
})


test_that("creating newdata fails on ungrouped data", {
	library(magrittr)
	library(dplyr)
	iris %<>% group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_error(make_newdata(iris, Sepal.length=5))
	expect_error(make_newdata(iris, expand="Sepal.length"))
	expect_error(make_newdata(iris, expand="Sepal.Length", length.out=-3))

})

