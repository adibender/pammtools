context("create newdata")


test_that("creating newdata works on ungrouped data", {
	library(checkmate)
	library(magrittr)
	library(dplyr)
	set.seed(123)
	iris2 <- iris %>%  group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_data_frame(make_newdata(iris2), any.missing=FALSE, nrows=1, ncols=5)
	expect_equal(colnames(make_newdata(iris2)), colnames(iris2))
	expect_data_frame(make_newdata(iris2, Sepal.Length=5), any.missing=FALSE, 
		nrows=1, ncols=5)
	expect_equal(make_newdata(iris2, Sepal.Length=5)$Sepal.Length, 5)
	expect_data_frame(make_newdata(iris2, Sepal.Length=c(5, 6)), 
		any.missing=FALSE, nrows=2, ncols=5)
	expect_data_frame(make_newdata(iris2, expand="Sepal.Length", length.out=2), 
		any.missing=FALSE, nrows=2, ncols=5)
	expect_equal(make_newdata(iris2, expand="Sepal.Length", length.out=2)$Sepal.Length, 
		c(4.4, 7.1))
})


test_that("creating newdata fails on ungrouped data", {
	library(magrittr)
	library(dplyr)
	set.seed(123)
	iris2 <- iris %>% group_by(Species) %>% sample_n(2) %>% ungroup()

	expect_error(make_newdata(iris2, Sepal.length=5))
	expect_error(make_newdata(iris2, expand="Sepal.length"))
	expect_error(make_newdata(iris2, expand="Sepal.Length", length.out=-3))

})

