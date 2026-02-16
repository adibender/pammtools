context("mgcv convenience functions")

test_that("mgcv convenience works", {
	
	library(mgcv)
	g <- gam(Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data=iris)
	expect_data_frame(s1d <- tidy_smooth(g), nrows=200, ncols=7)
  expect_equal(s1d$ci_upper - s1d$fit, 1.96 * s1d$se)

  g_re <- gam(Sepal.Length ~ s(Species, bs = "re"), data = iris)
  current_theme <- ggplot2::theme_get()
  expect_is(gg_re(g_re), "ggplot")
  expect_true(isTRUE(all.equal(ggplot2::theme_get(), current_theme)))

  g2 <- gam(Sepal.Length ~ te(Sepal.Width, Petal.Length), data = iris)
  expect_data_frame(s2d <- suppressWarnings(tidy_smooth2d(g2, ci = TRUE)), ncols = 9L)
  expect_equal(s2d$ci_upper - s2d$fit, 1.96 * s2d$se)

})
