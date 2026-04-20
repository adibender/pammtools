context("mgcv convenience functions")

test_that("mgcv convenience works", {

  library(mgcv)
  g <- gam(Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data = iris)
  z95 <- stats::qnorm(0.975)
  expect_data_frame(s1d <- tidy_smooth(g), nrows = 200, ncols = 7)
  expect_equal(s1d$ci_upper - s1d$fit, z95 * s1d$se)

  s1d_90 <- tidy_smooth(g, conf_level = 0.90)
  z90 <- stats::qnorm(0.95)
  expect_equal(s1d_90$ci_upper - s1d_90$fit, z90 * s1d_90$se)

  g_re <- gam(Sepal.Length ~ s(Species, bs = "re"), data = iris)
  current_theme <- ggplot2::theme_get()
  expect_is(gg_re(g_re), "ggplot")
  expect_true(isTRUE(all.equal(ggplot2::theme_get(), current_theme)))

  g2 <- gam(Sepal.Length ~ te(Sepal.Width, Petal.Length), data = iris)
  expect_data_frame(
    s2d <- suppressWarnings(tidy_smooth2d(g2, ci = TRUE)),
    ncols = 9L
  )
  expect_equal(s2d$ci_upper - s2d$fit, z95 * s2d$se)

  s2d_80 <- suppressWarnings(tidy_smooth2d(g2, ci = TRUE, conf_level = 0.80))
  z80 <- stats::qnorm(0.90)
  expect_equal(s2d_80$ci_upper - s2d_80$fit, z80 * s2d_80$se)

  expect_error(tidy_smooth(g, conf_level = 1.1))
  expect_error(tidy_smooth2d(g2, ci = TRUE, conf_level = 0))

})
