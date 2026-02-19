context("Plot and tidier consistency")

test_that("gg_tensor uses tidy_smooth2d output grid consistently", {
  fit <- mgcv::gam(Sepal.Length ~ te(Sepal.Width, Petal.Length), data = iris)

  tidy_grid <- suppressWarnings(tidy_smooth2d(fit, ci = FALSE))
  tensor_plot <- gg_tensor(fit, ci = FALSE)

  expect_equal(nrow(tidy_grid), nrow(tensor_plot$data))
  expect_equal(tidy_grid$fit, tensor_plot$data$fit, tolerance = 1e-12)
})

test_that("gg_re does not mutate global ggplot theme", {
  fit <- mgcv::gam(Sepal.Length ~ s(Species, bs = "re"), data = iris)
  old_theme <- ggplot2::theme_get()

  gg_re(fit)

  expect_equal(ggplot2::theme_get(), old_theme)
})

test_that("tidy smooth confidence intervals use standard normal scaling", {
  z95 <- stats::qnorm(0.975)

  fit_1d <- mgcv::gam(Sepal.Length ~ s(Sepal.Width), data = iris)
  smooth_1d <- tidy_smooth(fit_1d, ci = TRUE)
  expect_equal(
    smooth_1d$ci_upper - smooth_1d$fit,
    z95 * smooth_1d$se,
    tolerance = 1e-10
  )

  fit_2d <- mgcv::gam(Sepal.Length ~ te(Sepal.Width, Petal.Length), data = iris)
  smooth_2d <- suppressWarnings(tidy_smooth2d(fit_2d, ci = TRUE))
  expect_equal(
    smooth_2d$ci_upper - smooth_2d$fit,
    z95 * smooth_2d$se,
    tolerance = 1e-10
  )
})
