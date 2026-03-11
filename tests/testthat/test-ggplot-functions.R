context("ggplot functions")

test_that("ggplot2 stairstep helper handles all directions", {
  stair_df <- data.frame(
    x = c(1, 3, 5),
    ymin = c(2, 4, 6),
    ymax = c(3, 5, 7)
  )

  stair_hv <- pammtools:::ggplot2_stairstep(stair_df, direction = "hv")
  stair_vh <- pammtools:::ggplot2_stairstep(stair_df, direction = "vh")
  stair_mid <- pammtools:::ggplot2_stairstep(stair_df, direction = "mid")
  stair_single <- pammtools:::ggplot2_stairstep(
    stair_df[1, , drop = FALSE],
    direction = "hv"
  )

  expect_equal(stair_hv$x, c(1, 3, 3, 5, 5))
  expect_equal(stair_vh$x, c(1, 1, 3, 3, 5))
  expect_equal(stair_mid$x, c(1, 2, 2, 4, 4, 5))
  expect_identical(nrow(stair_single), 0L)
})

test_that("geom_stepribbon layer can be built", {
  ribbon_df <- data.frame(
    x = 1:4,
    ymin = c(0.1, 0.2, 0.3, 0.4),
    ymax = c(0.3, 0.4, 0.5, 0.6)
  )

  plt <- ggplot2::ggplot(
    ribbon_df,
    ggplot2::aes(x = x, ymin = ymin, ymax = ymax)
  ) +
    geom_stepribbon(direction = "hv")
  built <- ggplot2::ggplot_build(plt)
  grob <- ggplot2::ggplotGrob(plt)

  expect_identical(length(built$data), 1L)
  expect_true(all(c("x", "ymin", "ymax") %in% names(built$data[[1]])))
  expect_is(grob, "gtable")
})

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
