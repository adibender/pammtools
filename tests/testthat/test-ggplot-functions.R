context("ggplot functions")

data("prothr", package = "mstate")
prothr <- prothr[1:200, ] |>
  filter(Tstart != Tstop) |>
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  select(-trans)

ped_msm <- as_ped(
  data = prothr,
  formula = Surv(Tstart, Tstop, status) ~ .,
  transition = "transition",
  id = "id",
  timescale = "calendar"
)

pam_msm <- gam(
  ped_status ~ s(tend, by = transition, bs = "cr") + transition,
  data = ped_msm,
  family = poisson(),
  offset = offset
)

pam_msm_treat <- gam(
  ped_status ~ s(tend, by = transition, bs = "cr") + treat*transition,
  data = ped_msm,
  family = poisson(),
  offset = offset
)

ndf_msm <- make_newdata(ped_msm, 
  tend = unique(tend), 
  transition = unique(transition)) %>% 
  group_by(treat, transition) %>%
  add_trans_prob(pam_msm)

ndf_msm_treat <- make_newdata(ped_msm, 
  tend = unique(tend), 
  transition = unique(transition),
  treat = unique(treat)) %>% 
  group_by(treat, transition) %>%
  add_trans_prob(pam_msm_treat)

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

p_state_occupation = 
  gg_state_occupation(
    ndf_msm_treat, 
    group_var = "treat", 
    init_state = c(1,0,0)
  )

test_that("state occupation handles input correctly and returns ggplot object", {
  # 3-state model but init_state of length 4
  expect_error(
    gg_state_occupation(ndf_msm, init_state = c(1, 0, 0, 0))
  )
  p <- gg_state_occupation(ndf_msm, init_state = c(1, 0, 0))
  expect_s3_class(p, "ggplot")
  p_grouped <- gg_state_occupation(ndf_msm_treat, init_state = c(1, 0, 0), group_var = "treat")
  expect_s3_class(p_grouped, "ggplot")
})

test_that("state occupation plot has correct axis labels", {
  p <- gg_state_occupation(ndf_msm, init_state = c(1, 0, 0))
  expect_equal(p$labels$x, "Time")
  expect_equal(p$labels$y, "State occupation probability")
  expect_equal(p$labels$fill, "State")
})

test_that("state occupation applies facet wrap only when group_var is provided", {
  p <- gg_state_occupation(ndf_msm_treat, init_state = c(1, 0, 0), group_var = "treat")
  # Check that faceting is present
  expect_true(!is.null(p$facet))
  expect_false(inherits(p$facet, "FacetNull"))
  
  p <- gg_state_occupation(ndf_msm, init_state = c(1, 0, 0))
  expect_true(inherits(p$facet, "FacetNull"))
})

test_that("state occupation group labels from attributes are applied correctly", {
  p <- gg_state_occupation(ndf_msm_treat, init_state = c(1, 0, 0), group_var = "treat")
  # Factor levels should reflect the labels
  expect_true(is.factor(p$data[["treat"]]))
  expect_setequal(levels(p$data[["treat"]]), c("Placebo", "Prednisone"))
})

