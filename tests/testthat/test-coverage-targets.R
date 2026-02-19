context("Coverage increase targets")

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

test_that("seq_range validates n/by and rounds boundaries for pretty by sequences", {
  expect_error(seq_range(1:10, n = 5, by = 1), "May only specify one")

  seq_pretty <- seq_range(c(0.3, 2.2), by = 1, pretty = TRUE)
  expect_equal(seq_pretty, 0:3)
})

test_that("get_cut.list combines cuts in gap-timescale recurrent setting", {
  d1 <- data.frame(
    id = c(1, 1, 1),
    tstart = c(0, 1, 2),
    tstop = c(1, 2, 3),
    status = c(1, 0, 1),
    enum = c(1, 2, 3),
    age = c(50, 50, 50)
  )
  d2 <- data.frame(
    id = c(2, 2),
    tstart = c(0, 1.5),
    tstop = c(1.5, 2.5),
    status = c(0, 1),
    enum = c(1, 2),
    age = c(40, 40)
  )

  cuts <- pammtools:::get_cut.list(
    data = list(d1, d2),
    formula = survival::Surv(tstart, tstop, status) ~ enum + age,
    timescale = "gap"
  )

  expect_equal(cuts, c(1, 2.5, 3))
})

test_that("get_terms returns tidy partial effects for requested terms", {
  fit <- survival::coxph(
    survival::Surv(time, status) ~
      survival::pspline(karno) + survival::pspline(age),
    data = survival::veteran
  )

  terms_df <- get_terms(
    data = survival::veteran,
    fit = fit,
    terms = c("karno", "age")
  )

  expect_true(is.data.frame(terms_df))
  expect_identical(
    names(terms_df),
    c("term", "x", "eff", "se", "ci_lower", "ci_upper")
  )
  expect_equal(sort(unique(terms_df$term)), c("age", "karno"))
  expect_equal(as.integer(table(terms_df$term)[["age"]]), 100L)
  expect_equal(as.integer(table(terms_df$term)[["karno"]]), 100L)
  expect_true(all(terms_df$ci_lower <= terms_df$eff))
  expect_true(all(terms_df$ci_upper >= terms_df$eff))
})
