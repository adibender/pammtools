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

test_that("gg_smooth/get_terms handle TVE PAMMs, by-/fs-/sz-smooths and factors", {

  suppressPackageStartupMessages({
    library(mgcv)
    library(dplyr)
    library(survival)
  })
  data("tumor", package = "pammtools")
  set.seed(1)
  td <- tumor[1:350, ] %>% mutate(x2 = rnorm(n()))
  ped <- as_ped(
    Surv(days, status) ~ age + sex + charlson_score + complications + metastases + x2,
    data = td, cut = seq(0, 2500, by = 500)
  )
  pam <- pamm(
    ped_status ~ s(tend, k = 5) +
      s(tend, by = metastases, k = 5) +              # factor by-variable
      s(tend, by = x2, k = 5) +                      # numeric by-variable
      ti(age, k = 5) +                               # 1d ti -> included
      s(age, sex, bs = "fs", k = 5) +                # factor-smooth interaction
      s(charlson_score, complications, bs = "sz", k = 5) +
      te(age, tend) +                                # tensor -> excluded
      complications,                                 # parametric factor -> skipped
    data = ped
  )

  # (1) bare variable name "tend": multiple curves, no "size 800" error, factor-by
  #     collapsed into one facet with >= 2 levels, exactly 100 rows per (term, level)
  tt <- get_terms(ped, pam, terms = "tend")
  expect_s3_class(tt, "tbl_df")
  expect_identical(
    names(tt),
    c("term", "x", "level", "eff", "se", "ci_lower", "ci_upper")
  )
  counts <- tt %>% count(.data$term, .data$level)
  expect_true(all(counts$n == 100))
  meta_levels <- tt %>%
    filter(.data$term == "s(tend):metastases") %>%
    pull(.data$level) %>%
    unique()
  expect_setequal(meta_levels, c("no", "yes"))

  # (2) parametric factors warn-and-skip instead of erroring on range()
  expect_warning(
    get_terms(ped, pam, terms = c("complications", "metastases")),
    "not a univariate smooth"
  )

  # (3) terms omitted -> all 1d smooths; tensor/te excluded
  allt <- get_terms(ped, pam)
  facets <- unique(allt$term)
  expect_true(all(c("s(tend)", "ti(age)", "s(age,sex)") %in% facets))
  expect_false(any(grepl("te(age,tend)", facets, fixed = TRUE)))

  # (4) 1d ti picked up, te/2d not
  expect_true("ti(age)" %in% facets)

  # (5) bs="fs"/"sz": one facet, one level per category
  fs_levels <- allt %>% filter(.data$term == "s(age,sex)") %>%
    pull(.data$level) %>% unique()
  expect_setequal(fs_levels, c("male", "female"))
  sz_levels <- allt %>% filter(.data$term == "s(charlson_score,complications)") %>%
    pull(.data$level) %>% unique()
  expect_setequal(sz_levels, c("no", "yes"))

  # (6) re/mrf excluded
  ped_re <- ped %>% mutate(agf = as.factor(round(age / 10)))
  pam_re <- pamm(ped_status ~ s(tend, k = 5) + s(agf, bs = "re"), data = ped_re)
  expect_setequal(unique(get_terms(ped_re, pam_re)$term), "s(tend)")

  # (7) exact label selects a single smooth
  one <- get_terms(ped, pam, terms = "s(tend)")
  expect_equal(unique(one$term), "s(tend)")
  expect_equal(nrow(one), 100L)

  # (8) gg_smooth returns ggplot objects
  expect_s3_class(gg_smooth(ped, pam, terms = "tend"), "ggplot")
  expect_s3_class(gg_smooth(ped, pam), "ggplot")

  # (9) CI relationship
  z95 <- stats::qnorm(0.975)
  expect_equal(tt$ci_upper - tt$eff, z95 * tt$se)
  expect_true(all(tt$ci_lower <= tt$eff))

})
