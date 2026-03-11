context("Test cumulative coefficients functionality")

test_that("Cumulative coefficients work", {
  df <- tumor[1:30, c("days", "status", "age")]
  df$x1 <- as.factor(rep(letters[1:3], each = nrow(df) / 3L))

  ## pam
  ped <- as_ped(df, formula = Surv(days, status) ~ x1 + age)
  pam <- mgcv::gam(
    ped_status ~ s(tend) + x1 + age,
    data = ped,
    family = poisson(),
    offset = offset
  )
  cumu_coef_pam <- get_cumu_coef(pam, ped, terms = c("age", "x1"), nsim = 20L)
  expect_data_frame(cumu_coef_pam, nrows = 36L, ncols = 6L)
  expect_equal(unique(cumu_coef_pam$variable), c("age", "x1 (b)", "x1 (c)"))
  cumu_coef_pam <- get_cumu_coef(pam, ped, terms = c("(Intercept)", "age"))
  expect_data_frame(cumu_coef_pam, nrows = 24L, ncols = 6L)
})

test_that("Intercept cumulative coefficient honors nonstandard time_var", {
  df <- tumor[1:30, c("days", "status", "age")]
  df$x1 <- as.factor(rep(letters[1:3], each = nrow(df) / 3L))

  ped <- as_ped(df, formula = Surv(days, status) ~ x1 + age) |>
    rename(stop = tend)
  intvars <- attr(ped, "intvars")
  attr(ped, "intvars") <- replace(intvars, intvars == "tend", "stop")
  if (identical(attr(ped, "time_var"), "tend")) attr(ped, "time_var") <- "stop"
  pam <- mgcv::gam(
    ped_status ~ s(stop) + x1 + age,
    data = ped,
    family = poisson(),
    offset = offset
  )

  cumu_coef_pam <- get_cumu_coef(
    pam,
    ped,
    terms = "(Intercept)",
    time_var = "stop"
  )

  expect_true("time" %in% names(cumu_coef_pam))
  expect_false("stop" %in% names(cumu_coef_pam))
  expect_equal(unique(cumu_coef_pam$variable), "(Intercept)")
})
