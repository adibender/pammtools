context("Test cumulative coefficients functionality")

test_that("Cumulative coefficients work", {

  df <- tumor[1:30, c("days", "status", "age")]
  df$x1 <- as.factor(rep(letters[1:3], each = nrow(df) / 3L))

  ## aalen model
  library(timereg)
  mod       <- aalen(Surv(days, status) ~ x1 + age, data = df)
  cumu_coef_aalen <- get_cumu_coef(
    mod,
    df,
    terms = c("(Intercept)", "x1"))
  expect_data_frame(cumu_coef_aalen, nrows = 42L, ncols = 6L)
  expect_equal(unique(cumu_coef_aalen$variable), c("(Intercept)", "x1b", "x1c"))

  ## pam
  ped <- as_ped(df, formula = Surv(days, status)~ x1 + age)
  pam <- mgcv::gam(ped_status ~ s(tend) + x1 + age, data = ped,
    family = poisson(), offset = offset)
  cumu_coef_pam <- get_cumu_coef(pam, ped, terms = c("age", "x1"), nsim = 20L)
  expect_data_frame(cumu_coef_pam, nrows = 36L, ncols = 6L)
  expect_equal(unique(cumu_coef_pam$variable), c("age", "x1 (b)", "x1 (c)"))
  cumu_coef_pam <- get_cumu_coef(pam, ped, terms = c("(Intercept)", "age"))
  expect_data_frame(cumu_coef_pam, nrows = 24L, ncols = 6L)
  # ggplot(
  #   filter(cc2, variable == "(Intercept)"),
  #   aes(x = time, y = cumu_hazard)) +
  # geom_step() +
  # geom_step(data = filter(cumu_coef_aalen, variable == "(Intercept)"), col = 2)


})
