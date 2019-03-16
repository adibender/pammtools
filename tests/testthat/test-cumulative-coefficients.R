context("Test cumulative coefficients functionality")

test_that("Cumulative coefficients work", {

  data("tumor")

  # PAMs
  tumor <- tumor %>% dplyr::slice(1:100)
  tumor_ped <- tumor %>% as_ped(
    formula = Surv(days, status) ~ age + complications)
  pam <- mgcv::gam(
    formula = ped_status ~ s(tend, by = as.ordered(complications)) +
      s(tend, by = age),
      data = tumor_ped, family = poisson(), offset = offset)
  cumu_coef <- get_cumu_coef(pam, tumor_ped,
    terms = c("age", "complications"))
  expect_data_frame(cumu_coef, nrows = 82L, ncols = 6L)
  cumu_coef_pam <- get_cumu_coef(pam, tumor_ped,
    terms = c("(Intercept)", "age"))
  expect_data_frame(cumu_coef_pam, nrows = 82L, ncols = 6L)

  ## aalen model
  library(timereg)
  aalen <- aalen(Surv(days, status)~ age + complications, data = tumor)
  cumu_coef_aalen <- get_cumu_coef(aalen, terms = c("age", "complications"))
  expect_data_frame(cumu_coef_aalen, nrows = 86L, ncols = 6L)
  # cox aalen
  cox.aalen <- cox.aalen(Surv(days, status) ~ age + prop(complications),
    data = tumor)
  cumu_coef_cox.aalen <- get_cumu_coef(cox.aalen, terms = c("age"))
  expect_data_frame(cumu_coef_cox.aalen, nrows = 43, ncols = 6L)

})
