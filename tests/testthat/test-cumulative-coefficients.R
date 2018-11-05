context("Test cumulative coefficients functionality")

test_that("Cumulative coefficients work", {
  data("tumor")
  set.seed(123)
  tumor <- tumor %>% dplyr::sample_n(100)
  tumor_ped <- tumor %>% as_ped(
    formula = Surv(days, status) ~ age + complications)
  pam <- mgcv::gam(
    formula = ped_status ~ s(tend, by = as.ordered(complications)) +
      s(tend, by = age),
      data = tumor_ped, family = poisson(), offset = offset)
  expect_warning(bam <- mgcv::bam(
    formula = ped_status ~ s(tend, by = as.ordered(complications)) +
      s(tend, by = age),
      data   = tumor_ped, family = poisson(), offset = offset, method = "REML"))
  cumu_coef_bam <- get_cumu_coef(bam, tumor_ped,
    terms = c("age", "complications"))
  expect_data_frame(cumu_coef_bam, nrows = 100, ncols = 6)
  cumu_coef_pam <- get_cumu_coef(pam, tumor_ped,
    terms = c("(Intercept)", "age"))
  expect_data_frame(cumu_coef_pam, nrows = 100, ncols = 6)
  ## aalen model
  library(timereg)
  aalen <- aalen(Surv(days, status)~ age + complications, data = tumor)
  cumu_coef_aalen <- get_cumu_coef(aalen, terms = c("age", "complications"))
  expect_data_frame(cumu_coef_aalen, nrows = 104L, ncols = 6L)
  # cox aalen
  cox.aalen <- cox.aalen(Surv(days, status) ~ age + prop(complications),
    data = tumor)
  cumu_coef_cox.aalen <- get_cumu_coef(cox.aalen, terms = c("age"))
  expect_data_frame(cumu_coef_cox.aalen, nrows = 52, ncols = 6L)

})
