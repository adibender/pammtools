context("Test cumulative coefficients functionality")

test_that("Cumulative coefficients work for PAMMs", {
  data("patient")
  set.seed(123)
  patient <- patient %>% dplyr::sample_n(100)
  patient_ped <- patient %>% as_ped(
    formula = Surv(survhosp, PatientDied)~ Gender + ApacheIIScore,
    id = "CombinedID")
  pamm <- mgcv::gam(
    formula = ped_status ~ Gender + s(tend, by = as.ordered(Gender)) +
      s(tend, by = ApacheIIScore),
    data   = patient_ped,
    family = poisson(),
    offset = offset)
  # summary(pamm)
  # plot(pamm, page=1)
  cumu_coef_pam <- get_cumu_coef(pamm, patient_ped, terms=c("Gender", "ApacheIIScore"))
  expect_data_frame(cumu_coef_pam, nrows=36, ncols=6)
  # aalen <- timereg::aalen(survival::Surv(survhosp, PatientDied)~Gender + ApacheIIScore,
  #   data=patient)
  # cumu_coef_aalen <- get_cumu_coef(aalen, terms=c("Gender", "ApacheIIScore"))
  # expect_data_frame(cumu_coef_aalen, nrows=38, ncols=6)
  # bind_rows(cumu_coef_pam, cumu_coef_aalen) %>%
  #   ggplot(aes(x=time, y = cumu_hazard)) +
  #   geom_line(aes(col=method)) +
  #   geom_ribbon(aes(fill=method, ymin = cumu_lower, ymax=cumu_upper), alpha = 0.3) +
  #   facet_wrap(~variable)

})
