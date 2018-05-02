context("Test cumulative coefficients functionality")

test_that("Cumulative coefficients work for PAMMs", {
  data("tumor")
  set.seed(123)
  tumor<- tumor %>% dplyr::sample_n(100)
  tumor_ped <- tumor %>% as_ped(
    formula = Surv(days, status)~ age + complications)
  pamm <- mgcv::gam(
    formula = ped_status ~ s(tend, by = as.ordered(complications)) +
      s(tend, by =age),
    data   = tumor_ped,
    family = poisson(),
    offset = offset)
  # summary(pamm)
  # plot(pamm, page=1)
  cumu_coef_pam <- get_cumu_coef(pamm, tumor_ped, terms=c("age", "complications"))
  expect_data_frame(cumu_coef_pam, nrows=100, ncols=6)
  # aalen <- timereg::aalen(survival::Surv(survhosp, tumorDied)~Gender + ApacheIIScore,
  #   data=tumor)
  # cumu_coef_aalen <- get_cumu_coef(aalen, terms=c("Gender", "ApacheIIScore"))
  # expect_data_frame(cumu_coef_aalen, nrows=38, ncols=6)
  # bind_rows(cumu_coef_pam, cumu_coef_aalen) %>%
  #   ggplot(aes(x=time, y = cumu_hazard)) +
  #   geom_line(aes(col=method)) +
  #   geom_ribbon(aes(fill=method, ymin = cumu_lower, ymax=cumu_upper), alpha = 0.3) +
  #   facet_wrap(~variable)

})
