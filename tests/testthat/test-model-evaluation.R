context("Model evaluation helpers")

test_that("pec helpers work", {

  library(pec)
  data(tumor)
  ped <- tumor %>%
    as_ped(Surv(days, status) ~ complications, cut = seq(0, 500, by = 50))
  pam <- pamm(ped_status ~ complications, data = ped)
  suppressMessages({
    pec <- pec::pec(list(pam = pam), Surv(days, status) ~ 1, data = tumor,
      times = seq(.01, 500, by = 100), start = .01, exact = FALSE)
  })
  df_ibs <- as.data.frame(pec::crps(pec))
  expect_data_frame(df_ibs, nrow = 2, ncol = 3)
  expect_identical(colnames(df_ibs), c("method", "time", "IBS"))

})
