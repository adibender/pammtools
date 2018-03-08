context("Test formula special.")

test_that("Formula special func works as expected", {
  ## latency + covar (DLNM approach)
  f1 <- func(t-te, x)
  expect_list(f1, any.missing = FALSE, len = 7)
  expect_identical(f1$lgl_latency, c(TRUE, FALSE))
  expect_identical(f1$lgl_te, c(FALSE, FALSE))
  expect_identical(f1$lgl_time, c(FALSE, FALSE))
  expect_identical(f1$fun_covar, "x")
  expect_null(f1$by_var)

  ## latency + linear covar effect (WCE approach)
  f2 <- func(t-te, by=x)
  expect_list(f2, any.missing = FALSE, len = 7)
  expect_identical(f2$lgl_latency, TRUE)
  expect_identical(f2$lgl_te, FALSE)
  expect_identical(f2$lgl_time, FALSE)
  expect_null(f2$fun_covar)
  expect_identical(f2$by_var, "x")

  ## time-varying DLNM (with by term)
  f3 <- func(t, t-te, x, by = z)
  expect_list(f3, any.missing = FALSE, len = 7)
  expect_identical(f3$lgl_latency, c(FALSE, TRUE, FALSE))
  expect_identical(f3$lgl_te, c(FALSE, FALSE, FALSE))
  expect_identical(f3$lgl_time, c(TRUE, FALSE, FALSE))
  expect_identical(f3$fun_covar,"x")
  expect_identical(f3$by_var, "z")

  # General form
  f4 <- func(t, te, x, by = z)
  expect_list(f4, any.missing = FALSE, len = 7)
  expect_identical(f4$lgl_latency, c(FALSE, FALSE, FALSE))
  expect_identical(f4$lgl_te, c(FALSE, TRUE, FALSE))
  expect_identical(f4$lgl_time, c(TRUE, FALSE, FALSE))
  expect_identical(f4$fun_covar,"x")
  expect_identical(f4$by_var, "z")

})


context("Transformation of longitudinal covariates to functional covariates")


test_that("Covariate to matrix Transformation works", {
  event_df  <- filter(patient, CombinedID == 1116)
  tdc_df    <- filter(daily, CombinedID == 1116)
  ## check nesting
  nested_df <- nest_tdc(event_df, tdc_df, "Study_Day", "CombinedID",
    "survhosp", "PatientDied", 0:10, TRUE, 0, 0, "TDC")
  expect_tibble(nested_df, any.missing=FALSE, nrows=1, ncols=13)
  expect_identical(colnames(nested_df), c("Year", "CombinedicuID", "CombinedID", "Survdays",
      "PatientDied", "survhosp", "Gender", "Age", "AdmCatID", "ApacheIIScore",
      "BMI", "DiagID2", "TDC"))
  expect_identical(names(attributes(nested_df))[-c(1:3)],
    c("id_var", "time_var", "status_var", "te_var", "tdc_col", "cens_value",
      "breaks", "te", "id_n", "id_tseq", "id_teseq"))
  ## check data trafo
  expect_error(get_func(nested_df, ~func(t)))
  f1 <- get_func(nested_df, ~ func(t, t-te, caloriesPercentage))
  expect_list(f1, types=c("numeric", "numeric", "integer", "numeric"),
    any.missing=FALSE, len=4, names="named")
})
