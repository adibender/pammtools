context("Test formula special.")

test_that("Formula special func works as expected", {
  ## time + latency + covar (DLNM approach)
  f1 <- func(t, latency(te), x, te_var="te")
  expect_list(f1, any.missing = FALSE, len = 5)
  expect_identical(f1$latency_var, "te")
  expect_identical(f1$te_var, "te")
  expect_identical(f1$col_vars, c("t", "te", "x"))
  expect_function(f1$ll_fun, args=c("t", "te"))
  expect_identical(f1$suffix, NULL)

})


context("Transformation of longitudinal covariates to functional covariates")


test_that("Covariate to matrix Transformation works", {
  event_df  <- filter(patient, CombinedID == 1116)
  tdc_df    <- filter(daily, CombinedID == 1116)
  ## check nesting
  nested_df <- nest_tdc(
    data    = list(event_df, tdc_df),
    formula = Surv(survhosp, status)~.|func(Study_Day, caloriesPercentage) +
      func(proteinGproKG),
    cut     = 0:30,
    id  = "CombinedID")
  expect_tibble(nested_df, any.missing=FALSE, nrows=1, ncols=15)
  expect_identical(colnames(nested_df), c("Year", "CombinedicuID", "CombinedID", "Survdays",
      "PatientDied", "survhosp", "Gender", "Age", "AdmCatID", "ApacheIIScore",
      "BMI", "DiagID2", "Study_Day", "caloriesPercentage", "proteinGproKG"))
  expect_identical(names(attributes(nested_df))[-c(1:3)],
    c("id_var", "time_var", "status_var", "tdc_vars",
      "breaks", "id_n", "id_tseq", "id_teseq"))
  ## check data trafo
  expect_error(get_func(nested_df, ~func(t)))
  f1 <- get_func(nested_df, ~ func(survhosp,latency(Study_Day), caloriesPercentage,
    te_var = "Study_Day"))
  expect_list(f1, types=c("numeric", "numeric", "numeric", "integer"),
    any.missing=FALSE, len=4, names="named")
  f2 <- get_func(nested_df,
      ~func(survhosp,latency(Study_Day), caloriesPercentage, te_var = "Study_Day") +
      func(proteinGproKG, te_var = "Study_Day"))
  expect_list(f2, types=c(rep("numeric", 3), "integer", "numeric"),
    any.missing = FALSE, len=5, names="named")
})
