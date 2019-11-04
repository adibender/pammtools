context("Test formula special.")

test_that("Formula special 'func' works as expected", {
  ## time + latency + covar (DLNM approach)
  cumu1 <- eval_special(~ cumulative(t, latency(te), x, tz_var = "te"))[[1]]
  expect_list(cumu1, any.missing = TRUE, len = 5)
  expect_identical(cumu1$latency_var, "te")
  expect_identical(cumu1$tz_var, "te")
  expect_identical(cumu1$col_vars, c("t", "te", "x"))
  expect_function(cumu1$ll_fun, args = c("t", "tz"))
  expect_identical(cumu1$suffix, NULL)

})

test_that("Formula special 'concurrent' works as expected", {
  ## time + latency + covar (DLNM approach)
  ccr1 <- eval_special(~ concurrent(x1, x2, tz_var = "te"),
    special = "concurrent")[[1]]
  expect_list(ccr1, any.missing = TRUE, len = 5)
  expect_identical(ccr1$tz_var, "te")
  expect_identical(ccr1$col_vars, c("x1", "x2"))
  expect_function(ccr1$ll_fun, args = c("t"))
  expect_identical(ccr1$lag, 0)
  expect_identical(ccr1$suffix, NULL)

  data("pbc", package = "survival")
  event_df <- pbc %>%
    filter(id <= 5) %>%
    mutate(event = 1L*(status == 2)) %>%
    select(id, time, event, sex, bili, protime, albumin)
  tdc_df <- pbcseq %>%
    filter(id <= 5) %>%
    select(id, day, bili, protime, albumin)
  formula <- Surv(time, event)~ concurrent(bili, protime, tz_var = "day") +
    concurrent(albumin, tz_var  = "day")
  nested_fdf <- nest_tdc(list(event_df, tdc_df), formula, id = "id")
  ped_ccr <- as_ped(list(event_df, tdc_df), formula, id = "id")

})


context("Transformation of longitudinal covariates to functional covariates")


test_that("Covariate to matrix Transformation works", {
  event_df  <- filter(patient, CombinedID == 1116)
  tdc_df    <- filter(daily, CombinedID == 1116)
  ## check nesting
  nested_df <- nest_tdc(
    data    = list(event_df, tdc_df),
    formula = Surv(survhosp, status)~ . +
      cumulative(Study_Day, caloriesPercentage, tz_var="Study_Day") +
        cumulative(proteinGproKG, tz_var="Study_Day"),
    cut     = 0:30,
    id  = "CombinedID")
  expect_tibble(nested_df, any.missing=FALSE, nrows=1, ncols=15)
  expect_identical(colnames(nested_df), c("Year", "CombinedicuID", "CombinedID", "Survdays",
      "PatientDied", "survhosp", "Gender", "Age", "AdmCatID", "ApacheIIScore",
      "BMI", "DiagID2", "Study_Day", "caloriesPercentage", "proteinGproKG"))
  expect_identical(names(attributes(nested_df))[-c(1:3)],
    c("id_var", "time_var", "status_var", "tdc_vars",
      "breaks", "func_list", "id_n", "id_tseq", "id_tz_seq"))
  ## check data trafo
  expect_error(get_cumulative(nested_df, ~cumulative(t)))
  f1 <- get_cumulative(nested_df, ~ . +
      cumulative(survhosp, latency(Study_Day), caloriesPercentage, tz_var = "Study_Day"))
  expect_list(f1$func_mats, types=c("numeric", "numeric", "numeric", "integer"),
    any.missing=FALSE, len=4, names="named")
  f2 <- get_cumulative(nested_df,
      ~. + cumulative(survhosp, latency(Study_Day), caloriesPercentage, tz_var = "Study_Day") +
      cumulative(proteinGproKG, tz_var = "Study_Day"))
  expect_list(f2$func_mats, types=c(rep("numeric", 3), "integer", "numeric"),
    any.missing = FALSE, len=5, names="named")
})
