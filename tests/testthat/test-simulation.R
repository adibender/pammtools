context("Test simulation functions")

test_that("Simple simulation function works", {

  set.seed(24032018)
  df     <- cbind.data.frame(x1 = runif(3, -3, 3), x2=runif(3, 0, 6))
  sim_df <- sim_pexp(~ -3.5 -0.5*x1 + sqrt(x2), df, cut=0:10)
  expect_data_frame(sim_df, nrows = 3, ncols = 5)
  expect_identical(round(sim_df$time, 2), c(1.38, 7.14, 3.02))

})
