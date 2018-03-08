context("Transformation of longitudinal covariates to functional covariates")


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
