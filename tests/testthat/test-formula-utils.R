context("Formula utility functions")

test_that("Formula utilities work", {

  expect_identical(get_rhs_vars(~ x1 + sqrt(x2)), c("x1", "x2"))
  expect_identical(get_rhs_vars("~ x1 + sqrt(x2)"), c("x1", "x2"))
  expect_identical(get_tdc_vars( ~ x1 + cumulative(z.tz, tz_var = "tz")), "z.tz")
  expect_identical(get_ped_form(Surv(time, status) ~ x1 +
      cumulative(z.tz, tz_var = "tz")), Surv(time, status) ~ x1 )
  expect_true(has_lhs(Surv(time, status) ~ .))
  expect_identical(get_lhs_vars("Surv(time, status) ~ ."), c("time", "status"))
  expect_false(has_lhs( ~ .))

})
