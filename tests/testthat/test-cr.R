context("Test if as_ped_cr works")

test_that("For single risk result is identical to as_ped.", {
  data("veteran", package = "survival")
  veteran <- veteran[c(1:3, 135:137), ]
  ped <- as_ped(
    data    = veteran,
    formula = Surv(time, status) ~ trt + age,
    cut     = c(0, 100, 400))
  ped_cr <- as_ped_cr(
    data    = veteran,
    formula = Surv(time, status) ~ trt + age,
    cut     = c(0, 100, 400))
  expect_equal(sum(ped_cr != ped), 0)
})

test_that("Works for competing risks.", {
  data("veteran", package = "survival")
  veteran <- veteran[c(1:3, 135:137), ]
  ped <- as_ped(
    data    = veteran,
    formula = Surv(time, status) ~ trt + age,
    cut     = c(0, 100, 400))
  veteran$status[1] <- 2
  ped_cr <- as_ped_cr(
    data    = veteran,
    formula = Surv(time, status) ~ trt + age,
    cut     = c(0, 100, 400))
  expect_equal(sum(ped_cr != ped), 1)
})

context
