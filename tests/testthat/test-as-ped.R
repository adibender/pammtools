context("Test as_ped functions")

test_that("Attributes and class are appended when data already PED format", {
  # preparations
  data("veteran", package="survival")
  veteran <- veteran[c(1:3, 135:137), ]
  ped <- split_data(Surv(time, status)~ trt + age, data=veteran, cut=c(0, 100, 400))
  attributes(ped)[c("cut", "id_var", "intvars")] <- NULL
  class(ped) <- class(ped)[-1]
  # retransform to ped
  ped2 <- as_ped(ped)
  expect_data_frame(ped2, nrow = 10L, ncols = 8L)
  expect_is(ped2, "ped")
  expect_subset(c("ped_status", "tstart", "tend", "interval", "offset"), names(ped2))
  expect_is(attr(ped2, "cut"), "numeric")
  expect_is(attr(ped2, "intvars"), "character")
  expect_is(attr(ped2, "id_var"), "character")
  expect_equal(attr(ped2, "id_var"), "id")

})
