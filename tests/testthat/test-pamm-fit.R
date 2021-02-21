context("Test pamm wrapper function")


 test_that("pamm function works correctly", {

  data("veteran", package="survival")
  ped <- as_ped(Surv(time, status)~ trt + karno, data = veteran[1:20,])
  # gam engine
  pam <- pamm(ped_status ~ s(tend, k=3) + karno, data=ped)
  expect_is(pam, "pamm")
  expect_is(summary(pam), "summary.gam")
  expect_data_frame(int_info(pam), nrows = 18L, ncols = 5L)
  expect_identical(is.pamm(pam), TRUE)

  # check data trafo from pam object
  ped_new <- as_ped(pam, newdata = veteran[21:40, ])
  expect_data_frame(ped_new, nrows = 192L, ncols = 8L)
  expect_subset(ped_new$tend, ped$tend)

  # bam engine
  pam2 <- pamm(ped_status ~ s(tend, k = 3) + karno, data = ped, engine = "bam")
  expect_true(inherits(pam2, "bam"))
  expect_data_frame(int_info(pam2), nrows = 18L, ncols = 5L)
  expect_identical(is.pamm(pam2), TRUE)
  # pass arguments to bam
  pam3 <- pamm(ped_status ~ s(tend, k = 3) + karno, data = ped,
    engine = "bam", discrete = TRUE, method = "fREML")
  expect_true(inherits(pam3, "bam"))
  expect_data_frame(int_info(pam3), nrows = 18L, ncols = 5L)
  expect_identical(is.pamm(pam), TRUE)

 })
