context("Test pamm_cr function")

test_that("pamm_cr works for joint PED object.", {

  ped <- as_ped(
    data         = sir_adm2,
    formula      = Surv(time, status) ~ age + pneu
  )
  # gam engine
  pam <- pamm_cr(ped_status ~ s(tend, by = cause, k = 3L) + age : cause,
                 data = ped)
  expect_is(pam, "pamm")
  expect_is(summary(pam), "summary.gam")
  # expect_data_frame(int_info(pam), nrows = 43L, ncols = 5L)
  expect_identical(is.pamm(pam), TRUE)

  # bam engine
  pam2 <- pamm_cr(ped_status ~ s(tend, by = cause, k = 3L) + age : cause,
                  data = ped, engine = "bam")
  expect_true(inherits(pam2, "bam"))
  # expect_data_frame(int_info(pam2), nrows = 43L, ncols = 5L)
  expect_identical(is.pamm(pam2), TRUE)
  # pass arguments to bam
  pam3 <- pamm_cr(ped_status ~ s(tend, by = cause, k = 3L) + age : cause,
                  data = ped, engine = "bam", discrete = TRUE, method = "fREML")
  expect_true(inherits(pam3, "bam"))
  # expect_data_frame(int_info(pam3), nrows = 43L, ncols = 5L)
  expect_identical(is.pamm(pam), TRUE)

})

test_that("pamm_cr works for list PED object.", {
  # preparations
  ped <- as_ped_cr(
    data         = sir_adm2,
    formula      = Surv(time, status) ~ age + pneu,
    combine      = FALSE
  )
  # gam engine
  pam <- pamm_cr(ped_status ~ s(tend, k = 3L) + age, data = ped)
  expect_is(pam, "pamm_cr_list")
  expect_is(pam, "pamm_cr")
  expect_is(pam[[1L]], "pamm")
  expect_is(summary(pam)[[1L]], "summary.gam")
  expect_data_frame(int_info(pam[[1L]]), nrows = 36L, ncols = 5L)
  expect_data_frame(int_info(pam[[2L]]), nrows = 14L, ncols = 5L)
  expect_identical(is.pamm(pam[[1L]]), TRUE)

  # bam engine
  pam2 <- pamm_cr(ped_status ~ s(tend, k = 3L) + age, data = ped,
                  engine = "bam")
  expect_is(pam2, "pamm_cr_list")
  expect_is(pam2, "pamm_cr")
  expect_is(pam2[[1L]], "pamm")
  expect_is(pam2[[2L]], "bam")
  expect_is(summary(pam2)[[1L]], "summary.gam")
  expect_data_frame(int_info(pam2[[1L]]), nrows = 36L, ncols = 5L)
  expect_data_frame(int_info(pam2[[2L]]), nrows = 14L, ncols = 5L)
  expect_identical(is.pamm(pam2[[1L]]), TRUE)
  # pass arguments to bam
  pam3 <- pamm_cr(ped_status ~ s(tend, k = 3L) + age, data = ped,
                  engine = "bam", discrete = TRUE, method = "fREML")
  expect_is(pam3, "pamm_cr_list")
  expect_is(pam3, "pamm_cr")
  expect_is(pam3[[1L]], "pamm")
  expect_is(pam3[[2L]], "bam")
  expect_is(summary(pam3)[[1L]], "summary.gam")
  expect_data_frame(int_info(pam3[[1L]]), nrows = 36L, ncols = 5L)
  expect_data_frame(int_info(pam3[[2L]]), nrows = 14L, ncols = 5L)
  expect_identical(is.pamm(pam3[[1L]]), TRUE)

})
