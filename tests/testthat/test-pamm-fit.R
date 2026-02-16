context("Test pamm wrapper function")


 test_that("pamm function works correctly", {

  data("tumor")
  ped <- as_ped(Surv(days, status)~ complications + age, data = tumor[1:20,])
  # gam engine
  pam <- pamm(ped_status ~ s(tend, k=3) + age, data=ped)
  expect_is(pam, "pamm")
  expect_is(summary(pam), "summary.gam")
  expect_data_frame(int_info(pam), nrows = 9L, ncols = 5L)
  expect_identical(is.pamm(pam), TRUE)

  # check data trafo from pam object
  ped_new <- as_ped(pam, newdata = tumor[21:40, ])
  expect_data_frame(ped_new, nrows = 144L, ncols = 8L)
  expect_subset(ped_new$tend, ped$tend)

  # bam engine
  pam2 <- pamm(ped_status ~ s(tend, k = 3) + age, data = ped, engine = "bam")
  expect_true(inherits(pam2, "bam"))
  expect_data_frame(int_info(pam2), nrows = 9L, ncols = 5L)
  expect_identical(is.pamm(pam2), TRUE)
  # pass arguments to bam
  pam3 <- pamm(ped_status ~ s(tend, k = 3) + age, data = ped,
    engine = "bam", discrete = TRUE, method = "fREML")
  expect_true(inherits(pam3, "bam"))
  expect_data_frame(int_info(pam3), nrows = 9L, ncols = 5L)
  expect_identical(is.pamm(pam), TRUE)
  
  # warning if no offset in ped data
  ped_nooffset <- ped |> select(-offset)
  expect_warning(pamm(ped_status ~ s(tend, k=3) + age, data = ped_nooffset))

 })

test_that("pamm trafo_args uses as_ped for list-based data", {

  data("pbc", package = "survival")
  event_df <- pbc %>%
    filter(id <= 5) %>%
    mutate(status = 1L * (status == 2)) %>%
    select(id, time, status, sex)
  tdc_df <- pbcseq %>%
    filter(id <= 5) %>%
    select(id, day, bili)

  pam <- pamm(
    formula = ped_status ~ s(tend, k = 3),
    data = list(event_df, tdc_df),
    trafo_args = list(
      formula = Surv(time, status) ~ . + concurrent(bili, tz_var = "day"),
      id = "id"
    )
  )

  expect_is(pam, "pamm")
  expect_true(
    grepl(
      "concurrent",
      paste(deparse(pam[["trafo_args"]][["formula"]]), collapse = " ")
    )
  )

})
