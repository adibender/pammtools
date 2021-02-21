context("Test as_ped_cr functions")


test_that("Trafo works and attributes are appended.", {
  # preparations
  ped <- as_ped(
    data         = sir_adm,
    formula      = Surv(time, status) ~ age + pneu,
    cut          = c(0, 10, 100)
  )
  expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
  expect_is(ped, "ped_cr_union")
  expect_subset(c("ped_status", "tstart", "tend",
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(sum(as.numeric(ped$cause)), 36)

  # check that trafo can be recovered
  ped2 <- as_ped(ped, newdata = sir_adm)
  expect_equal(ped, ped2, check.attributes = FALSE)

  # check that list output identical for given cut points
  ped_list <- as_ped(
    data    = sir_adm,
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100),
    combine = FALSE)
  ped2 <- do.call(rbind, ped_list)
  expect_true(all.equal(do.call(rbind, ped_list), ped, check.attributes = FALSE))
  expect_identical(length(ped_list), 2L)
  expect_identical(class(ped_list), c("ped_cr_list", "ped_cr", "ped", "list"))
  expect_identical(names(attributes(ped_list)), c("class", "names", "trafo_args", "risks"))
  expect_identical(length(attr(ped_list, "trafo_args")$cut), 2L)

  # check that trafo can be recovered for ped list objects
  ped_list2 <- as_ped(ped_list, newdata = sir_adm)
  expect_equal(ped_list, ped_list2, check.attributes = FALSE)

  # test when split points not specified
  ped <- as_ped(data = sir_adm, formula = Surv(time, status) ~ .)
  expect_data_frame(ped, nrows = 56L, ncols = 10L)
  expect_equal(sum(as.numeric(ped$cause)), 84L)

  ped_list <- as_ped_cr(sir_adm, Surv(time, status) ~ ., combine = FALSE)
  expect_identical(attr(ped_list[[1]], "breaks"), c(4L, 10L, 24L, 37L, 101L))
  expect_identical(attr(ped_list[[2]], "breaks"), c(22L, 25L))
  ped_list2 <- as_ped(ped_list, newdata = sir_adm)
  expect_equal(ped_list, ped_list2, check.attributes = FALSE)

})

test_that("Trafo works for more than two risks.", {
  # preparations
  sir_adm$status[2] <- 3
  ped <- as_ped(
    data         = sir_adm,
    formula      = Surv(time, status) ~ age + pneu,
    cut          = c(0, 10, 100)
  )
  expect_data_frame(ped, nrow = 12L * 3L, ncols = 9L)
  expect_is(ped, "ped_cr_union")
  expect_subset(c("ped_status", "tstart", "tend",
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "breaks"), "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(sum(as.numeric(ped$cause)), 72)
  expect_equal(sum(ped$ped_status[ped$cause == 3L]), 1)

})
