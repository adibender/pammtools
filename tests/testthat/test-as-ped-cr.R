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
  expect_is(attr(ped, "breaks"), "list")
  expect_is(attr(ped, "breaks")[[1]], "numeric")
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
  ped <- as_ped_cr(
    data         = sir_adm,
    formula      = Surv(time, status) ~ .
  )
  expect_data_frame(ped, nrows = 56L, ncols = 10L)
  expect_equal(sum(as.numeric(ped$cause)), 84L)

  ped_list <- as_ped_cr(sir_adm, Surv(time, status) ~ ., combine = FALSE)
  expect_identical(attr(ped_list[[1]], "breaks"), c(4L, 10L, 24L, 37L, 101L))
  expect_identical(attr(ped_list[[2]], "breaks"), c(22L, 25L))
  ped_list2 <- as_ped(ped_list, newdata = sir_adm)
  expect_equal(ped_list, ped_list2, check.attributes = FALSE)


})

test_that("Default cut behaviour works.", {
  # preparations
  ped_union <- as_ped_cr(
    data    = sir_adm,
    formula = Surv(time, status) ~ age + pneu,
    combine = TRUE)
  ped_list <- as_ped_cr(
    data    = sir_adm,
    formula = Surv(time, status) ~ age + pneu,
    combine = FALSE)
  ped_union_1 <- ped_union[ped_union$cause == 1, 1:8]
  ped_union_2 <- ped_union[ped_union$cause == 2, 1:8]
  expect_equal(unique(ped_union$tend), c(4, 10, 22, 24, 25, 37, 101))
  expect_equal(unique(ped_union_1$tend), unique(ped_union$tend))
  expect_equal(unique(ped_union_1$tend), unique(ped_union$tend))
  expect_equal(unique(ped_list$`1`$tend), c(4, 10, 24, 37, 101))
  expect_equal(unique(ped_list$`2`$tend), c(22, 25))
})

test_that("Trafo works for list objects (with TDCs)", {
  data("patient")
  event_df  <- filter(patient, CombinedID %in% c(1110, 1116, 1316))
  event_df$PatientDied[3] <- 2L
  ped <- as_ped_cr(data = list(event_df),
                   formula = Surv(survhosp, PatientDied) ~ .,
                   cut = 0:30, id = "CombinedID")
  expect_data_frame(ped, nrows = 70 * 2, ncols = 16)
  tdc_df    <- filter(daily, CombinedID  %in% c(1110, 1116, 1316))
  ## check nesting
  ped <- as_ped_cr(
    data    = list(event_df, tdc_df),
    formula = Surv(survhosp, PatientDied) ~ . +
      cumulative(survhosp, Study_Day, caloriesPercentage, tz_var = "Study_Day") +
      cumulative(proteinGproKG, tz_var = "Study_Day"),
    cut     = 0:30,
    id  = "CombinedID")
  expect_subset("survhosp_Study_Day_mat", colnames(ped))
  expect_data_frame(ped, nrows = 140, ncols = 21)
  expect_identical(any(is.na(ped$caloriesPercentage_Study_Day)), FALSE)
  expect_identical(colnames(ped$Study_Day), paste0("Study_Day", 1:12))
  ped <- as_ped_cr(
    data    = list(event_df, tdc_df),
    formula = Surv(survhosp, PatientDied) ~ . +
      cumulative(Study_Day, caloriesPercentage, tz_var = "Study_Day") +
      cumulative(proteinGproKG, tz_var = "Study_Day"),
    id  = "CombinedID")
  expect_data_frame(ped, nrows = 2 * ((2 * 2) + 1), ncols = 20)
})

test_that("Single risk not accepted.", {
  data("patient")
  event_df  <- filter(patient, CombinedID %in% c(1110, 1116, 1316))
  expect_error(
    as_ped_cr(data = list(event_df), formula = Surv(survhosp, PatientDied)~ .,
              cut = 0:30, id = "CombinedID"),
    "There are no competing risks. Use as_ped() instead.", fixed = TRUE)
})

test_that("Tibbles are supported.", {
  # preparations
  ped <- as_ped_cr(
    data    = tibble(sir_adm),
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100)
    )
  # retransform to ped
  expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
  expect_is(ped, "ped_cr_union")
  expect_subset(c("ped_status", "tstart", "tend",
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "breaks"), "list")
  expect_is(attr(ped, "breaks")[[1]], "numeric")
  expect_is(attr(ped, "breaks")[[2]], "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(sum(as.numeric(ped$cause)), 36)
})

# if(requireNamespace("data.table")) {
# test_that("data.tables are supported.", {
#   # preparations
#   ped <- as_ped_cr(
#     data    = data.table::data.table(sir_adm),
#     formula = Surv(time, status) ~ age + pneu,
#     cut     = c(0, 10, 100)
#   )
#   # retransform to ped
#   expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
#   expect_is(ped, "ped_cr_union")
#   expect_subset(c("ped_status", "tstart", "tend",
#                   "interval", "offset", "cause"), names(ped))
#   expect_is(attr(ped, "breaks"), "list")
#   expect_is(attr(ped, "breaks")[[1]], "numeric")
#   expect_is(attr(ped, "breaks")[[2]], "numeric")
#   expect_is(attr(ped, "intvars"), "character")
#   expect_is(attr(ped, "id_var"), "character")
#   expect_equal(attr(ped, "id_var"), "id")
#   expect_equal(sum(as.numeric(ped$cause)), 36)
# })
# }
