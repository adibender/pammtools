context("Test as_ped_cr functions")

test_that("Trafo works and attributes are appended for joint PED object..", {
  # preparations
  data("sir.adm", package = "mvna")
  sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]
  ped <- as_ped_cr(
    data         = sir_adm,
    formula      = Surv(time, status) ~ age + pneu,
    cut          = c(0, 10, 100)
  )
  # retransform to ped
  expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
  expect_is(ped, "ped_cr_df")
  expect_subset(c("ped_status", "tstart", "tend", 
                  "interval", "offset", "cause"), names(ped))
  expect_is(attr(ped, "breaks"), "list")
  expect_is(attr(ped, "breaks")[[1]], "numeric")
  expect_is(attr(ped, "breaks")[[2]], "numeric")
  expect_is(attr(ped, "intvars"), "character")
  expect_is(attr(ped, "id_var"), "character")
  expect_equal(attr(ped, "id_var"), "id")
  expect_equal(sum(as.numeric(ped$cause)), 36)
  
  ped <- as_ped_cr(
    data         = sir_adm,
    formula      = Surv(time, status) ~ .
  )
  expect_data_frame(ped, nrows = 56L, ncols = 10L)
  expect_equal(sum(as.numeric(ped$cause)), 84L)
  
})

test_that("Trafo for list output is identical to joint output.", {
  # preparations
  data("sir.adm", package = "mvna")
  sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]
  ped_union <- as_ped_cr(
    data    = sir_adm,
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100),
    combine = TRUE)
  ped_list <- as_ped_cr(
    data    = sir_adm,
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100),
    combine = FALSE)
    ped_union_1 <- ped_union[ped_union$cause == 1, 1:8]
    ped_union_2 <- ped_union[ped_union$cause == 2, 1:8]
    expect_equal(sum(ped_union_1[, c(1:3, 5:8)] - ped_list[[1]][, c(1:3, 5:8)]), 0)
    expect_equal(sum(ped_union_2[, c(1:3, 5:8)] - ped_list[[2]][, c(1:3, 5:8)]), 0)
    expect_equal(ped_union_1[, 4], ped_list[[1]][, 4])
    expect_equal(ped_union_2[, 4], ped_list[[2]][, 4])
})

test_that("Default cut behaviour works.", {
  # preparations
  data("sir.adm", package = "mvna")
  sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]
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
  data("sir.adm", package = "mvna")
  sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]
  ped <- as_ped_cr(
    data    = tibble(sir_adm),
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100)
    )
  # retransform to ped
  expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
  expect_is(ped, "ped_cr_df")
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

test_that("data.tables are supported.", {
  # preparations
  data("sir.adm", package = "mvna")
  sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]
  ped <- as_ped_cr(
    data    = data.table(sir_adm),
    formula = Surv(time, status) ~ age + pneu,
    cut     = c(0, 10, 100)
  )
  # retransform to ped
  expect_data_frame(ped, nrow = 12L * 2L, ncols = 9L)
  expect_is(ped, "ped_cr_df")
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

