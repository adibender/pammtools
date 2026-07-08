# Regression tests for the safeguard that turns silently-wrong cumulative
# results (from a forgotten group_by()) into an explicit error. See
# stop_if_undergrouped_for_cumulation(). Uses only synthetic data.

make_strata_ped <- function(seed = 1L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    grp = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  as_ped(df, Surv(time, status) ~ grp)
}

make_cr_ped <- function(seed = 2L, n = 300L) {
  set.seed(seed)
  df <- data.frame(time = rexp(n, 0.2), status = sample(0:2, n, replace = TRUE))
  as_ped(df, Surv(time, status) ~ ., id = "id") %>%
    dplyr::mutate(cause = as.factor(cause))
}

test_that("add_surv_prob / add_cumu_hazard error on ungrouped multi-profile newdata", {
  ped <- make_strata_ped()
  pam <- pamm(ped_status ~ s(tend) + grp, data = ped)
  nd <- ped %>% make_newdata(tend = unique(tend), grp = unique(grp))

  expect_error(add_surv_prob(nd, pam), "repeated 'tend' values within a group")
  expect_error(
    add_cumu_hazard(nd, pam),
    "repeated 'tend' values within a group"
  )

  # correctly grouped -> no error
  expect_error(nd %>% dplyr::group_by(grp) %>% add_surv_prob(pam), NA)
  expect_error(nd %>% dplyr::group_by(grp) %>% add_cumu_hazard(pam), NA)
})

test_that("single-profile newdata needs no grouping (no false positive)", {
  ped <- make_strata_ped()
  pam0 <- pamm(ped_status ~ s(tend), data = ped)
  nd0 <- ped %>% make_newdata(tend = unique(tend))
  expect_error(add_surv_prob(nd0, pam0), NA)
  expect_error(add_cumu_hazard(nd0, pam0), NA)
})

test_that("check_grouping = FALSE opts out of the safeguard", {
  ped <- make_strata_ped()
  pam <- pamm(ped_status ~ s(tend) + grp, data = ped)
  nd <- ped %>% make_newdata(tend = unique(tend), grp = unique(grp))
  # opting out reproduces the old (unsafe) behaviour without erroring
  expect_error(add_surv_prob(nd, pam, check_grouping = FALSE), NA)
  expect_error(add_cumu_hazard(nd, pam, check_grouping = FALSE), NA)
})

test_that("add_cif errors on ungrouped newdata and works when grouped by cause", {
  pedc <- make_cr_ped()
  pamc <- pamm(ped_status ~ s(tend, by = cause), data = pedc)
  ndc <- pedc %>% make_newdata(tend = unique(tend), cause = unique(cause))

  expect_error(add_cif(ndc, pamc), "repeated 'tend' values within a group")
  expect_error(ndc %>% dplyr::group_by(cause) %>% add_cif(pamc, ci = FALSE), NA)
  # check_grouping = FALSE skips the profile guard, but get_cif()'s own
  # (pre-existing) cause check still catches ungrouped multi-cause input
  expect_error(
    add_cif(ndc, pamc, ci = FALSE, check_grouping = FALSE),
    "group by cause"
  )
})

test_that("add_trans_prob errors when covariate profiles are not separated", {
  skip_if_not_installed("mstate")
  data("prothr", package = "mstate")
  prothr <- prothr %>%
    dplyr::mutate(
      transition = as.factor(paste0(from, "->", to)),
      treat = as.factor(treat)
    ) %>%
    dplyr::filter(Tstart != Tstop, id <= 80) %>%
    dplyr::select(-trans)
  ped <- as_ped(
    prothr,
    Surv(Tstart, Tstop, status) ~ .,
    transition = "transition",
    id = "id",
    timescale = "calendar"
  )
  pam <- mgcv::bam(
    ped_status ~ s(tend, by = transition) + transition * treat,
    data = ped,
    family = poisson(),
    offset = offset,
    method = "fREML",
    discrete = TRUE
  )
  ndf <- make_newdata(
    ped,
    tend = unique(tend),
    treat = unique(treat),
    transition = unique(transition)
  )

  # grouped only by transition, but two `treat` profiles remain -> error
  expect_error(
    ndf %>% dplyr::group_by(transition) %>% add_trans_prob(pam),
    "repeated 'tend' values within a group"
  )
  # grouped by both profile columns -> no error
  expect_error(
    ndf %>% dplyr::group_by(treat, transition) %>% add_trans_prob(pam),
    NA
  )
  # check_grouping = FALSE opts out (old unsafe behaviour, no error)
  expect_error(
    ndf %>%
      dplyr::group_by(transition) %>%
      add_trans_prob(pam, check_grouping = FALSE),
    NA
  )
})
