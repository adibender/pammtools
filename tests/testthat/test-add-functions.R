context("Convenience functions for calculation of hazard and similar")

data("tumor")
ped <- tumor[1:200,] %>% as_ped(Surv(days, status)~ age + complications,
  cut = c(0, 50, 100, 200, 300, 400))

pam <- mgcv::gam(ped_status ~ s(tend, k = 5) + complications, data = ped,
  family = poisson(), offset = offset)
pam2 <- mgcv::gam(ped_status ~ s(tend, k = 5) + complications + s(age), data = ped,
  family = poisson(), offset = offset)
bam <- mgcv::bam(ped_status ~ s(tend, k = 5) + complications, data = ped,
  family = poisson(), offset = offset, method = "fREML", discrete = TRUE)
pem <- glm(ped_status ~ 0 + interval + complications, data = ped,
  family = poisson(), offset = offset)

pam3 <- mgcv::gam(ped_status ~ s(tend, k = 5, by = as.factor(complications)) + as.factor(complications),
  data = ped, family = poisson(), offset = offset)

data("prothr", package = "mstate")
prothr <- prothr[1:200,] |> 
  filter(Tstart != Tstop) |> 
  mutate(transition = as.factor(paste0(from, "->", to))) |> 
  select(-trans, -treat)
ped_msm <- as_ped(
  data       = prothr,
  formula    = Surv(Tstart, Tstop, status)~ .,
  # cut = seq(2000, 4000, length.out = 100),
  transition = "transition",
  id         = "id",
  timescale  = "calendar",
  tdc_specials="concurrent"
)

pam_msm <- gam(ped_status ~ s(tend, by=transition, bs="cr") + transition
               , data = ped_msm
               , family = poisson()
               , offset = offset)

test_that("hazard functions work for PAM", {

  expect_data_frame(haz <- add_hazard(ped_info(ped), bam), nrows = 5L,
    ncols = 11L)
  expect_data_frame(haz <- add_hazard(ped_info(ped), pam), nrows = 5L,
    ncols = 11L)
  expect_equal(all(haz$ci_lower < haz$hazard), TRUE)
  expect_equal(all(haz$ci_upper > haz$hazard), TRUE)
  expect_equal(round(haz$hazard, 3), c(0.001, 0.001, 0.001, 0.001, 0.001))
  expect_equal(round(haz$ci_lower, 3), c(0, 0, 0, 0, 0))
  expect_error(add_hazard(haz, pam))
  expect_data_frame(add_hazard(haz, pam, overwrite = TRUE),
    nrows = 5L, ncols = 11L)

  haz2 <- add_hazard(ped_info(ped), pam, type = "link")
  expect_equal(all(haz2$ci_lower < haz2$hazard), TRUE)
  expect_equal(all(haz2$ci_upper > haz2$hazard), TRUE)
  expect_equal(round(haz2$hazard, 2), c(-7.37, -7.39, -7.41, -7.43, -7.46))
  expect_equal(round(haz2$ci_lower, 2), c(-7.93, -7.86, -7.78, -7.83, -7.99))

  ## delta rule
  expect_data_frame(add_hazard(ped_info(ped), bam, ci_type = "delta"),
    nrows = 5L, ncols = 11L)
  haz3 <- add_hazard(ped_info(ped), pam, ci_type = "delta")
  expect_data_frame(haz3, nrows = 5L, ncols = 11L)
  expect_equal(round(haz3$hazard * 100, 2), c(.06, .06, .06, .06, .06))
  expect_equal(round(haz3$se * 100, 2), c(.02, .01, .01, .01, .02))
  expect_equal(round(haz3$ci_lower * 100, 2), c(.03, .03, .04, .04, .03))
  expect_equal(round(haz3$ci_upper * 100, 2), c(.10, .09, .08, .08, .09))

  ## simulation based ci (0.95)
  haz4 <- add_hazard(ped_info(ped), pam, ci_type = "sim")

  ## hazard with reference (i.e. hazard ratio)
  hr <- add_hazard(ped_info(ped), pam2, reference = list(age = c(30)))
  # hazard ratio is constant as age effect not time-varying
  expect_equal(round(hr$hazard, 3), rep(1.458, 5))
  # hr = 1 if reference = data
  hr2 <- ped_info(ped) %>%  add_hazard(pam2, reference = list(age = mean(.$age)))
  expect_equal(hr2$hazard, rep(1, 5))

  ## factor group variable
  ndf <- ped %>% make_newdata(tend = unique(tend), complications = unique(complications)) %>%
    group_by(complications)
  ndf1 <- ndf %>% add_cumu_hazard(pam3, ci = TRUE, ci_type = "default")
  ndf2 <- ndf %>% add_cumu_hazard(pam3, ci = TRUE, ci_type = "delta")
  ndf3 <- ndf %>% add_cumu_hazard(pam3, ci = TRUE, ci_type = "sim", nsim = 100L)
  expect_true(all(ndf1$cumu_hazard > ndf1$cumu_lower & ndf1$cumu_hazard < ndf1$cumu_upper))
  expect_true(all(ndf2$cumu_hazard > ndf2$cumu_lower & ndf2$cumu_hazard < ndf2$cumu_upper))
  expect_true(all(ndf3$cumu_hazard > ndf3$cumu_lower & ndf3$cumu_hazard < ndf3$cumu_upper))

})

test_that("hazard functions work for PEM", {

  expect_data_frame(haz <- add_hazard(ped_info(ped), pem),
    nrows = 5L, ncols = 11L)
  expect_error(add_hazard(haz, pem))
  expect_data_frame(add_hazard(haz, pem, overwrite = TRUE),
    nrows = 5L, ncols = 11L)

})


test_that("cumulative hazard functions work for PAM", {

  expect_data_frame(add_cumu_hazard(ped_info(ped), bam, ci = FALSE),
   nrows = 5L, ncols = 8L)
 expect_data_frame(haz <- add_cumu_hazard(ped_info(ped), pam, ci = FALSE),
   nrows = 5L, ncols = 8L)
  expect_data_frame(haz <- add_cumu_hazard(ped_info(ped), pam),
    nrows = 5L, ncols = 10L)
  expect_equal(round(haz$cumu_hazard, 2), c(.03, .06, .12, .18, .24))
  expect_equal(round(haz$cumu_lower, 2), c(.02, .04, .08, .12, .15))
  expect_equal(all(diff(haz$cumu_hazard) >= 0), TRUE)
  # overwrite works
  expect_data_frame(add_cumu_hazard(haz, pam, overwrite = TRUE),
    nrows = 5L, ncols = 10L)

  # error on wrong input
  expect_error(add_cumu_hazard(haz, pam))

  ## test that cumu_hazard works for grouped data
  grouped_haz <- ped %>% group_by(complications) %>%
    ped_info() %>%
    add_cumu_hazard(pam)
  expect_data_frame(grouped_haz, nrows = 10L, ncols = 10L)
  expect_equal(round(grouped_haz$cumu_hazard, 2),
    c(.03, .06, .12, .18, .24, .06, .13, .25, .37, .49))

  ## delta method
  haz2 <- ped_info(ped) %>% add_cumu_hazard(pam, ci_type = "delta")
  expect_equal(round(haz2$cumu_upper, 2), c(.05, .09, .18, .25, .33))
  expect_equal(round(haz2$cumu_lower, 2), c(.01, .03, .07, .11, .15))

  suppressWarnings(RNGversion("3.5.0"))
  ## sim CI (0.95)
  set.seed(123)
  haz3 <- ped_info(ped) %>% add_cumu_hazard(pam, ci_type = "sim")
  expect_equal(round(haz3$cumu_upper, 2), c(.06, .11, .19, .25, .34))
  expect_equal(round(haz3$cumu_lower, 2), c(.02, .04, .08, .13, .17))

  ## check that hazard columns are not deleted
  newdata <- ped_info(ped) %>% add_hazard(pam) %>%
    add_cumu_hazard(pam)
  expect_data_frame(newdata, nrows = 5L, ncols = 14L)
  newdata <- ped_info(ped) %>% add_hazard(pam, ci = FALSE) %>%
    add_cumu_hazard(pam)
  expect_data_frame(newdata, nrows = 5L, ncols = 11L)

})

test_that("cumulative hazard functions work for PEM", {

  expect_data_frame(haz <- add_cumu_hazard(ped_info(ped), pem),
    nrows = 5L, ncols = 10L)
  expect_error(add_cumu_hazard(haz, pem))
  expect_data_frame(add_cumu_hazard(haz, pem, overwrite = TRUE),
    nrows = 5L, ncols = 10L)

})

test_that("adding terms works for PAM", {

  # standard
  ndf2  <- make_newdata(ped, age = seq_range(age, 3))
  pred2 <- ndf2 %>% add_term(pam2, term = "age")
  expect_equal(round(pred2$fit, 3), c(-.604, -.236, .851))
  expect_data_frame(pred2, nrows = 3L, ncols = 12L)
  # with custom reference
  pred2 <- ndf2 %>%
    add_term(pam2, term = "age", reference = list(age = mean(.$age)))
  expect_equal(round(pred2$fit, 3), c(-.368, 0, 1.087))
  expect_data_frame(pred2, nrows = 3L, ncols = 12L)
  expect_equal(pred2$fit[2], 0)
  # with overall function application
  pred3 <- ndf2 %>% add_term(pam2, term = "age", reference = identity(.))
  expect_equal(pred3$fit, rep(0, 3))
  expect_data_frame(pred3, nrows = 3L, ncols = 12L)
  expect_equal(pred3$fit, rep(0, 3))
  # with separately created data frame
  df_mean <- sample_info(ndf2)
  pred4 <- ndf2 %>% add_term(pam2, term = "age", reference = df_mean)
  expect_equal(pred4$fit, pred2$fit)

})

test_that("adding terms works for PEM", {

  expect_data_frame(term <- add_term(ped_info(ped), pem, term = "complications"),
    nrows = 5L, ncols = 10L)
  expect_data_frame(ped_info(ped) %>%
      add_term(pem, term = "age", reference = list(age = mean(.$age))),
      nrows = 5L, ncols = 10L)
})

test_that("warns about / aborts for unknown intervals", {

  weird <- make_newdata(ped_info(ped), tend = c(150), interval = c("(1.4, 4]"))
  expect_warning(add_hazard(weird, pam), "not equivalent")
  expect_error(add_hazard(weird, pem), "not equivalent")

})

test_that("works for nonstandard baseline arguments", {

  pseudonymous <- ped %>% dplyr::rename(stop = tend, int = interval)
  pseudonymous <-  pseudonymous %>% dplyr::mutate(length = stop - tstart)
  ped <- ped %>% dplyr::mutate(intlen = tend - tstart)

  p_pam <- mgcv::gam(ped_status ~ s(stop, k = 5) + complications, data = pseudonymous,
    family = poisson(), offset = offset)
  p_pem <- glm(ped_status ~ 0 + int + complications, data = pseudonymous,
    family = poisson(), offset = offset)
  expect_equal(
    add_hazard(pseudonymous[1:5, ], p_pam, time_var = "stop")$hazard,
    add_hazard(ped[1:5, ], pam)$hazard)
  expect_equal(
    add_hazard(pseudonymous[1:5, ], p_pem, time_var = "int")$hazard,
    add_hazard(ped[1:5, ], pem)$hazard)

  expect_equal(
    add_cumu_hazard(pseudonymous[1:5, ], p_pam, time_var = "stop",
      interval_length = length)$cumu_hazard,
    add_cumu_hazard(ped[1:5, ], pam)$cumu_hazard)
  expect_equal(
    add_cumu_hazard(pseudonymous[1:5, ], p_pem, time_var = "int",
      interval_length = length)$cumu_hazard,
    add_cumu_hazard(ped[1:5, ], pem)$cumu_hazard)
  expect_equal(
    add_cumu_hazard(pseudonymous[1:5, ], p_pem, time_var = "int",
      interval_length = "length")$cumu_hazard,
    add_cumu_hazard(ped[1:5, ], pem)$cumu_hazard)

})


## test surv_prob
test_that("survival probabilities functions work for PAM", {

  suppressWarnings(RNGversion("3.5.0"))

  expect_data_frame(add_surv_prob(ped_info(ped), bam, ci = FALSE),
    nrows = 5L, ncols = 8L)
  expect_data_frame(surv <- add_surv_prob(ped_info(ped), pam, ci = FALSE),
    nrows = 5L, ncols = 8L)
  expect_data_frame(
    surv <- add_surv_prob(ped_info(ped), pam), nrows = 5L, ncols = 10L)
  stest <- sapply(surv[, c("surv_prob", "surv_lower", "surv_upper")],
    function(z) {
      all(z >= 0 & z <= 1)
    })
  expect_identical(all(stest), TRUE)
  expect_identical(round(surv$surv_prob, 2), c(0.97, 0.94, 0.88, 0.83, 0.79))
  expect_identical(round(surv$surv_lower, 2), c(0.95, 0.90, 0.83, 0.76, 0.68))
  expect_identical(round(surv$surv_upper, 2), c(0.98, 0.96, 0.92, 0.89, 0.86))
  # check that overwrite works
  expect_data_frame(add_surv_prob(surv, pam, overwrite = TRUE),
    nrows = 5L, ncols = 10L)
  # error on wrong input
  expect_error(add_surv_prob(surv, pam))

  ## test that cumu_hazard works for grouped data
  grouped_surv <- ped %>% group_by(complications) %>%
    ped_info() %>%
    add_surv_prob(pam)
  expect_data_frame(grouped_surv, nrows = 10L, ncols = 10L)
  expect_equal(round(grouped_surv$surv_prob, 2),
    c(0.97, 0.94, 0.88, .83, .79, .94, 0.88, .78, .69, .61))

  ## delta CI
  surv2 <- add_surv_prob(ped_info(ped), pam, ci_type = "delta")
  expect_equal(round(surv2$surv_lower, 2), c(.95, .91, .84, .78, .72))
  expect_equal(round(surv2$surv_upper, 2), c(.99, .97, .93, .89, .86))

  # sim CI
  set.seed(123)
  surv3 <- add_surv_prob(ped_info(ped), pam, ci_type = "sim")
  expect_equal(round(surv3$surv_lower, 2), c(.94, .90, .83, .78, .71))
  expect_equal(round(surv3$surv_upper, 2), c(.98, .96, .92, .88, .84))

})

test_that("CIF works", {

  set.seed(211758)
  df <- data.frame(time = rexp(20), status = sample(c(0,1, 2), 20, replace = TRUE))
  ped_cr <- as_ped(df, Surv(time, status)~., id = "id") %>%
    mutate(cause = as.factor(cause))
  pam <- pamm(ped_status ~ s(tend, by = cause), data = ped_cr)
  ndf <- ped_cr %>%
    make_newdata(tend = unique(tend), cause = unique(cause)) %>%
    group_by(cause) %>%
    add_cif(pam)
  expect_data_frame(ndf, nrows = 26L, ncols = 11L)
  expect_subset(c("cif", "cif_lower", "cif_upper"), colnames(ndf))
  expect_true(all(ndf$cif < ndf$cif_upper))
  expect_true(all(ndf$cif > ndf$cif_lower))
  expect_true(all(ndf$cif <= 1 & ndf$cif >= 0))
  expect_true(all(ndf$cif_lower <= 1 & ndf$cif_lower >= 0))
  expect_true(all(ndf$cif_upper <= 1 & ndf$cif_upper >= 0))

})

test_that("Transition Probability works", {
  
  ndf <- ped_msm %>%
    make_newdata(tend = unique(tend), transition = unique(transition)) %>%
    group_by(transition) %>%
    arrange(transition, tend) |>
    add_trans_prob(pam_msm, ci=T)
  expect_data_frame(ndf, nrows = 380L, ncols = 11L) 
  expect_subset(c("trans_prob", "trans_lower", "trans_upper"), colnames(ndf))
  expect_true(all(ndf$trans_prob < ndf$trans_upper))
  expect_true(all(ndf$trans_prob > ndf$trans_lower))
  expect_true(all(ndf$trans_prob <= 1 & ndf$trans_prob >= 0))
  expect_true(all(ndf$trans_lower <= 1 & ndf$trans_lower >= 0))
  expect_true(all(ndf$trans_upper <= 1 & ndf$trans_upper >= 0))
  
})
