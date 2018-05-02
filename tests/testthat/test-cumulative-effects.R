context("Cumulative effects (of time-dependent covariates)")

test_that("Cumulative effects are calculated correctly", {
  data("simdf_elra")
  ped <- as_ped(simdf_elra, Surv(time, status)~x1 + x2|
      cumulative(latency(tz2), z.tz2, tz_var="tz2"),
    cut = 0:10)
  pam <- mgcv::gam(ped_status ~ s(tend) + x1 + s(x2) + s(tz2_latency, by=z.tz2),
    data = ped, family=poisson(), offset=offset)
  # gg_slice(ped, pam, term="z.tz2", tz2_latency=seq(0, 12, by = 1), z.tz2=c(1))

  ndf <- make_newdata(ped, tz2_latency = seq(0, 6, by = 1), z.tz2=c(1))
  ndf <- ndf %>% add_hazard(pam)
  expect_equal(ndf$hazard, c(0.15, 0.17, 0.19, 0.20, 0.21, 0.21, 0.20), tolerance=10e-2)

  # tz grid with differences different than 1

  # generate exposures and add to data set
  n <- 250
  set.seed(123)
  # create data set with variables which will affect the hazard rate.
  df <- cbind.data.frame(x1 = runif(n, -3, 3), x2 = runif(n, 0, 6)) %>%
   tibble::as_tibble()
  # the formula which specifies how covariates affet the hazard rate
  f0 <- function(t) {
   dgamma(t, 8, 2) *6
  }
  tz3 <- c(-5, -3, 0, 3, 5)
  rng_z = function(nz) {
    as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
  }
  df <- df %>% add_tdc(tz3, rng_z)


  sim_df <- sim_pexp(
    formula = ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)|
    fcumu(t, tz3, z.tz3,
      f_xyz=function(t, tz, z) 5*(dnorm(t-tz,4,6)+dnorm(t-tz,25,4))*z,
      ll_fun=function(t, tz) {t - 2 >= tz}),
  data = df,
  cut = 0:10)

  ped <- as_ped(sim_df, Surv(time, status)~x1 + x2|
      cumulative(latency(tz3), z.tz3, tz_var="tz3"),
    cut = 0:10)
  expect_identical(ped$LL[1,], c(1:3, rep(0,2)))
  expect_identical(ped$LL[9,], c(1, 2, 3, 3, 2))
  expect_identical(ped$LL[10,], c(1, 2, 3, 3, 2))
  pam <- mgcv::gam(ped_status ~ s(tend) + x1 + s(x2) + s(tz3_latency, by=z.tz3),
    data = ped, family=poisson, offset=offset)
  # gg_slice(ped, pam, term="z.tz3", z.tz3=c(1), tz3_latency=seq(0, 12, by=1))
  ndf <- make_newdata(ped, tz3_latency = seq(0, 6, by = 1), z.tz3=c(1))
  ndf <- ndf %>% add_hazard(pam)
  expect_equal(ndf$hazard, c(0.11, 0.11, 0.12, 0.13, 0.14, 0.14, 0.14), tolerance=10e-2)
})
