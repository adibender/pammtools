context("Cumulative effects (of time-dependent covariates)")


test_that("Lag-lead is calculated correctly", {
  LL <- get_laglead(0:2, c(-2, -0.5, 0, 0.5, 2), ll_fun=function(t,tz) {t >= tz})
  expect_identical(LL$t, rep(0:2, each = 5))
  expect_identical(LL$tz, rep(c(-2, -0.5, 0, 0.5, 2), times = 3))
  expect_equal(LL$LL, c(rep(1, 3), rep(0,2), rep(1,4), 0, rep(1,5)))
})

test_that("Cumulative effects are calculated correctly", {
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
  # plot(survfit(Surv(time, status)~1, data=sim_df))

  ped <- as_ped(sim_df, Surv(time, status)~x1 + x2|
      cumulative(latency(tz3), z.tz3, tz_var="tz3"),
    cut = 0:10)
  # ped$LL[ped$LL!=0] <- 1
  ped5 <- subset(ped, id == 5)
  expect_identical(ped5$LL[1,], c(2.5, 2, 3, rep(0,2)))
  expect_identical(ped5$LL[9,], c(2.5, 2, 3, 3, 2))
  expect_identical(ped5$LL[10,], c(2.5, 2, 3, 3, 2))
  pam <- mgcv::gam(ped_status ~ s(tend) + x1 + s(x2) + s(tz3_latency, by=z.tz3),
    data = ped, family=poisson, offset=offset)
  # gg_slice(ped, pam, term="z.tz3", z.tz3=c(1), tz3_latency=unique(tz3_latency))
  ndf <- make_newdata(ped, tz3_latency = unique(tz3_latency), z.tz3=c(1))
  ndf <- ndf %>% add_term(pam, term = "z.tz3") %>% slice(1:7)
  expect_equal(ndf$fit, c(0.31, 0.30, 0.28, 0.25, 0.21, 0.17, 0.13),
    tolerance=10e-3)
})
