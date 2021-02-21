context("Test simulation functions")

test_that("Test that rpexp works", {
  expect_identical(length(rpexp(n = 1, rate = 1, t = 0)), 1L)
  expect_error(rpexp(n = 1, rate = 1, t = c(0, 1)))
  expect_error(rpexp(n = 1, rate = 1, t = 1))
  expect_error(rpexp(n = 1, rate = c(1, 1, 1), t = c(0, 2, 1)))
  expect_identical(rpexp(n=0, rate = 1, t = 0), numeric(0))
  expect_identical(length(rpexp(n=c(1, 3), rate = 1, t = 0)), 2L)
})

test_that("Simulation function works", {

  suppressWarnings(RNGversion("3.5.0"))

  set.seed(24032018)
  # standard data
  df     <- cbind.data.frame(x1 = runif (3, -3, 3), x2 = runif (3, 0, 6))
  sim_df <- sim_pexp(~ -3.5 - 0.5 * x1 + sqrt(x2), df, cut = 0:10)
  expect_data_frame(sim_df, nrows = 3, ncols = 5)
  expect_identical(round(sim_df$time, 2), c(1.38, 7.14, 3.02))

  # time-dependent covariates
  rng_z <- function(nz) {
    as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
  }
  tz1 <- 1:10
  df <- df %>% add_tdc(tz1, rng_z)

 # simulate data with cumulative effect
 sim_df <- sim_pexp(
   formula = ~ -3.5 - 0.5 * x1 + sqrt(x2) |
      fcumu(t, tz1, z.tz1,
        f_xyz = function(t, tz, z) {
          -1 * cos(t / 10 * pi) *
          0.8 * (dnorm(z, 1.5, 2) + 1.5 * dnorm(z, 7.5, 1)) *
          15 * dnorm(t - tz, 8, 10)
        },
        ll_fun = function(t, tz) t >= tz),
   data = df,
   cut = 0:10)

})
