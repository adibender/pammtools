context("Test simulation functions")

test_that("Test that rpexp works", {
  expect_identical(length(rpexp(n = 1, rate = 1, t = 0)), 1L)
  expect_error(rpexp(n = 1, rate = 1, t = c(0, 1)))
  expect_error(rpexp(n = 1, rate = 1, t = 1))
  expect_error(rpexp(n = 1, rate = c(1, 1, 1), t = c(0, 2, 1)))
  expect_identical(rpexp(n = 0, rate = 1, t = 0), numeric(0))
  expect_identical(length(rpexp(n = c(1, 3), rate = 1, t = 0)), 2L)
})

test_that("Simulation function works", {
  suppressWarnings(RNGversion("3.5.0"))

  set.seed(24032018)
  # standard data
  df <- cbind.data.frame(x1 = runif(3, -3, 3), x2 = runif(3, 0, 6))
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
    formula = ~ -3.5 -
      0.5 * x1 +
      sqrt(x2) |
      fcumu(
        t,
        tz1,
        z.tz1,
        f_xyz = function(t, tz, z) {
          -1 *
            cos(t / 10 * pi) *
            0.8 *
            (dnorm(z, 1.5, 2) + 1.5 * dnorm(z, 7.5, 1)) *
            15 *
            dnorm(t - tz, 8, 10)
        },
        ll_fun = function(t, tz) t >= tz
      ),
    data = df,
    cut = 0:10
  )
})

test_that("Competing risks simulation creates valid output", {
  set.seed(1202)
  sim_df <- pammtools:::sim_pexp_cr(
    formula = ~ -3 + 0.2 * x1 | -2 - 0.1 * x1,
    data = data.frame(x1 = c(-0.5, 0.2, 1.1, -1.2)),
    cut = 0:5
  )

  expect_data_frame(sim_df, nrows = 4, ncols = 8)
  expect_true(all(
    c(
      "id",
      "x1",
      "t",
      "hazard1",
      "hazard2",
      "time",
      "status",
      "type"
    ) %in%
      names(sim_df)
  ))
  expect_identical(sort(unique(sim_df$id)), seq_len(4L))
  expect_true(all(sim_df$t == 0))
  expect_true(all(sim_df$time >= 0 & sim_df$time <= 5))
  expect_true(all(sim_df$status %in% c(0L, 1L)))
  expect_true(all(sim_df$type %in% c(1L, 2L)))
  expect_true(all(sim_df$hazard1 > 0))
  expect_true(all(sim_df$hazard2 > 0))
})
