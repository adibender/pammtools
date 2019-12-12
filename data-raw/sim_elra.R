library(dplyr)
# set number of observations/subjects
n <- 250
set.seed(8042018)
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
 as_tibble()
rng_z <- function(nz) {
  as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
}
# two different exposure times  for two different exposures
tz1 <- 1:10
tz2 <- -5:5
# generate exposures and add to data set
df <- df %>%
  add_tdc(tz1, rng_z) %>%
  add_tdc(tz2, rng_z)
df

# baseline hazard
f0 <- function(t) {
 dgamma(t, 8, 2) * 6
}
# define tri-variate function of time, exposure time and exposure z
ft <- function(t, tmax) {
  -1 * cos(t / tmax * pi)
}

f_xyz1 <- function(t, tz, z) {
  ft(t, tmax = 10) * 0.8 * (dnorm(z, 1.5, 2) + 1.5 * dnorm(z, 7.5, 1)) * 15 *
    dnorm(t - tz, 8, 10)
}
f_xyz2 <- function(t, tz, z) {
  5 * (dnorm(t - tz, 4, 6) + dnorm(t - tz, 25, 4)) * z
}
# define lag-lead window function
ll_fun <- function(t, tz)  t >= tz
ll_fun2 <- function(t, tz) t >= tz + 2

# simulate data with cumulative effect
simdf_elra <- sim_pexp(
  formula = ~ -3.5 + f0(t) - 0.5 * x1 + sqrt(x2) |
     fcumu(t, tz1, z.tz1, f_xyz = f_xyz1, ll_fun = ll_fun) +
     fcumu(t, tz2, z.tz2, f_xyz = f_xyz2, ll_fun = ll_fun2),
  data = df,
  cut = 0:10)

usethis::use_data(simdf_elra, overwrite = TRUE)
