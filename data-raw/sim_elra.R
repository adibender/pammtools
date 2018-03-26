library(dplyr)
# set number of observations/subjects
n <- 250
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif(n, -3, 3), x2 = runif(n, 0, 6)) %>%
 as_tibble()
rng_z = function(nz) {
  as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
}
# two different exposure times  for two different exposures
te1 <- 1:10
te2 <- -5:5
# generate exposures and add to data set
df <- df %>%
  add_tdc(te1, rng_z) %>%
  add_tdc(te2, rng_z)
df

# define tri-variate function of time, exposure time and exposure z
ft <- function(t, tmax) {
  -1*cos(t/tmax*pi)
}
fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))
wpeak2 <- function(lag) 15*dnorm(lag,8,10)
wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))
f_ttez1 <- function(t, te, z) {
  ft(t, tmax=10) * 0.8*fdnorm(z)* wpeak2(t-te)
}
f_ttez2 <- function(t, te, z) {
  wdnorm(t-te)*z
}
# define lag-lead window function
ll_fun <- function(t, te) {t >= te}
# simulate data with cumulative effect
simdf_elra <- sim_pexp(
  formula = ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)|
     fcumu(t, te1, z.te1, f_xyz=f_ttez1, ll_fun=ll_fun) +
     fcumu(t, te2, z.te2, f_xyz=f_ttez2, ll_fun=ll_fun),
  data = df,
  cut = 0:10)

devtools::use_data(simdf_elra, overwrite=TRUE)
