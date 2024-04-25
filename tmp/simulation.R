# Pammtools
setwd("C:/Users/ra63liw/Documents/98_git/pammtools-multi-state/pammtools")
library(dplyr)
library(survival)
library(Formula)
devtools::load_all()
source("sandbox/helpers-msm.R")
# source("sandbox/sim-pexp-msm.R")
library(mstate)
# library(msm)
# library(dynpred)

# ------------------------------------------------------------------------------
#
# EXAMPLE FROM VIGNETTE / R HELPER
#
# ------------------------------------------------------------------------------
library(survival)
library(dplyr)
# library(pammtools)

# set number of observations/subjects
n <- 250
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
  as_tibble()
# the formula which specifies how covariates affect the hazard rate
f0 <- function(t) {
  dgamma(t, 8, 2) *6
}
form <- ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)
set.seed(24032018)
sim_df <- sim_pexp(form, df, 1:10)
head(sim_df)
plot(survfit(Surv(time, status)~1, data = sim_df ))

# for control, estimate with Cox PH
mod <- coxph(Surv(time, status) ~ x1 + pspline(x2), data=sim_df)
coef(mod)[1]
layout(matrix(1:2, nrow=1))
termplot(mod, se = TRUE)

# and using PAMs
layout(1)
ped <- sim_df %>% as_ped(Surv(time, status)~., max_time=10)
library(mgcv)
pam <- gam(ped_status ~ s(tend) + x1 + s(x2), data=ped, family=poisson, offset=offset)
coef(pam)[2]
plot(pam, page=1)

# ------------------------------------------------------------------------------
#
# MULTI-STATE EXAMPLE FROM VIGNETTE / R HELPER
#
# ------------------------------------------------------------------------------

t_mat <- matrix(data = NA, nrow = 3, ncol = 3)
t_mat[1,2] <- "log(0.7) + x1"
t_mat[1,3] <- "log(0.5) - 0.25 * x2"
t_mat[2,3] <- "log(0.9)"

n = 100
data <- cbind.data.frame(
 id = 1:n,
 x1 = runif(n, -3, 3),
 x2 = runif(n, 0, 6),
 from = 1,
 t = 0)

f1 <- function(x) log(0.7) + x
f2 <- function(x) log(0.5) - 0.25 * x

plot(data$x1, f1(data$x1))
plot(data$x2, f2(data$x2))

cut =  seq(0, 3, by = 0.01)

msm_df <- sim_pexp_msm_debug(
 t_mat = t_mat,
 data = data,
 cut = cut,
 keep_transitions_at_risk = TRUE)
head(msm_df)

# created the simulated data set. Transform data to be used to calculate transition probabilities

cal_msm_df <- as_ped_multistate(
  data       = msm_df,
  formula    = Surv(tstart, tstop, status)~ x1 + x2,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_msm_df)
head(cal_msm_df)

ctrl <- gam.control(trace = TRUE)

pam_msm <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                  + as.factor(transition)
                  + s(x1)
                  + x2
                  , data = cal_msm_df
                  , family=poisson()
                  , offset=offset
                  , control = ctrl)

summary(pam_msm)

plot(pam_msm, page = 1)

new_msm_pam <- make_newdata(cal_msm_df
                              , tend = unique(tend)
                              , transition=unique(transition)
                              #, age = quantile(age, probs=c(0.05, 0.5, 0.95))
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_msm) 

old_groups <- dplyr::groups(new_msm_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_msm_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 4, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 
