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
t_mat[1,2] <- "log(0.3)"
t_mat[1,3] <- "log(0.6)"
t_mat[2,3] <- "log(0.5)"

n = 100
data <- cbind.data.frame(
 id = 1:n,
 x1 = runif(n, -3, 3),
 x2 = runif(n, 0, 6),
 from = 1,
 t = 0)

head(data)
f1 <- function(x) log(0.3)
f2 <- function(x) log(0.6)

plot(data$x1, f1(data$x1))
plot(data$x2, f2(data$x2))

cut =  seq(0, 3, by = 0.01)

msm_df <- sim_pexp_msm(
 t_mat = t_mat,
 data = data,
 cut = cut,
 keep_transitions_at_risk = TRUE)
head(msm_df)

test <- msm_df %>% filter(status == 1, transition == "1->2") %>% arrange(tstop)
head(test)
test <- msm_df %>% filter(id == 12)
View(test)
# created the simulated data set. Transform data to be used to calculate transition probabilities

msm_df <- msm_df %>% add_counterfactual_transitions()

cal_msm_df <- as_ped_multistate(
  data       = msm_df,
  formula    = Surv(tstart, tstop, status)~ transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_msm_df)
head(cal_msm_df)

ctrl <- gam.control(trace = TRUE)

pam_msm <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                  + as.factor(transition)
                  , data = cal_msm_df
                  , family=poisson()
                  , offset=offset
                  , control = ctrl)

summary(pam_msm)

plot(pam_msm, page = 1)

new_msm_pam <- make_newdata(cal_msm_df
                              , tend = unique(tend)
                              , transition=unique(transition)
                              # , x1 = seq(-3,3, by = 3)
                              # , x2 = seq(0,6, by = 3)
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
ggplot(test_msm, aes(x=tend, y=cumu_hazard)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 3, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

ggplot(data=test_msm, aes(x=tend, y=trans_prob, fill=transition)) +
  geom_area()

test_msm_01 <- test_msm %>% 
  filter(transition == "1->2")

summary(test_msm_01$trans_prob)
hist(test_msm_01$trans_prob)
summary(test_msm_01$x1)
summary(test_msm_01$tend)


# ggplot(data=test_msm_01, aes(x=tend, y=x1, z=trans_prob)) +
#   geom_tile(aes(fill=trans_prob)) +
#   scale_fill_gradient2(
#     low  = "steelblue"
#     , high = "firebrick2"
#     , midpoint=0.5)+
#   stat_contour(col="grey30") +
#   xlim(c(-3,3))

n <- 250
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
  as_tibble()
# the formula which specifies how covariates affet the hazard rate
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
  

# competing risk
# ------------------------------------------------------------------------------
#
# MULTI-STATE SIMULATION USING sim_pexp_cr(form, data, cut)
#
# ------------------------------------------------------------------------------
# try and compare to beyersmann with alpha_01 = 0.3, alpha_02 = 0.6
f0 <- function(t) {
  dgamma(t, 4, 2) *6
}

# form <- formula(  ~ 0.5 + 0.25*x1**3| f0(t))

form <- formula(  ~ 0.3 | 0.6 )


n = 100
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  x2 = runif(n, 0, 6),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)

sim_df <- sim_pexp_cr(form, data, cut) %>%
  mutate(from = 1
         , to = type + 1
         , transition = case_when(
           type == 1 ~ "1->2",
           type == 2 ~ "1->3",
           .default = "err"
           )) %>%
  rename(tstart = t, tstop = time) %>%
  filter(status == 1) %>%
  add_counterfactual_transitions()

head(sim_df)

cal_sim_df <- as_ped_multistate(
  data       = sim_df,
  formula    = Surv(tstart, tstop, status)~ transition + x1,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_sim_df)
head(cal_sim_df)

ctrl <- gam.control(trace = TRUE)

pam_sim <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) 
                     + as.factor(transition)
                     + s(x1, by=as.factor(transition))
                     , data = cal_sim_df
                     , family=poisson()
                     , offset=offset
                     , method = "fREML"
                     , discrete = T
                     , control = ctrl)

summary(pam_sim)

plot(pam_sim, page = 1)

new_sim_pam <- make_newdata(cal_sim_df
                            , tend = unique(tend)
                            , transition=unique(transition)
                            , x1 = seq(-2,2, by = 0.2)
                            # , x2 = seq(0,6, by = 3)
) %>% 
  group_by(transition, x1) %>% 
  add_cumu_hazard(pam_sim) 

old_groups <- dplyr::groups(new_sim_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_sim_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot transitions
test_contour <- test_msm %>% mutate(tend = round(tend, 1), x1 = round(x1,2))
ggplot(test_contour, aes(x=tend, y=x1, z=trans_prob)) +
  geom_tile(aes(fill=trans_prob)) +
  scale_fill_gradient2(
    name = "probability"
    , low  = "steelblue"
    , high = "firebrick2"
    , midpoint=0.5)+
  stat_contour(col="grey30",lwd = 1.1) +
  # geom_vline(xintercept = c(25, 75), lty = 3) +
  # geom_hline(yintercept = c(8, 10, 12, 14), lty = 3) +
  facet_wrap(~ transition, ncol=3, labeller = label_both) +
  xlim(c(0,2)) +
  theme_bw() +
  theme(legend.position = "bottom"
        , strip.text = element_text(size = 14)
        , axis.text = element_text(size = 14))

ggplot(test_msm, aes(x=tend, y=cumu_hazard)) + 
  geom_line(aes(group = x1, col = x1)) + 
  facet_wrap(~transition, ncol = 2, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# competing risk
# ------------------------------------------------------------------------------
#
# MULTI-STATE SIMULATION USING sim_pexp_cr(form, data, cut) MULTIPLE TIMES
#
# ------------------------------------------------------------------------------
f0 <- function(t) {
  dgamma(t, 4, 2) *6
}

test <- runif(n, 0, 3)
plot(test, -f0(test))

# # define hazard for transitions 1->2 & 1->3
# form_stepone <- formula( ~ 0.5 + 0.05*x1**3 + f0(t) | 0.6 - 0.05*x1**2  )
# # define hazard for transition 2->3
# form_steptwo <- formula( ~ - 0.4 + sin(x2) + (t+x3 <= 3) * f0(t + x3))

# define hazard for transitions 1->2 & 1->3
form_stepone <- formula( ~ 0.5 + 0.25*x1**3| f0(t)  )
# define hazard for transition 2->3
form_steptwo <- formula( ~ 0.4 - x1**2)

n = 1000
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  x2 = runif(n, 0, 2*pi),
  x3 = rep(0, n),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)

sim_df_stepone <- sim_pexp_cr(form_stepone, data, cut) %>%
  mutate(from = 1
         , to = type + 1
         , transition = case_when(
           type == 1 ~ "1->2",
           type == 2 ~ "1->3",
           .default = "err"
         )) %>%
  rename(tstart = t, tstop = time) %>%
  filter(status == 1)

head(sim_df_stepone)

# second step


# take all relevant subjects, i.e. subjects with transition 1->2, i.e. currently in state 2
data_steptwo <- sim_df_stepone %>% 
  filter(to == 2, status == 1) %>%
  mutate(x3 = tstop) %>%
  select(id, x1, x2, x3, from, t = tstop)

dim(data_steptwo)
head(data_steptwo)

cut =  seq(0, 3, by = 0.01)

sim_df_steptwo <- sim_pexp_cr(form_steptwo, data_steptwo, cut) %>%
  mutate(from = 2
         , to = type + 2
         , transition = case_when(
           type == 1 ~ "2->3",
           .default = "err"
         )
         , time = time + t
         , hazard2 = 0) %>%
  rename(tstart = t, tstop = time) %>% 
  filter(status == 1)

head(sim_df_steptwo)

sim_df <- rbind(sim_df_stepone, sim_df_steptwo)
head(sim_df)

sim_df <- sim_df %>% add_counterfactual_transitions()
head(sim_df)

table(sim_df$transition)

# go on with pamm tools procedure
cal_sim_df <- as_ped_multistate(
  data       = sim_df,
  formula    = Surv(tstart, tstop, status)~ transition + x1,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_sim_df)
head(cal_sim_df)

ctrl <- gam.control(trace = TRUE)

# pam_sim <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
#                      + as.factor(transition)
#                      + s(x1, by=as.factor(transition))
#                      + s(x2, by=as.factor(transition))
#                      , data = cal_sim_df
#                      , family=poisson()
#                      , offset=offset
#                      , control = ctrl)
# 
# summary(pam_sim)
# plot(pam_sim, page = 1)
# plot(pam_sim, select = 7, ylim = c(-1,1))
# plot(pam_sim, select = 4, xlim = c(-2,2))

bam_sim <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) 
                     + as.factor(transition)
                     + s(x1, by=as.factor(transition))
                     , data = cal_sim_df
                     , family=poisson()
                     , offset=offset
                     , discrete = T
                     , method = "fREML"
                     , control = ctrl)

summary(bam_sim)

plot(bam_sim, page = 1)

par(mfrow=c(1,3))
plot(bam_sim, select = 4, xlim = c(-2,2), ylim=c(-2,2))
plot(bam_sim, select = 5, xlim = c(-2,2), ylim=c(-2,2))
plot(bam_sim, select = 6, xlim = c(-2,2), ylim=c(-2,2))
par(mfrow=c(1,1))

pam_sim <- bam_sim

new_sim_pam <- make_newdata(cal_sim_df
                            , tend = unique(tend)
                            , transition=unique(transition)
                            , x1 = seq(-2,2, by = 2)
) %>% 
  group_by(transition
           , x1) %>% 
  add_cumu_hazard(pam_sim) 

old_groups <- dplyr::groups(new_sim_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_sim_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot hazards
ggplot(test_msm, aes(x=tend, y=cumu_hazard)) + 
  geom_line() + 
  facet_wrap(~transition + x1
             , ncol = 3
             # , scales = "free_y"
             , labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylim(c(0,5)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition + x1
             , ncol = 3
             , scales = "free_y"
             , labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 2)) +
  ylim(c(0,1)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# ------------------------------------------------------------------------------
#
# simulate time only and check if it works
#
# ------------------------------------------------------------------------------
f0 <- function(t) {
  dgamma(t, 4, 2) *6
}

test <- runif(n, 0, 3)
plot(test, f0(test))

# define hazard for transitions 1->2 & 1->3
form_stepone <- formula( ~ f0(t) | 0.6)
# define hazard for transition 2->3
form_steptwo <- formula( ~ (t+x3 <= 3) * f0(t+x3))

n = 500
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  x2 = runif(n, 0, 2*pi),
  x3 = rep(0, n),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)

sim_df_stepone <- sim_pexp_cr(form_stepone, data, cut) %>%
  mutate(from = 1
         , to = type + 1
         , transition = case_when(
           type == 1 ~ "1->2",
           type == 2 ~ "1->3",
           .default = "err"
         )) %>%
  rename(tstart = t, tstop = time) %>%
  filter(status == 1)

head(sim_df_stepone)

# second step


# take all relevant subjects, i.e. subjects with transition 1->2, i.e. currently in state 2
data_steptwo <- sim_df_stepone %>% 
  filter(to == 2, status == 1) %>%
  mutate(x3 = tstop) %>%
  select(id, x1, x2, x3, from, t = tstop)

dim(data_steptwo)
head(data_steptwo)

cut =  seq(0, 3, by = 0.01)

sim_df_steptwo <- sim_pexp_cr(form_steptwo, data_steptwo, cut) %>%
  mutate(from = 2
         , to = type + 2
         , transition = case_when(
           type == 1 ~ "2->3",
           .default = "err"
         )
         , time = time + t
         , hazard2 = 0) %>%
  rename(tstart = t, tstop = time) %>% 
  filter(status == 1)

head(sim_df_steptwo)

sim_df <- rbind(sim_df_stepone, sim_df_steptwo)
head(sim_df)

sim_df <- sim_df %>% add_counterfactual_transitions()
head(sim_df)

table(sim_df$transition)

# go on with pamm tools procedure
cal_sim_df <- as_ped_multistate(
  data       = sim_df,
  formula    = Surv(tstart, tstop, status)~ transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_sim_df)
head(cal_sim_df)

ctrl <- gam.control(trace = TRUE)

# pam_sim <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
#                      + as.factor(transition)
#                      + s(x1, by=as.factor(transition))
#                      + s(x2, by=as.factor(transition))
#                      , data = cal_sim_df
#                      , family=poisson()
#                      , offset=offset
#                      , control = ctrl)
# 
# summary(pam_sim)
# plot(pam_sim, page = 1)
# plot(pam_sim, select = 7, ylim = c(-1,1))
# plot(pam_sim, select = 4, xlim = c(-2,2))

bam_sim <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) 
                     + as.factor(transition)
                     , data = cal_sim_df
                     , family=poisson()
                     , offset=offset
                     , discrete = T
                     , method = "fREML"
                     , control = ctrl)

summary(bam_sim)
# 
# plot(bam_sim, page = 1)
# plot(bam_sim, select = 7, ylim = c(-1,1))
# plot(bam_sim, select = 8, ylim = c(-1,1))
# plot(bam_sim, select = 9, ylim = c(-1,1))

pam_sim <- bam_sim

new_sim_pam <- make_newdata(cal_sim_df
                            , tend = unique(tend)
                            , transition=unique(transition)
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_sim) 

old_groups <- dplyr::groups(new_sim_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_sim_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot hazards
ggplot(test_msm, aes(x=tend, y=cumu_hazard)) + 
  geom_line() + 
  facet_wrap(~transition
             , ncol = 3
             # , scales = "free_y"
             , labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylim(c(0,5)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition
             , ncol = 3
             , scales = "free_y"
             , labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylim(c(0,1)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 




# test_msm_sub <- test_msm %>% filter(transition == "1->2") %>% mutate(tend = round(tend, 1)
#                                                                      , x1 = round(x1, 1))
# ggplot(test_msm_sub, aes(x=tend, y=x1, z=trans_prob)) + 
#   geom_contour() 
# 
# ggplot(test_msm_sub, aes(x=tend, y=x1, z=trans_prob)) +
#   geom_tile(aes(fill=trans_prob)) +
#   scale_fill_gradient2(
#     low  = "steelblue"
#     , high = "firebrick2"
#     , midpoint=0.5)+
#   stat_contour(col="grey30") +
#   xlim(c(0,3))




# ------------------------------------------------------------------------------
#
# MULTI-STATE SIMULATION USING MSTATE::MSSAMPLE()
#
# ------------------------------------------------------------------------------

library(mstate)
?mssample()

# ------------------------------------------------------------------------------
# copy vignette
# ------------------------------------------------------------------------------

# transition matrix for illness-death model
tmat <- trans.illdeath()
print(tmat)

# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (T&G)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
                 dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
                 x1=c(1,1,1,0,0,0),x2=c(6:1))
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 data=tg,keep=c("x1","x2"),trans=tmat)
# expanded covariates
tglong <- expand.covs(tglong,c("x1","x2"))

head(tglong)
dim(tglong)

# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
            data=tglong,method="breslow")

summary(cx)

# try with pammtools
pam_tglong <- tglong %>% 
  mutate(transition = paste0(from, "->", to)) %>%
  rename(tstart = Tstart, tstop = Tstop)

pam_tglong <- pam_tglong %>% add_counterfactual_transitions()

cal_tglong <- as_ped_multistate(
  data       = pam_tglong,
  formula    = Surv(tstart, tstop, status)~ x1.1 + x2.2,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

pam_tglong <- mgcv::gam(ped_status ~ tend*transition
                     + x1.1*transition
                     + x2.2*transition
                     , data = cal_tglong
                     , family=poisson()
                     , offset=offset)

summary(pam_tglong)

new_tglong <- make_newdata(cal_tglong
                            , tend = unique(tend)
                            , transition=unique(transition)) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_tglong) 

old_groups <- dplyr::groups(new_tglong)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_tglong <- new_tglong %>% ungroup(transition)
test_tglong <- group_split(res_tglong) |> 
  map(res_tglong, .f = ~ group_by(.x, transition)) |> 
  map(res_tglong, .f = ~ add_trans_prob(.x)) |>
  map(res_tglong, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot transitions
ggplot(test_tglong, aes(x=tend, y=trans_prob)) + 
  geom_line() +
  facet_wrap(~transition, ncol = 3, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0,10)) +
  ylim(c(0,1)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# new data, to check whether results are the same for transition 1 as T&G
newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
fit <- msfit(cx,newdata,trans=tmat)
tv <- unique(fit$Haz$time)
# mssample
set.seed(1234)
mssamp <- mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100)

# plot results
ggplot(data=mssamp, aes(x=time, y=pstate1)) + 
  geom_line()

ggplot(data=mssamp, aes(x=time, y=pstate2)) + 
  geom_line()

set.seed(1234)
paths(tmat)
mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100,output="path")
set.seed(1234)
mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100,output="data",do.trace=25)


# ------------------------------------------------------------------------------
#
# MULTI-STATE SIMULATION USING https://devinincerti.com/2019/01/01/sim-mstate.html
#
# ------------------------------------------------------------------------------

library(mstate)
library(survival)
data(ebmt4)

tmat <- mstate::transMat(x = list(c(2, 3, 5, 6), 
                                  c(4, 5, 6), 
                                  c(4, 5, 6), 
                                  c(5, 6),
                                  c(),
                                  c()),
                         names = c("Tx", "Rec", "AE", "Rec+AE", 
                                   "Rel", "Death"))
print(tmat)

msebmt <- msprep(data = ebmt4, trans = tmat, 
                 time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"), 
                 keep = c("match", "proph", "year", "agecl"))
msebmt[msebmt$id == 2, ]

# fitting
library(flexsurv)
library(data.table)
n_trans <- max(tmat, na.rm = TRUE)
fits_wei <- vector(mode = "list", length = n_trans)
msebmt$years <- msebmt$time/365.25
for (i in 1:n_trans){
  fits_wei[[i]] <- flexsurvreg(Surv(years, status) ~ match + proph + year + agecl ,
                               data = subset(msebmt, trans == i),
                               dist = "gompertz")
}

pat_2 <- data.frame(msebmt[msebmt$id == 2, 
                           c("match", "proph", "year", "agecl")][1, ])
head(pat_2)

yr_grid <- seq(0, 10, .1)

cumhaz_grid <- seq(0, max(msebmt$years), .01)
cumhaz_pat_2 <- msfit.flexsurvreg(fits_wei, trans = tmat, 
                                  t = cumhaz_grid,
                                  newdata = pat_2,
                                  variance = FALSE)
head(cumhaz_pat_2$Haz)
tail(cumhaz_pat_2$Haz)

# simulate with mstate
sim_stprobs_mstate_2 <- function(n_pats){
  mstate::mssample(Haz = cumhaz_pat_2$Haz, 
                   trans = tmat,
                   tvec = yr_grid,
                   clock = "reset",
                   M = n_pats) 
}


# compare simulations
n_pats <- seq(from = 200, to = 1000, by = 200)
stprobs2 <- vector(mode = "list", length = length(n_pats))
  # mstate
  mstate_stprobs2 <- sim_stprobs_mstate_2(n_pats[1])
  mstate_stprobs2 <- melt(mstate_stprobs2, id.vars = "time",
                          variable.name = "state_id",
                          value.name = "prob")
  mstate_stprobs2$state_id <- sub("pstate", "", mstate_stprobs2$state_id)
  mstate_stprobs2$state_id <- as.numeric(mstate_stprobs2$state_id)
  mstate_stprobs2$lab <- "mstate"
  
  stprobs2[[1]] <- rbind(mstate_stprobs2)
  stprobs2[[1]]$n_pats <- n_pats[1]
  print(1)
  
stprobs2 <- rbindlist(stprobs2)

# ------------------------------------------------------------------------------
# transform data to fit pammtools

head(msebmt)

msebmt_pamm <- msebmt %>% mutate(transition = paste0(from, "->", to)) %>%
  rename(tstart = Tstart,
         tstop = Tstop) %>%
  select(-trans, -years) %>%
  filter(status == 1) %>%
  add_counterfactual_transitions()
head(msebmt_pamm)

cal_msebmt <- as_ped_multistate(
  data       = msebmt_pamm,
  formula    = Surv(tstart, tstop, status) ~ transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_msebmt)
head(cal_msebmt)

ctrl <- gam.control(trace = TRUE)
pam_msebmt <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) + 
                   as.factor(transition)
                 , data = cal_msebmt
                 , family=poisson()
                 , method="fREML"
                 , discrete = T
                 , offset= offset
                 , control = ctrl)

summary(pam_msebmt)


new_pam <- make_newdata(cal_msebmt
                        , tend = unique(tend)
                        , transition=unique(transition)
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_msebmt) 

old_groups <- dplyr::groups(new_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

ggplot(test_msm, aes(x=tend, y=cumu_hazard)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 3, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 


# ------------------------------------------------------------------------------
#
# MULTI-STATE SIMULATION USING "simLexis":
# https://cran.r-project.org/web/packages/Epi/vignettes/simLexis.pdf
#
# ------------------------------------------------------------------------------
library(Epi)
library( survival )
library( splines )

data(DMlate)
dml <- Lexis( entry = list(Per=dodm, Age=dodm-dobth, DMdur=0 ),
                exit = list(Per=dox),
                exit.status = factor(!is.na(dodth),labels=c("DM","Dead")),
                data = DMlate )
dmi <- cutLexis( dml, cut = dml$doins,
                 pre = "DM",
                 new.state = "Ins",
                 new.scale = "t.Ins",
                 split.states = TRUE )

summary( dmi, timeScales=T )
boxes( dmi, boxpos = list(x=c(20,20,80,80),
                          y=c(80,20,80,20)),
       scale.R = 1000, show.BE = TRUE )

Si <- splitLexis( dmi, seq(0,20,1/4), "DMdur" )
summary( Si )

nk <- 5
ad.kn <- with( subset(Si,lex.Xst=="Dead"),
               quantile( Age+lex.dur , probs=(1:nk-0.5)/nk ) ) 
ai.kn <- with( subset(Si,lex.Xst=="Ins" & lex.Cst!=lex.Xst ),
               quantile( Age+lex.dur , probs=(1:nk-0.5)/nk ) )
di.kn <- with( subset(Si,lex.Xst=="Ins" & lex.Cst!=lex.Xst ),
               c(0,quantile( DMdur+lex.dur, probs=(1:(nk-1))/nk ) ))
dd.kn <- with( subset(Si,lex.Xst=="Dead"),
               c(0,quantile( DMdur+lex.dur, probs=(1:(nk-1))/nk ) ))
ti.kn <- with( subset(Si,lex.Xst=="Dead(Ins)"),
               c(0,quantile( t.Ins+lex.dur, probs=(1:(nk-1))/nk ) ))

DM.Ins <- glm( (lex.Xst=="Ins") ~ Ns( Age , knots=ai.kn ) +
                 Ns( DMdur, knots=di.kn ) +
                 I(Per-2000) + sex,
               family=poisson, offset=log(lex.dur),
               data = subset(Si,lex.Cst=="DM") )
DM.Dead <- glm.Lexis( Si, from = "DM", to = "Dead",
                      formula = ~ Ns( Age , knots=ad.kn ) +
                        Ns( DMdur, knots=dd.kn ) +
                        I(Per-2000) + sex )
Ins.Dead <- glm.Lexis( Si, from = "Ins",
                       formula = ~ Ns( Age , knots=ad.kn ) +
                         Ns( DMdur, knots=dd.kn ) +
                         Ns( t.Ins, knots=ti.kn ) +
                         I(Per-2000) + sex )


Tr <- list( "DM" = list( "Ins" = DM.Ins,
                         "Dead" = DM.Dead ),
            "Ins" = list( "Dead(Ins)" = Ins.Dead ) )
Cox.Dead <- coxph( Surv( DMdur, DMdur+lex.dur,
                         lex.Xst %in% c("Dead(Ins)","Dead")) ~
                     Ns( Age-DMdur, knots=ad.kn ) +
                     I(lex.Cst=="Ins") +
                     I(Per-2000) + sex,
                   data = Si )
round( ci.exp( Cox.Dead ), 3 )

Tr.c <- list( "DM" = list( "Ins" = Tr$DM$Ins,
                           "Dead" = Cox.Dead ),
              "Ins" = list( "Dead(Ins)" = Cox.Dead ) )

set.seed( 52381764 )
Nsim <- 5000
ini <- Si[NULL,1:9]
ini[1:2,"lex.id"] <- 1:2
ini[1:2,"lex.Cst"] <- "DM"
ini[1:2,"Per"] <- 1995
ini[1:2,"Age"] <- 60
ini[1:2,"DMdur"] <- 5
ini[1:2,"sex"] <- c("M","F")

head(Si)

system.time( simC <- simLexis( Tr.c,
                               ini,
                               t.range = 12,
                               N = Nsim ) )
head(simC)

dim(simC)

# follow the code of the simLexis vignette, count numbers in each state per time
CoxM <- pState( nState( subset(simC,sex=="M"),
                        at=seq(0,11,0.2),
                        from=60,
                        time.scale="Age" ),
                perm=c(1,2,4,3) )
plot( CoxM, border="black", col="transparent", lwd=3 )

# build pam data set on initial data set
pam_si <- Si %>% 
  mutate(id = lex.id
         , tstart = Per
         , tstop = tstart + lex.dur
         , from = case_when(
           lex.Cst == "DM" ~ 1,
           lex.Cst == "Ins" ~ 2,
           lex.Cst == "Dead" | lex.Cst == "Dead(Ins)" ~ 3,
           .default = 4
         )
         , to = case_when(
           lex.Xst == "DM" ~ 1,
           lex.Xst == "Ins" ~ 2,
           lex.Xst == "Dead" | lex.Xst == "Dead(Ins)" ~ 3,
           .default = 4
         )
         , status = case_when(
           from != to ~ 1,
           .default = 0
         )
         , transition = paste0(from, "->", to)) %>%
  filter(status == 1) %>%
  select(id, tstart, tstop, from, to, status, transition) %>%
  add_counterfactual_transitions() 
head(pam_si)
dim(pam_si)



Si %>% filter(lex.id == 6)

# estimate pam on simulated data
pam_simC <- simC %>% 
  mutate(id = lex.id
         , tstart = Per - 1995
         , tstop = tstart + lex.dur
         , from = case_when(
           lex.Cst == "DM" ~ 1,
           lex.Cst == "Ins" ~ 2,
           lex.Cst == "Dead" | lex.Cst == "Dead(Ins)" ~ 3,
           .default = 4
         )
         , to = case_when(
           lex.Xst == "DM" ~ 1,
           lex.Xst == "Ins" ~ 2,
           lex.Xst == "Dead" | lex.Xst == "Dead(Ins)" ~ 3,
           .default = 4
         )
         , status = case_when(
           from != to ~ 1,
           .default = 0
         )
         , transition = paste0(from, "->", to)) %>%
  filter(status == 1) %>%
  select(id, tstart, tstop, from, to, status, transition) %>%
  add_counterfactual_transitions()
head(pam_simC)
dim(pam_simC)

cal_simC <- as_ped_multistate(
  data       = pam_simC,
  formula    = Surv(tstart, tstop, status) ~ transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_simC)
head(cal_simC)

ctrl <- gam.control(trace = TRUE)
pam_simC <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) + 
                          as.factor(transition)
                        , data = cal_simC
                        , family=poisson()
                        , method="fREML"
                        , discrete = T
                        , offset= offset
                        , control = ctrl)

summary(pam_simC)


new_pam <- make_newdata(cal_simC
                        , tend = unique(tend)
                        , transition=unique(transition)
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_simC) 

old_groups <- dplyr::groups(new_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

ggplot(test_msm, aes(x=tend, y=cumu_hazard)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 3, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 



#-------------------------------------------------------------------------------
#
# TRY EFFECT OF TIME
#
#-------------------------------------------------------------------------------




# set number of observations/subjects
n <- 250
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
  as_tibble()
# the formula which specifies how covariates affet the hazard rate
f0 <- function(t) {
  dgamma(t, 8, 2) *6
}
form <- ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)

# understand what is done in the code
# lines 101 ff.
f1 <- formula(Formula(form), rhs = 1)
if (length(Formula(form))[2] > 1) {
  f2  <- formula(Formula(form), rhs = 2)
} else {
  f2 <- NULL
}

data <- df %>%
  mutate(
    id     = row_number(),
    time   = max(cut),
    status = 1)

ped  <- split_data(
  formula = Surv(time, status)~.,
  data    = select_if(data, is_atomic),
  cut     = cut,
  id      = "id")

ped <- ped %>%
  rename("t" = "tstart") %>%
  mutate(rate = exp(f_eval(f1, .)))

head(ped)

# time is used as a variable and set to tstart in split_data step. 
# idea: use this sampler to sample third transition

f0 <- function(t) {
  dgamma(t, 4, 2) *6
}

# form <- formula(  ~ 0.5 + 0.25*x1**3| f0(t))
form <- formula(  ~ 0.3 | 0.6 )

n = 100
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  x2 = runif(n, 0, 6),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)

sim_df <- sim_pexp_cr(form, data, cut) %>%
  mutate(from = 1
         , to = type + 1
         , transition = case_when(
           type == 1 ~ "1->2",
           type == 2 ~ "1->3",
           .default = "err"
         )) %>%
  rename(tstart = t, tstop = time) 

sim_df_2 <- sim_df %>%
  filter(to == 2) %>%
  rename(t = tstop) %>%
  mutate(from = 2
         , time = max(cut)
         , status = 1) %>%
  select(id, x1, x2, t, from, time, status)

# use again code for transition sampling
ped  <- split_data(
  formula = Surv(t, time, status)~.,
  data    = select_if(sim_df_2, is_atomic),
  cut     = cut,
  id      = "id")
  # filter(tstart >= t) %>% #exclude all time steps before transition is possible

ped <- ped %>%
  rename("tfirst" = "t"
         , "t" = "tstart") %>%
  mutate(rate = exp(f_eval(f1, .)))

ped <- ped %>%
  rename("t" = "tstart") %>%
  mutate(rate = exp(f_eval(f1, .)))

head(ped)

sim_df_step2 <- ped %>%
  group_by(id) %>%
  summarize(time = rpexp(rate = .data$rate, t = .data$t)) %>%
  mutate(
    status = 1L * (.data$time <= max(cut)),
    time   = pmin(.data$time, max(cut)))

# problem: t

head(sim_df_step2)


set.seed(24032018)
sim_df <- sim_pexp(form, df, 1:10)
head(sim_df)
plot(survfit(Surv(time, status)~1, data = sim_df ))

plot(1:10, f0(1:10))
