
# ---------------------------------------------------------------------------- #
# LIBRARIES
# ---------------------------------------------------------------------------- #
# Beyersmann
library(kmi)
library(mvna)
library(etm)

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

# ---------------------------------------------------------------------------- #
# DATAPREP
# ---------------------------------------------------------------------------- #

# load data used in BEYERSMANN
data(icu.pneu)

# BEYERSMANN
my.icu.pneu <- icu.pneu %>% 
  group_by(id) %>%
  mutate(event = event
         , to2  = ifelse(status == 1, as.character(event), ifelse((n() == 1 & status == 0) | (pneu == 1 & status == 0), "cens", "1"))
         , to   = ifelse((n() == 1 & status == 0) | (pneu == 1 & status == 0), "cens", as.character(ifelse(status == 1, min(event, 2), 1) ))
         , transition = paste0(as.character(pneu), "->", as.character(to))
  ) %>%
  rename(entry = start
         , exit = stop
         , from = pneu
  ) %>%
  select(id
         , entry
         , exit
         , from
         , to
         , to2
         , transition
         , age
         , sex
         # , adm.cens.exit
         # , status
         # , event
  )

# PAMMTOOLS
pamm.icu.pneu <- my.icu.pneu %>%
  mutate(status = ifelse(grepl("cens", transition), 0, 1)) %>%
  rename(tstart = entry
         , tstop = exit) %>%
  select(id
         , tstart
         , tstop
         , from
         , to
         , transition
         , status
         , age
         , sex)


my.pamm.icu.pneu_test <- pamm.icu.pneu %>% add_counterfactual_transitions()

cal_icu.pneu <- as_ped_multistate(
  data       = my.pamm.icu.pneu_test,
  formula    = Surv(tstart, tstop, status)~ age + sex + transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(pamm.icu.pneu)
dim(my.pamm.icu.pneu_test)
dim(cal_icu.pneu)

# ---------------------------------------------------------------------------- #
# Predict hazards and use add_trans_prob to calculate transition prob
# ---------------------------------------------------------------------------- #
# as factor: transition in cal_icu.pneu
# ggf in as_ped gleich als faktorvariable ausgeben.
pam <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) + as.factor(transition) + s(age) + sex, data = cal_icu.pneu, family=poisson(), offset=offset)

summary(pam)

plot(pam, xlim = c(0,100), ylim = c(-20, 20), page=1)
plot.gam(pam, select = 4, ylim=c(-1,1))

meeting_test <- make_newdata(cal_icu.pneu, tend = unique(tend), transition=unique(transition), age = quantile(age, probs=c(0.05, 0.95))) %>% 
  group_by(transition, age) %>% 
  add_cumu_hazard(pam)

# makenewdata anstelle von ped_info (setzt covariablen die nicht spezifiziert werden auf mittelwerte)


old_groups <- dplyr::groups(meeting_test)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- meeting_test %>% ungroup(transition)
test <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()
# # old code without grouping
# test <- meeting_test %>% add_trans_prob()
# ggplot(test, aes(x=tend, y=trans_prob)) + geom_line(aes(col=transition)) + xlim(c(0, 100)) + ylim(c(0, 1))

# plot transitions
ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=as.factor(age))) + 
  facet_wrap(~transition, ncol = 1) +
  xlim(c(0, 100)) + 
  ylim(c(0,1))

# ---------------------------------------------------------------------------- #
# RECALCULATE BEYERSMANN
# ---------------------------------------------------------------------------- #
# build transition matrix
tra.idm <- matrix(FALSE, 3, 3, dimnames = list(c(0, 1, 2), c(0, 1, 2)))
tra.idm[1, 2:3] <- TRUE
tra.idm[2, 3] <- TRUE

# print result of transition matrix
print(tra.idm)

mvna.idm <- mvna(my.icu.pneu, c("0", "1", "2"), tra.idm, "cens")

# re-create plots of book
if (require(lattice)){
  xyplot(mvna.idm, xlim=c(0, 100))
}
# same graphs as in Beyersmann --> correct

plot(mvna.idm$'0 1'$time, mvna.idm$'0 1'$na, xlim=c(0,100))

# estimate transition probabilites
etm.idm <- etm(my.icu.pneu, c("0", "1", "2"), tra.idm, "cens", s = 0)

par(mfrow=c(1,3))
plot(etm.idm, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.idm, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.idm, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")
par(mfrow=c(1,1))
