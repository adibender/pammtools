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


my.pamm.icu.pneu <- pamm.icu.pneu %>% add_counterfactual_transitions()

# # testing:
# pamm.icu.pneu %>% filter(to == "cens") # -> id = 30236
# 
# pamm.icu.pneu %>% filter(id == 30236)
# my.pamm.icu.pneu %>% filter(id == 30236)

cal_icu.pneu <- as_ped_multistate(
  data       = my.pamm.icu.pneu,
  formula    = Surv(tstart, tstop, status)~ age + sex + transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

# ---------------------------------------------------------------------------- #
# DATA CHECKS
# ---------------------------------------------------------------------------- #
# check if numbers are equal to Beyersmann page 6
table(my.icu.pneu$to, my.icu.pneu$to2)

# check censored data
dim(my.icu.pneu %>% 
      filter(to == "cens"))[1]
# = 21 -> correct

# check number of hospital acquired pneunomia
dim(my.icu.pneu %>% 
      filter(from == "1"))[1]
dim(my.icu.pneu %>% 
      filter(from == "1" & to != "cens"))[1]
# =108 without considering censoring, = 103 considering censoring -> correct
hosp.acq.pneu <- my.icu.pneu %>% 
  filter(from == "1")
table(hosp.acq.pneu$to2)
# 2: discharged alive = 21, 3: died = 82, censored = 5
non.acq.pneu <- my.icu.pneu %>% 
  filter(from == "0")
table(non.acq.pneu$to2)
# 2: discharged alive = 126, 3: died = 1063, censored = 16

# --> numbers are correct

# compare to Bayersmann p. 184
my.icu.pneu %>% filter(id %in% c(405, 410, 3163, 17743, 17776))
icu.pneu %>% filter(id %in% c(405, 410, 3163, 17743, 17776))

# compare to Bayersmann p. 184
my.icu.pneu %>% filter(id %in% c(405, 410, 3163, 17743, 17776))
icu.pneu %>% filter(id %in% c(405, 410, 3163, 17743, 17776))

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

plot(etm.idm, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.idm, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# ---------------------------------------------------------------------------- #
# calculate the transition probabilities using descrete time data
# ---------------------------------------------------------------------------- #
# calculate Aalen Johannsen
# example \delta A(1)^{01} = sum(ped_status) / count(*), data set filtered on transition "0->1"
test.pamm <- cal_icu.pneu %>% 
  # filter(transition == "0->2") %>%
  group_by(tend, transition) %>%
  summarize(n.risk.pamm = n()
            , n.event.pamm = sum(ped_status == 1)) %>%
  mutate(A_lj= ifelse(is.na(n.event.pamm / n.risk.pamm), 0, n.event.pamm / n.risk.pamm)) %>%
  rename(time = tend) %>%
  select(time, transition, A_lj, n.event.pamm, n.risk.pamm)
head(test.pamm)


t_k <- unique(test.pamm$time)
k <- length(t_k)
A <- array(rep(0, 3*3*k), dim=c(3,3,k))
cum_A <- array(rep(0, 3*3*k), dim=c(3,3,k))

M01 <- matrix(0, nrow = 3, ncol =3)
M01[1,1] <- -1
M01[1,2] <- 1
M02 <- matrix(0, nrow = 3, ncol =3)
M02[1,1] <- -1
M02[1,3] <- 1
M12 <- matrix(0, nrow = 3, ncol =3)
M12[2,2] <- -1
M12[2,3] <- 1

for (iter in 1:k){
  i01 <- which(test.pamm$time == t_k[iter] & test.pamm$transition == "0->1")
  i02 <- which(test.pamm$time == t_k[iter] & test.pamm$transition == "0->2")
  i12 <- which(test.pamm$time == t_k[iter] & test.pamm$transition == "1->2")
  
  alpha_01 <- if(length(i01) == 0) 0 else test.pamm[i01, 3]
  alpha_02 <- if(length(i02) == 0) 0 else test.pamm[i02, 3]
  alpha_12 <- if(length(i12) == 0) 0 else test.pamm[i12, 3]
  
  I <- diag(3)
  A[,,iter] <- I + M01 * as.numeric(alpha_01) +
    M02 * as.numeric(alpha_02) +
    M12 * as.numeric(alpha_12)
  
  # build Aalen Johannsen  
  if (iter == 1) {
    cum_A[,,iter] = A[,,iter]
  } else {
    cum_A[,,iter] = round(cum_A[,,iter-1] %*% A[,,iter],10) #use matrix multiplikation
  }
}

print(A)
print(cum_A)

A01 <- matrix(0, nrow = k, ncol =2)
A01[,2] <- cum_A[1,2,]
A01[,1] <- t_k
df.A01 <- as.data.frame(A01)
colnames(df.A01) <- c("time", "P.A")
print(df.A01)

A02 <- matrix(0, nrow = k, ncol =2)
A02[,2] <- cum_A[1,3,]
A02[,1] <- t_k
df.A02 <- as.data.frame(A02)
colnames(df.A02) <- c("time", "P.A")
print(df.A02)

A12 <- matrix(0, nrow = k, ncol =2)
A12[,2] <- cum_A[2,3,]
A12[,1] <- t_k
df.A12 <- as.data.frame(A12)
colnames(df.A12) <- c("time", "P.A")
print(df.A12)


# compare data frames from example with manually calculated AJ-Estimator
etm.idm.01 <- summary(etm.idm)$"0 1" %>% select(time, P, n.risk, n.event)
etm.idm.02 <- summary(etm.idm)$"0 2" %>% select(time, P, n.risk, n.event)
etm.idm.12 <- summary(etm.idm)$"1 2" %>% select(time, P, n.risk, n.event)

print(etm.idm.01)

comp.etm <- left_join(etm.idm.01, test.pamm, by ="time") %>%
  filter(transition == "0->1") %>%
  select(time, P, n.risk, n.risk.pamm, n.event, n.event.pamm)

print(comp.etm)


# compare probabilities
# 0->1
comp.etm <- left_join(etm.idm.01, df.A01, by ="time") %>%
  select(time, P, P.A, n.risk, n.event)
print(comp.etm)
# correct until time == 5 -> fixed, transition probabilities 1->2 have been wrong

# 0->2 
comp.etm <- left_join(etm.idm.02, df.A02, by ="time") %>%
  select(time, P, P.A, n.risk, n.event)
print(comp.etm)
# minor differences, probably due to different values for 0->1 -> fixed


# 1->2 
comp.etm <- left_join(etm.idm.12, df.A12, by ="time") %>%
  select(time, P, P.A, n.risk, n.event)
print(comp.etm)


# compare graphically with plots
par(mfrow=c(3,2))
# add transition 0->1
plot(df.A01, xlim = c(0,100), ylim = c(0, 0.1), type = "l", lwd = 2)
plot(etm.idm, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# add transition 0->2
plot(df.A02, ylim = c(0, 1), type = "l", lwd = 2, xlim = c(0, 100))
plot(etm.idm, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# add transition 1->2
plot(df.A12, ylim = c(0, 1), type = "l", lwd = 2, xlim = c(0, 100))
plot(etm.idm, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

par(mfrow=c(1,1))


# ---------------------------------------------------------------------------- #
# calculate the transition probabilities with pammtools
# ---------------------------------------------------------------------------- #

# as factor: transition in cal_icu.pneu
# ggf in as_ped gleich als faktorvariable ausgeben.
pam <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) + as.factor(transition) + s(age) + sex, data = cal_icu.pneu, family=poisson(), offset=offset)

summary(pam)

plot(pam, xlim = c(0,100), ylim = c(-20, 20), page=1)

meeting_test <- make_newdata(cal_icu.pneu, tend = unique(tend), transition=unique(transition)) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam) %>%
  mutate(delta_cumu_hazard = cumu_hazard - ifelse(is.na(lag(cumu_hazard)), 0, lag(cumu_hazard)),
         time = tend)


# plot cumu hazards
ggplot(meeting_test, aes(x=tend, y=cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 100)) + ylim(c(0, 20))
# makenewdata anstelle von ped_info (setzt covariablen die nicht spezifiziert werden auf mittelwerte)

# plot delta hazards
ggplot(meeting_test, aes(x=tend, y=delta_cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 100)) + ylim(c(0, 0.20))
# transition 0->2 look like containing negative values
# cumu hazards are relatively big -> high variation in delta hazard

# ---------------------------------------------------------------------------- #
# comparison cumu hazards pammtools vs beyersmann
# ---------------------------------------------------------------------------- #
# compare cumu hazards from beyersmann with pammtools. If correct, then bug is in calculation of transition prob.
meeting_test.01 <- meeting_test %>% filter(transition == "0->1")
plot.01 <- ggplot(meeting_test.01, aes(x=tend, y=cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 100)) + ylim(c(0, 0.5))
combined.plot.01 <- plot.01 + geom_line(data = mvna.idm$'0 1', aes(x=time, y=na), color="blue")
print(combined.plot.01)

# do for rest of the transitions
# transition: 0->2
meeting_test.02 <- meeting_test %>% filter(transition == "0->2")
plot.02 <- ggplot(meeting_test.02, aes(x=tend, y=cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 100)) + ylim(c(0, 8))
combined.plot.02 <- plot.02 + geom_line(data = mvna.idm$'0 2', aes(x=time, y=na), color="blue")
print(combined.plot.02)

# transition: 1->2
meeting_test.12 <- meeting_test %>% filter(transition == "1->2")
plot.12 <- ggplot(meeting_test.12, aes(x=tend, y=cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 100)) + ylim(c(0, 8))
combined.plot.12 <- plot.12 + geom_line(data = mvna.idm$'1 2', aes(x=time, y=na), color="blue")
print(combined.plot.12)

library(patchwork)
combined.plot <- combined.plot.01 + combined.plot.02 + combined.plot.12 + plot_layout(guides = "collect")
plot(combined.plot)
# graphical debugging result: no huge difference, looks ok

# analyze tables
joined.meeting_test.01 <- meeting_test.01 %>% 
  left_join(mvna.idm$'0 1', by = c("tend" = "time")) %>%
  mutate(diff = na - cumu_hazard
         , lag_na = na - lag(na)
         , lag_ch = cumu_hazard - lag(cumu_hazard)
         , lag_diff = lag_na - lag_ch) %>%
  select(tend
         , transition
         , na
         , cumu_hazard
         , diff
         , lag_na
         , lag_ch
         , lag_diff)
View(joined.meeting_test.01)

summary(joined.meeting_test.01)

joined.meeting_test.02 <- meeting_test.02 %>% 
  left_join(mvna.idm$'0 2', by = c("tend" = "time")) %>%
  mutate(diff = na - cumu_hazard
         , lag_na = na - lag(na)
         , lag_ch = cumu_hazard - lag(cumu_hazard)
         , lag_diff = lag_na - lag_ch) %>%
  select(tend
         , transition
         , na
         , cumu_hazard
         , diff
         , lag_na
         , lag_ch
         , lag_diff)

summary(joined.meeting_test.02)

joined.meeting_test.12 <- meeting_test.12 %>% 
  filter(tend < 100) %>%
  left_join(mvna.idm$'1 2', by = c("tend" = "time")) %>%
  mutate(diff = na - cumu_hazard
         , lag_na = na - lag(na)
         , lag_ch = cumu_hazard - lag(cumu_hazard)
         , lag_diff = lag_na - lag_ch) %>%
  select(tend
         , transition
         , na
         , cumu_hazard
         , diff
         , lag_na
         , lag_ch
         , lag_diff)

summary(joined.meeting_test.12)

# CONCLUSIO:
# Wrong transition probabilities do not arise from wrong cumu hazards.
# Bug must be within calculation of Aalen Johannsen.

# ---------------------------------------------------------------------------- #
# calculate transition probabilities
# ---------------------------------------------------------------------------- #
meeting_test.short <- meeting_test %>% filter(tend <= 100)

t_k <- unique(meeting_test$time)
k <- length(t_k)
A <- array(rep(0, 3*3*k), dim=c(3,3,k))
cum_A <- array(rep(0, 3*3*k), dim=c(3,3,k))

M01 <- matrix(0, nrow = 3, ncol =3)
M01[1,1] <- -1
M01[1,2] <- 1
M02 <- matrix(0, nrow = 3, ncol =3)
M02[1,1] <- -1
M02[1,3] <- 1
M12 <- matrix(0, nrow = 3, ncol =3)
M12[2,2] <- -1
M12[2,3] <- 1

for (iter in 1:k){
  i01 <- which(meeting_test$tend == t_k[iter] & meeting_test$transition == "0->1")
  i02 <- which(meeting_test$tend == t_k[iter] & meeting_test$transition == "0->2")
  i12 <- which(meeting_test$tend == t_k[iter] & meeting_test$transition == "1->2")
  
  # calculate delta A via difference in cumulative hazards.
  alpha_01 <- if(length(i01) == 0) 0 else meeting_test$delta_cumu_hazard[i01]
  alpha_02 <- if(length(i02) == 0) 0 else meeting_test$delta_cumu_hazard[i02]
  alpha_12 <- if(length(i12) == 0) 0 else meeting_test$delta_cumu_hazard[i12]
  
  I <- diag(3)
  A[,,iter] <- I + M01 * as.numeric(alpha_01) +
    M02 * as.numeric(alpha_02) +
    M12 * as.numeric(alpha_12)
  
  # build Aalen Johannsen  
  if (iter == 1) {
    cum_A[,,iter] = A[,,iter]
  } else {
    cum_A[,,iter] = round(cum_A[,,iter-1] %*% A[,,iter],10) #use matrix multiplikation
  }
}


A01 <- matrix(0, nrow = k, ncol =2)
A01[,2] <- cum_A[1,2,]
A01[,1] <- t_k
df.A01 <- as.data.frame(A01)
colnames(df.A01) <- c("time", "P.A")
print(df.A01)

A02 <- matrix(0, nrow = k, ncol =2)
A02[,2] <- cum_A[1,3,]
A02[,1] <- t_k
df.A02 <- as.data.frame(A02)
colnames(df.A02) <- c("time", "P.A")
print(df.A02)

A12 <- matrix(0, nrow = k, ncol =2)
A12[,2] <- cum_A[2,3,]
A12[,1] <- t_k
df.A12 <- as.data.frame(A12)
colnames(df.A12) <- c("time", "P.A")
print(df.A12)


# compare graphically with plots
par(mfrow=c(3,2))
# add transition 0->1
plot(df.A01, xlim = c(0,100), ylim = c(0, 0.1), type = "l", lwd = 2)
plot(etm.idm, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# add transition 0->2
plot(df.A02, ylim = c(0, 1), type = "l", lwd = 2, xlim = c(0, 100))
plot(etm.idm, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# add transition 1->2
plot(df.A12, ylim = c(0, 1), type = "l", lwd = 2, xlim = c(0, 100))
plot(etm.idm, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

par(mfrow=c(1,1))


# ---------------------------------------------------------------------------- #
# calculate the transition probabilities with pammtools, using hazards
# ---------------------------------------------------------------------------- #

# ped.test <- cal_icu.pneu %>%
#   mutate(ped_status = case_when(
#     transition == "0->1" & ped_status != 0 ~ 1
#     , transition == "0->2" & ped_status != 0 ~ 2
#     , transition == "1->2" & ped_status != 0 ~ 3
#     , TRUE ~ 0
#   ))

ped.test01 <- cal_icu.pneu %>% filter(transition == "0->1") %>% select(-transition) %>% mutate(interval = factor(interval, levels = unique(interval)))
ped.test02 <- cal_icu.pneu %>% filter(transition == "0->2") %>% select(-transition) %>% mutate(interval = factor(interval, levels = unique(interval)))
ped.test12 <- cal_icu.pneu %>% filter(transition == "1->2") %>% select(-transition) %>% mutate(interval = factor(interval, levels = unique(interval)))

ped.test01 <- cal_icu.pneu %>% filter(transition == "0->1")  
ped.test02 <- cal_icu.pneu %>% filter(transition == "0->2") 
ped.test12 <- cal_icu.pneu %>% filter(transition == "1->2") 

print(ped.test02 %>% filter(ped_status == 1))
print(ped.test12 %>% filter(ped_status == 1))


# calculate transition probabilities by predicting hazards
pam01 <- mgcv::gam(ped_status ~ s(tend), data = ped.test01, family=poisson(), offset=offset)
pam02 <- mgcv::gam(ped_status ~ s(tend), data = ped.test02, family=poisson(), offset=offset)
pam12 <- mgcv::gam(ped_status ~ s(tend), data = ped.test12, family=poisson(), offset=offset)

output01 <- ped_info(ped.test01) %>% add_cumu_hazard(pam01)
output02 <- ped_info(ped.test02) %>% add_cumu_hazard(pam02)
output12 <- ped_info(ped.test12) %>% add_cumu_hazard(pam12)


output01 <- ped_info(ped.test01) %>% add_cumu_hazard(pam01)
output02 <- ped_info(ped.test02) %>% add_cumu_hazard(pam02)
output12 <- ped_info(ped.test12) %>% add_cumu_hazard(pam12)
# warning because unique intervals are not the same for ped.test12 and output12

# calculate difference in hazards, i.e. \delta A
for (iter in 1:k){
  i01 <- which(output01$tend == t_k[iter])
  i02 <- which(output02$tend == t_k[iter])
  i12 <- which(output12$tend == t_k[iter])
  
  # calculate delta A via difference in cumulative hazards.
  alpha_01 <- if(length(i01) == 0) 0 else output01[i01, 7] - ifelse(is.na(lag(output01[i01, 7])),0, lag(output01[i01, 7]))
  alpha_02 <- if(length(i02) == 0) 0 else output02[i02, 7] - ifelse(is.na(lag(output02[i02, 7])),0, lag(output02[i02, 7]))
  alpha_12 <- if(length(i12) == 0) 0 else output12[i12, 7] - ifelse(is.na(lag(output12[i12, 7])),0, lag(output12[i12, 7]))
  
  print(alpha_02)
  
  I <- diag(3)
  A[,,iter] <- I + M01 * as.numeric(alpha_01) +
    M02 * as.numeric(alpha_02) +
    M12 * as.numeric(alpha_12)
  
  # build Aalen Johannsen  
  if (iter == 1) {
    cum_A[,,iter] = A[,,iter]
  } else {
    cum_A[,,iter] = round(cum_A[,,iter-1] %*% A[,,iter],10) #use matrix multiplikation
  }
}


A01 <- matrix(0, nrow = k, ncol =2)
A01[,2] <- cum_A[1,2,]
A01[,1] <- t_k
df.A01 <- as.data.frame(A01)
colnames(df.A01) <- c("time", "P.A")
print(df.A01)

A02 <- matrix(0, nrow = k, ncol =2)
A02[,2] <- cum_A[1,3,]
A02[,1] <- t_k
df.A02 <- as.data.frame(A02)
colnames(df.A02) <- c("time", "P.A")
print(df.A02)

A12 <- matrix(0, nrow = k, ncol =2)
A12[,2] <- cum_A[2,3,]
A12[,1] <- t_k
df.A12 <- as.data.frame(A12)
colnames(df.A12) <- c("time", "P.A")
print(df.A12)

# compare graphically with plots
par(mfrow=c(3,2))
# add transition 0->1
plot(df.A01, ylim = c(0, 0.1), xlim = c(0,100), type = "l", lwd = 2)
plot(etm.idm, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# add transition 0->2
plot(df.A02, ylim = c(0, 1), type = "l", lwd = 2, xlim = c(0, 100))
plot(etm.idm, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

# add transition 1->2
plot(df.A12, ylim = c(0, 1), type = "l", lwd = 2, xlim = c(0, 100))
plot(etm.idm, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")

par(mfrow=c(1,1))

# transition probabilities 0->2 and 1->2 look alike.
# check if the values are the same
comp.02.12 <- left_join(df.A02, df.A12, by = "time") %>%
  mutate(diff = round(P.A.x - P.A.y, 2))

print(comp.02.12)
# values are indeed equal.

# Debug transitions 0->2, 1->2 for iter = 1
iter <- 1

# check if row ids are different
i01 <- which(output01$tend == t_k[iter])
i02 <- which(output02$tend == t_k[iter])
i12 <- which(output12$tend == t_k[iter])

print(paste0(i01, ", ", i02, ", ", i12))
head(output01)
head(output02)
head(output12)
# => same ids are rectified. Tables contain different content
# => used wrong names

# Debug warning message about different time intervals
original_intervals <- unique(model.frame(ped.test12)[["tend"]])
prediction_intervals <- unique(output12[["tend"]])

print(paste0(length(original_intervals), ", ", length(prediction_intervals)))
# different lengths of original and prediction intervals




# calculate delta A via difference in cumulative hazards.
alpha_01 <- if(length(i01) == 0) 0 else output01[i01, 7] - ifelse(is.na(lag(output01[i01, 7])),0, lag(output01[i01, 7]))
alpha_02 <- if(length(i02) == 0) 0 else output02[i02, 7] - ifelse(is.na(lag(output02[i02, 7])),0, lag(output02[i02, 7]))
alpha_12 <- if(length(i12) == 0) 0 else output12[i12, 7] - ifelse(is.na(lag(output12[i12, 7])),0, lag(output12[i12, 7]))

print(paste0(alpha_01, ", ", alpha_02, ", ", alpha_12))
# alphas are all the same!
# outputij all have the same entries

# compare graphically with plots
# take a closer look into 0->1
par(mfrow=c(2,1))
# add transition 0->1
plot(df.A01, ylim = c(0, 0.3), xlim = c(0,100), type = "l", lwd = 2)
plot(etm.idm, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.3),
     xlim = c(0, 100), xlab = "Days",
     ci.fun = "cloglog")
par(mfrow=c(1,1))
# over-estimates the transition probability within t = c(5,20)

# ---------------------------------------------------------------------------- #
# calculate the transition probabilities with pammtools, using CIF
# ---------------------------------------------------------------------------- #

head(my.pamm.icu.pneu)
# only contains status %in% c(0, 1),
# status should take the following values
# 0: in case nothing happens / censoring
# 1: 0->1
# 2: 0->2
# 3: 1->2
# with changed status, treat multistate as competing risk

cr.pamm.icu.pneu <- my.pamm.icu.pneu %>%
  filter(transition != "1->2") %>%  
  mutate(status = case_when(
    status == 0 ~ 0
    , status == 1 & transition == "0->1" ~ 1
    , status == 1 & transition == "0->2" ~ 2
    , status == 1 & transition == "1->2" ~ 3
    , TRUE ~ -1
  ))

# if you use status > 1 then invalid status value. 
# information is transported via transition

cal_icu.pneu <- as_ped_multistate(
  data       = my.pamm.icu.pneu,
  formula    = Surv(tstart, tstop, status)~ 1,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

pam <- mgcv::gam(ped_status ~ tend*transition, data = cal_icu.pneu, family=poisson(), offset=offset)
ped <- ped_info(cal_icu.pneu) %>% add_hazard(pam)

ped <- cal_icu.pneu %>% 
  add_hazard(pam) %>%
  filter(transition == "1->2") %>%
  group_by(tend, transition) %>%
  select(-id, -ped_status) %>%
  distinct()

ped <- ped_info(cal_icu.pneu) %>% 
  add_hazard(pam)

ped

int_info <- int_info(cal_icu.pneu)
ped <- cal_icu.pneu %>% 
  left_join(int_info, by = join_by(tstart, tend, interval)) %>%
  add_hazard(pam) %>%
  filter(transition == "0->1") %>%
  group_by(tend, transition) %>%
  select(-id, -ped_status) %>%
  distinct()

print(ped)

plot(ped$hazard)

head(ped)
summary(ped)
summary(pam)

head(cr.pamm.icu.pneu)
unique(cr.pamm.icu.pneu$status)

cr.icu.pneu <- my.icu.pneu %>%
  filter(transition != "1->2") %>% 
  mutate(status = as.numeric(ifelse(to2 == "cens", 0, to2))
         , tstart = entry
         , tstop = exit) %>%
  select(id, tstart, tstop, status, age, sex)

# example from website 
# https://adibender.github.io/pammtools/articles/competing-risks.html
data("fourD", package = "etm")

head(fourD)

unique(fourD$status)
unique(fourD$time)
head(my.pamm.icu.pneu)
cut <- sample(fourD$time, 100)
ped <- fourD %>%
  select(-medication, - treated) %>%
  as_ped(Surv(time, status)~., id = "id", cut = cut, combine = FALSE)
str(ped, 1)


# ---------------------------------------------------------------------------- #
# Debug ped_info
# ---------------------------------------------------------------------------- #

# #code from https://github.com/adibender/pammtools/blob/multi-state/R/interval-information.R
# ped_info.ped <- function(ped) {
#   
#   int_df <- int_info(ped)
#   sdf    <- sample_info(ped)
#   if (is.null(sdf)) {
#     return(int_df)
#   } else {
#     bind_cols(
#       int_df %>% slice(rep(seq_len(nrow(int_df)), times = nrow(sdf))),
#       sdf %>% slice(rep(seq_len(nrow(sdf)), each = nrow(int_df)))) %>%
#       grouped_df(vars = group_vars(sdf))
#   }
#   
# }

int_df <- int_info(cal_icu.pneu)
print(int_df)
# works correctly

sdf <- sample_info(cal_icu.pneu)
print(sdf)
# contains only transition "0->1"

# debug
# # code from https://github.com/adibender/pammtools/blob/multi-state/R/make-newdata.R
# sample_info.data.frame <- function(x) {
# 
# cn  <- colnames(x)
# num <- summarize_if (x, .predicate = is.numeric, ~mean(., na.rm = TRUE))
# fac <- summarize_if (x, .predicate = compose("!", is.numeric), modus)
# 
# nnames <- intersect(names(num), names(fac))
# 
# if (length(nnames) != 0) {
#   suppressMessages(
#     x <- left_join(num, fac) %>% group_by(!!!lapply(nnames, as.name))
#   )
# } else {
#   x <- bind_cols(num, fac)
# }
# 
# return(select(x, one_of(cn)))
# 
# }
# 
# sample_info.ped <- function(x) {
#   # is.grouped_df
#   # remove "noise" information on interval variables
#   grps <- group_vars(x)
#   iv <- attr(x, "intvars")
#   id_var <- attr(x, "id_var")
#   x <- x %>%
#     group_by(!!sym(id_var)) %>%
#     slice(1) %>%
#     ungroup() %>%
#     grouped_df(grps) %>%
#     select(-one_of(iv))
#   if (test_data_frame(x, min.rows = 1, min.cols = 1)) {
#     sample_info.data.frame(x)
#   } else {
#     NULL
#   }
#   
# }

cn  <- colnames(cal_icu.pneu)
num <- summarize_if (cal_icu.pneu, .predicate = is.numeric, ~mean(., na.rm = TRUE))
fac <- summarize_if (cal_icu.pneu, .predicate = compose("!", is.numeric), modus)

nnames <- intersect(names(num), names(fac))

x <- bind_cols(num, fac)

print(x)

select(x, one_of(cn))

x <- cal_icu.pneu
x <- x %>%
 group_by(!!sym(id_var)) %>%
 slice(1) %>%
 ungroup() %>%
 grouped_df(grps) %>%
 select(-one_of(iv))

x <- cal_icu.pneu
x <- x %>%
  group_by(!!sym(id_var), transition) %>%
  slice(1) %>%
  ungroup() %>%
  grouped_df(grps) %>%
  select(-one_of(iv))

unique(x$transition)

# 
# head(ped.test)
# ped <- pamm.icu.pneu %>% as_ped_multistate(
#   data       = my.pamm.icu.pneu,
#   formula    = Surv(tstart, tstop, status)~ 1 ,
#   transition = "transition",
#   id         = "id",
#   censor_code = 0,
#   timescale  = "calendar")
# 
# # try using etm
# my.pamm.pneu <- cal_icu.pneu %>% 
#   rename(entry = tstart, exit = tend) %>% 
#   select(id
#         , entry
#         , exit
#         , from
#         , to
#         , age
#         , sex)
# etm.pamm <- etm(my.pamm.pneu, c("0", "1", "2"), tra.idm, "cens", s = 0)
# plot(etm.pamm, tr.choice = "0 1", conf.int = TRUE,
#      lwd = 2, legend = FALSE, ylim = c(0, 0.1),
#      xlim = c(0, 100), xlab = "Days",
#      ci.fun = "cloglog")
# 
# summary(etm.idm)$"0 1"
# 
# 
# table(cal_icu.pneu$from, cal_icu.pneu$to)
# table(cal_icu.pneu$transition)
# 
# 
# 
# # Simulate multi-state survival data (replace this with your own data)
# set.seed(123)
# n <- 100
# time <- rexp(n, rate = 0.1)
# state <- sample(1:3, n, replace = TRUE)
# 
# ped <- tumor[1:50,] %>% as_ped(Surv(days, status)~ age)
# pam <- mgcv::gam(ped_status ~ s(tend)+age, data = ped, family=poisson(), offset=offset)
# ped_info(ped) %>% add_hazard(pam, type="link")
# ped_info(ped) %>% add_hazard(pam, type = "response")
# ped_info(ped) %>% add_cumu_hazard(pam)


# ------------------------------------------------------------------------------
#
# SIMULATION EXAMPLE sim.msm
# https://rdrr.io/cran/msm/man/simmulti.msm.html
#
# ------------------------------------------------------------------------------

library(msm)
?sim.msm

### Simulate 100 individuals with common observation times
sim.df <- data.frame(subject = rep(1:100, rep(13,100)), time = rep(seq(0, 24, 2), 100))
qmatrix <- rbind(c(-0.11,   0.1,  0.01 ),
                 c(0.05,   -0.15,  0.1 ),
                 c(0.02,   0.07, -0.09))
data <- simmulti.msm(sim.df, qmatrix)
head(data)

data %>% group_by(subject) %>%
  mutate(from = state,
                to = ifelse(is.na(lead(state)), state, lead(state)),
                tstart = time,
                tstop = ifelse(is.na(lead(time)), time, lead(time)),
         transition = paste0(state, "->",ifelse(is.na(lead(state)), state, lead(state)))) %>% 
  filter(from != to)
