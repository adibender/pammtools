# ------------------------------------------------------------------------------
#
# SIMULATION EXAMPLE SECTION 3.2. + 4.2 
# => EXTEND TO MULTISTATE
# => EXTEND TO TIME DEPENDENT THIRD STATE TRANSITION
# => BUT INDEPENDENT OF TRANSITION
#
# CODE INCLUDES:
# 1 COMPARISON OF MVNA AND TWO-STEP SIM_PEXP_CR RESULTS
# 2 COMPARISON OF SIM_PEXP_MSM AND TWO-STEP SIM_PEXP_CR RESULTS
#
# ------------------------------------------------------------------------------

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
library(mstate)

# ---------------------------------------------------------------------------- #
# CODE
#
# 1 COMPARISON OF MVNA AND TWO-STEP SIM_PEXP_CR RESULTS
#
# ---------------------------------------------------------------------------- #

# first part remains equal to above
# here: alpha_01 = 0.3, alpha_02 = 0.6, alpha_12 = 0.5
n = 750

# from section 3.2
event.times <- rexp(n, 0.9)

# time dependent hazard for transition from 1->2, i.e. depending on how long 
# subject stays in state 1, it transitions.
f12 <- function(t) {
  dgamma(t, 4, 2) * 6
}

# meeting 11.06.2024
# Möglichkeit 1 --> funktioniert und ändert die transition probabilities
f12 <- function(t, t0) {
  dgamma(t, 4, 2) *6 - dgamma(t0, 1, 1) *6
}

# meeting 11.06.2024
# # Möglichkeit 2
# f12 <- function(t, t0) {
#   dgamma(t, 4, 2) *6 - 1.5*t0
# }

# alpha_01 / (alpha_01 + alpha_02) = 1/3
f.cause <- rbinom(n, size = 1, prob = 1/3)
f.cause <- ifelse(f.cause == 0, 2, 1)

# censoring times uniformly between 0 and 5
cens.times <- runif(n,0,5)
obs.times <- pmin(event.times, cens.times)
obs.cause <- c(event.times <= cens.times) * f.cause

# now section 4.2
# define illess death
tra <- matrix(FALSE, ncol = 3, nrow = 3)
dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
tra[1, 2:3] <- TRUE
tra[2, 3] <- TRUE

# print transition matrix
tra

# build base data set containing all necessary structure
id <- seq_along(obs.cause)
from <- rep(0, length(obs.cause))
to <- as.factor(ifelse(obs.cause == 0, "cens", obs.cause))
my.data <- data.frame(id, from, to, tstart = 0, time = obs.times, cens.times, x3=0)
head(my.data)
summary(my.data$time)

# second competing risk scenario
my.data.steptwo <- my.data %>% filter(to == 1) %>% mutate(x3 = time) 
head(my.data.steptwo)

# define new sample size based on subjects in risk set, i.e. in state 1 (ill)
n.steptwo <- dim(my.data.steptwo)[1]

# define cutpoints to be able to let pammtools sampler run
# make sure it has same length as censoring times
cut =  seq(0, 5, by = 0.05)

# sample data with competing risk sampler
my.data.steptwo <- sim_pexp_cr(formula( ~ f12(tend, x3)), my.data.steptwo, cut) %>%
  rename(tstart = t, tstop = time) %>%
  mutate(tstart = x3
         , tstop = pmin(x3+tstop, cens.times)
         , from = as.factor(1)
         , to = as.factor(ifelse(tstop >= cens.times, "cens", 2))
  ) %>%
  select(id, from, to, tstart, tstop, cens.times, x3)

head(my.data.steptwo)
summary(my.data.steptwo$tstop)

# rename columns to bind rows
my.data <- my.data %>% rename(tstop = time)

# no need for rbinom since only one transition is possible
my.data.complete <- rbind(my.data, my.data.steptwo) %>% rename(entry = tstart, exit = tstop) %>% arrange(id, entry, exit)
dim(my.data.complete)
head(my.data.complete)

# calculate cause specific hazards
my.nelaal <- mvna(my.data.complete, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal
         , strip=strip.custom(bg="white")
         , ylab="Cumulative Hazard"
         , lwd=2
         , xlim=c(0,5)
         , ylim=c(0,10))
}
# same graphs as in Beyersmann --> correct

plot(my.nelaal$'0 1'$time, my.nelaal$'0 1'$na)

# estimate transition probabilites
etm.my.data<- etm(my.data.complete, c("0", "1", "2"), tra, "cens", s = 0)

par(mfrow = c(1,3))
plot(etm.my.data, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.my.data, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")

plot(etm.my.data, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")
par(mfrow = c(1,1))



# try with pammtools
pamm.my.data <- my.data.complete %>%
  rename(tstart = entry
         , tstop = exit) %>%
  mutate(status = ifelse(to == "cens", 0, 1),
         transition = paste0(from, "->", to),
         dummy_12 = ifelse(from == 1, 1, 0)) %>%
  select(id
         , tstart
         , tstop
         , from
         , to
         , transition
         , status
         , x3
         , dummy_12)

head(pamm.my.data)

pamm.my.data %>% filter(to == "cens" & from == 1) # -> z.B. id = 103

pamm.my.data <- pamm.my.data %>% add_counterfactual_transitions() 
head(pamm.my.data)

pamm.my.data %>% filter(id == 103) # no potential transition from 1->2 after add_counterfactual_transitions()

table(pamm.my.data$from, pamm.my.data$to)

cal.my.data <- as_ped_multistate(
  data       = pamm.my.data,
  formula    = Surv(tstart, tstop, status)~ transition + x3,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar") 

cal.my.data.twostep <- cal.my.data

head(cal.my.data)

ggplot(data=cal.my.data, aes(x=transition, y=x3)) + geom_boxplot()

# anpassung nach Meeting vom 11.06.2024
#
pam <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) + 
                   as.factor(transition)
                 , data = cal.my.data
                 , family=poisson()
                 , offset=offset
                 , method = "fREML"
                 , discrete = T)

# # Modell mit welchem korrekte Daten erzeugt werden
# pam <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) + s(t0, by=)
#                    as.factor(transition)
#                  , data = cal.my.data
#                  , family=poisson()
#                  , offset=offset
#                  , method = "fREML"
#                  , discrete = T)

summary(pam)
plot(pam, page = 1, ylim = c(-4,4), xlim=c(0,3))


new_pam <- make_newdata(cal.my.data
                        , tend = unique(tend)
                        , transition=unique(transition)
                        # , x3 = c(seq(0,3,by=1), mean(x3))
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam) 

head(new_pam)

old_groups <- dplyr::groups(new_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

par(mfrow = c(1,3))
test_msm.01 <- test_msm %>% filter(transition == "0->1")
plot.01 <- ggplot(test_msm.01, aes(x=tend, y=cumu_hazard)) + 
  geom_line(aes(col=transition)) + 
  xlim(c(0, 3)) + 
  ylim(c(0, 1)) + 
  facet_wrap(~ x3)
combined.plot.01 <- plot.01 + geom_line(data = my.nelaal$'0 1', aes(x=time, y=na), color="blue")
print(combined.plot.01)

test_msm.02 <- test_msm %>% filter(transition == "0->2")
plot.02 <- ggplot(test_msm.02, aes(x=tend, y=cumu_hazard)) + 
  geom_line(aes(col=transition)) + 
  xlim(c(0, 3)) + 
  ylim(c(0, 1)) + 
  facet_wrap(~ x3)
combined.plot.02 <- plot.02 + geom_line(data = my.nelaal$'0 2', aes(x=time, y=na), color="blue")
print(combined.plot.02)

test_msm.12 <- test_msm %>% filter(transition == "1->2")
plot.12 <- ggplot(test_msm.12, aes(x=tend, y=cumu_hazard)) + 
  geom_line(aes(col=transition)) + 
  xlim(c(0, 3)) + 
  ylim(c(0, 1)) + 
  facet_wrap(~ x3)
combined.plot.12 <- plot.12 + geom_line(data = my.nelaal$'1 2', aes(x=time, y=na), color="blue")
print(combined.plot.12)
par(mfrow = c(1,1))

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 3, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylim(c(0,0.8)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

par(mfrow = c(1,3))
plot(etm.my.data, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE,  ylim = c(0, 1), xlab = "Days",
     ci.fun = "cloglog")
title("transition: 0->1")
lines(test_msm.01$tend, test_msm.01$trans_prob, lwd=2, col="firebrick2")

plot(etm.my.data, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")
title("transition: 0->2")
lines(test_msm.02$tend, test_msm.02$trans_prob, lwd=2, col="firebrick2")

plot(etm.my.data, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")
title("transition: 1->2")
lines(test_msm.12$tend, test_msm.12$trans_prob, lwd=2, col="firebrick2")
par(mfrow = c(1,1))

# ---------------------------------------------------------------------------- #
# CODE
# 
# 2 COMPARISON OF SIM_PEXP_MSM AND TWO-STEP SIM_PEXP_CR RESULTS
#
# ---------------------------------------------------------------------------- #

t_mat <- matrix(data = NA, nrow = 3, ncol = 3)
t_mat[1,2] <- "log(0.3)"
t_mat[1,3] <- "log(0.6)"
t_mat[2,3] <- "f12(t)"

n = 1000
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  x2 = runif(n, 0, 6),
  from = 1,
  t = 0)

cut =  seq(0, 5, by = 0.05)

msm_tmp <- sim_pexp_msm(
  t_mat = t_mat,
  data = data,
  cut = cut,
  keep_transitions_at_risk = TRUE)
head(msm_tmp)

# include censoring as before
cens_df <- data.frame(cens.times)
cens_df$id <- seq.int(nrow(cens_df))

head(cens_df)
msm_df <- left_join(msm_tmp, cens_df, by = "id") %>%
  filter(cens.times > tstart)%>%
  mutate(status = ifelse(cens.times <= tstop, 0, status) 
         , to = ifelse(cens.times <= tstop, "cens", as.character(to - 1))
         , from = from - 1
         , transition = paste0(from, "->", to))
head(msm_df)
tail(msm_df)



# define illess death for mvna
tra <- matrix(FALSE, ncol = 3, nrow = 3)
dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
tra[1, 2:3] <- TRUE
tra[2, 3] <- TRUE

# print transition matrix
tra

# no need for rbinom since only one transition is possible
my.data.complete <- msm_df %>% rename(entry = tstart, exit = tstop) %>% arrange(id, entry, exit) %>% distinct()
dim(my.data.complete)
head(my.data.complete)

# calculate cause specific hazards
my.nelaal <- mvna(my.data.complete, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal
         , strip=strip.custom(bg="white")
         , ylab="Cumulative Hazard"
         , lwd=2
         , xlim=c(0,5)
         , ylim=c(0,10))
}


etm.my.data<- etm(my.data.complete, c("0", "1", "2"), tra, "cens", s = 0)

par(mfrow = c(1,3))
plot(etm.my.data, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.my.data, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")

plot(etm.my.data, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")
par(mfrow = c(1,1))



# try with pammtools
pamm.my.data <- my.data.complete %>%
  rename(tstart = entry
         , tstop = exit) %>%
  mutate(status = ifelse(to == "cens", 0, 1),
         transition = paste0(from, "->", to)) %>%
  select(id
         , tstart
         , tstop
         , from
         , to
         , transition
         , status)

head(pamm.my.data)

pamm.my.data %>% filter(to == "cens" & from == 1) # -> z.B. id = 103

pamm.my.data <- pamm.my.data %>% add_counterfactual_transitions()
head(pamm.my.data)

pamm.my.data %>% filter(id == 103) # no potential transition from 1->2 after add_counterfactual_transitions()

table(pamm.my.data$from, pamm.my.data$to)

cal.my.data <- as_ped_multistate(
  data       = pamm.my.data,
  formula    = Surv(tstart, tstop, status)~ transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

head(cal.my.data)

pam <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) + 
                   as.factor(transition)
                 , data = cal.my.data
                 , family=poisson()
                 , offset=offset
                 , method = "fREML"
                 , discrete = T)

summary(pam)
plot(pam, page = 1)


new_pam <- make_newdata(cal.my.data
                        , tend = unique(tend)
                        , transition=unique(transition)
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam) 

old_groups <- dplyr::groups(new_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

par(mfrow = c(1,3))
test_msm.01 <- test_msm %>% filter(transition == "0->1")
plot.01 <- ggplot(test_msm.01, aes(x=tend, y=cumu_hazard)) + 
  geom_line(aes(col=transition)) + 
  xlim(c(0, 3)) + 
  ylim(c(0, 1))
combined.plot.01 <- plot.01 + geom_line(data = my.nelaal$'0 1', aes(x=time, y=na), color="blue")
print(combined.plot.01)

test_msm.02 <- test_msm %>% filter(transition == "0->2")
plot.02 <- ggplot(test_msm.02, aes(x=tend, y=cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 3)) + ylim(c(0, 1))
combined.plot.02 <- plot.02 + geom_line(data = my.nelaal$'0 2', aes(x=time, y=na), color="blue")
print(combined.plot.02)

test_msm.12 <- test_msm %>% filter(transition == "1->2")
plot.12 <- ggplot(test_msm.12, aes(x=tend, y=cumu_hazard)) + 
  geom_line(aes(col=transition)) + 
  xlim(c(0, 3)) + 
  ylim(c(0, 1))
combined.plot.12 <- plot.12 + geom_line(data = my.nelaal$'1 2', aes(x=time, y=na), color="blue")
print(combined.plot.12)
par(mfrow = c(1,1))

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 3, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylim(c(0,0.8)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

par(mfrow = c(1,3))
plot(etm.my.data, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE,  ylim = c(0, 1), xlab = "Days",
     ci.fun = "cloglog")
title("transition: 0->1")
lines(test_msm.01$tend, test_msm.01$trans_prob, lwd=2, col="firebrick2")

plot(etm.my.data, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")
title("transition: 0->2")
lines(test_msm.02$tend, test_msm.02$trans_prob, lwd=2, col="firebrick2")

plot(etm.my.data, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlab = "Days",
     ci.fun = "cloglog")
title("transition: 1->2")
lines(test_msm.12$tend, test_msm.12$trans_prob, lwd=2, col="firebrick2")
par(mfrow = c(1,1))

# ANALYZE DIFFERENCES
diff_12 <- cal.my.data %>% filter(transition == "1->2" & ped_status == 1)
summary(diff_12$tstart)

diff_12_twostep <- cal.my.data.twostep %>% filter(transition == "1->2" & ped_status == 1)
summary(diff_12_twostep$tstart)

# comparison of start end end times
msm_df_12 <- msm_df %>% filter(transition == "1->2")
msm_df_01 <- msm_df %>% filter(transition == "0->1")
summary(msm_df_12$tstop)
summary(msm_df_01$tstop)

my.data.steptwo.12 <- my.data.steptwo %>% filter(to != "cens")
summary(my.data.steptwo$tstart)

