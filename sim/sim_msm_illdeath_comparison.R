
# ------------------------------------------------------------------------------
#
# SIMULATION EXAMPLE SECTION 3.2. + 4.2 => EXTEND TO MULTISTATE
# => COMPARE RESULTS FROM SIMPLE SAMPLING METHOD WITH PAMMTOOLS SIM_PEXP_MSM
#
# SUCCESSFUL
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
# ---------------------------------------------------------------------------- #

# here: alpha_01 = 0.3, alpha_02 = 0.6, alpha_12 = 0.5
n = 750

event.times <- rexp(n, 0.9)

# alpha_01 / (alpha_01 + alpha_02) = 1/3

# alpha01 = h(t,age,sex) = x1 + t^2
# simpexpcr
# bei zweitem datensatz müssen covariates aktualisiert werden

f.cause <- rbinom(n, size = 1, prob = 1/3)
f.cause <- ifelse(f.cause == 0, 2, 1)
cens.times <- runif(n,0,5)
obs.times <- pmin(event.times, cens.times)
obs.cause <- c(event.times <= cens.times) * f.cause

# now section 4.2
tra <- matrix(FALSE, ncol = 3, nrow = 3)
dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
tra[1, 2:3] <- TRUE
tra[2, 3] <- TRUE
tra
id <- seq_along(obs.cause)
from <- rep(0, length(obs.cause))
to <- as.factor(ifelse(obs.cause == 0, "cens", obs.cause))
my.data <- data.frame(id, from, to, tstart = 0, time = obs.times, cens.times)
head(my.data)

my.data.steptwo <- my.data %>% filter(to == 1)
head(my.data.steptwo)

# define new sample size based on subjects in risk set
n.steptwo <- dim(my.data.steptwo)[1]

event.times.steptwo <- rexp(n.steptwo, 0.5)
# no need for rbinom since only one transition is possible

my.data.steptwo <- data.frame(cbind(my.data.steptwo, event.times.steptwo)) %>%
  mutate(tstart = time
         , time = pmin(event.times.steptwo + time, cens.times)
         , from = as.factor(1)
         , to = as.factor(ifelse((event.times.steptwo + time) >= cens.times, "cens", 2))
  ) %>%
  select(id, from, to, tstart, time, cens.times)

head(my.data.steptwo)
my.data.steptwo %>% filter(to == 2)

my.data.complete <- rbind(my.data, my.data.steptwo)
dim(my.data.complete)
head(my.data.complete)

#View(my.data)

# try to sample with sim_pexp
t_mat <- matrix(data = NA, nrow = 3, ncol = 3)
t_mat[1,2] <- "log(0.3)"
t_mat[1,3] <- "log(0.6)"
t_mat[2,3] <- "log(0.5)"

n = 1000
data <- cbind.data.frame(
  id = 1:n,
  #x1 = runif(n, -3, 3),
  #x2 = runif(n, 0, 6),
  from = 1,
  t = 0)

head(data)
# f1 <- function(x) log(0.3)
# f2 <- function(x) log(0.6)
# 
# plot(data$x1, f1(data$x1))
# plot(data$x2, f2(data$x2))

cut = seq(0, 3, by = 0.01)

msm_df <- sim_pexp_msm(
  t_mat = t_mat,
  data = data,
  cut = cut,
  keep_transitions_at_risk = TRUE)
head(msm_df)

my.data <- msm_df %>% filter(status == 1) %>% mutate(from = from - 1, to = as.factor(to-1), time = tstop)
table(my.data$from, my.data$to)

head(my.data)


# calculate cause specific hazards
library(mvna)
my.nelaal <- mvna(my.data.complete, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal
         , strip=strip.custom(bg="white")
         , ylab="Cumulative Hazard"
         , lwd=2
         , xlim=c(0,3)
         , ylim=c(0,1.5))
}
# same graphs as in Beyersmann --> correct

plot(my.nelaal$'0 1'$time, my.nelaal$'0 1'$na)

# estimate transition probabilites
library(etm)
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
# my.data <- my.data %>% mutate(status = ifelse(to == "cens", 0, 1),
#                               tstart = 0,
#                               tstop = time,
#                               transition = paste0(from, "->", to))

pamm.my.data <- my.data.complete %>%
  mutate(status = ifelse(to == "cens", 0, 1),
         tstop = time,
         transition = paste0(from, "->", to)) %>%
  select(id
         , tstart
         , tstop
         , from
         , to
         , transition
         , status)

head(pamm.my.data)


pamm.my.data <- pamm.my.data %>% add_counterfactual_transitions()
head(pamm.my.data)

cal.my.data <- as_ped_multistate(
  data       = pamm.my.data,
  formula    = Surv(tstart, tstop, status)~ transition,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

head(cal.my.data)

pam <- mgcv::bam(ped_status ~ tend*as.factor(transition) + 
                   as.factor(transition)
                 , data = cal.my.data
                 , family=poisson()
                 , offset=offset
                 , method = "fREML"
                 , discrete = T)

summary(pam)


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

pdf("tmp/example/simulation_aj_equality.pdf", width = 8, height = 2.5)
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
dev.off()
