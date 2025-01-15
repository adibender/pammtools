# ------------------------------------------------------------------------------
#
# SIMULATION EXAMPLE SECTION 3.2. + 4.2
#
# => BASIC SIMULATION WITH CONSTANT HAZARDS
# => COMPARISON OF PAMMTOOLS AND MVNA
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

# here: alpha_01 = 0.3, alpha_02 = 0.6
n = 100

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
tra
id <- seq_along(obs.cause)
from <- rep(0, length(obs.cause))
to <- as.factor(ifelse(obs.cause == 0, "cens", obs.cause))
my.data <- data.frame(id, from, to, time = obs.times)
head(my.data)

#View(my.data)

# calculate cause specific hazards
library(mvna)
my.nelaal <- mvna(my.data, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal, strip=strip.custom(bg="white"), ylab="Cumulative Hazard", lwd=2)
}
# same graphs as in Beyersmann --> correct

plot(my.nelaal$'0 1'$time, my.nelaal$'0 1'$na)

# estimate transition probabilites
library(etm)
etm.my.data<- etm(my.data, c("0", "1", "2"), tra, "cens", s = 0)

par(mfrow = c(1,2))
plot(etm.my.data, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.8), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.my.data, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 0.8),
     xlab = "Days",
     ci.fun = "cloglog")
par(mfrow = c(1,1))



# try with pammtools
# my.data <- my.data %>% mutate(status = ifelse(to == "cens", 0, 1),
#                               tstart = 0,
#                               tstop = time,
#                               transition = paste0(from, "->", to))

pamm.my.data <- my.data %>%
  mutate(status = ifelse(to == "cens", 0, 1),
         tstart = 0,
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

pam <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) + 
                   as.factor(transition)
                 , data = cal.my.data
                 , family=poisson()
                 , offset=offset)

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

test_msm.01 <- test_msm %>% filter(transition == "0->1")
plot.01 <- ggplot(test_msm.01, aes(x=tend, y=cumu_hazard)) + geom_line(aes(col=transition)) + xlim(c(0, 5)) + ylim(c(0, 1))
combined.plot.01 <- plot.01 + geom_line(data = my.nelaal$'0 1', aes(x=time, y=na), color="blue")
print(combined.plot.01)

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 4, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylim(c(0,0.8)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# --> correct sampling, same results
