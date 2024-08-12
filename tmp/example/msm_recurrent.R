# ---------------------------------------------------------------------------- #
# LIBRARIES
# ---------------------------------------------------------------------------- #
# Pammtools
setwd("C:/Users/ra63liw/Documents/98_git/pammtools-multi-state/pammtools")
library(dplyr)
library(survival)
library(Formula)
devtools::load_all()
source("sandbox/helpers-msm.R")
# source("sandbox/sim-pexp-msm.R")
library(mstate)
library(msm)
library(mvna)
# library(dynpred)
library(effects)
library(tidyr)
library(etm)

devtools::install("C:/Users/ra63liw/Documents/98_git/pammtools-multi-state/pammtools")

# try mstate::prothr data set
# LINK: https://cran.r-project.org/web/packages/mstate/mstate.pdf

# code from mstate documentation
data(prothr)
tmat <- attr(prothr, "trans")
pr0 <- subset(prothr, treat=="Placebo") |> filter(Tstart != Tstop)
attr(pr0, "trans") <- tmat
pr1 <- subset(prothr, treat=="Prednisone") |> filter(Tstart != Tstop)
attr(pr1, "trans") <- tmat
c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=pr0)
c1 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=pr1)


msf0 <- msfit(c0, trans=tmat)
msf1 <- msfit(c1, trans=tmat)

# Comparison as in Figure 2 of Titman (2015)
# Aalen-Johansen
pt0 <- probtrans(msf0, predt=0)[[1]] # changed predt from 1000 to 0
pt1 <- probtrans(msf1, predt=0)[[2]] # changed predt from 1000 to 0
par(mfrow=c(1,2))
plot(pt0$time, pt0$pstate2, type="s", lwd=2, xlim=c(0,4000), ylim=c(0,1),
     xlab="Time since randomisation (days)", ylab="Probability")
lines(pt1$time, pt1$pstate1, type="s", lwd=2, lty=3)
legend("topright", c("Placebo", "Prednisone"), lwd=2, lty=1:2, bty="n")
title(main="Aalen-Johansen")
# Landmark Aalen-Johansen
LMpt0 <- LMAJ(msdata=pr0, s=0, from=2) # changed predt from 1000 to 0
LMpt1 <- LMAJ(msdata=pr1, s=0, from=2) # changed predt from 1000 to 0
plot(LMpt0$time, LMpt0$pstate1, type="s", lwd=2, xlim=c(0,4000), ylim=c(0,1),
     xlab="Time since randomisation (days)", ylab="Probability")
lines(LMpt1$time, LMpt1$pstate1, type="s", lwd=2, lty=3)
legend("topright", c("Placebo", "Prednisone"), lwd=2, lty=1:2, bty="n")
title(main="Landmark Aalen-Johansen")
par(mfrow=c(1,1))

# plot hazards
mstate_dat0 <- msf0$Haz %>% mutate(transition = case_when(
  trans == 1 ~ "1->2",
  trans == 2 ~ "1->3",
  trans == 3 ~ "2->1",
  trans == 4 ~ "2->3",
  .default = "cens"
  )
  , treat = "Placebo")
mstate_dat1 <- msf1$Haz %>% mutate(transition = case_when(
  trans == 1 ~ "1->2",
  trans == 2 ~ "1->3",
  trans == 3 ~ "2->1",
  trans == 4 ~ "2->3",
  .default = "cens"
  )
, treat = "Prednisone")

mstate_dat <- rbind(mstate_dat0, mstate_dat1)


ggplot(data = mstate_dat, aes(x=time, y=Haz)) + 
  geom_line(aes(col=treat)) +
  facet_wrap(~transition, ncol = 4, labeller = label_both)

# data analysis / preparation
head(prothr)
table(prothr$from, prothr$to)
# transitions 1->2, 1->3, 2->1, 2->3 possible.
# classical multi-state setup with recurrent events. 

# try analysis with etm
# transform days in fractions of year
data <- prothr %>% mutate(Tstart = Tstart / 365.25, 
                          Tstop = Tstop /365.25,
                          from = from - 1,
                          to = to - 1)
head(data)
dim(data)


# build transition matrix for mvna
tra <- matrix(FALSE, ncol = 3, nrow = 3)
dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
tra[1, 2:3] <- TRUE
tra[2, c(1,3)] <- TRUE

# print transition matrix
tra

my.data <- prothr %>% 
  mutate(from = from -1, to = to-1) %>%
  rename(entry = Tstart, exit = Tstop) %>% 
  arrange(id, entry, exit) %>% 
  filter(!(entry == exit)) %>%
  distinct()
dim(my.data)
head(my.data)

my.data %>% filter(tstart== tstop)

table(my.data$from, my.data$to)

# calculate cause specific hazards
my.nelaal <- mvna(my.data, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal
         , strip=strip.custom(bg="white")
         , ylab="Cumulative Hazard"
         , lwd=2
         # , xlim=c(0,5)
         # , ylim=c(0,10)
         )
}
# differentiate between placebo and not
my.data.placebo <- my.data %>% filter(treat=="Placebo")
my.data.prednisone <- my.data %>% filter(treat=="Prednisone")
table(my.data.placebo$from, my.data.placebo$to)

my.nelaal <- mvna(my.data.placebo, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal
         , strip=strip.custom(bg="white")
         , ylab="Cumulative Hazard"
         , lwd=2
         # , xlim=c(0,5)
         # , ylim=c(0,10)
  )
}
my.nelaal <- mvna(my.data.prednisone, c("0", "1", "2"), tra, "cens")
if (require(lattice)){
  xyplot(my.nelaal
         , strip=strip.custom(bg="white")
         , ylab="Cumulative Hazard"
         , lwd=2
         # , xlim=c(0,5)
         # , ylim=c(0,10)
  )
}


# estimate transition probabilites
etm.prothr.placebo <- etm(my.data.placebo, c("0", "1", "2"), tra, "cens", s = 0)
etm.prothr.prednisone <- etm(my.data.prednisone, c("0", "1", "2"), tra, "cens", s = 0)

etm.prothr <- etm.prothr.placebo
par(mfrow=c(2,2))
plot(etm.prothr, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.prothr, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.prothr, tr.choice = "1 0", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.prothr, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")
par(mfrow=c(1,1))

etm.prothr <- etm.prothr.prednisone
par(mfrow=c(2,2))
plot(etm.prothr, tr.choice = "0 1", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.prothr, tr.choice = "0 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.prothr, tr.choice = "1 0", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")

plot(etm.prothr, tr.choice = "1 2", conf.int = TRUE,
     lwd = 2, legend = FALSE, ylim = c(0, 1),
     xlim = c(0, 4000), xlab = "Days",
     ci.fun = "cloglog")
par(mfrow=c(1,1))

# data already contains all possible transitions,
# no need for add_counterfactual_transition
data <- prothr %>% 
  select(-trans) %>%
  rename(tstart = Tstart, tstop = Tstop) %>%
  mutate(transition = paste0(from, "->", to),
         time = (tstop - tstart) /365.25) |>
  filter(tstart != tstop)

# # count number of transitions
# test <- data %>% filter(status == 1)
# table(test$transition)

# my.data <- data %>% add_counterfactual_transitions() %>% distinct()
dim(my.data)

head(my.data)

table(my.data$transition)

cal.my.data <- as_ped_multistate(
  data       = data,
  formula    = Surv(tstart, tstop, status)~ .,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar", 
  tdc_specials="concurrent"
  )


ctrl <- gam.control(trace = TRUE)
pam <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                 + as.factor(transition) * as.factor(treat)
                 , data = cal.my.data
                 , family=poisson()
                 , offset=offset
                 , control = ctrl)

summary(pam)
plot(pam, ylim = c(-1, 1), page=1)


new_pam <- make_newdata(cal.my.data
                        , tend = unique(tend)
                        , transition=unique(transition)
                        , treat = unique(treat)
) %>% 
  group_by(transition, treat) %>% 
  add_cumu_hazard(pam) 

%>%
  add_trans_prob(pam)

# new_haz <- make_newdata(cal.my.data
#                         , tend = unique(tend)
#                         , transition=unique(transition)
#                         , treat = unique(treat)
# ) %>% 
#   group_by(transition, treat) %>% 
#   add_hazard(pam)

old_groups <- dplyr::groups(new_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_msm <- new_pam %>% ungroup(transition)
test_msm <- group_split(res_msm) |> 
  map(res_msm, .f = ~ group_by(.x, transition)) |> 
  map(res_msm, .f = ~ add_trans_prob(.x)) |>
  map(res_msm, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()


# alternative without trans prob
test_msm <- new_pam

head(test_msm)
head(mstate_dat)

joined_test_mstate <- left_join(test_msm, mstate_dat, by=c('tend' = 'time', 'transition', 'treat')) %>%
  select(tend
         , treat
         , transition
         , cumu_hazard
         , cumu_lower
         , cumu_upper
         # , trans_prob
         , Haz) %>%
  filter(!is.na(Haz))

# joined_haz_mstate <- left_join(new_haz, mstate_dat, by=c('tend' = 'time', 'transition', 'treat')) %>%
#   select(tend, treat, transition, hazard, Haz) %>%
#   filter(!is.na(Haz))

dim(joined_test_mstate)
head(joined_test_mstate)
table(joined_test_mstate$transition)

# # 
# dim(joined_haz_mstate)
# head(joined_haz_mstate)

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=as.factor(treat))) + 
  facet_wrap(~transition, ncol = 2, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylim(c(0,0.8)) +
  xlim(c(0, 4000)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# plot hazards
long_mstate <- mstate_dat %>% 
  rename(tend = time, cumu_hazard = Haz) %>% 
  mutate(package = "mstate") %>%
  select(tend, treat, transition, cumu_hazard, package)
long_msm <- test_msm %>% 
  mutate(package = "pammtools") %>% 
  select(tend, treat, transition, cumu_hazard, package)

long_haz_df <- rbind(long_mstate, long_msm) %>%
  mutate(transition = case_when(transition == "1->2" ~ "0->1"
                                , transition == "1->3" ~ "0->2"
                                , transition == "2->1" ~ "1->0"
                                , transition == "2->3" ~ "1->2")
         )

comparison_nelaal <- ggplot(long_haz_df, aes(x=tend, y=cumu_hazard, col=treat, linetype = package)) + 
  geom_line(linewidth=1) +
  facet_wrap(~transition, ncol = 4, labeller = label_both) +
  scale_color_manual(values = c("firebrick2"
                                , "steelblue")
                     )+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +18:
  ylab("Cumulative Hazards") +
  xlab("time in days") +
  ylim(c(0,6)) +
  scale_linetype_manual(values=c("dotted", "solid")) +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
  )

comparison_nelaal
# paper figure
ggsave("tmp/example/hazards_comparison_cox_mstate_recurrent.pdf"
       , plot = comparison_nelaal
       , width = 16
       , height = 10)

# poster figure
ggsave("tmp/example/hazards_comparison_cox_mstate_recurrent.png"
       , plot = comparison_nelaal
       , width = 35
       , height = 10
       , dpi = 300
       , units = "cm"
)
dev.off()

# ---------------------------------------------------------------------------- #
#
# extend example
#
# ---------------------------------------------------------------------------- #
  
# anzahl visits in Datensatz hinzufügen
prothr_event <- prothr %>% filter(status == 1)
head(prothr_event) #--> good example id==2

prothr_event <- prothr_event %>%
  group_by(id, trans) %>%      # Step 1: Group by id and trans
  arrange(Tstart, .by_group = TRUE) %>% # Step 2: Arrange by Tstart within each group
  mutate(visits = row_number(),
         visits_cap = ifelse(visits > 3, 3, visits)) %>% # Step 3: Add a column for row count within each group
  select(id, trans, Tstart, Tstop, visits, visits_cap) # Step 4: select unique identifier and number of visits

head(prothr_event)
table(prothr_event$visits_cap)

test <- left_join(prothr, prothr_event, by=c('id', 'trans', 'Tstart', 'Tstop'))
head(test)

test <- test %>% 
  group_by(id, trans) %>% 
  arrange(id, trans, Tstart) %>%
  #fill(color, age, gender) %>% #default direction down
  fill(visits, visits_cap, .direction = "downup") %>%
  mutate(visits_cap = ifelse(is.na(visits_cap), 1, visits_cap)
         , visits = ifelse(is.na(visits), 0, visits-1))
head(test)

test %>% filter(visits_cap == 3)

table(test$visits_cap)

# pipeline to calculate gam
data <- test %>% 
  rename(tstart = Tstart, tstop = Tstop) %>%
  mutate(transition = paste0(from, "->", to),
         time = (tstop - tstart) /365.25)

# # count number of transitions
# test <- data %>% filter(status == 1)
# table(test$transition)

my.data <- data %>% add_counterfactual_transitions()
dim(my.data)

head(my.data)

table(my.data$transition)

cal.my.data <- as_ped_multistate(
  data       = data,
  formula    = Surv(tstart, tstop, status)~ .,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar", 
  tdc_specials="concurrent"
) %>%
  mutate(transition = as.factor(transition)
         , treat = as.factor(treat))

ctrl <- gam.control(trace = TRUE)
pam <- mgcv::bam(ped_status ~ s(tend, by=transition) 
                 + transition * treat
                 + visits * transition
                 , data = cal.my.data
                 , family=poisson()
                 , offset=offset
                 , method = "fREML"
                 , control = ctrl
                 )

summary(pam)
plot(pam, ylim = c(-1, 1), page=1)

new_prothr <- make_newdata(cal.my.data
                        , tend = unique(tend)
                        , transition=unique(transition)
                        , visits = 0:3
                        , treat = unique(treat)
) %>% 
  group_by(transition, visits, treat) %>% 
  add_cumu_hazard(pam)

# new_haz <- make_newdata(cal.my.data
#                         , tend = unique(tend)
#                         , transition=unique(transition)
#                         , treat = unique(treat)
# ) %>% 
#   group_by(transition, treat) %>% 
#   add_hazard(pam)

old_groups <- dplyr::groups(new_prothr)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_prothr <- new_prothr %>% ungroup(transition)
test_prothr <- group_split(res_prothr) |> 
  map(res_prothr, .f = ~ group_by(.x, transition)) |> 
  map(res_prothr, .f = ~ add_trans_prob(.x)) |>
  map(res_prothr, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

ggplot(test_prothr, aes(x=tend, y=trans_prob, col = treat)) + 
  geom_line() + 
  facet_wrap(~transition + visits, ncol = 4, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylim(c(0,1)) +
  ylab("Transition Probability") +
  xlab("time in days") +
  theme_bw()


# ----

# ----
# cav 
# results not equal. needs more work.
# ----

cav_preprocessed <- cav %>% group_by(PTNUM) %>%
  mutate(from = state,
         to = ifelse(is.na(lead(state)), state, lead(state)),
         tstart = years,
         tstop = ifelse(is.na(lead(years)), years, lead(years)),
         transition = paste0(state , "->",ifelse(is.na(lead(state)), state, lead(state))),
         status = ifelse(from == to, 0, 1)) %>% 
  filter(from != to) %>%
  rename(id = PTNUM) %>%
  select(-state, -statemax, -years, -firstobs)

cav_preprocessed %>% filter(id == 100002)

table(cav_preprocessed$from, cav_preprocessed$to)

cav_preprocessed <- cav_preprocessed %>% add_counterfactual_transitions()
table(cav_preprocessed$from, cav_preprocessed$to)

cav_cal <- as_ped_multistate(
  data       = cav_preprocessed,
  formula    = Surv(tstart, tstop, status)~ .,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar", 
  tdc_specials="concurrent"
) %>%
  mutate()

cav_cal %>% filter(id == 100002)
table(cav_cal$from, cav_cal$to)
head(cav_cal)

cav_cal <- cav_cal %>% mutate(from = from -1, to = to-1, transition = as.factor(paste0(from, "->", to)))
head(cav_cal)

ctrl <- gam.control(trace = TRUE)
pam_cav <- mgcv::bam(ped_status ~ s(tend, by=transition) 
                 + transition
                 , data = cav_cal
                 , family=poisson()
                 , offset=offset
                 , method = "fREML"
                 , descrete = T
                 , control = ctrl
)

summary(pam_cav)
plot(pam_cav, page=1)

new_cav <- make_newdata(cav_cal
                        , tend = unique(tend)
                        , transition=unique(transition)
) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_cav)

# new_haz <- make_newdata(cal.my.data
#                         , tend = unique(tend)
#                         , transition=unique(transition)
#                         , treat = unique(treat)
# ) %>% 
#   group_by(transition, treat) %>% 
#   add_hazard(pam)

old_groups <- dplyr::groups(new_cav)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_cav <- new_cav %>% ungroup(transition)
test_cav <- group_split(res_cav) |> 
  map(res_cav, .f = ~ group_by(.x, transition)) |> 
  map(res_cav, .f = ~ add_trans_prob(.x)) |>
  map(res_cav, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

ggplot(test_cav, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 3, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  ylim(c(0,1)) +
  ylab("Transition Probability") +
  xlab("time in days") +
  theme_bw() 

# continue with vignette as comparison
# Link: https://cran.r-project.org/web/packages/msm/vignettes/msm-manual.pdf

Q <- rbind ( c(0, 0.25, 0, 0.25),
             c(0.166, 0, 0.166, 0.166),
             c(0, 0.25, 0, 0.25),
             c(0, 0, 0, 0) )
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                qmatrix = Q, deathexact = 4)
cav.msm

# 10 years probability matrix of model cav.msm
toi = 4
pmatrix.msm(cav.msm, t=toi)

test_cav %>% filter(transition=="0->1", between(tend, toi - 0.01, toi + 0.01)) %>% select(tend, trans_prob)
test_cav %>% filter(transition=="0->2", between(tend, toi - 0.01, toi + 0.01)) %>% select(tend, trans_prob)
test_cav %>% filter(transition=="0->3", between(tend, toi - 0.01, toi + 0.01)) %>% select(tend, trans_prob)

hazard.msm(cav.msm)
