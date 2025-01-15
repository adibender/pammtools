
library(RColorBrewer)
library(gridExtra)
library(scam)
library(dplyr)

data(cancer, package="survival")

head(mgus2)

# get overview over data
# censored observations
nrow(mgus2 %>% filter(pstat == 0 & death == 0))

# transition from pcm to death
nrow(mgus2 %>% filter(pstat == 1 & death == 1))

# transition to pcm
nrow(mgus2 %>% filter(pstat == 1))

# transition from non-pcm to death
nrow(mgus2 %>% filter(pstat == 0 & death == 1))

# define status column
my.mgus2 <- mgus2 %>%
  mutate(status = ifelse(pstat == 1 | death == 1, 1, 0)
         , ptime = ifelse(ptime == futime, ptime-1, ptime)
         # , transition = case_when(
         #   pstat == 0 & death == 0 ~ "cens"
         #   , pstat == 1 & death == 0 ~ "0->1"
         #   , pstat == 1 & death == 1 ~ "1->2"
         #   , pstat == 0 & death == 1 ~ "0->2"
         # )
         )

# build data sets for each transition to have "long format"
my.mgus2.pcm <- my.mgus2 %>% 
  filter(pstat == 1) %>%
  mutate(tstart = 0
         , tstop = ptime
         , from = 0
         , to = 1
         , transition = "0->1") %>%
  select(-pstat, -ptime, -death, -futime)

my.mgus2.death <- my.mgus2 %>% 
  filter(death == 1 & pstat == 0) %>%
  mutate(tstart = 0
         , tstop = futime
         , from = 0
         , to = 2
         , transition = "0->2") %>%
  select(-pstat, -ptime, -death, -futime)

my.mgus2.pcmdeath <- my.mgus2 %>% 
  filter(death == 1 & pstat == 1) %>%
  mutate(tstart = ptime
         , tstop = futime
         , from = 1
         , to = 2
         , transition = "1->2") %>%
  select(-pstat, -ptime, -death, -futime)

my.mgus2.cens <- my.mgus2 %>% 
  filter(status == 0) %>%
  mutate(tstart = 0
         , tstop = futime
         , from = 0
         , to = 0
         , transition = "cens") %>%
  select(-pstat, -ptime, -death, -futime)

# merge transitions
my.mgus2.pam <- bind_rows(my.mgus2.pcm, my.mgus2.death, my.mgus2.pcmdeath, my.mgus2.cens)

# pamm pipeline
my.mgus2.pam <- my.mgus2.pam %>% add_counterfactual_transitions()

cal.my.mgus2.pam <- as_ped_multistate(
  data       = my.mgus2.pam,
  formula    = Surv(tstart, tstop, status)~ .,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

# dim check for run time
dim(my.mgus2.pam)
dim(cal.my.mgus2.pam)

# check progress, runs approx 30-45min
ctrl <- gam.control(trace = TRUE)
pam <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                 + as.factor(transition)
                 + s(mspike, by=as.factor(transition))
                 + s(hgb, by=as.factor(transition))
                 + s(age, by=as.factor(transition))
                 + sex
                 , data = cal.my.mgus2.pam
                 , family=poisson()
                 , offset=offset
                 , control = ctrl)

summary(pam)
plot(pam, ylim = c(-1, 1), page=1)

# fit gam with forced monotonic decline at the end of the data

mpam <- scam::scam(ped_status ~ s(tend, by=as.factor(transition)) 
             + as.factor(transition)
             + s(hgb, bs = "mpd", by=as.factor(transition))
             + s(age, by= as.factor(transition))
             , data = cal.my.mgus2.pam
             , family = poisson()
             , offset = offset)

# visualize splines
plot(mpam, xlim = c(0,20), ylim = c(-1, 1), page=1)

summary(mpam)
class(mpam)

plot(mpam, xlim = c(0,20), ylim = c(-1, 1), page=1)
plot(mpam, select=4, ylim = c(-0.75, 2))
plot(mpam, select=6, ylim = c(-0.75, 2))
# 
# plot.gam(pam, select=4, xlim=c(0.5,2), ylim=c(-1,1))
# plot.gam(pam, select=5, xlim=c(0.5,2), ylim=c(-1,1))
plot.gam(pam, select=6, xlim=c(0.5,2), ylim=c(-1,2))
plot.gam(pam, select=7, xlim=c(10,16), ylim=c(-1,2))
plot.gam(pam, select=8, xlim=c(6,15), ylim=c(-1,2))

plot.gam(pam, select=10)
# plot.gam(pam, select=9, xlim=c(10,16), ylim=c(-1,1))
# plot.gam(pam, select=10)
# plot.gam(pam, select=11)
# plot.gam(pam, select=12)


# prepare add_trans_prob pipeline
prep.cal.my.mgus2.pam <- make_newdata(cal.my.mgus2.pam
                                      , tend = unique(tend)
                                      , transition=unique(transition)
                                      , age = quantile(age, probs=seq(0, 0.75, by=0.25), na.rm = T)
                                      , hgb = c(11,11.5,14.5,15)
                                      ) %>% 
  group_by(transition
           , hgb
           , age
           ) %>% 
  add_cumu_hazard(pam) 
traceback()
unique(prep.cal.my.mgus2.pam$hgb)

#-------------------------------------------------------------------------------
# Example
# hgb groups (high hgb, low hgb differences) plots 
# grouped by age
#-------------------------------------------------------------------------------

# workaround for grouped data -> include in add_trans_prob() when time
old_groups <- dplyr::groups(prep.cal.my.mgus2.pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- prep.cal.my.mgus2.pam %>% ungroup(transition)
test <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

test$hgb <- as.factor(test$hgb)

# plot transitions

transition_ggplot <- ggplot(test, aes(x=tend, y=trans_prob, linetype=hgb)) + 
  geom_line(aes(col=hgb)) + 
  facet_wrap(~transition + as.factor(age), ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 100)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw()
ggsave("tmp/example/transition_probabilities.pdf", plot = transition_ggplot, width = 10)
  

ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=as.factor(age))) + 
  facet_wrap(~transition + hgb, ncol = 4, scales = "free_y") +
  xlim(c(0, 100)) +
  scale_color_brewer(palette = "Blues")


#-------------------------------------------------------------------------------
# Example
# current example in abstract - hgb contour plots for linear and non-linear fit
# no age effect
#-------------------------------------------------------------------------------

# linear hgb
ctrl <- gam.control(trace = TRUE)
pam_lin_hgb <- mgcv::gam(ped_status ~ tend*as.factor(transition)
                         + as.factor(transition)
                         + sex
                         + hgb * as.factor(transition)
                         , data = cal.my.mgus2.pam
                         , family=poisson()
                         , offset=offset
                         , control = ctrl)


summary(pam_lin_hgb)



# PAM

ctrl <- gam.control(trace = TRUE)
pam_hgb <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                     + as.factor(transition)
                     + sex
                     + s(hgb, by=as.factor(transition))
                     , data = cal.my.mgus2.pam
                     , family=poisson()
                     , offset=offset
                     , control = ctrl)
summary(pam_hgb)

# test tensor spline for tend and hgb
cal.my.mgus2.pam$transition <- as.factor(cal.my.mgus2.pam$transition)
ctrl <- gam.control(trace = TRUE)
pam_hgbtend <- mgcv::bam(ped_status ~ s(tend, by=transition) 
                     + transition
                     + sex
                     + te(hgb, age, by=transition)
                     , data = cal.my.mgus2.pam
                     , family=poisson()
                     , offset=offset
                     , method = "fREML"
                     , discrete = TRUE
                     , control = ctrl)

summary(pam_hgbtend)

plot(pam_hgbtend, page = 1)
plot(pam_hgbtend, select=4)

gg_tensor(pam_hgbtend)

test_cal_pam <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , hgb = seq(6, 16, by = 0.1)) %>% 
  group_by(transition, hgb) %>% 
  add_trans_prob(pam_hgb) 


# test new wrapper
test_cal_pam <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , hgb = seq(6, 16, by = 0.1)) %>% 
  group_by(transition, hgb)
debug <- test_cal_pam %>% add_trans_prob(pam_hgb)
debug <- test_cal_pam %>% add_cumu_hazard(pam_hgb)

combined_contour <- ggplot(debug, aes(x=tend, y=hgb, z=trans_prob)) +
  geom_tile(aes(fill=trans_prob)) +
  scale_fill_gradient2(
    name = "transition\nprobability"
    , low  = "steelblue"
    , high = "firebrick2"
    , midpoint=0.5)+
  stat_contour(col="grey30",lwd = 1.1) +
  facet_wrap(~ transition
             , ncol=3
             , labeller = label_both
  ) +
  xlim(c(0,100)) +
  xlab("time in days") +
  ylab("hemoglobin") +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
  ) 

combined_contour

# workaround for grouped data -> include in add_trans_prob() when time
old_groups <- dplyr::groups(test_cal_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- test_cal_pam %>% ungroup(transition)
test <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

test_sub <- test %>% filter(transition == "0->2")

# PAM LINEAR
test_cal_lin <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , hgb = seq(6, 16, by = 0.1)) %>% 
  group_by(transition, hgb) %>% 
  add_cumu_hazard(pam_lin_hgb) 
# workaround for grouped data -> include in add_trans_prob() when time
old_groups <- dplyr::groups(test_cal_lin)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- test_cal_lin %>% ungroup(transition)
test_lin <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition))|> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

test_lin_sub <- test_lin %>% 
  filter(transition == "0->2")

# combine plots
test_lin_sub <- test_lin_sub %>% mutate(model = "linear")
test_sub <- test_sub %>% mutate(model = "non-linear")

test <- rbind(test_lin_sub, test_sub)

test$model <- as.factor(test$model)

table(test$model)

combined_contour <- ggplot(test, aes(x=tend, y=hgb, z=trans_prob)) +
  geom_tile(aes(fill=trans_prob)) +
  scale_fill_gradient2(
    name = "transition\nprobability"
    , low  = "steelblue"
    , high = "firebrick2"
    , midpoint=0.5)+
  stat_contour(col="grey30",lwd = 1.1) +
  facet_wrap(transition ~ model
             , ncol=2
             , labeller = label_both
             ) +
  xlim(c(0,100)) +
  xlab("time in days") +
  ylab("hemoglobin") +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
  ) 

combined_contour

# combine spline plots
hgb_df <- cal.my.mgus2.pam %>%
  make_newdata(hgb = seq(6, 16, by = 0.1), transition = unique(transition)) %>%
  add_term(pam_hgb, term = "hgb")
time_df <- cal.my.mgus2.pam %>%
  make_newdata(tend = unique(tend), transition = unique(transition)) %>%
  add_term(pam_hgb, term = "tend")

hgb_df <- hgb_df %>% filter(transition == "0->2") %>% mutate(model = "non-linear")
time_df <- time_df %>% filter(transition == "0->2") %>% mutate(model = "non-linear")

hgb_pp <- ggplot(hgb_df, aes(x = hgb)) +
  geom_line(aes(y = fit), lwd = 1.1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = .2) +
  facet_wrap(transition ~ model, labeller = label_both) +
  ylab("s(hgb, 4.65):transition0->2") +
  xlab("hgb") + coord_cartesian(ylim = c(-0.5, 2.5)) +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
         ) 

hgb_pp

time_pp <- ggplot(time_df, aes(x = tend)) +
  geom_line(aes(y = fit), lwd = 1.1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = .2) +
  facet_wrap(transition ~ model, labeller = label_both) +
  ylab("s(tend, 7.7):transition0->2") +
  xlab("time in days") + coord_cartesian(ylim = c(-0.5, 1)) +
  xlim(c(0,100)) +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
         )
  

time_pp

combined_spline_contour <- grid.arrange(hgb_pp, time_pp, combined_contour, ncol=3, widths=c(1/4, 1/4, 2/4))
combined_spline_contour
# save for paper
ggsave("tmp/example/transition_probabilities_hgb_contour.pdf", plot = combined_spline_contour, width = 16) 

# save for poster
ggsave("tmp/example/transition_probabilities_hgb_contour.png"
       , plot = combined_spline_contour
       , width = 35
       , height = 10
       , dpi = 300
       , units = "cm"
       )
dev.off()
# 
# combined_ggplot <- ggplot(test, aes(x=tend, y=trans_prob)) + 
#   geom_line(aes(group=hgb, col=hgb)) + 
#   facet_wrap(transition~ model, ncol = 2, labeller = label_both) +
#   # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
#   # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
#   xlim(c(0, 200)) +
#   ylab("Transition Probability") +
#   xlab("time") +
#   theme_bw() +
#   theme(legend.position = "bottom")
# combined_ggplot
# ggsave("tmp/example/transition_probabilities_hgb.pdf", plot = combined_ggplot, width = 10)
# 
# pdf("tmp/example/hgb_spline.pdf", width = 10)
# plot(pam_hgb, select=5, ylim=c(-0.5,1.5), xlim=c(6,16))
# dev.off()

#-------------------------------------------------------------------------------
# SIMULATION STUDY
# model graph with different states
#-------------------------------------------------------------------------------

# define true functions for plotting
f0 <- function(t) {
  dgamma(t, 8, 2) * 12 - 0.5 
}

fx <- function(x, transition) {
  (transition == "1->2")*(0.5 + 0.25*x**3) + (transition == "2->3")*(0.4 - x**2)
}

ft <- function(t, transition) {
  (transition == "1->3")*f0(t)
}


# # define hazard for transitions 1->2 & 1->3
# form_stepone <- formula( ~ 0.5 + 0.05*x1**3 + f0(t) | 0.6 - 0.05*x1**2  )
# # define hazard for transition 2->3
# form_steptwo <- formula( ~ - 0.4 + sin(x2) + (t+x3 <= 3) * f0(t + x3))

# define hazard for transitions 1->2 & 1->3
form_stepone <- formula( ~ 0.5 + 0.25*x1**3| f0(t)  )
# define hazard for transition 2->3
form_steptwo <- formula( ~ 0.4 - x1**2)

nSim = 20

n = 2500
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)
seq_x = seq(-2, 2, by = 0.1)
seq_t = seq(0, 3, by = 0.1)

library(foreach)
library(doParallel)

for (i in 1:nSim) {
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

# second step
# take all relevant subjects, i.e. subjects with transition 1->2, i.e. currently in state 2
data_steptwo <- sim_df_stepone %>% 
  filter(to == 2, status == 1) %>%
  select(id, x1, from, t = tstop)

sim_df_steptwo <- sim_pexp_cr(form_steptwo, data_steptwo, cut) %>%
  mutate(from = 2
         , to = type + 2
         , transition = case_when(
           type == 1 ~ "2->3",
           .default = "err"
         )
         , status = ifelse(time + t > 3, 0, status)
         , time = min(time + t, 3)
         , hazard2 = 0) %>%
  rename(tstart = t, tstop = time) %>% 
  filter(status == 1)

sim_df <- rbind(sim_df_stepone, sim_df_steptwo)
sim_df <- sim_df %>% add_counterfactual_transitions()
# 
# test <- sim_df %>% mutate(time = ifelse(tstop > 2, 1, 0)) %>% filter(status == 1)
# table(test$time, test$transition)

# go on with pamm tools procedure
cal_sim_df <- as_ped_multistate(
  data       = sim_df,
  formula    = Surv(tstart, tstop, status)~ transition + x1,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")
# 
# dim(cal_sim_df)
# head(cal_sim_df)

ctrl <- gam.control(trace = TRUE)
bam_sim <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) 
                     + as.factor(transition)
                     + s(x1, by=as.factor(transition))
                     , data = cal_sim_df
                     , family=poisson()
                     , offset=offset
                     , discrete = T
                     , method = "fREML"
                     , control = ctrl)

# summary(bam_sim)

plot(bam_sim, select = 6, xlim = c(-2,2), ylim = c(-1,1))

if(i == 1) {
x1_df <- cal_sim_df %>%
  make_newdata(x1 = seq_x, transition = unique(transition)) %>%
  add_term(bam_sim, term = "x1")%>%
  mutate(true.value = fx(x=x1, transition=transition)
         , shift = case_when(
           transition == "1->2" ~ coef(bam_sim)[1]
           , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
           , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
           , TRUE ~ 0
         )
         , fit = fit + shift
         , ci_lower = ci_lower + shift
         , ci_upper = ci_upper + shift
         , iter = i)

tend_df <- cal_sim_df %>%
  make_newdata(tend = unique(tend), transition = unique(transition)) %>%
  add_term(bam_sim, term = "tend") %>%
  mutate(true.value = ft(tend, transition)
         , shift = case_when(
           transition == "1->2" ~ coef(bam_sim)[1]
           , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
           , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
           , TRUE ~ 0
         )
         , fit = fit + shift
         , ci_lower = ci_lower + shift
         , ci_upper = ci_upper + shift
         , iter = i
  )
} else {
  x1_df_temp <- cal_sim_df %>%
    make_newdata(x1 = seq_x, transition = unique(transition)) %>%
    add_term(bam_sim, term = "x1") %>%
    mutate(true.value = fx(x=x1, transition=transition)
           , shift = case_when(
             transition == "1->2" ~ coef(bam_sim)[1]
             , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
             , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
             , TRUE ~ 0
           )
           , fit = fit + shift
           , ci_lower = ci_lower + shift
           , ci_upper = ci_upper + shift
           , iter = i)
    x1_df <- rbind(x1_df, x1_df_temp)
    
    tend_df_temp <- cal_sim_df %>%
      make_newdata(tend = unique(tend), transition = unique(transition)) %>%
      add_term(bam_sim, term = "tend") %>%
      mutate(true.value = ft(tend, transition)
             , shift = case_when(
               transition == "1->2" ~ coef(bam_sim)[1]
               , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
               , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
               , TRUE ~ 0
             )
             , fit = fit + shift
             , ci_lower = ci_lower + shift
             , ci_upper = ci_upper + shift
             , iter = i
      )
    tend_df <- rbind(tend_df, tend_df_temp)
}

}
# 
# simulation_x <- x1_df
# simulation_t <- tend_df

# sim_x_compl <- rbind(simulation_x, x1_df)
# sim_t_compl <- rbind(simulation_t, tend_df)

# x1_df_new <- x1_df %>% mutate(iter = iter + 15)
# tend_df_new <- tend_df %>% mutate(iter = iter + 15)

sim_x_compl <- rbind(simulation_x, x1_df_new)
sim_t_compl <- rbind(simulation_t, tend_df_new)

sim_x_compl_new <- sim_x_compl %>%
  mutate(transition = case_when(
    transition == "1->2" ~"0->1"
    , transition == "1->3" ~"0->2"
    , transition == "2->3" ~"1->2"
  )) %>%
  group_by(transition, x1) %>% 
  summarise(fit = fit
            , iter = iter
            , true.value = true.value
            , mean_fit = mean(fit)
            , mean_ci_up = quantile(fit, 0.95)
            , mean_ci_lo = quantile(fit, 0.05)
            )

sim_t_compl_new <- sim_t_compl %>%
  mutate(tend = round(tend, 2)
         , transition = case_when(
           transition == "1->2" ~"0->1"
           , transition == "1->3" ~"0->2"
           , transition == "2->3" ~"1->2"
         )) %>%
  group_by(transition, tend) %>% 
  summarise(fit = fit
            , iter = iter
            , true.value = true.value
            , mean_fit = mean(fit)
            , mean_ci_up = mean(ci_upper)
            , mean_ci_lo = mean(ci_lower)
  )

test <- sim_x_compl %>% 
  mutate(cover = as.numeric((true.value >= ci_lower & true.value <= ci_upper))) %>% 
  select(x1, fit, ci_lower, ci_upper, true.value, transition, iter, cover)

test1 <- test %>% 
  select(transition, cover) %>% 
  group_by(transition) %>% 
  summarise(mean = mean(cover))
test1

ggplot(test1, aes(x=x1, y=mean, col=transition)) +
  geom_line()

gg_x1 <- ggplot(sim_x_compl_new, aes(x = x1)) +
  geom_line(aes(y = (1+1.5*(transition == "1->2"))*fit - (transition == "1->2")*0.75, group = iter), lwd = 0.5, color = "grey") +
  geom_line(aes(y = (1+1.5*(transition == "1->2"))*mean_fit - (transition == "1->2")*0.75), color="firebrick2", lwd = 1.1) +
  geom_line(aes(y = (1+1.5*(transition == "1->2"))*mean_ci_up - (transition == "1->2")*0.75), color="firebrick2", lwd = 1.1 , linetype="dotted") +
  geom_line(aes(y = (1+1.5*(transition == "1->2"))*mean_ci_lo - (transition == "1->2")*0.75), color="firebrick2", lwd = 1.1, linetype="dotted") +
  geom_line(aes(y = true.value + (transition == "1->2")*(0.1)), lwd = 1.1) +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = iter), alpha = .2) +
  ylab("f(x)") +
  xlab("x") + coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme( strip.text = element_text(size = 14)
         , axis.text = element_text(size = 14)) +
  facet_wrap(~ transition, labeller = label_both)

gg_x1

sim_t_compl_sub <- sim_t_compl_new %>% filter(transition == "0->2")

gg_tend <- ggplot(sim_t_compl_sub, aes(x = tend)) +
  geom_line(aes(y = fit, group = iter), lwd = 0.5, color = "grey") +
  stat_smooth(aes(y = mean_fit), color="firebrick2", lwd = 1.1) +
  geom_line(aes(y = true.value), lwd=1.1) +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .2) +
  ylab("f(t)") +
  xlab("t") + coord_cartesian(ylim = c(-4, 4)) +
  theme_bw() +
  theme( strip.text = element_text(size = 14)
         , axis.text = element_text(size = 14)) +
  facet_wrap(~ transition, labeller = label_both)


gg_tend 


combined_simulation_contour <- grid.arrange(gg_x1, gg_tend, ncol=2, widths = c(3/4, 1/4))
combined_simulation_contour
ggsave("tmp/example/simulation_splines.pdf", plot = combined_simulation_contour, width = 16, height = 5)
dev.off()

plot(bam_sim, page = 1)
plot(bam_sim, select = 2, xlim = c(0,3), ylim=c(-4,4))

par(mfrow=c(1,3))
plot(bam_sim, select = 4, xlim = c(-2,2), ylim=c(-2,2))
plot(bam_sim, select = 5, xlim = c(-2,2), ylim=c(-2,2))
plot(bam_sim, select = 6, xlim = c(-2,2), ylim=c(-2,2))
par(mfrow=c(1,1))

pam_sim <- bam_sim

new_sim_pam <- make_newdata(cal_sim_df
                            , tend = unique(tend)
                            , transition=unique(transition)
                            , x1 = seq(-2,2, by = 0.2)
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
             , scales = "free_y"
             , labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 3)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

# plot transitions
ggplot(test_msm, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(group = x1, col = x1)) + 
  facet_wrap(~transition
             , ncol = 3
             , scales = "free_y"
             , labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0,3)) +
  ylim(c(0,1)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 

test_contour <- test_msm %>% mutate(tend = round(tend, 1), x1 = round(x1,2))

simulation_contour <- ggplot(test_contour, aes(x=tend, y=x1, z=trans_prob)) +
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
  xlim(c(0,3)) +
  theme_bw() +
  theme(legend.position = "bottom"
        , strip.text = element_text(size = 14)
        , axis.text = element_text(size = 14))

simulation_contour


# combine spline plots
x1_df <- cal_sim_df %>%
  make_newdata(x1 = seq(-2, 2, by = 0.1), transition = unique(transition)) %>%
  add_term(pam_sim, term = "x1")

x1_pp <- ggplot(x1_df, aes(x = x1)) +
  geom_line(aes(y = fit), lwd = 1.1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = .2) +
  ylab("f(x1)") +
  xlab("x1") + coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_bw() +
  theme( strip.text = element_text(size = 14)
         , axis.text = element_text(size = 14)) +
  facet_wrap(~ transition, labeller = label_both)

x1_pp


combined_simulation_contour <- grid.arrange(x1_pp, simulation_contour, ncol=1, heights = c(0.8, 1))
combined_simulation_contour
ggsave("tmp/example/transition_probabilities_simulation_contour.pdf", plot = combined_simulation_contour, width = 16, height = 10)
dev.off()

#-------------------------------------------------------------------------------
# Example
# real contour plots
# no age effect
#-------------------------------------------------------------------------------
ctrl <- gam.control(trace = TRUE)
pam_te_hgb <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                         + as.factor(transition)
                         + te(tend, hgb, by=as.factor(transition))
                         , data = cal.my.mgus2.pam
                         , family=poisson()
                         , offset=offset
                         , control = ctrl)


summary(pam_te_hgb)

plot(pam_te_hgb, select=6)

#-------------------------------------------------------------------------------
# Example
# https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf
# p. 57
#-------------------------------------------------------------------------------
library(survival)
data(nafld, package="survival")
# recalculate example
# create transition figure
state5 <- c("0MC", "1MC", "2MC", "3MC", "death")
tmat <- matrix(0L, 5, 5, dimnames=list(state5, state5))
tmat[1,2] <- tmat[2,3] <- tmat[3,4] <- 1
tmat[-5,5] <- 1
statefig(rbind(4,1), tmat)

# rebuild dataset
ndata <- tmerge(nafld1[,1:8], nafld1, id=id, death= event(futime, status))
ndata <- tmerge(ndata, subset(nafld3, event=="nafld"), id,
                  nafld= tdc(days))
ndata <- tmerge(ndata, subset(nafld3, event=="diabetes"), id = id,
                  diabetes = tdc(days), e1= cumevent(days))
ndata <- tmerge(ndata, subset(nafld3, event=="htn"), id = id,
                  htn = tdc(days), e2 = cumevent(days))
ndata <- tmerge(ndata, subset(nafld3, event=="dyslipidemia"), id=id,
                  lipid = tdc(days), e3= cumevent(days))
ndata <- tmerge(ndata, subset(nafld3, event %in% c("diabetes", "htn",
                                                     "dyslipidemia")),
                  id=id, comorbid= cumevent(days))
summary(ndata)
head(ndata)

with(ndata, if (any(e1>1 | e2>1 | e3>1)) stop("multiple events"))
ndata$cstate <- with(ndata, factor(diabetes + htn + lipid, 0:3,
                                     c("0mc", "1mc", "2mc", "3mc")))
temp <- with(ndata, ifelse(death, 4, comorbid))
ndata$event <- factor(temp, 0:4,
                      c("censored", "1mc", "2mc", "3mc", "death"))
ndata$age1 <- ndata$age + ndata$tstart/365.25 # analysis on age scale
ndata$age2 <- ndata$age + ndata$tstop/365.25
check1 <- survcheck(Surv(age1, age2, event) ~ nafld + male, data=ndata,
                      id=id, istate=cstate)
check1


nfit1 <- coxph(list(Surv(age1, age2, event) ~ nafld + male,
                    "0mc":state("1mc", "2mc", "3mc") ~ nafld+ male / common,
                    2:3 + 2:4 ~ nafld + male / common,
                    0:"death" ~ male / common),
               data=ndata, id=id, istate=cstate)

summary(nfit1)
str(nfit1)
nfit1$transitions

#-------------------------------------------------------------------------------
# PAMMTOOLS
# prepare data set for as_ped function
#-------------------------------------------------------------------------------

head(ndata)
str(ndata)

# censoring is coded with 0
# instead of string names, name states 0,1,...,4, where 4 = death
ndata$from <- as.numeric(with(ndata, factor(diabetes + htn + lipid, 0:3,
                                   c(0, 1, 2, 3))))

ndata$to <- as.numeric(factor(temp, 0:4,
                      c(0, 1, 2, 3, 4)))
ndata <- ndata %>% 
  mutate(status = case_when(
    as.character(event) == "censored" ~ 0,
    TRUE ~ 1)
    , from = from - 1
    , to = to - 1
    )
head(ndata)

str(ndata)

perc <- 0.1
ndata_sample <- sample_n(ndata, round(dim(ndata)[1]*perc,0))
dim(ndata_sample)


table(ndata$status)
ndata_pamm <- ndata_sample %>% add_counterfactual_transitions()
head(ndata_pamm)
dim(ndata_pamm)

length(unique(ndata_pamm$tstop))
length(unique(ndata_pamm$transition))

cal_ndata_pamm <- as_ped_multistate(
  data       = ndata_pamm,
  formula    = Surv(tstart, tstop, status)~ age + height + bmi + weight + male,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(cal_ndata_pamm)
ctrl <- gam.control(trace = TRUE)
# dauert sehr lange
# pam_ndata <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
#                  + as.factor(transition)
#                  + s(bmi)
#                  + age
#                  , data = cal_ndata_pamm
#                  , family=poisson()
#                  , offset=offset
#                  , control = ctrl)

pam_ndata <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) 
                       + as.factor(transition)
                       , data = cal_ndata_pamm
                       , family=poisson()
                       , offset=offset
                       , discrete = T #include to make it faster
                       , method = "fREML" #include to make it faster
                       , control = ctrl)

summary(pam_ndata)

head(cal_ndata_pamm)


new_ndata_pam <- make_newdata(cal_ndata_pamm
                              , tend = unique(tend)
                              , transition=unique(transition)
                              #, age = quantile(age, probs=c(0.05, 0.5, 0.95))
                              ) %>% 
  group_by(transition) %>% 
  add_cumu_hazard(pam_ndata) 

old_groups <- dplyr::groups(new_ndata_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_ndata <- new_ndata_pam %>% ungroup(transition)
test_ndata <- group_split(res_ndata) |> 
  map(res_ndata, .f = ~ group_by(.x, transition)) |> 
  map(res_ndata, .f = ~ add_trans_prob(.x)) |>
  map(res_ndata, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot transitions
ggplot(test_ndata, aes(x=tend, y=trans_prob)) + 
  geom_line() + 
  facet_wrap(~transition, ncol = 4, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 6700)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 






























