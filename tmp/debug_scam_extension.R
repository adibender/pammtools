data(cancer, package="survival")
head(mgus2)

library(mgcv)
library(scam)

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

# fit gam for benchmark
ctrl <- gam.control(trace = TRUE)
pam_debug <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                 + as.factor(transition)
                 + s(age, by=as.factor(transition))
                 , data = cal.my.mgus2.pam
                 , family=poisson()
                 , offset=offset
                 , control = ctrl)

summary(pam_debug)

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

par(mfrow=c(1,3))
plot(pam_hgb, select=4, ylim=c(-2,2), xlim=c(6,16))
plot(pam_hgb, select=5, ylim=c(-0.5,2), xlim=c(6,16))
plot(pam_hgb, select=6, ylim=c(-2,2), xlim=c(6,16))
par(mfrow=c(1,1))

plot(pam_hgb, ylim = c(-1, 1), page=1)




pam_lin_hgb <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                           + as.factor(transition)
                         + sex
                           + hgb * as.factor(transition)
                           , data = cal.my.mgus2.pam
                           , family=poisson()
                           , offset=offset
                           , control = ctrl)


summary(pam_lin_hgb)

ctrl <- gam.control(trace = TRUE)
pam_debug_linear <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                         + as.factor(transition)
                         + age * as.factor(transition)
                         , data = cal.my.mgus2.pam
                         , family=poisson()
                         , offset=offset
                         , control = ctrl)

summary(pam_debug_linear)

plot(pam_debug_linear, ylim = c(-1, 1), page=1)

plot(pam_debug, select=5,  ylim = c(-1, 1), xlim=c(50,80))

# fit gam with forced monotonic decline at the end of the data
scam_debug <- scam::scam(ped_status ~ s(tend, by=as.factor(transition)) 
                   + as.factor(transition)
                   + s(hgb, bs = "mpd", by=as.factor(transition))
                   , data = cal.my.mgus2.pam
                   , family = poisson()
                   , offset = offset)

# visualize splines
summary(scam_debug)
plot(scam_debug, xlim = c(0,20), ylim = c(-1, 1), page=1)
plot(scam_debug, select=5, ylim = c(-0.75, 2))

# plot s(tend, by=transition)
par(mfrow=c(2,3))
plot(pam_debug, select=1, ylim = c(-0.75, 2))
plot(pam_debug, select=2, ylim = c(-0.75, 2))
plot(pam_debug, select=3, ylim = c(-0.75, 2))
plot(scam_debug, select=1, ylim = c(-0.75, 2))
plot(scam_debug, select=2, ylim = c(-0.75, 2))
plot(scam_debug, select=3, ylim = c(-0.75, 2))
par(mfrow=c(1,1))
# kaum unterschiede.
# transition 0->2: das scam modell hat hat ein geringeres minimum bei t=320 als das gam modell
# transition 1->2: das scam modell ist leicht gecurved, allerdings sehr breite Konfidenzintervalle

# plot s(hgb, by=transition)
par(mfrow=c(2,3))
plot(pam_debug, select=4, ylim = c(-0.75, 2), xlim = c(6, 16))
plot(pam_debug, select=5, ylim = c(-0.75, 2), xlim = c(6, 16))
plot(pam_debug, select=6, ylim = c(-0.75, 2), xlim = c(6, 16))
plot(scam_debug, select=4, ylim = c(-0.75, 2), xlim = c(6, 16))
plot(scam_debug, select=5, ylim = c(-0.75, 2), xlim = c(6, 16))
plot(scam_debug, select=6, ylim = c(-0.75, 2), xlim = c(6, 16))
par(mfrow=c(1,1))
# probleme beim Einschränken des Definitionsbereichs von (6,16) bei SCAM.

par(mfrow=c(1,3))
plot(pam_debug, select=7, xlim=c(50,90), ylim=c(-2,2))
plot(pam_debug, select=8, xlim=c(50,90), ylim=c(-2,2))
plot(pam_debug, select=9, xlim=c(50,90), ylim=c(-2,2))
par(mfrow=c(1,1))

# compare design matrices
X_scam <- predict.scam(scam_debug, newdata = cal.my.mgus2.pam, type = "lpmatrix")
X_pam <- predict.gam(pam_debug, newdata = cal.my.mgus2.pam, type = "lpmatrix")

# prepare add_trans_prob pipeline
test_cal_pam_linear <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , age = seq(50, 90, by = 1)) %>% 
  group_by(transition, age) %>% 
  add_cumu_hazard(pam_debug_linear) 

test_cal_pam <- make_newdata(cal.my.mgus2.pam
                                    , tend = unique(tend)
                                    , transition=unique(transition)
                                    #, hgb = seq(7, 16, by = 0.5)
                                    , age = seq(50, 90, by = 1)) %>% 
  group_by(transition, age) %>% 
  add_cumu_hazard(pam_debug) 


test_cal_scam <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             , hgb = seq(8, 16, by = 1)) %>% 
  group_by(transition, hgb) %>% 
  add_cumu_hazard(scam_debug)

test_cal_comp <- test_cal_pam %>% 
  left_join(test_cal_scam, by=c("tstart", "tend", "transition", "hgb"), suffix=c("", "_scam")) %>%
  select(tend, transition, hgb, cumu_hazard, cumu_lower, cumu_upper, cumu_hazard_scam)
View(test_cal_comp)

ggplot(data=test_cal_comp, aes(x=tend, col=transition)) + 
  geom_line(aes(y=cumu_hazard), linetype = "solid") +
  geom_line(aes(y=cumu_hazard_scam), linetype = "dashed") +
  facet_wrap(~as.factor(hgb))

diff <- data.frame(cbind(exp(X_scam %*% coef(scam_debug)), exp(X_pam %*% coef(pam_debug))))
diff <- diff %>% mutate(diff_scam_pam = X1 - X2) %>% arrange(diff_scam_pam)
View(diff)
summary(diff$diff_scam_pam)
# differences come from different identifiability constraints



# SCAM

# workaround for grouped data -> include in add_trans_prob() when time
old_groups <- dplyr::groups(test_cal_scam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- test_cal_scam %>% ungroup(transition)
test <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

test$hgb <- as.factor(test$hgb)

# plot transitions

transition_ggplot <- ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=hgb)) + 
  facet_wrap(~transition + as.factor(age), ncol = 4, scales = "free_y") +
  #scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  #scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 100)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 
transition_ggplot
# ggsave("tmp/example/transition_probabilities.pdf", plot = transition_ggplot, width = 10)

# PAM
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

test$hgb <- as.factor(test$hgb)
test$age <- as.factor(test$age)

# PAM LINEAR
# workaround for grouped data -> include in add_trans_prob() when time
old_groups <- dplyr::groups(test_cal_pam_linear)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- test_cal_pam_linear %>% ungroup(transition)
test_linear <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

test$hgb <- as.factor(test$hgb)
test$age <- as.factor(test$age)

# plot transitions

transition_ggplot <- ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=hgb)) + 
  facet_wrap(transition~age, ncol = 4, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 100)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw()
transition_ggplot

transition_linear_ggplot <- ggplot(test_linear, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(group=age, col = age)) + 
  facet_wrap(~transition, ncol = 4, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 400)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw()
transition_linear_ggplot

transition_ggplot <- ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(group=age, col = age)) + 
  facet_wrap(~transition, ncol = 4, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 400)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw()
transition_ggplot

# ------------------------------------------------------------------------------
# HGB
# ------------------------------------------------------------------------------

# PAM

test_cal_pam <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , hgb = seq(6, 16, by = 1)) %>% 
  group_by(transition, hgb) %>% 
  add_cumu_hazard(pam_hgb) 

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

test$hgb <- as.factor(test$hgb)
test$age <- as.factor(test$age)

test_sub <- test %>% filter(transition == "0->2")

transition_ggplot <- ggplot(test_sub, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(group=hgb, col=hgb)) + 
  facet_wrap(~transition, ncol = 2, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 200)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw()
transition_ggplot

# PAM LINEAR

test_cal_lin <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , hgb = seq(6, 16, by = 1)) %>% 
  group_by(transition, hgb) %>% 
  add_cumu_hazard(pam_lin_hgb) 
# workaround for grouped data -> include in add_trans_prob() when time
old_groups <- dplyr::groups(test_cal_pam_linear)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- test_cal_lin %>% ungroup(transition)
test_lin <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

test_lin_sub <- test_lin %>% filter(transition == "0->2")

transition_ggplot_lin <- ggplot(test_lin_sub, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(group=hgb, col=hgb)) + 
  facet_wrap(~transition, ncol = 3, scales = "free_y", labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 200)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() 
transition_ggplot_lin

# combine plots
test_lin_sub <- test_lin_sub %>% mutate(model = "linear")
test_sub <- test_sub %>% mutate(model = "non-linear")

test <- rbind(test_lin_sub, test_sub)

test$model <- as.factor(test$model)

table(test$model)

combined_ggplot <- ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(group=hgb, col=hgb)) + 
  facet_wrap(transition~ model, ncol = 2, labeller = label_both) +
  # scale_color_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#33a02c"))+
  # scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  xlim(c(0, 200)) +
  ylab("Transition Probability") +
  xlab("time") +
  theme_bw() +
  theme(legend.position = "bottom")
combined_ggplot
ggsave("tmp/example/transition_probabilities_hgb.pdf", plot = combined_ggplot, width = 10)

pdf("tmp/example/hgb_spline.pdf", width = 10)
plot(pam_hgb, select=5, ylim=c(-0.5,1.5), xlim=c(6,16))
dev.off()
