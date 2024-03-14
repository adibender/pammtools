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
library(scam)
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


test_cal_pam <- make_newdata(cal.my.mgus2.pam
                             , tend = unique(tend)
                             , transition=unique(transition)
                             #, hgb = seq(7, 16, by = 0.5)
                             , hgb = seq(6, 16, by = 0.1)) %>% 
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

plot(pam_hgb, select=5, ylim=c(-0.5,1.5), xlim=c(6,16))

ggplot(test_sub, aes(x=tend, y=hgb, z=trans_prob)) +
  geom_contour_filled() +
  xlim(c(0,100))


ggplot(test_sub, aes(x=tend, y=hgb, z=trans_prob)) +
  geom_tile(aes(fill=trans_prob)) +
  scale_fill_gradient2(low  = "steelblue", high = "firebrick2", midpoint=0.5)+
  stat_contour(col="grey30") +
  xlim(c(0,100))

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

ggplot(test_lin_sub, aes(x=tend, y=hgb, z=trans_prob)) +
  geom_tile(aes(fill=trans_prob)) +
  scale_fill_gradient2(
      low  = "steelblue"
    , high = "firebrick2"
    , midpoint=0.5)+
  stat_contour(col="grey30") +
  xlim(c(0,100))

# combine plots
test_lin_sub <- test_lin_sub %>% mutate(model = "linear")
test_sub <- test_sub %>% mutate(model = "non-linear")

test <- rbind(test_lin_sub, test_sub)

test$model <- as.factor(test$model)

table(test$model)

combined_contour <- ggplot(test, aes(x=tend, y=hgb, z=trans_prob)) +
  geom_tile(aes(fill=trans_prob)) +
  scale_fill_gradient2(
    name = "probability"
    , low  = "steelblue"
    , high = "firebrick2"
    , midpoint=0.5)+
  stat_contour(col="grey30",lwd = 1.1) +
  # geom_vline(xintercept = c(25, 75), lty = 3) +
  # geom_hline(yintercept = c(8, 10, 12, 14), lty = 3) +
  facet_wrap(transition ~ model, ncol=2, labeller = label_both) +
  xlim(c(0,100)) +
  theme_bw() +
  theme(legend.position = "right"
        , strip.text = element_text(size = 14)
        , axis.text = element_text(size = 14))

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

library(RColorBrewer)
library(gridExtra)
Greens  <- RColorBrewer::brewer.pal(9, "Greens")
Purples <- RColorBrewer::brewer.pal(9, "Purples")

hgb_pp <- ggplot(hgb_df, aes(x = hgb)) +
  geom_line(aes(y = fit), lwd = 1.1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = .2) +
  ylab("s(hgb, 4.65):transition0->2") +
  xlab("hgb") + coord_cartesian(ylim = c(-0.5, 2.5)) +
  theme_bw() +
  theme( strip.text = element_text(size = 14)
         , axis.text = element_text(size = 14)) +
  facet_wrap(transition ~ model, labeller = label_both)

hgb_pp

time_pp <- ggplot(time_df, aes(x = tend)) +
  geom_line(aes(y = fit), lwd = 1.1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = .2) +
  ylab("s(tend, 7.7):transition0->2") +
  xlab("tend") + coord_cartesian(ylim = c(-0.5, 1)) +
  xlim(c(0,100)) +
  theme_bw() +
  theme( strip.text = element_text(size = 14)
         , axis.text = element_text(size = 14)) +
  facet_wrap(transition ~ model, labeller = label_both)

time_pp

combined_spline_contour <- grid.arrange(hgb_pp, time_pp, combined_contour, ncol=3, widths=c(1/4, 1/4, 2/4))
combined_spline_contour
ggsave("tmp/example/transition_probabilities_hgb_contour.pdf", plot = combined_spline_contour, width = 16)

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
    )
head(ndata)

str(ndata)

perc <- 0.5
ndata_sample <- sample_n(ndata, round(dim(ndata)[1]*perc,0))
dim(ndata_sample)


table(ndata$status)
ndata_pamm <- ndata %>% add_counterfactual_transitions()
head(ndata_pamm)
dim(ndata_pamm)

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
                       + male
                       + s(bmi, by=as.factor(transition))
                       + height*male
                       + height
                       + age
                       , data = cal_ndata_pamm
                       , family=poisson()
                       , offset=offset
                       , control = ctrl)

summary(pam_ndata)


new_ndata_pam <- make_newdata(cal_ndata_pamm
                              , tend = unique(tend)
                              , transition=unique(transition)
                              #, age = quantile(age, probs=c(0.05, 0.5, 0.95))
                              ) %>% 
  group_by(transition, age) %>% 
  add_cumu_hazard(pam_ndata) 

old_groups <- dplyr::groups(new_ndata_pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_ndata <- new_ndata_pam %>% ungroup(transition)
test <- group_split(res_ndata) |> 
  map(res_ndata, .f = ~ group_by(.x, transition)) |> 
  map(res_ndata, .f = ~ add_trans_prob(.x)) |>
  map(res_ndata, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot transitions
ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=as.factor(age))) + 
  facet_wrap(~transition, ncol = 1) +
  xlim(c(30, 80)) +
  ylim(c(0,0.05))





























