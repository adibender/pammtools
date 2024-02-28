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


my.mgus2.pam <- my.mgus2.pam %>% add_counterfactual_transitions()

cal.my.mgus2.pam <- as_ped_multistate(
  data       = my.mgus2.pam,
  formula    = Surv(tstart, tstop, status)~ .,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

pam <- mgcv::gam(ped_status ~ s(tend, by=as.factor(transition)) 
                 + as.factor(transition)
                 + s(mspike, by=as.factor(transition))
                 + s(hgb, by=as.factor(transition))
                 + s(age, by=as.factor(transition))
                 + sex
                 , data = cal.my.mgus2.pam
                 , family=poisson()
                 , offset=offset)

summary(pam)
plot(pam, xlim = c(0,20), ylim = c(-1, 1), page=1)
# 
# plot.gam(pam, select=4, xlim=c(0.5,2), ylim=c(-1,1))
# plot.gam(pam, select=5, xlim=c(0.5,2), ylim=c(-1,1))
# plot.gam(pam, select=6, xlim=c(0.5,2), ylim=c(-1,1))
# plot.gam(pam, select=7, xlim=c(10,16), ylim=c(-1,1))
# plot.gam(pam, select=8, xlim=c(10,16), ylim=c(-1,1))
# plot.gam(pam, select=9, xlim=c(10,16), ylim=c(-1,1))
# plot.gam(pam, select=10)
# plot.gam(pam, select=11)
# plot.gam(pam, select=12)


# hgb = quantile(....) --> Vermutung: aktuell noch nicht möglich, so wie funktion definiert ist

prep.cal.my.mgus2.pam <- make_newdata(cal.my.mgus2.pam
                                      , tend = unique(tend)
                                      , transition=unique(transition)
                                      , age = quantile(age, probs=c(0.05, 0.5, 0.95))
                                      ) %>% 
  group_by(transition, age) %>% 
  add_cumu_hazard(pam) 

old_groups <- dplyr::groups(prep.cal.my.mgus2.pam)
# transition is needed in the add_trans_prob because the transitions probabilities
# depend on each other
res_data <- prep.cal.my.mgus2.pam %>% ungroup(transition)
test <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition)) |> 
  map(res_data, .f = ~ add_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()

# plot transitions
ggplot(test, aes(x=tend, y=trans_prob)) + 
  geom_line(aes(col=as.factor(age))) + 
  facet_wrap(~transition, ncol = 1, scales="free_y") +
  xlim(c(0, 100))

#-------------------------------------------------------------------------------
# cox multistate
#-------------------------------------------------------------------------------
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

cfit1 <- coxph(Surv(etime, event=="pcm") ~ age + sex, mgus2)
cfit2 <- coxph(Surv(etime, event=="death") ~ age + sex, mgus2)

# predicted competing risk curves for a 72 year old with mspike of 1.2
# (median values), male and female.
# The survfit call is a bit faster without standard errors.
newdata <- expand.grid(sex=c("F", "M"), age=72, mspike=1.2)

AJmat <- matrix(list(), 3,3)
AJmat[1,2] <- list(survfit(cfit1, newdata, std.err=FALSE))
AJmat[1,3] <- list(survfit(cfit2, newdata, std.err=FALSE))
csurv  <- survfit(AJmat, p0 =c(entry=1, PCM=0, death=0))



