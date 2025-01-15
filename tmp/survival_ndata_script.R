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


# instead of string names, name states 0,1,...,4, where 4 = death
ndata$from <- as.numeric(with(ndata, factor(diabetes + htn + lipid, 0:3,
                                            c(0, 1, 2, 3))))

ndata$to <- as.numeric(factor(temp, 0:4,
                              c(0, 1, 2, 3, 4)))
# censoring is coded with 0 in status
ndata <- ndata %>% 
  mutate(status = case_when(
    as.character(event) == "censored" ~ 0,
    TRUE ~ 1)
  )
head(ndata)

# perc <- 0.5
# ndata_sample <- sample_n(ndata, round(dim(ndata)[1]*perc,0))
# dim(ndata_sample)

#check number of events
table(ndata$status)

# add all possible transitions
ndata_pamm <- ndata %>% add_counterfactual_transitions()
head(ndata_pamm)
dim(ndata_pamm)

# build as ped data set in calendar format
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

# excessive runtime
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