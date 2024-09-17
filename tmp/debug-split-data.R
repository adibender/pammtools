
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
theme_set(theme_bw())
# devtools::install_github("thomasp85/patchwork")
# library(patchwork)
library(survival)
library(mgcv)
setwd("C:/Users/ra63liw/Documents/98_git/pammtools-multi-state/pammtools")
devtools::load_all()
source("sandbox/helpers-msm.R")


library(checkmate)

data("prothr", package = "mstate")
prothr_sub <- prothr %>% filter(Tstart != Tstop) %>% select(1:8) %>% group_by(trans, treat) %>% slice(1:2)
prothr_sub
ped_prothr <- as_ped_multistate(data = prothr_sub,
                     Surv(Tstart, Tstop, status) ~.,
                     max_time = 1000,
                     transition = "trans",
                     timescale  = "calendar"
) %>%
  select(1:9)
ped_prothr

# Bug?
ped_prothr %>% filter(id == 2) %>% arrange(tstart, tend, trans)

# investigate survival::survSplit

# Load necessary library
library(tibble)

# Create the dataset as a tibble
data <- tibble(
  id = c(1, 2, 7, 8),
  from = c(2, 2, 2, 2),
  to = c(3, 3, 3, 3),
  trans = c(4, 4, 4, 4),
  Tstart = c(0, 0, 0, 0),
  Tstop = c(151, 251, 202, 211),
  status = c(1, 0, 1, 0),
  treat = factor(c("Placebo", "Placebo", "Prednisone", "Prednisone"))
)

# Print the dataset
data

survSplit(formula = Surv(Tstart, Tstop, status) ~ . + trans, 
          data = data,
          cut = sort(unique(c(seq(0,1000, by = 100), 151, 202))))

get_cut.default(data, Surv(Tstart, Tstop, status) ~ . + trans)

outcome_vars <- get_lhs_vars(Surv(Tstart, Tstop, status) ~ . + trans)
event = 1L
unique(data[[outcome_vars[2]]][1L * (data[[outcome_vars[3]]]) == event])

# now fix max_time, throws error message if as_ped is used for big data set
devtools::load_all()

data("prothr", package = "mstate")
prothr <- prothr %>% filter(Tstart != Tstop)

ped_prothr <- as_ped(data = prothr_sub,
                                Surv(Tstart, Tstop, status) ~.,
                                max_time = 1000,
                                transition = "trans",
                                timescale  = "calendar"
) %>%
  select(1:9)
ped_prothr
# ==> works on small data set
ped_prothr <- as_ped(data = prothr,
                     Surv(Tstart, Tstop, status) ~.,
                     max_time = 1000,
                     transition = "trans",
                     timescale  = "calendar"
) %>%
  select(1:9)
ped_prothr

max(unique(ped_prothr$tend))


