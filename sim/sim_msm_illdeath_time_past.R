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
library(gridExtra)
library(latex2exp)

#-------------------------------------------------------------------------------
# SIMULATION STUDY
#
# illness-death
#
# x dependency for 1->2.      1->3,   non-linear
# t dependency for      2->3,         non-linear
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1
# define true functions for plotting
#-------------------------------------------------------------------------------

# t dependency for      2->3,         non-linear
# meeting 11.06.2024
# Möglichkeit 1 --> funktioniert und ändert die transition probabilities
f12 <- function(t, t0) {
  (dgamma(t, 8, 2) *8 - 0.2*(t0-3)**2) #dgamma(t0, 1, 1) *6
}

plot(seq(0,3,by=0.1), f12(seq(0,3,by=0.1), 0))

plot(seq(0,3,by=0.1), f12(0, seq(0,3,by=0.1)))

ft <- function(t, t0, transition) {
  (transition == "2->3")*f12(t, t0)
}

# x dependency for 1->2.      1->3,   non-linear
fx <- function(x, transition) {
  (transition == "1->2")*(0.5 + 0.25*x**3) + (transition == "1->3")*(0.4 - x**2)
}

form_stepone <- formula( ~ 0.5 + 0.25*x1**3| 0.4 - x1**2  )
# define hazard for transition 2->3
form_steptwo <- formula( ~ f12(t, x3))

#-------------------------------------------------------------------------------
# 2
# Initialize Simulation parameters and build simulation data frame
#-------------------------------------------------------------------------------

nSim = 50

n = 1000
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)
seq_x = seq(-2, 2, by = 0.1)
seq_t = seq(0, 3, by = 0.1)

for (i in 1:nSim) {
  sim_df_stepone <- sim_pexp_cr(form_stepone, data, cut) %>%
    mutate(from = 1
           , to = type + 1
           , transition = case_when(
             type == 1 ~ "1->2",
             type == 2 ~ "1->3",
             .default = "err"
           )
           , x3 = 0) %>%
    rename(tstart = t, tstop = time) %>%
    filter(status == 1)
  
  # second step
  # take all relevant subjects, i.e. subjects with transition 1->2, i.e. currently in state 2
  data_steptwo <- sim_df_stepone %>%
    filter(to == 2, status == 1) %>%
      mutate(x3 = tstop) %>%
    select(id
           , x1
           , x3
           , from
           , t = tstop)
  
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
  
  #-------------------------------------------------------------------------------
  # 3
  # plot intermediate results for comparison
  #-------------------------------------------------------------------------------
  if(nSim == 1) {
    sim_df_mvna <- sim_df %>% mutate(from = from -1
                                     , to = to-1
                                     , transition = paste0(from, "->", to)
                                     , entry = tstart
                                     , exit = tstop)
    tra <- matrix(FALSE, ncol = 3, nrow = 3)
    dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
    tra[1, 2:3] <- TRUE
    tra[2, 3] <- TRUE
    
    # print transition matrix
    tra
    
    my.nelaal <- mvna(sim_df_mvna, c("0", "1", "2"), tra, "cens")
    if (require(lattice)){
      xyplot(my.nelaal
             , strip=strip.custom(bg="white")
             , ylab="Cumulative Hazard"
             , lwd=2
             , xlim=c(0,3)
             , ylim=c(0,10))
    }
    
    # estimate transition probabilites
    etm.my.data<- etm(sim_df_mvna, c("0", "1", "2"), tra, "cens", s = 0)
    
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
  }
  
  sim_df <- sim_df %>% add_counterfactual_transitions()
  # 
  # test <- sim_df %>% mutate(time = ifelse(tstop > 2, 1, 0)) %>% filter(status == 1)
  # table(test$time, test$transition)
  
  # go on with pamm tools procedure
  cal_sim_df <- as_ped_multistate(
    data       = sim_df,
    formula    = Surv(tstart, tstop, status)~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")
  
  
  ctrl <- gam.control(trace = TRUE)
  bam_sim <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition)) 
                       + as.factor(transition)
                       + s(x1, by=as.factor(transition))
                       + s(x3, by=as.factor(transition))
                       , data = cal_sim_df
                       , family=poisson()
                       , offset=offset
                       , discrete = T
                       , method = "fREML"
                       , control = ctrl)
  
  if(i == 1) {
    x1_df <- cal_sim_df %>%
      make_newdata(x1 = seq_x
                   , transition = unique(transition)) %>%
      add_term(bam_sim, term = c("x1", "tend", "x3")) %>%
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
  
    tend_df <- make_newdata(cal_sim_df
                            , tend = seq_range(tend, 500)
                            , x3 = seq(0, 3, by=0.1)
                            , transition = c("2->3")
                            )%>%
      add_term(bam_sim, term = c("x1", "tend", "x3")) %>%
      mutate(true.value = ft(t=tend, t0=x3, transition=transition)
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
      make_newdata(x1 = seq_x
                   , transition = unique(transition)) %>%
      add_term(bam_sim, term = c("x1", "tend", "x3")) %>%
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
    
    tend_df_temp <- make_newdata(cal_sim_df
                   , tend = seq_range(tend, 500)
                   , x3 = seq(0, 3, by=0.1)
                   , transition = c("2->3")
      )%>%
      add_term(bam_sim, term = c("x1", "tend", "x3")) %>%
      mutate(true.value = ft(t=tend, t0=x3, transition=transition)
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

x1_df_new <- x1_df %>% mutate(iter = iter + 12+16)
tend_df_new <- tend_df %>% mutate(iter = iter + 12+16)

sim_x_compl <- rbind(sim_x_compl, x1_df_new)
sim_t_compl <- rbind(sim_t_compl, tend_df_new)

# sim_x_compl <- rbind(simulation_x, x1_df)
# sim_t_compl <- rbind(simulation_t, tend_df)



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
  mutate(transition = case_when(
           transition == "1->2" ~"0->1"
           , transition == "1->3" ~"0->2"
           , transition == "2->3" ~"1->2"
         )) %>%
  group_by(transition, tend, x3) %>% 
  summarise(fit = fit
            , iter = iter
            , true.value = true.value
            , mean_fit = mean(fit)
            , mean_ci_up = mean(ci_upper)
            , mean_ci_lo = mean(ci_lower)
  )

sim_x_compl_sub <- sim_x_compl_new %>% filter(transition != "1->2")
gg_x1 <- ggplot(sim_x_compl_sub, aes(x = x1)) +
  geom_line(aes(y = fit+0.2*(transition == "0->1"), group = iter), lwd = 0.5, color = "grey") +
  geom_line(aes(y = true.value), lwd = 1.1) +
  geom_line(aes(y = mean_fit+0.2*(transition == "0->1")), color="firebrick2", lwd = 1.1) +
  # geom_line(aes(y = mean_ci_up), color="firebrick2", lwd = 1.1 , linetype="dotted") +
  # geom_line(aes(y = mean_ci_lo), color="firebrick2", lwd = 1.1, linetype="dotted") +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = iter), alpha = .2) +
  facet_wrap(~ transition, labeller = label_both) +
  ylab("f(x)") +
  xlab("x") + coord_cartesian(ylim = c(-2,2)) +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
  ) 

gg_x1

sim_t_compl_sub <- sim_t_compl_new %>% filter(transition == "1->2", x3 %in% c(0.5, 1, 1.5)) %>% rename(t0 = x3)

gg_tend <- ggplot(sim_t_compl_sub, aes(x = tend)) +
  geom_line(aes(y = fit, group = iter, color = "fit"), lwd = 0.5) +
  stat_smooth(aes(y = mean_fit, color="mean"), lwd = 1.1) +
  geom_line(aes(y = true.value, color="true"), lwd=1.1) +
  # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .2) +
  coord_cartesian(ylim = c(-3, 3)
                  , xlim=c(0,2.5)) +
  facet_wrap(~ transition + t0, labeller = label_both) +
  ylab("f(t, t0)") +
  xlab("t") +
  scale_color_manual(name = "Log-Hazards",
                     values = c("fit" = "grey"
                                , "mean" = "firebrick2"
                                , "true" = "black"
                    ),
                    labels = unname(c("simulation fit"
                                      , "average fit" 
                                      , "f(t, t_0)"
                                      ))
  ) +
  theme_bw() +
  theme( strip.text = element_text(size = 20)
         , axis.text = element_text(size = 14)
         , axis.title = element_text(size = 20)
         , legend.text = element_text(size = 20)
         , legend.title = element_text(size = 20)
         ) 

  
  


gg_tend 

combined.plot.simulation <- grid.arrange(gg_x1, gg_tend, ncol=2, widths=c(2/6, 4/6))

ggsave("tmp/example/simulation_past_time.png"
       , plot = combined.plot.simulation
       , width = 35
       , height = 10
       , dpi = 300
       , units = "cm")
dev.off()

class(combined.plot.simulation)

# 
# test <- sim_t_compl_sub %>% filter(iter == 1)
# dim(test)
# 
# 
# gg_tend <- ggplot(sim_t_compl_sub, aes(x = x3)) +
#   geom_line(aes(y = fit, group = iter), lwd = 0.5, color = "grey") +
#   stat_smooth(aes(y = mean_fit), color="firebrick2", lwd = 1.1) +
#   geom_line(aes(y = true.value), lwd=1.1) +
#   geom_ribbon(aes(ymin = mean_ci_lo, ymax = mean_ci_up), alpha = .2) +
#   ylab("f(t)") +
#   xlab("t") + coord_cartesian(ylim = c(-4, 4), xlim=c(0,3)) +
#   theme_bw() +
#   theme( strip.text = element_text(size = 14)
#          , axis.text = element_text(size = 14)) +
#   facet_wrap(~ transition, labeller = label_both)
# 
# 
# gg_tend 
