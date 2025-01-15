
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
