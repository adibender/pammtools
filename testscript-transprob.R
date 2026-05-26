data(prothr, package = "mstate")

setwd("C:/Users/ra63liw/Documents/98_git/pammtools")
devtools::load_all("C:/Users/ra63liw/Documents/98_git/pammtools")

# add time spent in abnormal prothrombin state
prothr_with_cov <- prothr %>%
  arrange(id, Tstart) %>%        # ensure chronological order
  group_by(id) %>%
  mutate(
    duration = Tstop - Tstart,
    time_apl =
      cumsum(duration * (from == 2)) - duration * (from == 2)
  ) %>%
  ungroup()

prothr |> head()
my.prothr <- prothr_with_cov |> 
  mutate(transition = as.factor(paste0(from, "->", to))) |>
  filter(Tstart != Tstop) # filter instantaneous transitions

ped <- pammtools::as_ped(
  data       = my.prothr,
  formula    = Surv(Tstart, Tstop, status)~ .,
  transition = "trans",
  id         = "id",
  timescale  = "calendar"
)

pam <- mgcv::bam(
  ped_status ~ s(tend, by=transition) + treat*transition,
  data = ped,
  offset = offset,
  family = poisson(),
  method = "fREML"
  )
summary(pam)

pam_int <- mgcv::bam(
  ped_status ~ s(tend, by=interaction(treat, transition)) + treat*transition,
  data = ped,
  offset = offset,
  family = poisson(),
  method = "fREML",
  discrete = TRUE
)
summary(pam_int)

ndf = make_newdata(
  ped, tend=unique(tend), transition=levels(transition), treat=levels(treat)
  ) |>
  group_by(treat, transition) |>
  add_cumu_hazard(pam_int)
  
ggplot(ndf, aes(x=tend, y=cumu_hazard)) + 
  geom_stephazard(aes(color=treat)) + 
  geom_stepribbon(aes(ymin = cumu_lower, ymax = cumu_upper, fill=treat), alpha = 0.1) +
  facet_wrap(~transition)

ndf <- pammtools::make_newdata(ped
                               , tend  = unique(tend)
                               , treat  = unique(treat)
                               , transition = unique(transition)) |>
  # filter(tend > 200) |>
  dplyr::group_by(treat, transition)

ndf_tp <- ndf |>
  dplyr::arrange(treat, transition, tend) |>
  add_trans_prob(pam, ci = FALSE)

ggplot(ndf_tp, aes(x=tend, y=trans_prob, col = treat)) + geom_line() + facet_wrap(~transition)



p = attributes(ndf_tp)$matrix |> dplyr::filter(treat == "Placebo") |> pull(trans_prob_matrix) # matrix for every time point

# state occupation
v <- c(0,1,0,0)
x <- p[[1]]  # 4 x 4 x 741
res <- apply(x, 3, function(mat) v %*% mat)

df <- as.data.frame(t(res))
colnames(df) <- paste0("state_", seq_len(ncol(df)))
df$time <- seq_len(nrow(df))

df_long <- df |>
  dplyr::mutate(time = seq_len(n())) |>
  tidyr::pivot_longer(
    cols = starts_with("state_"),
    names_to = "state",
    values_to = "prob"
  )

ggplot(df_long, aes(x=time, y = prob, fill=state)) + 
  geom_area(color = "black")

gg_state_occupation(ndf_tp, group_var = "treat", init_state = c(1,0,0))

## error testing no groups:
ndf <- pammtools::make_newdata(ped
                               , tend  = unique(tend)
                               , transition = unique(transition)) |>
  dplyr::group_by(transition)

ndf_tp <- ndf |>
  dplyr::arrange(transition, tend) |>
  add_trans_prob(pam, ci = FALSE)
attributes(ndf_tp)

gg_state_occupation(ndf_tp, init_state = c(1,0,0))

Rprof()
ndf_msm_treat <- make_newdata(
  ped, 
  tend       = unique(tend), 
  transition = unique(transition),
  treat      = unique(treat)
) |>
  group_by(treat, transition) |>
  #arrange(treat, transition, tend) |>
  add_trans_prob(pam_int, ci = TRUE)
Rprof(NULL)
summaryRprof()
ggplot(ndf_msm_treat, aes(x=tend, y=trans_prob, col = treat)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper, fill = treat), alpha = 0.1) +
  facet_wrap(~transition)
gg_state_occupation(ndf_msm_treat, group_var = "treat", init_state = c(1,0,0))
