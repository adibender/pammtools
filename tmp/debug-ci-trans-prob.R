# test simulation

newdata <- test_cal_pam
object <- pam_hgb
n_sim = 10
alpha = 0.05

X             <- predict.gam(object, newdata = newdata, type = "lpmatrix")
coefs         <- coef(object)

# newdata$cumu_hazard_man <- cumsum(newdata$intlen * exp(X %*% coefs))

newdata$hazard_man <- exp(X %*% coefs)
newdata <- newdata |> add_hazard(object, ci = FALSE)

View(newdata %>% filter(hazard != hazard_man)) # --> no difference

V             <- object$Vp
sim_coef_mat  <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
  exp(X %*% z))

newdata <- split(newdata, group_indices(newdata)) %>%
  map_dfr(get_sim_ci_transprob, object = object, nsim = n_sim)

sim_data <- split(newdata, group_indices(newdata)) %>%
  map_dfr(get_sim_ci_transprob, object = object, nsim = n_sim)

# copied from repo. adjust to fit transition prob purpose
get_sim_ci_transprob <- function(newdata, object, alpha = 0.05, nsim = 100L, ...) {
  
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp
  coefs <- coef(object)
  
  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    cumsum(newdata$intlen * exp(X %*% z)))
  
  # newdata$cumu_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  # newdata$cumu_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)
  
  sim_fit_mat
  
}

# X <- prep_X(object, newdata)
# coefs <- coef(object)
# newdata$hazard <- unname(drop(X %*% coefs))
# 
# V             <- object$Vp
# sim_coef_mat  <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
# sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
#   cumsum(newdata$intlen * exp(X %*% z)))

dim(sim_fit_mat)
dim(newdata)
head(sim_fit_mat) # per column one iteration.

test_add_trans_prob_ci <- apply(sim_fit_mat, 2, function(z) {
  newdata$cumu_hazard <- z
  old_groups <- dplyr::groups(newdata)
  res_data <- newdata %>% ungroup(transition)
  sim_trans_prob <- group_split(res_data) |> 
    map(res_data, .f = ~ group_by(.x, transition))|> 
    map(res_data, .f = ~ get_trans_prob(.x)) |>
    map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
    bind_rows() |>
    select(trans_prob)
  
  return(sim_trans_prob)
}) # per row, get one value.


# need different cumu hazards in newdata to calculate differnt transition probabilities


newdata <- newdata %>% add_cumu_hazard(object)

old_groups <- dplyr::groups(newdata)
res_data <- newdata %>% ungroup(transition)
out_data <- group_split(res_data) |> 
  map(res_data, .f = ~ group_by(.x, transition))|> 
  map(res_data, .f = ~ get_trans_prob(.x)) |>
  map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  bind_rows()


# for comparison of results
combined_contour <- ggplot(out_data, aes(x=tend, y=hgb, z=trans_prob)) +
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
