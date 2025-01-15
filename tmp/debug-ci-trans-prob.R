# test simulation

newdata <- test_cal_pam
object <- pam_hgb
n_sim = 10
alpha = 0.1

newdata <- newdata %>% add_trans_ci(pam_hgb, alpha = 0.1)



# copied from repo. adjust to fit transition prob purpose
get_sim_cumu <- function(newdata, ...) {
  
  newdata$cumu_hazard <- cumsum(newdata$intlen * newdata$hazard)
  
  newdata
  
}

add_trans_ci <- function(newdata, object, n_sim=100L, alpha=0.05, ...) {
  
  X             <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  coefs         <- coef(object)
  V             <- object$Vp
  
  sim_coef_mat <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    exp(X %*% z))
  
  # create list with replicated newdata
  nlst <- as.list(replicate(n_sim, newdata, simplify=F))

  # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
  nlst <- lapply(1:n_sim, function(i) {
    nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
    # split by group and calculate cumu hazard
    nlst[[i]] <- split(nlst[[i]], group_indices(nlst[[i]]))%>%
      map_dfr(get_sim_cumu)
    
    old_groups <- dplyr::groups(nlst[[i]])
    res_data <- nlst[[i]] %>% ungroup(transition)
    nlst[[i]] <- group_split(res_data) |> 
      map(res_data, .f = ~ group_by(.x, transition))|> 
      map(res_data, .f = ~ get_trans_prob(.x)) |>
      map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
      bind_rows()
    
    nlst[[i]]
  })
  
  sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
  newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
  newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)
  
  newdata
}

X             <- predict.gam(object, newdata = newdata, type = "lpmatrix")
coefs         <- coef(object)
V             <- object$Vp


sim_coef_mat <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
  exp(X %*% z))

# create list with replicated newdata
nlst <- as.list(replicate(n_sim, newdata, simplify=F))

# i = 1
# newdata <- nlst[[i]]
# newdata <- cbind(newdata, hazard = sim_fit_mat[, i]) # add hazard
# # split by group and calculate cumu hazard
# newdata <- split(newdata, group_indices(newdata))%>%
#   map_dfr(get_sim_cumu)
# 
# newdata <- newdata %>% rename(cumu_hazard_man = cumu_hazard)
# newdata <- newdata %>% add_cumu_hazard(pam_hgb)

# copied from repo. adjust to fit transition prob purpose
get_sim_cumu <- function(newdata, ...) {

  newdata$cumu_hazard <- cumsum(newdata$intlen * newdata$hazard)

  newdata

}

# # copied from repo. adjust to fit transition prob purpose
# get_sim_cumu <- function(newdata, alpha = 0.05, nsim = 1L, ...) {
#   
#   X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
#   V     <- object$Vp
#   coefs <- coef(object)
#   
#   sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
#   # sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
#   #   cumsum(newdata$intlen * exp(X %*% z)))
#   
#   # newdata$cumu_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
#   # newdata$cumu_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)
#   
#   newdata$cumu_hazard <- apply(sim_coef_mat, 1, function(z)
#     cumsum(newdata$intlen * exp(X %*% z)))
#   
#   newdata
#   
# }

# add cumu-hazard in each element and calculate trans_prob with perturbed hazards
nlst <- lapply(1:n_sim, function(i) {
  nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
  # split by group and calculate cumu hazard
  nlst[[i]] <- split(nlst[[i]], group_indices(nlst[[i]]))%>%
    map_dfr(get_sim_cumu)
  
  old_groups <- dplyr::groups(nlst[[i]])
  res_data <- nlst[[i]] %>% ungroup(transition)
  nlst[[i]] <- group_split(res_data) |> 
    map(res_data, .f = ~ group_by(.x, transition))|> 
    map(res_data, .f = ~ get_trans_prob(.x)) |>
    map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
    bind_rows()
  
  nlst[[i]]
})

sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)

View(newdata)

newdata <- newdata %>% add_trans_prob(pam_hgb)

# for comparison of results
combined_contour <- ggplot(newdata, aes(x=tend, y=hgb, z=trans_lower)) +
  geom_tile(aes(fill=trans_lower)) +
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

## testing

## CIs around tend spline
testdata <- newdata %>% filter(transition == "0->2"
                               , hgb == 10)



View(testdata)


ggplot(testdata, aes(x=tend, y=trans_prob)) + geom_line() + 
  geom_line(aes(y = trans_lower), linetype = "dashed") +
  geom_line(aes(y = trans_upper), linetype = "dashed")


## CIs around hgb spline
testdata <- newdata %>% filter(transition == "0->2"
                               , tend == 100)

ggplot(testdata, aes(x=hgb, y=trans_prob)) + geom_line() + 
  geom_line(aes(y = trans_lower), linetype = "dashed") +
  geom_line(aes(y = trans_upper), linetype = "dashed")

combined_contour <- ggplot(newdata, aes(x=tend, y=hgb, z=trans_lower)) +
  geom_tile(aes(fill=trans_lower)) +
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

## CIs for prothr data (website)

# test simulation

newdata <- ndf |> group_by(treat, transition) |> arrange(treat, transition, tend)
object <- pam
n_sim = 100
alpha = 0.1

newdata <- newdata %>% add_trans_ci(pam, alpha = 0.1)



# copied from repo. adjust to fit transition prob purpose
get_sim_cumu <- function(newdata, ...) {
  
  newdata$cumu_hazard <- cumsum(newdata$intlen * newdata$hazard)
  
  newdata
  
}

add_trans_ci <- function(newdata, object, n_sim=100L, alpha=0.05, ...) {
  
  X             <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  coefs         <- coef(object)
  V             <- object$Vp
  
  sim_coef_mat <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    exp(X %*% z))
  
  # create list with replicated newdata
  nlst <- as.list(replicate(n_sim, newdata, simplify=F))
  
  # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
  nlst <- lapply(1:n_sim, function(i) {
    nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
    # split by group and calculate cumu hazard
    nlst[[i]] <- split(nlst[[i]], group_indices(nlst[[i]]))%>%
      map_dfr(get_sim_cumu)
    
    old_groups <- dplyr::groups(nlst[[i]])
    res_data <- nlst[[i]] %>% ungroup(transition)
    nlst[[i]] <- group_split(res_data) |> 
      map(res_data, .f = ~ group_by(.x, transition))|> 
      map(res_data, .f = ~ get_trans_prob(.x)) |>
      map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
      bind_rows()
    
    nlst[[i]]
  })
  
  sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
  newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
  newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)
  
  newdata
}

X             <- predict.gam(object, newdata = newdata, type = "lpmatrix")
coefs         <- coef(object)
V             <- object$Vp

groups_array <- group_indices(newdata)

groups_trans <- newdata %>% ungroup(transition) %>% group_indices()


sim_coef_mat <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
  exp(X %*% z))

# create list with replicated newdata
nlst <- as.list(replicate(n_sim, newdata[,c("tend", "transition", "intlen")], simplify=F))

nlst <- lapply(1:n_sim, function(i) {
  nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
  # split by group and calculate cumu hazard
  nlst[[i]] <- split(nlst[[i]], groups_array) %>% map_dfr(get_sim_cumu) 
  
  nlst[[i]] <- split(nlst[[i]], groups_trans) %>% map_dfr(get_trans_prob)
  
  
  # 
  # old_groups <- dplyr::groups(nlst[[i]])
  # res_data <- nlst[[i]] %>% ungroup(transition)
  # nlst[[i]] <- group_split(res_data) |>
  #   map(res_data, .f = ~ group_by(.x, transition))|>
  #   map(res_data, .f = ~ get_trans_prob(.x)) |>
  #   map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
  #   bind_rows()
  
  nlst[[i]]
})

View(nlst[[1]])

# i = 1
# newdata <- nlst[[i]]
# newdata <- cbind(newdata, hazard = sim_fit_mat[, i]) # add hazard
# # split by group and calculate cumu hazard
# newdata <- split(newdata, group_indices(newdata))%>%
#   map_dfr(get_sim_cumu)
# 
# newdata <- newdata %>% rename(cumu_hazard_man = cumu_hazard)
# newdata <- newdata %>% add_cumu_hazard(pam_hgb)

# copied from repo. adjust to fit transition prob purpose
get_sim_cumu <- function(newdata, ...) {
  
  newdata$cumu_hazard <- cumsum(newdata$intlen * newdata$hazard)
  
  newdata
  
}

# # copied from repo. adjust to fit transition prob purpose
# get_sim_cumu <- function(newdata, alpha = 0.05, nsim = 1L, ...) {
#   
#   X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
#   V     <- object$Vp
#   coefs <- coef(object)
#   
#   sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
#   # sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
#   #   cumsum(newdata$intlen * exp(X %*% z)))
#   
#   # newdata$cumu_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
#   # newdata$cumu_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)
#   
#   newdata$cumu_hazard <- apply(sim_coef_mat, 1, function(z)
#     cumsum(newdata$intlen * exp(X %*% z)))
#   
#   newdata
#   
# }

# add cumu-hazard in each element and calculate trans_prob with perturbed hazards
nlst <- lapply(1:n_sim, function(i) {
  nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
  # split by group and calculate cumu hazard
  nlst[[i]] <- split(nlst[[i]], group_indices(nlst[[i]]))%>%
    map_dfr(get_sim_cumu)
  
  old_groups <- dplyr::groups(nlst[[i]])
  res_data <- nlst[[i]] %>% ungroup(transition)
  nlst[[i]] <- group_split(res_data) |>
    map(res_data, .f = ~ group_by(.x, transition))|>
    map(res_data, .f = ~ get_trans_prob(.x)) |>
    map(res_data, .f = ~ group_by(.x, !!!old_groups)) |>
    bind_rows()
  
  nlst[[i]]
})

sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)

View(newdata)

newdata <- newdata |> add_trans_prob(pam)

ggplot(newdata, aes(x=tend)) + 
  geom_line(aes(y=trans_prob, col=treat), linetype= "solid") +
  geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper, fill=treat), alpha = .3) +
  facet_wrap(~transition) +
  labs(y = "Transition Probability", x = "time", color = "Treatment", fill= "Treatment")
