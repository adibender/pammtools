add_hazard_cr <- function(newdata, object, ci = TRUE, se_mult = 2,
                          ci_type = c("default", "delta", "sim"),
                          overwrite = FALSE, time_var = NULL, ...) {
  UseMethod("add_hazard_cr", object)
}

add_hazard_cr.pamm_cr <- function(newdata, object, 
                                  ci = TRUE, se_mult = 2, 
                                  ci_type = c("default", "delta", "sim"),
                                  overwrite = FALSE, time_var = NULL, ...) {
  dots <- list(...)
  dots$newdata <- newdata
  dots$object <- object
  dots$ci <- ci
  dots$se_mult <- se_mult
  dots$ci_type <- ci_type
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  do.call(add_hazard, dots)
}

add_hazard_cr.pamm_cr_list <- function(newdata, object, ci = TRUE, se_mult = 2, 
                                       ci_type = c("default", "delta", "sim"),
                                       overwrite = FALSE, time_var = NULL, ...) {
  dots <- list(...)
  dots$ci <- ci
  dots$se_mult <- se_mult
  dots$ci_type <- ci_type
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  data_list <- newdata
  res <- vector(mode = "list", length = length(data_list))
  for (i in 1:length(data_list)) {
    dots$object <- object[[i]]
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_hazard, dots)
    res[[i]]$cause <- names(object)[i]
  }
  res <- do.call(rbind, res)
  res$cause <- as.factor(res$cause)
  res
}

add_cumu_hazard_cr <- function(newdata, object, ci = TRUE, se_mult = 2,
                          ci_type = c("default", "delta", "sim"),
                          overwrite = FALSE, time_var = NULL, ...) {
  UseMethod("add_cumu_hazard_cr", object)
}

add_cumu_hazard_cr.pamm_cr <- function(newdata, object, ci = TRUE, se_mult = 2, 
                                       ci_type = c("default", "delta", "sim"),
                                       overwrite = FALSE, time_var = NULL, 
                                       ...) {
  dots <- list(...)
  dots$object <- object
  dots$ci <- ci
  dots$se_mult <- se_mult
  dots$ci_type <- ci_type
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  
  data_list <- newdata %>% split(newdata[["cause"]])
  res <- vector(mode = "list", length = length(data_list))
  for (i in 1:length(data_list)) {
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_cumu_hazard, dots)
  }
  do.call(rbind, res)
}

add_cumu_hazard_cr.pamm_cr_list <- function(newdata, object, ci = TRUE, se_mult = 2, 
                                            ci_type = c("default", "delta", "sim"),
                                            overwrite = FALSE, time_var = NULL, ...) {
  dots <- list(...)
  dots$ci <- ci
  dots$se_mult <- se_mult
  dots$ci_type <- ci_type
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  data_list <- newdata
  res <- vector(mode = "list", length = length(data_list))
  for (i in 1:length(data_list)) {
    dots$object <- object[[i]]
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_cumu_hazard, dots)
    res[[i]]$cause <- names(object)[i]
  }
  res <- do.call(rbind, res)
  res$cause <- as.factor(res$cause)
  res
}

add_surv_prob_cr <- function(newdata, object, ci = TRUE, se_mult = 2,
                               ci_type = c("default", "delta", "sim"),
                               overwrite = FALSE, time_var = NULL, ...) {
  UseMethod("add_surv_prob_cr", object)
}

add_surv_prob_cr.pamm_cr <- function(newdata, object, 
                                       ci = TRUE, se_mult = 2, 
                                       ci_type = c("default", "delta", "sim"),
                                       overwrite = FALSE, time_var = NULL, ...) {
  dots <- list(...)
  dots$object <- object
  dots$ci <- ci
  dots$se_mult <- se_mult
  dots$ci_type <- ci_type
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  
  data_list <- newdata %>% split(newdata[["cause"]])
  res <- vector(mode = "list", length = length(data_list))
  for (i in 1:length(data_list)) {
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_surv_prob, dots)
  }
  do.call(rbind, res)
}

add_surv_prob_cr.pamm_cr_list <- function(newdata, object, ci = TRUE, se_mult = 2, 
                                            ci_type = c("default", "delta", "sim"),
                                            overwrite = FALSE, time_var = NULL, ...) {
  dots <- list(...)
  dots$ci <- ci
  dots$se_mult <- se_mult
  dots$ci_type <- ci_type
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  data_list <- newdata
  res <- vector(mode = "list", length = length(data_list))
  for (i in 1:length(data_list)) {
    dots$object <- object[[i]]
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_surv_prob, dots)
    res[[i]]$cause <- names(object)[i]
  }
  res <- do.call(rbind, res)
  res$cause <- as.factor(res$cause)
  res
}

add_cif <- function(newdata, object, ci = TRUE,
                    overwrite = FALSE, time_var = NULL, alpha = 0.05,
                    n_sim = 500L, ...) {
  UseMethod("add_cif", object)
}

add_cif.pamm_cr <- function(newdata, object, 
                            ci = TRUE, overwrite = FALSE, time_var = NULL, 
                            alpha = 0.05, n_sim = 500L, ...) {
  dots <- list(...)
  dots$object <- object
  dots$ci <- ci
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  
  data_list <- newdata %>% split(newdata[["cause"]])
  cumu_hazards <- vector(mode = "list", length = length(data_list))
  hazards <- cumu_hazards
  cif_lower <- hazards
  cif_upper <- hazards
  cif_median <- hazards
  res <- hazards
  for (i in 1:length(data_list)) {
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_cumu_hazard, dots)
    X        <- predict.gam(object, newdata = data_list[[i]], type = "lpmatrix")
    V        <- object$Vp
    coefs    <- coef(object)
    sim_coef_mat <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
    hazards[[i]] <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))
    cumu_hazards[[i]] <- apply(hazards[[i]], 2, function(z) cumsum(z) * 
                                 (newdata$tend - newdata$tstart))
  }
  overall_survival <- exp( - Reduce(`+`, cumu_hazards[[i]]))
  lagged_overall_survival <- lag(overall_survival)
  lagged_overall_survival[1] <- 1
  averaged_overall_survival <- (overall_survival + lagged_overall_survival) / 2
  for (i in 1:length(hazards)) {
    current_cif <- hazards[[i]] * (newdata$tend - newdata$tstart) * averaged_overall_survival
    current_cif <- apply(current_cif, 2, cumsum)
    cif_median[[i]] <- apply(current_cif, 1, quantile, probs = 0.5)
    cif_lower[[i]] <- apply(current_cif, 1, quantile, probs = alpha / 2)
    cif_upper[[i]] <- apply(current_cif, 1, quantile, probs = 1 - alpha / 2)
  }
  message(cat("Prediction for the CIF is based on simulated values."))
  res <- do.call(rbind, res)
  cif <- do.call(c, cif_median)
  res[["cumu_hazard"]] <- cif
  if (!ci) {
    res %>% rename(cif = cumu_hazard)
  } else {
    cif_lower <- do.call(c, cif_lower)
    cif_upper <- do.call(c, cif_upper)
    res[, c("cumu_hazard", "cumu_lower", "cumu_upper")] <- 
      cbind(cif, cif_lower, cif_upper)
    res %>% rename(cif = cumu_hazard, 
                   cif_lower = cumu_lower, cif_upper = cumu_upper)
  }
}

add_cif.pamm_cr_list <- function(newdata, object, 
                                 ci = TRUE, overwrite = FALSE, time_var = NULL, 
                                 alpha = 0.05, n_sim = 500L, ...) {
  dots <- list(...)
  dots$ci <- ci
  dots$overwrite <- overwrite
  dots$time_var <- time_var
  
  data_list <- newdata
  cumu_hazards <- vector(mode = "list", length = length(data_list))
  hazards <- cumu_hazards
  cif_lower <- hazards
  cif_upper <- hazards
  cif_median <- hazards
  res <- hazards
  for (i in 1:length(data_list)) {
    dots$object <- object[[i]]
    dots$newdata <- data_list[[i]]
    res[[i]] <- do.call(add_cumu_hazard, dots)
    X        <- predict.gam(dots$object, newdata = data_list[[i]], type = "lpmatrix")
    V        <- dots$object$Vp
    coefs    <- coef(dots$object)
    sim_coef_mat <- mvtnorm::rmvnorm(n_sim, mean = coefs, sigma = V)
    hazards[[i]] <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))
    cumu_hazards[[i]] <- apply(hazards[[i]], 2, function(z) cumsum(z) * 
                                 (dots$newdata$tend - dots$newdata$tstart))
  }
  overall_survival <- exp( - Reduce(`+`, cumu_hazards[[i]]))
  lagged_overall_survival <- lag(overall_survival)
  lagged_overall_survival[1] <- 1
  averaged_overall_survival <- (overall_survival + lagged_overall_survival) / 2
  for (i in 1:length(hazards)) {
    current_cif <- hazards[[i]] * 
      (newdata[[i]]$tend - newdata[[i]]$tstart) * averaged_overall_survival
    current_cif <- apply(current_cif, 2, cumsum)
    cif_median[[i]] <- apply(current_cif, 1, quantile, probs = 0.5)
    cif_lower[[i]] <- apply(current_cif, 1, quantile, probs = alpha / 2)
    cif_upper[[i]] <- apply(current_cif, 1, quantile, probs = 1 - alpha / 2)
    res[[i]]$cause <- names(object)[i]
  }
  message(cat("Prediction for the CIF is based on simulated values."))
  res <- do.call(rbind, res)
  cif <- do.call(c, cif_median)
  res[["cumu_hazard"]] <- cif
  res$cause <- as.factor(res$cause)
  if (!ci) {
    res %>% rename(cif = cumu_hazard)
  } else {
    cif_lower <- do.call(c, cif_lower)
    cif_upper <- do.call(c, cif_upper)
    res[, c("cumu_hazard", "cumu_lower", "cumu_upper")] <- 
      cbind(cif, cif_lower, cif_upper)
    res %>% rename(cif = cumu_hazard, 
                   cif_lower = cumu_lower, cif_upper = cumu_upper)
  }
}

int_info.ped_cr_union <- function(x, ...) {

  ii <- int_info.ped(x, ...)
  ii <- ii[rep(seq_len(nrow(ii)), each = length(unique(x[["cause"]]))), ]
  ii$cause <- rep(unique(ped[["cause"]]), nrow(ii) / length(unique(x[["cause"]])))
  ii[order(ii[["cause"]], ii[["tend"]]), ]
  
}
