#' Add (cumulative) hazards / survival probabilities for competing risks
#' 
#' This function is the main subroutine of add_cumu_hazard_cr(),
#' add_hazard_cr() and add_surv_prob().
#' For details see add_cumu_hazard_cr() and add_hazard_cr().
#' @author Philipp Kopper
hazard_adder_cr <- function(newdata, object, hazard_function, ci, se_mult, 
                            ci_type = c("default", "delta", "sim"), 
                            overwrite, time_var, result = "df", invert = FALSE, 
                            ...) {
  measure <- vector(mode = "list", length = length(object))
  for (i in 1:length(object)) {
    measure[[i]] <- hazard_function(newdata, object[[i]], ci = ci, 
                                    se_mult = se_mult, 
                                    ci_type = ci_type, 
                                    overwrite = overwrite, 
                                    time_var = time_var, ...)
    added_cols <- !(colnames(measure[[i]]) %in% colnames(newdata))
    measure[[i]] <- measure[[i]][ , added_cols, drop = FALSE]
    attr(measure, "risks")[i] <- attr(object, "risks")[i]
    colnames(measure[[i]]) <- paste(attr(object, "risks")[i], 
                                    colnames(measure[[i]]), sep = "_")
  }
  if (result == "df") {
    new_cols <- Reduce(cbind, measure)
    if (invert) {
      new_cols <- 1 - new_cols
      colnames(new_cols) <- gsub("surv", "cif", colnames(new_cols))
      colnames(new_cols) <- gsub("_prob", "", colnames(new_cols))
    }
    return(cbind(newdata, new_cols))
  } else {
    return(list(newdata = newdata, measure = measure))
  }
}

#' Add predicted (cumulative) hazard to data set for competing risks
#' 
#' Add (cumulative) hazard based on the provided data set and model for 
#' competing risks. 
#' If ci=TRUE confidence intervals are also added. Their width can be 
#' controlled via the se_mult argument. This is a wrapper around predict.gam.
#' Additionally, it is a wrapper for the non cr functions predicting the 
#' (cumulative) hazard rate.
#' @param newdata A data frame or list containing the values of the model
#'  covariates at which predictions are required. If this is not provided then 
#'  predictions corresponding to the original data are returned. If newdata is
#'  provided then it should contain all the variables needed for prediction: a
#'  warning is generated if not. See details for use with 
#'  link{linear.functional.terms}.
#' @param object a fitted pem_cr or pam_cr object as produced by gam_cr() or 
#' glm_cr().
#' @param ci Logical indicating whether to include confidence intervals. 
#'  Defaults to TRUE.
#' @param se_mult Factor by which standard errors are multiplied for 
#' calculating the confidence intervals.
#' @ci_type The method by which standard errors/confidence intervals will be
#'  calculated. Default transforms the linear predictor at respective intervals.
#'  "delta" calculates CIs based on the standard error calculated by the Delta 
#'  method. "sim" draws the property of interest from its posterior based on 
#'  the normal distribution of the estimated coefficients. CIs are given by 
#'  respective quantiles.
#' @param overwrite Should hazard columns be overwritten if already present in
#'  the data set? Defaults to FALSE. If TRUE, columns with names 
#'  c("hazard", "se", "lower", "upper") will be overwritten.
#' @param time_var Name of the variable used for the baseline hazard. If not
#'  given, defaults to "tend" for gam fits, else "interval". The latter is 
#'  assumed to be a factor, the former numeric.
#' @param ... Further arguments passed to add_hazard / add_cumu_hazard, 
#' predict.gam and get_hazard.
#' @param interval_length The variable in newdata containing the interval 
#' lengths. Can be either bare unquoted variable name or character. Defaults 
#' to "intlen".
#' @return a data.frame (or tibble) containing the original data and new columns
#' featuring the (cumulative) hazard (and if ci = TRUE the respective 
#' confidence intervals) for all risks.
#' @export
#' @author Philipp Kopper
add_cumu_hazard_cr <- function(newdata, object, 
                               ci = TRUE, se_mult = 2, 
                               ci_type = c("default", "delta", "sim"),
                               overwrite = FALSE, time_var = NULL, ...) {
  hazard_adder_cr(newdata, object, hazard_function = add_cumu_hazard, ci, 
                  se_mult, ci_type, overwrite, time_var, ...)
}

#' Add predicted (cumulative) hazard to data set for competing risks
#' 
#' Add (cumulative) hazard based on the provided data set and model for 
#' competing risks. 
#' If ci=TRUE confidence intervals are also added. Their width can be 
#' controlled via the se_mult argument. This is a wrapper around predict.gam.
#' Additionally, it is a wrapper for the non cr functions predicting the 
#' (cumulative) hazard rate.
#' @param newdata A data frame or list containing the values of the model
#'  covariates at which predictions are required. If this is not provided then 
#'  predictions corresponding to the original data are returned. If newdata is
#'  provided then it should contain all the variables needed for prediction: a
#'  warning is generated if not. See details for use with 
#'  link{linear.functional.terms}.
#' @param object a fitted pem_cr or pam_cr object as produced by gam_cr() or 
#' glm_cr().
#' @param type When this has the value "link" (default) the linear predictor
#'  (possibly with associated standard errors) is returned. When type="terms"
#'  each component of the linear predictor is returned seperately 
#'  (possibly with standard errors): this includes parametric model components, 
#'  followed by each smooth component, but excludes any offset and any 
#'  intercept. type="iterms" is the same, except that any standard errors 
#'  returned for smooth components will include the uncertainty about the 
#'  intercept/overall mean. When type="response" predictions on the scale of 
#'  the response are returned (possibly with approximate standard errors). 
#'  When type="lpmatrix" then a matrix is returned which yields the values of 
#'  the linear predictor (minus any offset) when postmultiplied by the 
#'  parameter vector (in this case se.fit is ignored). The latter option is 
#'  most useful for getting variance estimates for quantities derived from 
#'  the model: for example integrated quantities, or derivatives of smooths. A 
#'  linear predictor matrix can also be used to implement approximate 
#'  prediction outside R (see example code, below).
#' @param ci Logical indicating whether to include confidence intervals. 
#'  Defaults to TRUE.
#' @param se_mult Factor by which standard errors are multiplied for 
#' calculating the confidence intervals.
#' @ci_type The method by which standard errors/confidence intervals will be
#'  calculated. Default transforms the linear predictor at respective intervals.
#'  "delta" calculates CIs based on the standard error calculated by the Delta 
#'  method. "sim" draws the property of interest from its posterior based on 
#'  the normal distribution of the estimated coefficients. CIs are given by 
#'  respective quantiles.
#' @param overwrite Should hazard columns be overwritten if already present in
#'  the data set? Defaults to FALSE. If TRUE, columns with names 
#'  c("hazard", "se", "lower", "upper") will be overwritten.
#' @param time_var Name of the variable used for the baseline hazard. If not
#'  given, defaults to "tend" for gam fits, else "interval". The latter is 
#'  assumed to be a factor, the former numeric.
#' @param ... Further arguments passed to add_hazard / add_cumu_hazard, 
#' predict.gam and get_hazard.
#' @param interval_length The variable in newdata containing the interval 
#' lengths. Can be either bare unquoted variable name or character. Defaults 
#' to "intlen".
#' @return a data.frame (or tibble) containing the original data and new columns
#' featuring the (cumulative) hazard (and if ci = TRUE the respective 
#' confidence intervals) for all risks.
#' @export
#' @author Philipp Kopper
add_hazard_cr <- function(newdata, object, type = c("link", "response"), 
                          ci = TRUE, se_mult = 2, 
                          ci_type = c("default", "delta", "sim"),
                          overwrite = FALSE, time_var = NULL, ...) {
  hazard_adder_cr(newdata, object, hazard_function = add_hazard, type, ci, 
                  se_mult, ci_type, overwrite, time_var, ...)
}

#' Add survival probability estimates for competing risks
#' 
#' Given suitable data (i.e. data with all columns used for estimation of 
#' the model), this functions adds a column surv_prob containing survival 
#' probabilities for the specified covariate and follow-up information (and 
#' CIs surv_lower, surv_upper if ci=TRUE).
#' @param newdata A data frame or list containing the values of the model 
#' covariates at which predictions are required. If this is not provided 
#' then predictions corresponding to the original data are returned. If 
#' newdata is provided then it should contain all the variables needed 
#' for prediction: a warning is generated if not. See details for use 
#' with link{linear.functional.terms}.
#' @param object 	a fitted pem_cr / pam_cr object as produced by gam_cr() 
#' or glm_cr().
#' @param ci Logical indicating whether to include confidence intervals. 
#' Defaults to TRUE.
#' @param se_mult Factor by which standard errors are multiplied for 
#' calculating the confidence intervals.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to FALSE. If TRUE, columns with names c("hazard", 
#' "se", "lower", "upper") will be overwritten.
#' @param time_var Name of the variable used for the baseline hazard. If not
#' given, defaults to "tend" for gam fits, else "interval". The latter is 
#' assumed to be a factor, the former numeric.
#' @param interval_length The variable in newdata containing the interval 
#' lengths. Can be either bare unquoted variable name or character. Defaults 
#' to "intlen".
#' @param ... Further arguments passed to add_surv_prob, predict.gam and 
#' get_hazard
#' @return a data.frame (or tibble) containing the original data and new columns
#' featuring the predicted survival probability (and if ci = TRUE the 
#' respective confidence intervals) for all risks.
#' @export
#' @author Philipp Kopper
add_surv_prob_cr <- function(newdata, object, 
                             ci = TRUE, se_mult = 2, 
                             ci_type = c("default", "delta", "sim"),
                             overwrite = FALSE, time_var = NULL, ...) {
  hazard_adder_cr(newdata, object, hazard_function = add_surv_prob, ci, 
                  se_mult, ci_type, overwrite, time_var, ...)
}

#' Add cumulative incidence estimates for competing risks
#' 
#' Given suitable data (i.e. data with all columns used for estimation of 
#' the model), this functions adds a column cif containing 
#' cumulative subdistribution hazards for the specified covariate and follow-up 
#' information (and CIs cif_lower, cif_upper if ci=TRUE).
#' @param newdata A data frame or list containing the values of the model 
#' covariates at which predictions are required. If this is not provided 
#' then predictions corresponding to the original data are returned. If 
#' newdata is provided then it should contain all the variables needed 
#' for prediction: a warning is generated if not. See details for use 
#' with link{linear.functional.terms}.
#' @param object 	a fitted pem_cr / pam_cr object as produced by gam_cr() 
#' or glm_cr().
#' @param ci Logical indicating whether to include confidence intervals. 
#' Defaults to TRUE.
#' @param se_mult Factor by which standard errors are multiplied for 
#' calculating the confidence intervals.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to FALSE. If TRUE, columns with names c("hazard", 
#' "se", "lower", "upper") will be overwritten.
#' @param time_var Name of the variable used for the baseline hazard. If not
#' given, defaults to "tend" for gam fits, else "interval". The latter is 
#' assumed to be a factor, the former numeric.
#' @param interval_length The variable in newdata containing the interval 
#' lengths. Can be either bare unquoted variable name or character. Defaults 
#' to "intlen".
#' @param ... Further arguments passed to add_cumu_hazard, predict.gam and 
#' get_hazard
#' @return a data.frame (or tibble) containing the original data and new columns
#' featuring the predicted cumulatice incidence for all risks. 
#' (and if ci = TRUE the respective confidence intervals).
#' @export
#' @author Philipp Kopper
#' add_cif(pam_cs, ped_cs, ci = FALSE)
add_cif <- function(newdata, model, ped, ci = TRUE, ...) {
  UseMethod("add_cif", object = model)
}

add_cif.default <- function(newdata, model, ped, ci = TRUE, ...) {
  add_cif.cs(newdata, model, ped, ci, ...)
}

add_cif.sh <- function(newdata, model, ped = NULL, ci = TRUE, se_mult = 2, 
                       ci_type = "default", overwrite = FALSE, 
                       time_var = NULL, ...) {
  res <- hazard_adder_cr(newdata, model, hazard_function = add_surv_prob, ci, 
                         se_mult, ci_type, overwrite, time_var, invert = TRUE, ...)
}

add_cif.cs <- function(newdata, model, ped, 
                    ci = TRUE, alpha = 0.05,
                    overwrite = FALSE, time_var = NULL, keep = TRUE,
                    ...) {
  hazards <- hazard_adder_cr(newdata, model, hazard_function = add_hazard, 
                             type = "response", ci = FALSE, overwrite = overwrite, 
                             time_var = time_var, ...)
  hazards <- hazards[, (ncol(newdata) + 1): (ncol(hazards))]
  cumu_hazards <- hazard_adder_cr(newdata, model, hazard_function = 
                                    add_cumu_hazard, ci = FALSE, 
                                  overwrite = overwrite, 
                                  time_var = time_var, ...) 
  cumu_hazards <- cumu_hazards[, (ncol(newdata) + 1): (ncol(cumu_hazards))]
  overall_survival <- exp( - apply(cumu_hazards, 1, sum))
  lagged_overall_survival <- lag(overall_survival)
  lagged_overall_survival[1] <- 1
  cif <- vector(mode = "list", length = length(model))
  for (i in 1:length(cif)) {
    cif[[i]] <- cumsum(hazards[, i]  * newdata$intlen * 
                         apply(cbind(overall_survival, lagged_overall_survival),
                               1, mean))
  }
    add_this <- add_cif_ci(cif = cif, objects = model, ci = ci, 
                           newdata = newdata, ...)
  if (keep) {
    cbind(newdata, Reduce(cbind, add_this))
  } else {
    Reduce(cbind, add_this)
  }
}

add_cif_ci <- function(cif, objects, newdata, ci, alpha = 0.05, nsim = 500L) {
  if (ci) {
    hazards <- vector(mode = "list", length = length(objects))
    cumu_hazards <- hazards
    ci_lower <- hazards
    ci_upper <- hazards
    cif_median <- hazards
    for (i in 1:length(objects)) {
      X     <- predict.gam(objects[[i]], newdata = newdata, type = "lpmatrix")
      V     <- objects[[i]]$Vp
      coefs <- coef(objects[[i]])
      sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
      hazards[[i]] <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))
      cumu_hazards[[i]] <- apply(hazards[[i]], 2, function(z) cumsum(z) * 
                                   (newdata$tend - newdata$tstart))
    }
    overall_survival <- exp( - Reduce("+", cumu_hazards))
    lagged_overall_survival <- apply(overall_survival, 2, lag)
    lagged_overall_survival[1, ] <- 1
    averaged_overall_survival <- (overall_survival + lagged_overall_survival) / 2
    for (i in 1:length(objects)) {
      current_cif <- hazards[[i]] * (newdata$tend - newdata$tstart) * averaged_overall_survival
      current_cif <- apply(current_cif, 2, cumsum)
      cif_median[[i]] <- apply(current_cif, 1, quantile, probs = 0.5)
      ci_lower[[i]] <- apply(current_cif, 1, quantile, probs = alpha / 2)
      ci_upper[[i]] <- apply(current_cif, 1, quantile, probs = 1 - alpha / 2)
    }
  }
  add_this <- vector(mode = "list", length = length(objects))
  for (i in 1:length(objects)) {
    if (ci) {
      add_this[[i]] <- cbind(cif_median[[i]], 
                             ci_lower[[i]],
                             ci_upper[[i]])
      add_name <- paste("cif", c("", "_lower", "_upper"), sep = "")
      colnames(add_this[[i]]) <- paste(attr(objects, "risks")[i], 
                                       add_name, sep = "_")
    } else {
      add_this[[i]] <- as.data.frame(cif[[i]])
      colnames(add_this[[i]]) <- paste(attr(objects, "risks")[i], 
                                       "cif", sep = "_")
    }
  }
  return(add_this)
}




#' Helper function for the final output of the CFI columns
#' 
#' This function outputs the CFIs for all risks and if necessary the CIs.
#' @author Philipp Kopper 
add_cif_columns <- function(object, cif, ci, se_mult, cif_var = NULL) {
  add_this <- vector(mode = "list", length = length(object))
  for (i in 1:length(object)) {
    if (ci) {
      add_this[[i]] <- cbind(cif[[i]], 
                             cif[[i]] - se_mult * sqrt(cif_var[[i]]),
                             cif[[i]] + se_mult * sqrt(cif_var[[i]]))
      add_name <- paste("cif", c("", "_lower", "_upper"), sep = "")
      colnames(add_this[[i]]) <- paste(attr(object, "risks")[i], 
                                       add_name, sep = "_")
    } else {
      add_this[[i]] <- as.data.frame(cif[[i]])
      colnames(add_this[[i]]) <- paste(attr(object, "risks")[i], 
                                       "cif", sep = "_")
    }
  }
  return(add_this)
}


add_overall_surv <- function(newdata, object, data, 
                             ci = TRUE, se_mult = 2,
                             overwrite = FALSE, time_var = NULL, ...) {
  cif <- add_cif(newdata, object, data, ci, se_mult, overwrite, time_var, 
                 keep = FALSE, ...)
  
  if (ci) {
    add_this <- as.data.frame(matrix(0, ncol = 3, nrow = nrow(newdata)))
    colnames(add_this) <- c("overall_surv", 
                            "overall_surv_lower", "overall_surv_upper")
    n_risks <- length(attr(object, "risks"))
    for (i in 1:n_risks) {
      add_this[, 1] <-  rowSums(cif[, seq(1, n_risks * 3, by = 3)])
      add_this[, 2] <-  rowSums(cif[, seq(2, n_risks * 3, by = 3)])
      add_this[, 3] <-  rowSums(cif[, seq(3, n_risks * 3, by = 3)])
    }
    data.frame(newdata, 1 - add_this)
  } else {
    data.frame(newdata, overall_surv = 1 - rowSums(cif))
  }
}

