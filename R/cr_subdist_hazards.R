#' This is a wrapper for the general data transformation function provided by 
#' the pammtools package (as_ped). The usage is identical to as_ped (except of 
#' one issue: the status can be a character or factor. It will be memorised and 
#' used later on in the methods. This is helpful when interpreting result.
#' For more details use ?as_ped. 
#' @param data Either an object inheriting from data frame or in case of 
#' time-dependent covariates a list of data frames, where the first data frame
#' contains the time-to-event information and static covariates while the 
#' second (and potentially further data frames) contain information on 
#' time-dependent covariates and the times at which they have been observed.
#' @param formula A two sided formula with a Surv object on the left-hand-side
#' and covariate specification on the right-hand-side (RHS). The RHS can be an
#' extended formula, which specifies how TDCs should be transformed using
#' specials concurrent and cumulative.
#' @param ... Further arguments passed to the data.frame method and possibly 
#' to survSplit
#' @param cut Break points, used to partition the follow up into intervals. 
#' If unspecified, all unique event times will be used.
#' @max_time If cut is unspecified, this will be the last possible event time.
#' All event times after max_time will be administratively censored at max_time.
#' @return 
#' @export
#' @import pammtools randomForest
#' @author Philipp Kopper
#' library(pammtools)
#' set.seed(31072019)
#' df <- cbind.data.frame(
#' x1 = runif (n, -3, 3),
#' x2 = runif (n, 0, 6))
#' # two component formula specifying cause specific hazards
#' form <- ~ log(0.5) + x1 + x2 | log(0.9) + (-1) * x1 + (-1 * x2)
#' df <- sim_pexp_cr(form, df, seq(0, 3, by = 0.25)) %>%
#'  mutate(
#'    cens_time = runif(n(), 0.5, 3),
#'    status = if_else(cens_time < time, 0, 1),
#'    time = pmin(time, cens_time),
#'    type = status * type)
#' df <- df[, - c(5, 6)]
#' colnames(df)[7] <- "obs_times"
#' ped_cr <- as_ped_cr(data = df, Surv(obs_times, status) ~ ., 
#' id = "id", cut = seq(0, max(df$obs_times), 0.25))
as_ped_cr_subdist <- function(data, formula, keep_status = TRUE, 
                              censor_code = 0L, imputation_set = NULL, 
                              ntree = 500L, ...) {
  assert_data_frame(data)
  assert_formula(formula)
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  id_str <- id ## to be improved
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  status <- unique(true_status)
  status <- status[status != censor_code]
  ped <- vector(mode = "list", length = length(status))
  if (is.null(imputation_set)) {
    predictors <- !(colnames(data) %in% c(id_str, status_str, time_str))
  } else {
    predictors <- imputation_set
  }
  for (i in 1:length(status)) {
    current_data <- data
    to_be_predicted <- !(data[[status_str]] %in% c(status[i], censor_code))
    predict_set <- current_data[to_be_predicted, ]  
    model_set <- current_data[(current_data[[status_str]] == censor_code), ]
    rf <- randomForest(x = model_set[, predictors], 
                       y = model_set[[time_str]],
                       ntree = ntree)
    prediction <- predict(rf, newdata = predict_set)
    current_data[to_be_predicted, time_str] <- max(current_data[to_be_predicted, time_str])
      #pmax(current_data[to_be_predicted, time_str], prediction)
    current_data[[status_str]][data[[status_str]] != status[i]] <- 0L
    current_data[[status_str]][data[[status_str]] == status[i]] <- 1L
    ped[[i]] <- as_ped(current_data, formula, ...)
    class(ped[[i]]) <- c("ped", "data.frame")
  }
  class(ped) <- "ped_cr_subdist"
  attr(ped, "risks") <- status
  ped
}

impute_censorship <- function(data, status, time) {
  
}

pem_cr.default <- function(formula, family = poisson, data, offset, ...) {
  pem_cr.ped_cr(formula, family = poisson, data, offset, ...)
}

#' Wrapping function for PEMs with competing risks using stats::glm
#' 
#' This function serves as a wrapper for the glm function (stats) to be 
#' applicable to competing risks. Competing risks are separately modelled 
#' in their own glm. The results can be interpreted synoptically, 
#' using the corresponding methods (summary etc.) or separately using the glm
#' methods for each entry alone.
#' The input is entered as for a single risk PEM.
#' 
#' @param formula an object of class formula or a string convertible to it.
#' The model formula corresponding to a glm (from the stats package).
#' @param family the family of the glm to be modelled. Only poisson is
#'  reasonable. Hence, any other input will not be accepted.
#' @param data A data.frame of class ped_cr that features time-to-event data
#' and convariates. (see https://adibender.github.io/pammtools/)
#' @param offset The offset for each observation. Contained in data.
#' @param ... additional arguments passed to the glm function.
#' @return a list of glms - one entry for a single competing risk.
#' @import checkmate
#' @importFrom stats glm
#' @export
#' @examples 
#' library(pammtools)
#' set.seed(31072019)
#' df <- cbind.data.frame(
#' x1 = runif (n, -3, 3),
#' x2 = runif (n, 0, 6))
#' # two component formula specifying cause specific hazards
#' form <- ~ log(0.5) + x1 + x2 | log(0.9) + (-1) * x1 + (-1 * x2)
#' df <- sim_pexp_cr(form, df, seq(0, 3, by = 0.25)) %>%
#'  mutate(
#'    cens_time = runif(n(), 0.5, 3),
#'    status = if_else(cens_time < time, 0, 1),
#'    time = pmin(time, cens_time),
#'    type = status * type)
#' df <- df[, - c(5, 6)]
#' colnames(df)[7] <- "obs_times"
#' ped_cr <- as_ped_cr(data = df, Surv(obs_times, status) ~ ., id = "id",
#' cut = seq(0, max(df$obs_times), 0.25))
#' pem_cr <- glm_cr(ped_status ~ interval + x1 + x2, 
#'                  data = ped_cr, offset = offset, family = poisson())
#' @author Philipp Kopper
pem_cr.ped_cr_subdist <- function(formula, family = poisson, data, offset, ...) {
  check_input(formula, data, offset)
  res <- fit_cr_subdist(formula, family, data, offset, m_type = "glm", ...)
  class(res) <- "pem_cr"
  attr(res, "risks") <- attr(data, "risks")
  #for methods
  return(res)
}

#' Fitting of cr PEMs / PAMMs
#' 
#' This function fits a gams/glms depending on the input for a cr framework.
#' It is the operative sub-routine of glm_cr and gam_cr.
#' @param formula an object of class formula or a string convertible to it.
#' The model formula corresponding to a gam or glm.
#' @param family the family of the gam to be modelled. Only poisson is
#'  reasonable. Hence, any other input will not be accepted.
#' @param data A data.frame of class ped_cr that features time-to-event data
#' and convariates. (see https://adibender.github.io/pammtools/)
#' @param offset The offset for each observation. Contained in data.
#' @param m_type A character indicating the type of model to be fit: either 
#' "glm" or "gam".
#' @param ... additional arguments passed to the gam/glm function.
#' @return a list of gams or glms - one entry for a single competing risk.
#' @importFrom mgcv gam
#' @importFrom stats glm
#' @author Philipp Kopper
fit_cr_subdist <- function(formula, family, data, offset, m_type, ...) {
  crs <- attr(data, "risks")
  n_crs <- length(crs)
  res <- vector(mode = "list", length = n_crs)
  names(res) <- crs
  for (i in 1:n_crs) {
    # this function is supposed to make a ped_cr object to a ped object
    # where we only investiagte one of the competing risks
    current_data <- data[[i]]
    command <- paste(m_type, "(formula = formula, family = family, ", 
                     "data = current_data, offset = offset, ...)", sep = "")
    res[[i]] <- eval(parse(text = command)) # verpönt
  }
  return(res)
}

modify_ped_cr <- function(ped, formula, data, type, censorship, censor_independent, censor_formula, ...) {
  if (type == "cs") {
    res <- ped
    attr(res, "type") <- "cs"
    res
  } else {
    ### here do the assertions
    res <- modify_ped_cr_sh(ped, formula, data, censorship, censor_independent, censor_formula, ...)
    attr(res, "type") <- "sh"
    res
  }
}

modify_ped_cr_sh <- function(ped, formula, data, censorship, censor_independent, censor_formula, ...) {
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  if (censorship %in% c("none", "admin")) {
    modify_ped_cr_sh_nocensor(ped, formula, data, censor_formula, ...)
  } else {
    modify_ped_cr_sh_censor(ped, formula, data, censor_independent, censor_formula, ...)
  }
}

modify_ped_cr_sh_nocensor <- function(ped, formula, data, censor_formula = "glm", ...) {
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  status <- unique(true_status)
  status <- status[status != censor_code]
  ped_sets <- vector(mode = "list", length = length(status))
  for (i in 1:length(status)) {
    current_data <- data
    current_data[!(data[[status_str]] %in% c(0, status[i])), time_str] <- 
      max(data[[time_str]])
    current_data[[status_str]][data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][data[[status_str]] == status[i]] <- 1
    ped_sets[[i]] <- as_ped(current_data, formula, ...)
  }
  class(ped_sets) <- c("ped_cr", "ped", "data.frame")
  attr(ped_sets, "risks") <- attr(data, "risks")
  attr(ped_sets, "type") <- "sh_noncensor"
  ped_sets
}

modify_ped_cr_sh_censor <- function(ped, formula, data, censor_independent, 
                                    censor_formula, ...) {
  if (censor_independent) {
    modify_ped_cr_sh_censor_independent(ped, formula, data, censor_formula, ...)
  } else {
    modify_ped_cr_sh_censor_dependent(ped, formula, data, censor_formula, ...)
  }
}

modify_ped_cr_sh_censor_independent <- function(ped, formula, data, censor_formula,
                                                threshold = 0.0001, ...) {
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  status <- unique(true_status)
  status <- as.integer(status[status != censor_code])
  status <- status[order(status)]
  ped_status <- vector(mode = "list", length = length(status))
  censor_data <- data
  censor_data[data[[status_str]] == 0, status_str] <- 1
  censor_data[data[[status_str]] != 0, status_str] <- 0
  censor_ped <- as_ped(censor_data, formula, ...)
  if (censor_formula == "glm") {
    censor_model <- glm(ped_status ~ interval, data = censor_ped, family = "poisson", offset = offset)
  } else {
    censor_model <- gam(ped_status ~ s(tend), data = censor_ped, family = "poisson", offset = offset)
  }
  int_frame <- int_info(censor_ped)
  new_frame <- int_frame[int_frame$interval %in% unique(censor_ped$interval), ]
  new_frame <- add_surv_prob(new_frame, censor_model, ci = FALSE)
  int_frame <- rbind(new_frame, 
                     cbind(int_frame[nrow(new_frame):nrow(int_frame),],
                           surv_prob = rep(new_frame$surv_prob[nrow(new_frame)],
                                           nrow(int_frame) - nrow(new_frame) + 1)))
  ped_sets <- vector(mode = "list", length = length(status))
  for (i in 1:length(status)) {
    modified_data <- data
    modified_data[!(data[[status_str]] %in% c(0, status[i])), time_str] <- max(data[[time_str]])
    current_data <- modified_data
    current_data[[status_str]][modified_data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][modified_data[[status_str]] == status[i]] <- 1
    ped_sets[[i]] <- as_ped(current_data, formula, ...)
    actual_event <- as.data.frame(cbind(id = 1:nrow(data), 
                                        actual_event = data[[time_str]]))
    actual_event$actual_event <- pmin(actual_event$actual_event, 
                                      max(ped$tend))
    intervals <- unique(cbind(ped$tstart, ped$tend))
    actual_event$actual_tend <- rep(0, nrow(actual_event))
    for (j in 1:nrow(actual_event)) {
      actual_event$actual_tend[j] <- intervals[(actual_event$actual_event[j] > intervals[, 1]) &
                                            (actual_event$actual_event[j] <= intervals[, 2]), 2]
    }
    actual_event <- merge(actual_event, int_frame[, c("tend", "surv_prob")], by.x = "actual_tend", 
                          by.y = "tend", all = FALSE)
    colnames(actual_event)[4] <- "km_actual" 
    ped_sets[[i]] <- merge(ped_sets[[i]], int_frame[, c("tend", "surv_prob")], by = "tend")
    ped_sets[[i]] <- merge(ped_sets[[i]], actual_event[, c("id", "km_actual")] , by = "id")
    ped_sets[[i]] <- ped_sets[[i]][order(ped_sets[[i]][["id"]], 
                                         ped_sets[[i]][["tstart"]]), ]
    ped_sets[[i]]$weight <- pmin(1L, ped_sets[[i]]$surv_prob / ped_sets[[i]]$km_actual)
    #ped_sets[[i]]$weight <- ifelse(ped_sets[[i]]$weight != 1L, 
     #                              pmax(0, log(ped_sets[[i]]$km_actual) / log(ped_sets[[i]]$surv_prob)),
      #                             ped_sets[[i]]$weight)
    #ped_sets[[i]]$weight <- pmax(0, pmin(1, exp(ped_sets[[i]]$km_actual - ped_sets[[i]]$surv_prob)))
    ped_sets[[i]] <- ped_sets[[i]][ped_sets[[i]][["weight"]] > threshold, ]
    ped_sets[[i]] <- ped_sets[[i]][ , !(colnames(ped_sets[[i]]) %in% c("surv_prob", "km_actual"))]
    print(ped_sets[[i]]$weight)
  }
  ped_sets
}


modify_ped_cr_sh_censor_independent <- function(ped, formula, data, censor_formula = "gam", ...) {
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  status <- unique(true_status)
  status <- as.integer(status[status != censor_code])
  status <- status[order(status)]
  ped_status <- vector(mode = "list", length = length(status))
  censor_data <- data
  censor_data[data[[status_str]] == 0, status_str] <- 1
  censor_data[data[[status_str]] != 0, status_str] <- 0
  #censor_ped <- as_ped(censor_data, formula, cut = cut, id = "id")#,...)
  censor_ped <- as_ped(censor_data, formula,...)
  if (censor_formula == "gam") {
    censor_model <- gam(ped_status ~ s(tend), data = censor_ped, family = "poisson", offset = offset)
    #censor_model <- glm(ped_status ~ interval, data = censor_ped, family = "poisson", offset = offset)
  } else {
    censor_model <- glm(ped_status ~ interval, data = censor_ped, family = "poisson", offset = offset)
    #censor_model <- gam(censor_formula, data = censor_ped, family = "poisson", offset = offset)
  }
  #predicted_hazards <- add_hazard(censor_ped, censor_model, ci = FALSE)
  #predicted_hazards$increment <- predicted_hazards$tend - predicted_hazards$tstart
  #predicted_hazards$predicted_add_time <- rexp(nrow(predicted_hazards), 
   #                                            predicted_hazards$hazard * 
    #                                             (1 / predicted_hazards$increment))
  int_frame <- int_info(censor_ped)
  new_frame <- int_frame[int_frame$interval %in% unique(censor_ped$interval), ]
  new_frame <- add_hazard(new_frame, censor_model, ci = FALSE)
  int_frame <- rbind(new_frame, 
                     cbind(int_frame[nrow(new_frame):nrow(int_frame),],
                           hazard = rep(new_frame$hazard[nrow(new_frame)],
                                           nrow(int_frame) - nrow(new_frame) + 1)))
  int_frame$increment <- int_frame$tend - int_frame$tstart
  ped_sets <- vector(mode = "list", length = length(status))
  actual_event <- as.data.frame(cbind(id = 1:nrow(data), 
                                      actual_event = data[[time_str]],
                                      status = data[[status_str]]))
  #actual_event$actual_event <- pmin(actual_event$actual_event, max(ped$tend))
  #actual_event$actual_tend <- rep(0, nrow(actual_event))
  #intervals <- unique(cbind(ped$tstart, ped$tend))
  #for (j in 1:nrow(actual_event)) {
  #  actual_event$actual_tend[j] <- intervals[(actual_event$actual_event[j] > intervals[, 1]) &
  #                                             (actual_event$actual_event[j] <= intervals[, 2]), 2]
  #}
  #actual_event <- merge(actual_event, int_frame[, c("tend", "hazard", "increment")], 
  #                      by.x = "actual_tend", by.y = "tend", all.y = FALSE)
  replacement <- max(data[[time_str]])
    #actual_event$actual_event + rexp(actual_event$hazard / actual_event$increment)
  print(actual_event$id[order(actual_event$id)])
  actual_event$new_time <- ifelse(actual_event[[status_str]] != 0, 
                                  replacement, actual_event$actual_event)
  actual_event <- actual_event[order(actual_event$id), ]
  cd <- ped_sets
  for (i in 1:length(status)) {
    modified_data <- data
    modified_data[, time_str] <- ifelse( 
      !(data[[status_str]] %in% c(0, status[i])),
      actual_event$new_time, modified_data[[time_str]])
    current_data <- modified_data
    #duplicates <- current_data
    current_data[[status_str]][modified_data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][modified_data[[status_str]] == status[i]] <- 1
    ped_sets[[i]] <- as_ped(current_data, formula, cut = cut, id = "id")#, ...)
    #ped_sets[[i]] <- as_ped(current_data, formula, ...)
    cd[[i]] <- current_data
  }
  attr(ped_sets, "data_base") <- cd
  ped_sets
}


