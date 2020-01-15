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
