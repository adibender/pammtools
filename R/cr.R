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
#' @author Philipp Kopper
glm_cr <- function(formula, family = poisson, data, offset, ...) {
  check_input(formula, data, offset)
  res <- fit_cr(formula, family, data, offset, type = "glm", ...)
  class(res) <- "pem_cr"
  attr(res, "risks") <- attr(data, "risks")
  #for methods
  return(res)
}

#' Summary method for competing risk PEMs (piece-wise exp. models)
#' 
#' This function summarises the underlying models of a pem_cr object.
#' The summaries are returned list-wise with each element belonging to
#' one competing risk.
#' @param pem_cr An object of class pem_cr where all elements are one 
#' glm object. Each element should correspond to one partial competing risks
#' model of a PEM.
#' @return A list of summaries.
#' @author Philipp Kopper
summary.pem_cr <- function(pem_cr) {
  summary_list <- vector(mode = "list", length = length(pem_cr))
  names(summary_list) <- names(pem_cr)
  for (i in 1:length(pem_cr)) {
    summary_list[[i]] <- summary(pem_cr[[i]])
  }
  names(summary_list) <- attr(pem_cr, "risks")
  summary_list
}

#' Print method for competing risk PEMs (piece-wise exp. models)
#' @param summary_list a list of summaries where each element is one summary
#' for a glm. Each element should correspond to one partial competing risks
#' model of a PEM.
#' @return A (printed) list of summaries.
#' @author Philipp Kopper
print.pem_cr <- function(summary_list) {
  for (i in 1:length(summary_list)) {
    print(summary_list[[i]])
  }
}

#' Wrapping function for PAMs with competing risks using mgcv::gam
#' 
#' This function serves as a wrapper for the gam function (mgcv) to be 
#' applicable to competing risks. Competing risks are separately modelled 
#' in their own gam. The results can be interpreted synoptically, 
#' using the corresponding methods (summary etc.) or separately using the glm
#' methods for each entry alone.
#' The input is entered as for a single risk PAM.
#' 
#' @param formula an object of class formula or a string convertible to it.
#' The model formula corresponding to a gam (from the mgcv package).
#' @param family the family of the gam to be modelled. Only poisson is
#'  reasonable. Hence, any other input will not be accepted.
#' @param data A data.frame of class ped_cr that features time-to-event data
#' and convariates. (see https://adibender.github.io/pammtools/)
#' @param offset The offset for each observation. Contained in data.
#' @param ... additional arguments passed to the gam function.
#' @return a list of gams - one entry for a single competing risk.
#' @import checkmate
#' @importFrom mgcv gam
#' @export
#' @author Philipp Kopper
gam_cr <- function(formula, family = gaussian(), 
                   data = list(), offset = NULL, ...) {
  check_input(formula, data, offset)
  res <- fit_cr(formula, family, data, offset, type = "gam", ...)
  class(res) <- "pem_cr"
  attr(res, "risks") <- attr(data, "risks")
  #for methods
  return(res)
}

#' Summary method for competing risk PAMs (piece-wise additive models)
#' 
#' This function summarises the underlying models of a pem_cr object.
#' The summaries are returned list-wise with each element belonging to
#' one competing risk.
#' @param pem_cr An object of class pem_cr where all elements are one 
#' glm object. Each element should correspond to one partial competing risks
#' model of a PEM.
#' @return A list of summaries.
#' @author Philipp Kopper
summary.pam_cr <- function(pam_cr) {
  summary.pem_cr(pam_cr)
}

#' Print method for competing risk PAMs (piece-wise additive models)
#' @param summary_list a list of summaries where each element is one summary
#' for a gam. Each element should correspond to one partial competing risks
#' model of a PAM.
#' @return A (printed) list of summaries.
#' @author Philipp Kopper
print.pam_cr <- function(summary_list) {
  print.pem_cr(summary_list)
}

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
#' @author Philipp Kopper
as_ped_cr <- function(data, formula, keep_status = TRUE, censor_code = 0, ...) {
  assert_data_frame(data)
  assert_formula(formula)
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  stati <- unique(true_status)
  stati <- stati[stati != censor_code]
  ped_stati <- vector(mode = "list", length = length(stati))
  for (i in 1:length(stati)) {
    current_data <- data
    current_data[[status_str]][data[[status_str]] != stati[i]] <- 0
    current_data[[status_str]][data[[status_str]] == stati[i]] <- 1
    current_status <- as_ped(current_data, formula, ...)$ped_status
    ped_stati[[i]] <- current_status * stati[i]
  }
  ped <- as_ped(current_data, formula, ...)
  ped$ped_status <- Reduce("+", ped_stati)
  class(ped) <- c("ped_cr", "ped", "data.frame")
  attr(ped, "risks") <- attr(data, "risks")
  ped
}

#' Helper function which creates the CR data frame for one specific risk.
#' 
#' @param data a ped_cr data.frame. 
#' @param cr The competing risk which we aim a data frame for. 
#' Needs to have the same name / number as in data.
#' @return a ped data.frame object for a single competing risk. 
#' (All other risks are treated as censoring.)
#' @author Philipp Kopper
modify_cr_data <- function(data, cr) {
  data$ped_status <- ifelse(data$ped_status == cr, 1, 0)
  data
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
fit_cr <- function(formula, family, data, offset, m_type, ...) {
  crs <- unique(data$ped_status)
  crs <- crs[!(crs == 0)]
  n_crs <- length(crs)
  res <- vector(mode = "list", length = n_crs)
  names(res) <- crs
  for (i in 1:n_crs) {
    # this function is supposed to make a ped_cr object to a ped object
    # where we only investiagte one of the competing risks
    current_data <- modify_cr_data(data, cr = crs[i])
    command <- paste(m_type, "(formula = formula, family = family, ", 
                    "data = current_data, offset = offset, ...)", sep = "")
    res[[i]] <- eval(parse(text = command))
  }
  return(res)
}

#' Input checking for cr PEMs / PAMMs
#' 
#' This function fits a gams/glms depending on the input for a cr framework.
#' It is the operative sub-routine of glm_cr and gam_cr.
#' @param formula an object of class formula or a string convertible to it.
#' The model formula corresponding to a gam or glm.
#' @param data A data.frame of class ped_cr that features time-to-event data
#' and convariates. (see https://adibender.github.io/pammtools/)
#' @param offset The offset for each observation. Contained in data.
#' @return An assertion if there is false input.
#' @import checkmate
#' @author Philipp Kopper
check_input <- function(formula, data, offset) {
  assert_formula(formula)
  assert_data_frame(data)
  if (is.null(offset)) stop("You need to specifiy an offset.")
}

#' Conversion to numeric risks
#' 
#' This function converts character or factor status into numeric ones.
#' This is necessary for further processing.
#' @param data A data.frame of class ped_cr that features time-to-event data
#' and convariates. (see https://adibender.github.io/pammtools/)
#' @param status_str A character indicating the colname of the status variable.
#' @param censor_code A value (character/numeric/factor) indicating as what
#' censoring has been coded.
#' @return data with (converted) numeric status column
make_numeric <- function(data, status_str, censor_code) {
  risks <- unique(data[[status_str]])
  risks <- risks[risks != censor_code]
  data[data[[status_str]] == censor_code, status_str] <- 0
  for (i in 1:length(risks)) {
    data[data[[status_str]] == risks[i], status_str] <- i
  }
  data[[status_str]] <- as.numeric(data[[status_str]])
  attr(data, "risks") <- risks
  data
}

#' Add (cumulative) hazards / survival probabilities for competing risks
#' 
#' This function is the main subroutine of add_cumu_hazard_cr(),
#' add_hazard_cr() and add_surv_prob().
#' For details see add_cumu_hazard_cr() and add_hazard_cr().
#' @author Philipp Kopper
hazard_adder_cr <- function(newdata, object, hazard_function, type, ci, se_mult, 
                            ci_type, overwrite, time_var, ...) {
  new_cols <- as.data.frame(matrix(0, nrow = nrow(newdata), 
                                   ncol = length(object) + 
                                     length(object) * ci * 2))
  j <- 1
  for (i in 1:length(object)) {
    if (ci) {
      new_object <- hazard_function(newdata, object[[i]], ci = ci, 
                                    se_mult = se_mult, 
                                    ci_type = ci_type, 
                                    overwrite = overwrite, 
                                    time_var = time_var, ...)
      where <- which(colnames(new_object) == "cumu_hazard")
      new_cols[, j:(j + 2)] <- new_object[, where:(where + 2)]
      colnames(new_cols)[j:(j + 2)] <- paste(attr(object, "risks")[i], 
                                             colnames(new_object)
                                             [where:(where + 2)], sep = "_")
      j <- j + 3
    } else {
      new_object <- hazard_function(newdata, object[[i]], ci = ci, 
                                    se_mult = se_mult, 
                                    ci_type = ci_type, 
                                    overwrite = overwrite, 
                                    time_var = time_var, ...)
      where <- which(colnames(new_object) == "cumu_hazard")
      new_cols[, i] <- new_object[, where]
      colnames(new_cols)[i] <- paste(attr(object, "risks")[i], 
                                     colnames(new_object)[where], sep = "_")
    }
  }
  return(cbind(newdata, new_cols))
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
add_cumu_hazard_cr <- function(newdata, object, type = c("link", "response"), 
                               ci = TRUE, se_mult = 2, 
                               ci_type = c("default", "delta", "sim"),
                               overwrite = FALSE, time_var = NULL, ...) {
  hazard_adder_cr(newdata, object, hazard_function = add_cumu_hazard, type, ci, 
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
add_surv_prob_cr <- function(newdata, object, type = c("link", "response"), 
                             ci = TRUE, se_mult = 2, 
                             ci_type = c("default", "delta", "sim"),
                             overwrite = FALSE, time_var = NULL, ...) {
  hazard_adder_cr(newdata, object, hazard_function = add_surv_prob, type, ci, 
                  se_mult, ci_type, overwrite, time_var, ...)
}
