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
glm_cr <- function(formula, family = poisson, data, offset, ...) {
  check_input(formula, data, offset)
  res <- fit_cr(formula, family, data, offset, m_type = "glm", ...)
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
    pem_cr[[i]]$call <- ""
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
    cat(paste("Risk:", names(summary_list)[i]))
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
#' #' library(pammtools)
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
#' pem_cr <- gam_cr(ped_status ~ interval + x1 + x2, 
#'                  data = ped_cr, offset = offset, family = poisson())
gam_cr <- function(formula, family = gaussian(), 
                   data = list(), offset = NULL, ...) {
  check_input(formula, data, offset)
  res <- fit_cr(formula, family, data, offset, m_type = "gam", ...)
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
#' @import pammtools
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
as_ped_cr <- function(data, formula, keep_status = TRUE, censor_code = 0, ...) {
  assert_data_frame(data)
  assert_formula(formula)
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  status <- unique(true_status)
  status <- status[status != censor_code]
  ped_status <- vector(mode = "list", length = length(status))
  for (i in 1:length(status)) {
    current_data <- data
    current_data[[status_str]][data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][data[[status_str]] == status[i]] <- 1
    current_status <- as_ped(current_data, formula, ...)$ped_status
    ped_status[[i]] <- current_status * status[i]
  }
  ped <- as_ped(current_data, formula, ...)
  ped$ped_status <- Reduce("+", ped_status)
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
  crs <- attr(data, "risks")
  n_crs <- length(crs)
  res <- vector(mode = "list", length = n_crs)
  names(res) <- crs
  for (i in 1:n_crs) {
    # this function is supposed to make a ped_cr object to a ped object
    # where we only investiagte one of the competing risks
    current_data <- modify_cr_data(data, cr = crs[i])
    command <- paste(m_type, "(formula = formula, family = family, ", 
                    "data = current_data, offset = offset, ...)", sep = "")
    res[[i]] <- eval(parse(text = command)) # verpönt
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
  if (is.numeric(risks)) {
    #otherwise not intuitive.
    risks <- risks[order(risks)]
  }
  attr(data, "risks") <- risks
  data[data[[status_str]] == censor_code, status_str] <- 0
  for (i in 1:length(risks)) {
    data[data[[status_str]] == risks[i], status_str] <- i
  }
  data[[status_str]] <- as.numeric(data[[status_str]])
  data
}
