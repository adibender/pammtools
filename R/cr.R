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
pem_cr <- function(formula, family = poisson, ped, offset, ...) {
  #check_input(formula, ped, offset)
  res <- vector(mode = "list", length = length(ped))
  for (i in 1:length(res)) {
    res[[i]] <- glm(formula = formula, family = family, 
                    data = ped[[i]], offset = offset, ...)
  }
  names(res) <- attr(ped, "risks")
  which_type <- c("sh", "cs")[c("sh", "cs") %in% class(ped)]
  class(res) <- c("pem_cr", which_type)
  attr(res, "risks") <- attr(ped, "risks")
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
pam_cr <- function(formula, family = gaussian(), 
                   ped = list(), offset = NULL, bam = FALSE, ...) {
  #check_input(formula, ped, offset)
  res <- vector(mode = "list", length = length(ped))
  for (i in 1:length(res)) {
    if (bam) {
      res[[i]] <- bam(formula = formula, family = family, 
                      data = ped[[i]], offset = offset, ...)
    } else {
      res[[i]] <- gam(formula = formula, family = family, 
                      data = ped[[i]], offset = offset, ...)
    }
  }
  names(res) <- attr(ped, "risks")
  which_type <- c("sh", "cs")[c("sh", "cs") %in% class(ped)]
  class(res) <- c("pam_cr", "pem_cr", which_type)
  attr(res, "risks") <- attr(ped, "risks")
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


as_ped_cr_cs <- function(data, formula, censor_code = 0, cut = NULL, ...) {
  assert_data_frame(data)
  assert_formula(formula)
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
    current_data[[status_str]][data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][data[[status_str]] == status[i]] <- 1
    ped_sets[[i]] <- as_ped(current_data, formula, cut = cut)
    class(ped_sets[[i]]) <- c("ped", "data.frame")
  }
  names(ped_sets) <- status
  class(ped_sets) <- c("cs", "ped_cr")
  attr(ped_sets, "risks") <- status
  ped_sets
}


as_ped_cr_sh <- function(data, formula, censor_code = 0, cut = NULL, 
                         max_time = NULL, ...) {
  time_str <- all.vars(formula)[1]
  status_str <- all.vars(formula)[2]
  if (is.null(max_time)) max_time <- max(data[[time_str]])
  true_time <- data[[time_str]]
  data <- make_numeric(data, status_str, censor_code)
  true_status <- data[[status_str]]
  status <- unique(true_status)
  status <- status[status != censor_code]
  ped_sets <- vector(mode = "list", length = length(status))
  for (i in 1:length(status)) {
    current_data <- data
    current_data[!(data[[status_str]] %in% c(0, status[i])), time_str] <- 
      max_time
    current_data[[status_str]][data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][data[[status_str]] == status[i]] <- 1
    ped_sets[[i]] <- as_ped(current_data, formula, cut = cut, ...)
    class(ped_sets[[i]]) <- c("ped", "data.frame")
  }
  names(ped_sets) <- status
  class(ped_sets) <- c("sh", "ped_cr")
  attr(ped_sets, "risks") <- status
  ped_sets
}

as_ped_cr_cens <- function(data, formula, censor_code = 0, cut = NULL, 
                           censor_formula = NULL, ...) {
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
  censor_ped <- as_ped(censor_data, formula, cut = cut, ...)
  if (is.null(censor_formula)) {
    censor_formula <- as.formula(ped_status ~ s(tend))
  }
  censor_model <- mgcv::gam(censor_formula, data = censor_ped, family = "poisson", 
                            offset = offset)
    #censor_model <- glm(ped_status ~ interval, data = censor_ped, family = "poisson", offset = offset)
  predicted_hazards <- add_hazard(censor_ped, censor_model, ci = FALSE)
  predicted_hazards$increment <- predicted_hazards$tend - predicted_hazards$tstart
  predicted_hazards$predicted_add_time <- rexp(nrow(predicted_hazards), 
                                               predicted_hazards$hazard)
  add_time <- rep(0, nrow(data))
  tend <- rep(0, nrow(data))
  j <- 1
  for (i in 1:(nrow(predicted_hazards) - 1)) {
    if (predicted_hazards$id[i] != predicted_hazards$id[i + 1]) {
      add_time[j] <- predicted_hazards$predicted_add_time[i]
      tend[j] <- predicted_hazards$tend[i]
      j <- j + 1
    }
  }
  add_time[j] <- 
    predicted_hazards$predicted_add_time[nrow(predicted_hazards)]
  actual_event <- as.data.frame(cbind(id = 1:nrow(data), 
                                      actual_event = data[[time_str]],
                                      status = data[[status_str]]))
  replacement <- pmin(max(data[[time_str]]), actual_event$actual_event + add_time)
  actual_event$new_time <- ifelse(actual_event[[status_str]] != 0, 
                                  replacement, actual_event$actual_event)
  actual_event <- actual_event[order(actual_event$id), ]
  ped_sets <- vector(mode = "list", length = length(status))
  cd <- ped_sets
  for (i in 1:length(status)) {
    modified_data <- data
    modified_data[, time_str] <- ifelse( 
      !(data[[status_str]] %in% c(0, status[i])),
      actual_event$new_time, modified_data[[time_str]])
    current_data <- modified_data
    current_data[[status_str]][modified_data[[status_str]] != status[i]] <- 0
    current_data[[status_str]][modified_data[[status_str]] == status[i]] <- 1
    #ped_sets[[i]] <- as_ped(current_data, formula, cut = cut, id = "id")#, ...)
    ped_sets[[i]] <- as_ped(current_data, formula, cut = cut, ...)
    class(ped_sets[[i]]) <- c("ped", "data.frame")
    cd[[i]] <- current_data
  }
  attr(ped_sets, "data_base") <- cd
  names(ped_sets) <- status
  class(ped_sets) <- c("sh", "ped_cr")
  attr(ped_sets, "risks") <- status
  ped_sets
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
    res[[i]] <- eval(parse(text = command)) # verpönt: invoke!
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
  if (!("ped_cr_subdist" %in% class(data))) {
    assert_data_frame(data)
  } else {
    for (i in 1:length(data)) {
      assert_data_frame(data[[i]])
    }
  }
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

