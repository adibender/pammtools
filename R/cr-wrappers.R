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

pem_cr <- function(formula, family = poisson, ped, offset, ...) {
  UseMethod("pem_cr", ped)
}

pem_cr.ped_cr_df <- function(formula, family = poisson, ped, offset, ...) {
  res <- glm(formula = formula, family = family, data = ped, 
             offset = offset, ...)
  class(res) <- c("pem_cr", "gam", "glm")
  attr(res, "risks") <- attr(ped, "risks")
  return(res)
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
pem_cr.ped_cr_list <- function(formula, family = poisson, ped, offset, ...) {
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

pam_cr <- function(formula, family = poisson(), 
                   ped = list(), offset = NULL, bam = FALSE, ...) {
  UseMethod("pam_cr", ped)
}

pam_cr.ped_cr_df <- function(formula, family = poisson(), 
                             ped = list(), offset = NULL, bam = FALSE, ...) {
  if (bam) {
    res <- bam(formula = formula, family = family, 
               data = ped[[i]], offset = offset, ...)
  } else {
    res <- gam(formula = formula, family = family, 
               data = ped[[i]], offset = offset, ...)
  }
  class(res) <- c("pam_cr", "pem_cr", "gam", "glm")
  attr(res, "risks") <- attr(ped, "risks")
  return(res)
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
pam_cr.ped_cr_list <- function(formula, family = poisson(), 
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
  class(res) <- c("pam_cr_list", "pem_cr_list")
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
summary.pam_cr_list <- function(pam_cr) {
  summary.pem_cr(pam_cr)
}

#' Print method for competing risk PAMs (piece-wise additive models)
#' @param summary_list a list of summaries where each element is one summary
#' for a gam. Each element should correspond to one partial competing risks
#' model of a PAM.
#' @return A (printed) list of summaries.
#' @author Philipp Kopper
print.pam_cr_list <- function(summary_list) {
  print.pem_cr(summary_list)
}