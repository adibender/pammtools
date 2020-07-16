
pamm_cr <- function(formula, data = list(), method = "REML", ...,
                    trafo_args = NULL, engine = "gam") {
  UseMethod("pamm_cr", data)
}

pamm_cr.ped_cr_list <- function(formula, data = list(), method = "REML", ...,
                                trafo_args = NULL, engine = "gam") {
  res <- vector(mode = "list", length(data))
  for (i in 1:res) {
    res[[i]] <- pamm(formula, data[[i]], method, ..., trafo_args, engine)
  }
  names(res) <- names(data)
  class(res) <- c("pamm_cr")
  attr(res, "risks") <- attr(data, "risks")
  attr(res, "attr_ped") <- 
    list(
      breaks = attr(data, "breaks"),
      id_var = attr(data, "id_var"),
      intvars = attr(data, "int_vars")
    )
  attr(res, "trafo_args") <- attr(data, "trafo_args")
  return(res)
}

pamm_cr.ped_cr_union <- function(formula, data = list(), method = "REML", ...,
                                 trafo_args = NULL, engine = "gam") {
  pamm(formula, data, method, ..., trafo_args, engine)
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
pamm_cr.ped_cr_list <- function(formula, family = poisson(), 
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
summary.pamm_cr_list <- function(pam_cr) {
  summary_list <- vector(mode = "list", length = length(pem_cr))
  names(summary_list) <- names(pem_cr)
  for (i in 1:length(pem_cr)) {
    pem_cr[[i]]$call <- ""
    summary_list[[i]] <- summary(pem_cr[[i]])
  }
  names(summary_list) <- attr(pem_cr, "risks")
  summary_list
}

#' Print method for competing risk PAMs (piece-wise additive models)
#' @param summary_list a list of summaries where each element is one summary
#' for a gam. Each element should correspond to one partial competing risks
#' model of a PAM.
#' @return A (printed) list of summaries.
#' @author Philipp Kopper
print.pamm_cr_list <- function(summary_list) {
  for (i in 1:length(summary_list)) {
    cat(paste("Risk:", names(summary_list)[i]))
    print(summary_list[[i]])
  }
}
