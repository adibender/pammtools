#' Fit a competing risk piece-wise exponential additive model
#'
#' A thin wrapper around \code{\link[mgcv]{gam}},
#' however, some arguments are prespecified:
#' \code{family=poisson()}, \code{offset=data$offset} and \code{method="REML"}.
#' The first two can not be overwritten. The \code{method} argument
#' can be specified as usual, but defaults to \code{GCV.cp} in
#' \code{\link[mgcv]{gam}}.
#' The wrapper works conditionally on the argument \code{data}.
#' If the competing risks data is supplied in a list, two separate single risk
#' models are estimated w.r.t. the provided \code{formula}.
#' If the data is supplied in a joint PED object (\code{ped_cr_union}),
#' \code{pamm} is straight-forwardly called. Hence, \code{formula} must reflect
#' how the competing risks are supposed to be modeled. (E.g. as an interaction
#' term.)
#' Unlike \code{pamm} this function does not support on-the-spot transformation
#' to a PED object (via \code{trafo_args}) if no PED object is supplied.
#'
#' @param formula A GAM formula, or a list of formulae (see \code{formula.gam}
#' and also \code{gam.models}). These are exactly like the formula for a GLM
#' except that smooth terms, \code{s}, \code{te}, \code{ti} and \code{t2}, can
#' be added to the right hand side to specify that the linear predictor depends
#' on smooth functions of predictors (or linear functionals of these).
#' @param data A PED object (either a \code{ped_cr_union} or \code{ped_cr_list}
#' object) containing the model status and event time variable such as
#' covariates required by the formula. By default the variables are taken from
#' \code{environment(formula)}: typically the environment from which \code{gam}
#' is called.
#' @param method Optimisation method directly handed over to \code{gam}.
#' @param ... Further arguments passed to \code{engine}.
#' @param engine Character name of the function that will be called to fit the
#' model. The intended entries are either \code{"gam"} or \code{"bam"}
#' (both from package \code{mgcv}).
#' @import mgcv
#' @importFrom stats poisson
#' @rdname pamm_cr
#' @seealso \code{\link[mgcv]{gam}}
#' @export
pamm_cr <- function(
  formula,
  data = list(),
  method = "REML",
  ...,
  engine = "gam") {
  UseMethod("pamm_cr", data)
}

#' @rdname pamm_cr
#' @export
pamm_cr.default <- function(
  formula,
  data = list(),
  method = "REML",
  ...,
  engine = "gam") {

  res <- try(
    pamm_cr.ped_cr_union(
      formula,
      data,
      method,
      ...,
      trafo_args = NULL,
      engine = engine),
    silent = TRUE)
  message("Using pamm_cr.ped_cr_union().")
  if (inherits(res, "try-error")) {
    stop(
      cat(
        "Using pamm_cr.ped_cr_union() for the provided object of class",
        class(data), "fails."))
  } else {
    warning(
      cat(
        "There is no pamm_cr method for class", class(data),
        ". pamm_cr.ped_cr_union() is tried instead.",
        "Please double check your results."))
    return(res)
  }
}

#' @rdname pamm_cr
#' @export
pamm_cr.ped_cr_list <- function(
  formula,
  data = list(),
  method = "REML",
  ...,
  engine = "gam") {
    res <- vector(mode = "list", length(data))
    for (i in 1:length(res)) {
      dots <- list(...)
      dots$formula <- formula
      dots$family  <- poisson()
      dots$data   <- data[[i]]
      dots$offset <- data[[i]]$offset
      pamm_fit        <- do.call(engine, dots)
      class(pamm_fit) <- c("pamm", class(pamm_fit))
      pamm_fit        <- append_ped_attr(pamm_fit, data[[i]])
      pamm_fit[["trafo_args"]] <- attr(data[[i]], "trafo_args")
      pamm_fit
      res[[i]] <- pamm_fit
    }
    names(res) <- names(data)
    class(res) <- c("pamm_cr_list", "pamm_cr")
    attr(res, "risks") <- attr(data, "risks")
    attr(res, "attr_ped") <- list(
      breaks = attr(data, "breaks"),
      id_var = attr(data, "id_var"),
      intvars = attr(data, "int_vars"))
  attr(res, "trafo_args") <- attr(data, "trafo_args")
  return(res)

}

#' @rdname pamm_cr
#' @export
pamm_cr.ped_cr_union <- function(
  formula,
  data = list(),
  method = "REML",
  ...,
  engine = "gam") {
    dots <- list(...)
    dots$formula <- formula
    dots$family  <- poisson()
    dots$data   <- data
    dots$offset <- data$offset
    pamm_fit        <- do.call(engine, dots)
    class(pamm_fit) <- c("pamm_cr", "pamm", class(pamm_fit))
    pamm_fit        <- append_ped_attr(pamm_fit, data)
    pamm_fit[["trafo_args"]] <- attr(data, "trafo_args")
    pamm_fit
}

#' Summary method for competing risk PAMs (piece-wise additive models)
#'
#' This function summarises the underlying models of a pamm_cr_list object.
#' The summaries are returned list-wise with each element belonging to
#' one competing risk.
#' @param pam_cr An object of class pamm_cr_list where all elements are one
#' pamm object. Each element should correspond to one partial competing risks
#' model of a PAM(M).
#' @return A list of summaries.
#' @author Philipp Kopper
summary.pamm_cr_list <- function(pam_cr) {
  summary_list <- vector(mode = "list", length = length(pam_cr))
  names(summary_list) <- names(pam_cr)
  for (i in 1:length(pam_cr)) {
    pam_cr[[i]]$call <- ""
    summary_list[[i]] <- summary(pam_cr[[i]])
  }
  names(summary_list) <- attr(pam_cr, "risks")
  summary_list
}

#' Print method for competing risk PAMs (piece-wise additive models)
#'
#' @param summary_list a list of summaries where each element is one summary
#' for a \code{pamm}. Each element should correspond to one partial competing
#' risks model of a PAM(M).
#' @return A (printed) list of summaries.
#' @author Philipp Kopper
print.pamm_cr_list <- function(summary_list) {
  for (i in 1:length(summary_list)) {
    cat(paste("Risk:", names(summary_list)[i]))
    print(summary_list[[i]])
  }
}
