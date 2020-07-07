#' Transform data to Piece-wise Exponential Data (PED) for a Cause-Specific
#' hazards competing risks model
#'
#' This is a wrapper for the general data transformation function provided by 
#' the pammtools package (as_ped). The use is identical to as_ped.
#' However, the function creates one single PED object for each risk.
#' The results are stored in a list
#' For more details use ?as_ped. 
#' For each risk, all competing risks are assumed to be right-censored.
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
#' @param censor_code Either a string or integer (depending on data) with out
#' outlines which level in the status variable is associated with censorship.
#' @return A list of PED data.frames.
#' @export
#' @author Philipp Kopper
as_ped_cr <- function(data, formula, censor_code = 0L, 
                      output = c("data.frame", "list"), cut = NULL, ...) {
  if (!(is.data.frame(data))) {
    f_data <- data[[1L]]
    assert_data_frame(data[[1L]])
    if (length(data) == 2L) {
      assert_data_frame(data[[2L]])
    }
  } else {
    f_data <- data
    assert_data_frame(f_data)
  }
  assert_formula(formula)
  output <- match.arg(output, c("data.frame", "list"))
  time_str <- all.vars(formula)[1L]
  status_str <- all.vars(formula)[2L]
  true_time <- f_data[[time_str]]
  f_data <- make_numeric(f_data, status_str, censor_code)
  true_status <- f_data[[status_str]]
  status <- unique(true_status)
  status <- status[status != censor_code]
  if (length(status) < 2) {
    stop("There are no competing risks. Use as_ped() instead.")
  } 
  cut <- check_cuts(cut, status)
  ped_sets <- vector(mode = "list", length = length(status))
  for (i in 1:length(status)) {
    current_data <- f_data
    current_data[[status_str]][f_data[[status_str]] != status[i]] <- 0L
    current_data[[status_str]][f_data[[status_str]] == status[i]] <- 1L
    if (!is.data.frame(data) & length(data) == 2L) {
      current_data <- list(current_data, data[[2L]])
    }
    ped_sets[[i]] <- as_ped(data = current_data, formula = formula, 
                            cut = cut[[i]], ...)
    if (output != "list") {
      ped_sets[[i]]$cause <- i
    }
    class(ped_sets[[i]]) <- c("ped", "data.frame")
  }
  if (output == "list") {
    ped <- ped_sets
    class(ped) <- c("ped_cr_list", "ped_cr")
    attr(ped, "risks") <- status
    return(ped)
  } else {
    ped <- as.data.frame(Reduce("rbind", ped_sets))
    class(ped) <- c("ped_cr_df", "ped_cr", "data.frame")
    attr(ped, "intvars") <- c(attr(ped, "intvars"), "cause")
    attr(ped, "breaks") <- cut
    attr(ped, "trafo_args")[["cut"]] <- cut
    attr(ped, "risks") <- status
    return(ped)
  }
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

check_cuts <- function(cut, status) {
  if (is.null(cut)) {
    return(rep(list(cut), length(status)))
  }
  if ((!is.list(cut)) & (!is.numeric(cut))) {
    stop("cut must be either a numeric vector or a list of numeric vectors.")
  }
  if (is.list(cut)) {
    if (any((!unlist(lapply(cut, is.numeric))) & (!unlist(lapply(cut, is.null))))) {
      stop("cut must be either a numeric vector or a list of numeric vectors.")
    }
  }
  if ((length(cut) != length(status)) & (!is.numeric(cut))) {
    stop("cut must have length 1 or identical length as the number of status in the data.")
  }
  if (is.numeric(cut)) {
    return(rep(list(cut), length(status)))
  } else {
    return(cut)
  }
}
