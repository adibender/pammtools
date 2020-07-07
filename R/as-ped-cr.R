#' Transform competing risks data to Piece-wise Exponential Data (PED)
#'
#' This is the competing risks wrapper for the general data transformation 
#' function provided by the
#' \code{pammtools} package. 
#' The data transformation works identically compared to \code{as_ped}.
#' Like for the single risk case, two main applications must be distinguished:
#' \enumerate{
#'  \item Transformation of standard time-to-event data.
#'  \item Transformation of time-to-event data with time-dependent covariates (TDC).
#' }
#' For the latter, the type of effect one wants to estimate is also
#' important for the data transformation step.
#' In any case, the data transformation is specified by a two sided formula.
#' In case of TDCs, the right-hand-side of the formula can contain formula specials
#' \code{concurrent} and \code{cumulative}.
#' See the \href{https://adibender.github.io/pammtools//articles/data-transformation.html}{data-transformation}
#' vignette for details.
#'
#' This function wraps the \code{as_ped} function and repeats the data 
#' transfomation for each single risk. The data transformation procedure can be
#' different for all associated risks, e.g. the user can supply distinct cut 
#' points for each single risk. This function returns either a data.frame 
#' which complies with the data.frames proposed 
#' \href{https://arxiv.org/abs/2006.15442}{here} or a list of single risk 
#' data.frames. (This depends on the argument "output".)
#' Currently our function only supports cause-specific competing 
#' risks. That means that for each single risk a competing event is treated 
#' (and recoded) as censoring event.
#'
#' @rdname as_ped_cr
#' @param data Either an object inheriting from data frame or in case of
#' time-dependent covariates a list of data frames (of length 2), where the 
#' first data frame contains the time-to-event information and static 
#' covariates while the second contains information on time-dependent
#' covariates and the times at which they have been observed.
#' @param formula A two sided formula with a \code{\link[survival]{Surv}} object
#' on the left-hand-side and covariate specification on the right-hand-side (RHS).
#' The RHS can be an extended formula, which specifies how TDCs should be transformed
#' using specials \code{concurrent} and \code{cumulative}.
#' @param cut Either a numeric vector or a list of numeric vectors.
#' If provided as list, the length of the list must equal the number of 
#' competing risks.
#' Break points, used to partition the follow up into intervals.
#' If unspecified, all unique event times will be used.
#' @param max_time If \code{cut} is unspecified, this will be the last
#' possible event time. All event times after \code{max_time}
#' will be administratively censored at \code{max_time}.
#' @param censor_code Either a string or integer (depending on data) with out
#' outlines which level in the status variable is associated with censorship.
#' @param output A chracter value: either "data.frame" or "list". 
#' Will be matched via \code{match.arg}.
#' @param ... Further arguments passed to the \code{as_ped} function or its 
#' methods (\code{data.frame} method) and eventually to 
#' \code{\link[survival]{survSplit}}
#' @importFrom Formula Formula
#' @examples
#' library(mvna)
#' data("sir.adm", package = "mvna")
#' sir_adm <- sir.adm[c(1, 36, 85), ]
#' sir_adm %>% as_ped_cr(Surv(time, status) ~ age + sex, 
#'                       cut = list(c(0, 5, 20), c(0, 10, 25)))
#' sir_adm %>% as_ped_cr(Surv(time, status) ~ age + sex, output = "list")
#' sir_adm$status <- c("death", "discharge", "cens")
#' sir_adm %>% as_ped_cr(Surv(time, status) ~ age + sex, censor_code = "cens")
#' sir_adm %>% as_ped_cr(Surv(time, status) ~ age + sex, output = "list", 
#'                       censor_code = "cens")
#' @return A data frame class \code{ped} in piece-wise exponential data format
#' with an additional column indicating which cause is observed (if output == 
#' "data.frame".) or a named list of single risk data.frames (if output == 
#' "list).
#' @export
#' @author Philipp Kopper
as_ped_cr <- function(data, formula, cut = NULL, max_time, censor_code = 0L,
                      output = c("data.frame", "list"), ...) {
  df <- check_data(data)
  data <- df[[1]]
  f_data <- df[[2]]
  assert_formula(formula)
  output <- match.arg(output, c("data.frame", "list"))
  time_str <- all.vars(formula)[1L]
  status_str <- all.vars(formula)[2L]
  true_time <- f_data[[time_str]]
  f_data <- make_numeric(f_data, status_str, censor_code)
  true_status <- f_data[[status_str]]
  if (is.data.frame(data)) {
    status <- unique(data[[status_str]])
  } else {
    status <- unique(data[[1L]][[status_str]])
  }
  status_num <- unique(f_data[[status_str]])
  status <- status[status != censor_code]
  if (length(status) < 2) {
    stop("There are no competing risks. Use as_ped() instead.")
  } 
  cut <- check_cuts(cut, status)
  ped_sets <- vector(mode = "list", length = length(status))
  for (i in 1:length(status)) {
    current_data <- f_data
    current_data[[status_str]][f_data[[status_str]] != status_num[i]] <- 0L
    current_data[[status_str]][f_data[[status_str]] == status_num[i]] <- 1L
    if (!is.data.frame(data) & length(data) == 2L) {
      current_data <- list(current_data, data[[2L]])
    }
    ped_sets[[i]] <- as_ped(data = current_data, formula = formula, 
                            cut = cut[[i]], ...)
    if (output != "list") {
      ped_sets[[i]]$cause <- as.factor(as.character(status[i]))
    }
    class(ped_sets[[i]]) <- c("ped", "data.frame")
  }
  if (output == "list") {
    ped <- ped_sets
    names(ped) <- status
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

check_data <- function(data) {
  if (!(is.data.frame(data)) & !(is.list(data))) {
    stop("data must be either a data.frame or a list.")
  }
  if (!(is.data.frame(data))) {
    f_data <- data[[1L]]
    assert_data_frame(f_data)
    if (length(data) == 1L) {
      data <- data[[1L]]
    }
    if (length(data) == 2L) {
      assert_data_frame(data[[2L]])
    }
  } else {
    f_data <- data
    assert_data_frame(f_data)
  }
  return(list(data, f_data))
}
