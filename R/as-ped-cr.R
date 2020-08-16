#' Transform competing risks data to Piece-wise Exponential Data (PED)
#'
#' Transformation of competing risks data to Piece-wise Exponantial Data (PED).
#' The data transformation works equivantly to the single event case (see
#' \code{\link{as_ped}} for details). More details about data transformation are
#' also provided in the \href{https://adibender.github.io/pammtools//articles/data-transformation.html}{data-transformation
#' vignette}.
#'
#' This function wraps the \code{as_ped} function and repeats the data
#' transfomation for each single risk. The data transformation procedure can be
#' different for all associated risks, e.g. the user can supply distinct cut
#' points for each single risk. This function returns either a data.frame
#' which complies with the data.frames proposed
#' \href{https://arxiv.org/abs/2006.15442}{here} or a list of single risk
#' PED objects. (This depends on the argument "output".)
#' Currently our function only supports cause-specific competing
#' risks. That means that for each single risk a competing event is treated
#' (and recoded) as censoring event.
#'
#' @rdname as_ped_cr
#' @inheritParams as_ped
#' @param cut Either a numeric vector or a list of numeric vectors.
#' If provided as list, the length of the list must equal the number of
#' competing risks.
#' Specifies split points, used to partition the follow up into intervals.
#' If unspecified, all unique event times will be used.
#' This argument interacts with argument \code{combine}. If \code{cut} is unspecified
#' and \code{combine} is \code{FALSE}, each risk has a different set of cuts.
#' If set to \code{combine = TRUE}, the cut points are a union of unique event
#' times across all competing events.
#' @param censor_code Either a string or integer (depending on data) with out
#' outlines which level in the status variable is associated with censorship.
#' @param combine A logical value. Defaults to \code{TRUE}. Indicates which
#' output type should be used. \code{TRUE} results in a joint PED object for
#' all risks, i.e., PED objects for each cause are stacked and returned as one data frame.
#' \code{FALSE} results in a named list of PED objects.
#' @param ... Further arguments passed to the \code{as_ped} function or its
#' methods (\code{data.frame} method) and eventually to
#' \code{\link[survival]{survSplit}}
#' @importFrom Formula Formula
#' @examples
#' data("pbc", package = "survival")
#' pbc_tmp <- pbc[c(1, 2, 5), ]
#' # creates cause specific data sets with cause specific, custom interval split points
#' # and stacks them
#' pbc_tmp %>% as_ped_cr(Surv(time, status) ~ age + sex,
#'                       cut = list(c(0, 200, 500), c(0, 500, 1550)))
#' # returns cause specific data sets as a list
#' pbc_tmp %>% as_ped_cr(Surv(time, status) ~ age + sex, combine = FALSE)
#' # provide censoring code if status is not numeric
#' pbc_tmp$status <- c("death", "discharge", "censoring")
#' # by default, a union of event times for all causes will be used as split points
#' pbc_tmp %>% as_ped_cr(Surv(time, status) ~ age + sex, censor_code = "censoring")
#' pbc_tmp %>% as_ped_cr(Surv(time, status) ~ age + sex, combine = FALSE,
#'                       censor_code = "censoring")
#' @return A data frame class \code{ped} in piece-wise exponential data format
#' with an additional column indicating which cause is observed (if output ==
#' "data.frame".) or a named list of single risk data.frames (if output ==
#' "list).
#' @export
#' @author Philipp Kopper
as_ped_cr <- function(
  data,
  formula,
  cut = NULL,
  max_time,
  censor_code = 0L,
  combine = TRUE,
  ...) {

  assert_logical(combine)
  df <- check_data(data)
  data <- df[[1L]]
  event_data <- df[[2L]]
  assert_formula(formula)
  time_str <- all.vars(formula)[1L]
  status_str <- all.vars(formula)[2L]
  true_time <- event_data[[time_str]]
  event_data <- make_numeric(event_data, status_str, censor_code)
  true_status <- event_data[[status_str]]
  if (is.data.frame(data)) {
    status <- unique(data[[status_str]])
  } else {
    status <- unique(data[[1L]][[status_str]])
  }
  status_num <- unique(event_data[[status_str]])
  status <- status[status != censor_code]
  if (length(status) < 2L) {
    stop("There are no competing risks. Use as_ped() instead.")
  }
  cut <- check_cuts(cut, status, combine, true_time)
  ped_sets <- vector(mode = "list", length = length(status))
  for (i in 1L:length(status)) {
    current_data <- event_data
    current_data[[status_str]][event_data[[status_str]] != status_num[i]] <- 0L
    current_data[[status_str]][event_data[[status_str]] == status_num[i]] <- 1L
    if (!is.data.frame(data) & length(data) == 2L) {
      current_data <- append(list(current_data), list(data[[2L]]))
    } else if (!is.data.frame(data) & length(data) > 2L) {
      current_data <- append(list(current_data), data[2L:length(data)])
    }
    ped_sets[[i]] <- as_ped(data = current_data, formula = formula,
                            cut = cut[[i]], ...)
    if (combine) {
      ped_sets[[i]]$cause <- as.factor(as.character(status[i]))
    }
    class(ped_sets[[i]]) <- c("ped", "data.frame")
  }
  if (!combine) {
    ped <- ped_sets
    names(ped) <- status
    class(ped) <- c("ped_cr_list", "ped_cr", "list")
    attr(ped, "risks") <- status
    return(ped)
  } else {
    ped <- do.call(rbind, ped_sets)
    class(ped) <- c("ped_cr_union", "ped_cr", "ped", "data.frame")
    attr(ped, "intvars") <- c(attr(ped, "intvars"), "cause")
    attr(ped, "breaks") <- unique(unlist(cut))
    attr(ped, "trafo_args")[["cut"]] <- unique(unlist(cut))
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
  data[data[[status_str]] == censor_code, status_str] <- 0L
  for (i in 1L:length(risks)) {
    data[data[[status_str]] == risks[i], status_str] <- i
  }
  data[[status_str]] <- as.numeric(data[[status_str]])
  data
}

check_cuts <- function(cut, status, combine, times) {
  if (is.null(cut)) {
    if (!combine) {
      return(rep(list(cut), length(status)))
    } else {
      return(rep(list(unique(times)), length(status)))
    }
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
    stop("data must be either a data.frame (tibble / data.table etc.) or a list.")
  }
  if (!(is.data.frame(data))) {
    event_data <- data[[1L]]
    assert_data_frame(event_data)
    if (length(data) == 1L) {
      data <- data[[1L]]
    }
    if (length(data) == 2L) {
      assert_data_frame(data[[2L]])
    }
  } else {
    event_data <- data
    assert_data_frame(event_data)
  }
  return(list(data, event_data))
}
