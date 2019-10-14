#' Extract time-dependent covariates from data set
#'
#' For all covariates in the data set, this functions checks if the values
#' of the covariate changes per ID or other grouping variable. Returns the names
#' of variables that change over time.
#'
#' @param data A data frame (potentially) containing time-dependent covariates.
#' @param id_var A character indicating the grouping variable. For each covariate
#' it will be checked if their values change within a group specified by
#' \code{id_var}.
#' @import dplyr
#' @return A character vector containing names of variables that are not constant
#' in each group (\code{id_var}).
#' @keywords internal
get_tdc <- function(data, id_var) {

  data %>% group_by_(.dots = list(id_var)) %>%
    summarize_all(.funs = ~any(length(unique(.)) > 1)) %>%
    select_if (any) %>%
    names()

}

#' @inherit get_tdc
has_tdc <- function(data, id_var) {

  data %>% group_by(!!sym(id_var)) %>%
    summarize_all(.funs = ~any(length(unique(.)) > 1)) %>%
    select(-one_of(id_var)) %>%
    summarize_all(any) %>% unlist() %>% any()

}


#' Extract unique cut points when time-dependent covariates present
#'
#' Given a data frame with one row per subject containing event times and
#' a data frame containing time points at which a time-dependent covariate changes
#' its value, returns the unique time at which either event occurs or a
#' time-dependent covariate changes its value.
#'
#' @inheritParams get_tdc
#' @param tdc_df A data frame containing information on time-dependent variables
#' in long format. Needs to contain a "time" variable indicating when the
#' TDCs change their value. Must have same name as variable indicating event
#' times in \code{event_df}.
#' @param time_var A character, specifies the column of the event or
#' censoring time in \code{event_df} and the time of measurement for
#' the time-dependent covariates in \code{tdc_df}.
#' @param status_var As \code{time_var}, but specifies column containing the
#' event indicator. Can be missing in the \code{tdc_df}.
#' @param cens_value The value that indicates censoring in the
#' \code{status_var} column.
#' @import dplyr
#' @keywords internal
combine_cut <- function(
  event_df,
  tdc_df,
  time_var,
  status_var,
  tz_var,
  cens_value = 0) {

  tdc_time   <- tdc_df %>% pull(tz_var) %>% unlist() %>% unique()
  event_time <- event_df %>% select(one_of(time_var)) %>% unlist()
  event_time <- event_time[event_df[[status_var]] != cens_value] %>% unique()

  union(tdc_time, event_time) %>% sort()

}
