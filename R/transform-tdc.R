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

#' Create piece-wise exponential data in case of time-dependent covariates
#'
#' Given to data frames, one containing event time information (one row per subject)
#' and one containing information on time-dependent covariates, creates
#' piece-wise exponential data (with one split per event time and time at
#' which a TDC changes its value).
#'
#' @inherit combine_cut
#' @inheritParams split_data
#' @param event_df Data frame (or similar) containing survival information.
#' @param tdc_df Data frame (or similar) containing information on time-dependent
#' covariates
#' @param id_var The ID variable name, identifying subjects.
#' @param tz_var The time variable in \code{tdc_df} indicating time points at
#' which time-dependent covariate (tdc) was observed.
#' Needs to be the same name in both data sets.
#' @param entry_time If scalar, the time-point at which the follow up for each
#' observation unit begins. (Eventually, support for subject specific
#' entry time could be supported through this argument).
#' @importFrom tidyr fill
#' @importFrom rlang sym
#' @export
split_tdc <- function(
  formula,
  event_df,
  tdc_df,
  tz_var,
  id_var     = "id",
  time_var   = "time",
  status_var = "status",
  cens_value = 0,
  entry_time = 0,
  ...) {

  warning("'split_tdc' is deprecated and will be removed in the near future.
    See 'as_ped' and 'concurrent' for the same functionality.")

  assert_data_frame(event_df)
  assert_data_frame(tdc_df)
  assert_subset(c(id_var, time_var, status_var), colnames(event_df))
  assert_subset(c(id_var, tz_var), colnames(tdc_df))

  common_id <- warn_partial_overlap(event_df[[id_var]], tdc_df[[id_var]])
  event_df  <- filter(event_df, !!sym(id_var) %in% common_id)
  tdc_df    <- filter(tdc_df, !!sym(id_var) %in% common_id)
  # intervals must be split at each event time and time at which the TDC
  # changes its value
  utime <- union(entry_time, combine_cut(event_df, tdc_df, time_var, status_var,
    tz_var, cens_value = cens_value))
  # for joining, we remove baseline information of variables that are present
  # as TDC variables in tdc_df
  tdc_vars <- setdiff(get_tdc(tdc_df, id_var),
    c(id_var, time_var, tz_var, status_var))
  tdc_in_event_df <- intersect(tdc_vars, names(event_df))
  if (!(length(tdc_in_event_df) == 0)) {
    event_df <- event_df %>%  select(-one_of(tdc_vars))
  }
  ped <- split_data(formula, data = event_df, cut = utime, id = id_var,
    zero = entry_time, ...)

  #
  tdc_df <- tdc_df %>% select(one_of(c(id_var, tz_var, tdc_vars)))


  ped <- left_join(ped, tdc_df, by = c(id_var, "tstart" = tz_var)) %>%
    group_by(!!sym(id_var)) %>%
    fill(setdiff(tdc_vars, c(id_var, time_var, status_var)))

  attr(ped, "id_var") <- id_var
  attr(ped, "time_var") <- time_var
  attr(ped, "status_var") <- status_var
  attr(ped, "tz_var") <- tz_var
  attr(ped, "cens_value") <- cens_value
  attr(ped, "id_n") <- ped %>% group_by(!!sym(id_var)) %>%
    summarize(id_n = n()) %>% pull(.data$id_n)
  attr(ped, "id_tseq")    <- attr(ped, "id_n") %>% map(seq_len) %>% unlist()
  attr(ped, "id_tz_seq") <- rep(seq_along(common_id), times = attr(ped, "id_n"))
  attr(ped, "tz") <- tdc_df %>% pull(tz_var) %>% unique() %>% sort()

  class(ped) <- c("ped", class(ped))

  return(ped)

}
