#' Extract time-dependent covariates from data set 
#' 
#' For all covariates in the data set, this functions checks if the values 
#' of the covariate changes per ID or other grouping variable. Returns the names
#' of variables that change over time. 
#' 
#' @param event_df A data frame (potentially) containing time-dependent covariates. 
#' @param id_var A character indicating the grouping variable. For each covariate
#' it will be checked if their values change within a group specified by 
#' \code{id_var}.
#' @import dplyr 
#' @return A character vector containing names of variables that are not constant 
#' in each group (\code{id_var}).
get_tdc <- function(event_df, id_var) {

	event_df %>% group_by_(.dots=list(id_var)) %>% 
		summarize_all(function(.x) any(length(unique(.x)) > 1)) %>% 
		select_if(any) %>% 
		names()

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
combine_cut <- function(
	event_df, 
	tdc_df, 
	time_var, 
	status_var, 
	cens_value=0) {

	
	tdc_time   <- tdc_df %>% select(one_of(time_var)) %>% unlist() %>% unique()
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
#' @param id_var The ID variable name, identifying subjects. 
#' Needs to be the same name in both data sets.
#' @importFrom tidyr fill_ fill
#' @examples
#' data("pbc", package="survival")# loads both, pbc and pbcseq
#' pbc$status = 1*(pbc$status == 2)
#' pbcseq$time <- pbcseq$day # time of event/measurement must be equal in both data 
#' pbcseq$day  <- NULL
#' pbc_ped     <- split_tdc(Surv(time, status)~., pbc, pbcseq, "id", "time", "status")
#' @export 
split_tdc <- function(
	formula,
	event_df, 
	tdc_df, 
	id_var, 
	time_var, 
	status_var,
	cens_value=0) {

	# intervals must be split at each event time and time at which the TDC 
	# changes its value
	utime <- combine_cut(event_df, tdc_df, time_var, status_var, 
		cens_value=cens_value)
	# for joining, we remove baseline information of variables that are present 
	# as TDC variables in tdc_df
	tdc <- setdiff(get_tdc(tdc_df, id_var), c(id_var, time_var, status_var))
	event_df %<>% select(-one_of(tdc))
	ped <- split_data(formula, data = event_df, cut=utime, id=id_var)

	# 
 	tdc_df %<>% select(one_of(c(id_var, time_var, tdc)))
 	
	ped %>% left_join(tdc_df, by=c(id_var, "tstart"=time_var)) %>% 
		group_by_(.dots=list(id_var)) %>% 
		fill_(setdiff(tdc, c(id_var, time_var, status_var)))

}