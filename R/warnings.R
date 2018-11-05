#' Warn if new t_j are used
#'
#' @keywords internal
warn_about_new_time_points <- function(newdata, object, time_var) {

  is_pam <- inherits(object, "gam")

  original_intervals <- if (is_pam) {
    unique(model.frame(object)[[time_var]])
  } else levels(model.frame(object)[[time_var]])
  prediction_intervals <- if (is_pam) {
    unique(newdata[[time_var]])
  } else levels(factor(newdata[[time_var]]))
  new_ints <- which(!(prediction_intervals %in% original_intervals))
  n_out <- pmin(10, length(new_ints))
  if (length(new_ints)) {
   message <- paste0("Intervals in <newdata> contain values (",
     paste(prediction_intervals[new_ints[1:n_out]], collapse = ","),
     " ...) not used in original fit.",
     " Setting intervals to values not used for original fit in <object>",
     "can invalidate the PEM assumption and yield incorrect predictions.")
   if (is_pam) warning(message) else stop(message)
  }
}


#' @keywords internal
#' @importFrom dplyr intersect union setequal
warn_partial_overlap <- function(event_id, tdc_id) {
  common_id <- intersect(event_id, tdc_id)
  union_id  <- union(event_id, tdc_id)
  if (!setequal(common_id, union_id)) {
    warning("Not all IDs are present in both data sets.
      IDs that do not appear in both data sets will be removed.")
  }

  invisible(common_id)

}


status_error <- function(data, formula) {

  outcome_vars <- get_lhs_vars(formula)
  if (!any(1L * (unique(data[[outcome_vars[2]]])) == 1)) {
    stop(paste("No events in data! Check your", outcome_vars[2], "variable."))
  }

}
