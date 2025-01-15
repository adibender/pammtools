#' Warn if new t_j are used
#'
#' @keywords internal
warn_about_new_time_points <- function(object, newdata, ...) {

  UseMethod("warn_about_new_time_points", object)

}

#' @inherit warn_about_new_time_points
#' @keywords internal
warn_about_new_time_points.glm <- function(object, newdata, time_var, ...) {

  is_pam <- (inherits(object, "gam" ) | inherits( object, "scam"))

  if(is_pam & is.null(object$model)){
    return(invisible())
  }

  original_intervals <- if (is_pam) {
    unique(model.frame(object)[[time_var]])
  } else levels(model.frame(object)[[time_var]])
  prediction_intervals <- if (is_pam) {
    unique(newdata[[time_var]])
  } else levels(factor(newdata[[time_var]]))
  new_ints <- which(!(prediction_intervals %in% original_intervals))
  n_out <- pmin(10, length(new_ints))
  if (length(new_ints)) {
   message <- paste0(
    "Time points/intervals in new data not equivalent to time points/intervals during model fit.",
    " Setting intervals to values not used for original fit",
    "can invalidate the PEM assumption and yield incorrect predictions.")
   if (is_pam) warning(message) else stop(message)
  }
}


#' @rdname warn_about_new_time_points
warn_about_new_time_points.pamm <- function(object, newdata, ...) {

  if (inherits(object, "pamm")) {
    int_original <- int_info(object)
    if ("interval" %in% colnames(newdata)) {
      int_new <- unique(newdata[["interval"]])
      if(!all(int_new %in% int_original[["interval"]])) {
        warning(
          paste0(
            "Time points/intervals in new data not equivalent to time points/intervals during model fit.",
            " Setting intervals to values not used for original fit",
            "can invalidate the PEM assumption and yield incorrect predictions."
          )
        )
      }

    }
  }

}

# #' @keywords internal
# #' @importFrom dplyr intersect union setequal
# warn_partial_overlap <- function(event_id, tdc_id) {
#   common_id <- intersect(event_id, tdc_id)
#   union_id  <- union(event_id, tdc_id)
#   if (!setequal(common_id, union_id)) {
#     warning("Not all IDs are present in both data sets.
#       IDs that do not appear in both data sets will be removed.")
#   }

#   invisible(common_id)

# }


status_error <- function(data, formula, censor_code = 0L) {

  outcome_vars <- get_lhs_vars(formula)
  if (!any(unique(data[[outcome_vars[length(outcome_vars)]]]) != censor_code)) {
    stop(paste(
      "No events in data! Check your",
      outcome_vars[length(outcome_vars)],
      "variable."))
  }

}
