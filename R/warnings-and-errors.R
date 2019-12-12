#' Warn if new t_j are used
#'
#' @importFrom stats model.frame
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
