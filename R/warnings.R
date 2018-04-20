#' @keywords internal
#' @importFrom dplyr intersect union setequal
warn_partial_overlap <- function(event_id, tdc_id) {
  common_id <- intersect(event_id, tdc_id)
  union_id  <- union(event_id, tdc_id)
  if(!setequal(common_id, union_id)) {
    warning("Not all IDs are present in both data sets.
      IDs that do not appear in both data sets will be removed.")
  }

  invisible(common_id)

}


status_error <- function(data, formula) {

  outcome_vars <- get_lhs_vars(formula)
  if (!any(1L*(unique(data[[outcome_vars[2]]])) == 1)) {
    stop(paste0("No events in data! Check your ", outcome_vars[2], " variable."))
  }

}
