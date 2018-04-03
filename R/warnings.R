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
