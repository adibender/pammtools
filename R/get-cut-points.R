#' Obtain interval break points
#'
#' Default method words for data frames.
#' The list method applies the default method to each data set within the list.
#'
#'
#' @import Formula
#' @keywords internal
get_cut <- function(data, formula, cut = NULL, ...) {
  UseMethod("get_cut", data)
}

#' @rdname get_cut
#' @inherit get_cut
get_cut.default <- function(
  data,
  formula,
  cut      = NULL,
  max_time = NULL,
  event    = 1L,
  ...) {

  if (is.null(cut)) {
    outcome_vars <- get_lhs_vars(formula)
    if (length(outcome_vars) == 2) {
      cut <- unique(data[[outcome_vars[1]]][1L * (data[[outcome_vars[2]]]) == event])
    } else {
      cut_start <- unique(data[[outcome_vars[1]]])
      cut_end <- unique(data[[outcome_vars[2]]])
      cut <- union(cut_start, cut_end)
    }
    if (!is.null(max_time)) {
      cut <- cut[cut < max_time]
      cut <- c(cut, max_time)
    }
  }
  # sort interval cut points in case they are not (so that interval factor
  # variables will be in correct ordering)
  sort(cut)

}
