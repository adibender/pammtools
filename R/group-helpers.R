#' Return ungrouped data frame without grouping variables
#'
#' @param data A data frame (potentially `grouped_df`) from which grouping
#' variables will be removed.
#' @import dplyr
#' @return Returns \code{data} without grouping variables (and group properties).
#' @keywords internal
rm_grpvars <- function(data) {

  if(is.grouped_df(data)) {
    grp.vars <- group_vars(data)
    data     <- data %>%
      ungroup() %>%
      select(-one_of(grp.vars))
  }

  return(data)

}
