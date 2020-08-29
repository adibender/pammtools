#' Checks if data contains timd-dependent covariates
#'
#' @param data A data frame (potentially) containing time-dependent covariates.
#' @param id_var A character indicating the grouping variable. For each covariate
#' it will be checked if their values change within a group specified by
#' \code{id_var}.
#' @import dplyr
#' @return Logical. \code{TRUE} if data contains time-dependent covariates, else \code{FALSE}.
#' @keywords internal
has_tdc <- function(data, id_var) {

  data %>% group_by(!!sym(id_var)) %>%
    summarize_all(.funs = ~any(length(unique(.)) > 1)) %>%
    select(-one_of(id_var)) %>%
    summarize_all(any) %>% unlist() %>% any()

}
