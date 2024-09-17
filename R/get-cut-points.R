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

  outcome_vars <- get_lhs_vars(formula)
  if (is.null(cut)) {
    if (length(outcome_vars) == 2) {
      cut <- unique(data[[outcome_vars[1]]][1L * (data[[outcome_vars[2]]]) == event])
    } else {
      cut_start <- unique(data[[outcome_vars[1]]])
      cut_end   <- unique(data[[outcome_vars[2]]])
      cut       <- union(cut_start, cut_end)
    }
    if (!is.null(max_time)) {
      cut <- cut[cut < max_time]
      cut <- c(cut, max_time)
    }
  } else {
    if (length(outcome_vars) == 2) {
      # sort interval cut points in case they are not (so that interval factor
      # variables will be in correct ordering)
      cut <- sort(unique(cut))
    } else {
      # sort interval cut points in case they are not (so that interval factor
      # variables will be in correct ordering)
      # add transitions within interval cut points to not lose information
      cut <- sort(unique(union(cut, data[[outcome_vars[2]]][1L * (data[[outcome_vars[3]]]) == event &
                         1L * (data[[outcome_vars[2]]]) < max(cut)])
                         )
                  )
    }
    if (!is.null(max_time)) {
      cut <- cut[cut < max_time]
      cut <- c(cut, max_time)
    }
  }
  
  return(sort(unique(cut)))

}

#' @rdname get_cut
#' @inherit get_cut
get_cut.list <- function (
  data,
  formula,
  cut       = NULL,
  max_time  = NULL,
  event     = 1L,
  timescale = "gap",
  ...) {

  lhs_vars <- get_lhs_vars(formula)
  if (length(lhs_vars) == 3 & timescale == "gap") {
    rhs_vars <- get_rhs_vars(formula)
    formula_cuts <- as.formula(
      paste0("Surv(", lhs_vars[2], ",", lhs_vars[3], ") ~ ",
        paste(rhs_vars, collapse = "+")))
  } else {
    formula_cuts <- formula
  }
  cuts <- map(
    .x = data,
    .f = ~get_cut.default(
      data     = .x,
      formula  = formula_cuts,
      cut      = cut,
      max_time = max_time,
      event    = event,
      ...)
  )

  cuts <- Reduce(union, cuts)
  sort(unique(cuts))

}
