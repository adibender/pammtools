#' Extract variables from the left-hand-side of a formula
#'
#' @rdname formula_helpers
#' @param formula A \code{\link{formula}} object.
#' @import Formula
#' @keywords internal
get_lhs_vars <- function(formula) {


  if (is.character(formula) ) formula <- as.formula(formula)
  formula(Formula(formula), lhs = TRUE, rhs = FALSE) %>% all.vars()

}

#' Extract variables from the right-hand side of a formula
#'
#' @rdname formula_helpers
#' @inherit get_lhs_vars
#' @keywords internal
get_rhs_vars <- function(formula) {

  if (is.character(formula) ) formula <- as.formula(formula)
  formula(Formula(formula), lhs = FALSE, rhs = TRUE) %>% all.vars()

}


#' @inherit get_lhs_vars
#' @keywords internal
get_tdc_vars <- function(formula, specials = "cumulative") {

  f2      <- formula(Formula(formula), lhs = FALSE, rhs = 2)
  terms_f <- terms(f2, specials = specials)
  all.vars(terms_f)

}

#' @inherit get_lhs_vars
#' @keywords internal
get_tdc_form <- function(formula) {
  formula(Formula(formula), lhs = FALSE, rhs = 2)
}

#' @inherit get_lhs_vars
#' @keywords internal
get_ped_form <- function(formula) {
  formula(Formula(formula), lhs = 1, rhs = 1)
}



#' @keywords internal
has_tdc_form <- function(formula) {

  formula <- Formula(formula)
  length_form <- length(formula)
  length_form[2] > 1

}

has_lhs <- function(formula) {
  length(Formula(formula))[1] > 0
}
