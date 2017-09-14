#' Extract variables from the left-hand-side of a formula
#'
#' @rdname formula_helpers 
#' @param formula A \code{\link{formula}} object. 
#' @importFrom stats as.formula
get_lhs_vars <- function(formula) {

	if (is.character(formula) ) formula <- as.formula(formula)
	all.vars(formula[-3])

}

#' Extract variables from the right-hand side of a formula
#' 
#' @rdname formula_helpers 
#' @inheritParams get_lhs_vars
get_rhs_vars <- function(formula) {

	if (is.character(formula) ) formula <- as.formula(formula)
	all.vars(formula[-2])

}