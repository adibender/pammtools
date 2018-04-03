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
	formula(Formula(formula), lhs = FALSE, rhs=TRUE) %>% all.vars()

}


#' @inherit get_lhs_vars
get_tdc_vars <- function(formula, specials = "func") {

  f2      <- formula(Formula(formula), lhs=FALSE, rhs = 2)
  terms_f <- terms(f2, specials = specials)
  all.vars(terms_f)

}


#' A formula special for defining functional covariates
#'
#' @rdname func
#' @importFrom purrr map
#' @export
#' @keywords internal
func <- function(...,
  te_var,
  ll_fun = function(t, te) {t >= te},
  suffix = NULL) {

  vars        <- as.list(substitute(list(...)))[-1]
  vars_chr    <- vars %>% map(~as.character(.))
  lgl_latency <- map_lgl(vars_chr, ~any(. %in% "latency"))

  if (any(lgl_latency)) {
    latency_var <- unlist(vars_chr)[unlist(vars_chr) != "latency"][lgl_latency]
    col_vars    <- unlist(vars_chr)[unlist(vars_chr) != "latency"]
  } else {
    latency_var <- ""
    col_vars    <- unlist(vars_chr)
  }

  list(
    col_vars    = col_vars,
    latency_var = latency_var,
    te_var      = te_var,
    suffix      = suffix,
    ll_fun      = ll_fun)

}
