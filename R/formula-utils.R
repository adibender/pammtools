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
#' @keywords internal
get_rhs_vars <- function(formula) {

  if (is.character(formula) ) formula <- as.formula(formula)
  formula(Formula(formula), lhs = FALSE, rhs = TRUE) %>% all.vars()

}


#' @inherit get_lhs_vars
#' @keywords internal
get_tdc_vars <- function(
  formula,
  specials = "cumulative",
  data     = NULL) {

  f_specials <- get_tdc_form(formula, data = data, tdc_specials = specials)
  terms_f    <- terms(f_specials, specials = specials)
  all.vars(terms_f)

}

#' @inherit get_lhs_vars
#' @keywords internal
get_tdc_form <- function(
  formula,
  data         = NULL,
  tdc_specials = c("concurrent", "cumulative"),
  invert       = FALSE) {

  terms   <- terms(formula, data = data, specials = tdc_specials)
  labels <- attr(terms, "term.labels")
  ind_tdc <- map(tdc_specials, ~grep(.x, labels)) %>% unlist()

  if(invert) {
    if(length(ind_tdc) > 0) {
      formula(terms[ind_tdc * -1])
    } else {
      formula
    }
  } else {
    formula(terms[ind_tdc])
  }


}



#' @inherit get_lhs_vars
#' @keywords internal
get_ped_form <- function(
  formula,
  data         = NULL,
  tdc_specials = c("concurrent", "cumulative")) {

  get_tdc_form(formula, data = data, tdc_specials = tdc_specials, invert = TRUE)

}



#' @keywords internal
has_tdc_form <- function(
  formula,
  tdc_specials = c("concurrent", "cumulative")) {

  form_chr <- as.character(formula) %>% paste0(collapse = "")


  any(map_lgl(tdc_specials, ~grepl(.x, form_chr)))

}

has_lhs <- function(formula) {

  length(Formula(formula))[1] > 0

}
