#' Create nested data frame from data with time-dependent covariates
#'
#' Provides methods to nest data with time-dependent covariates (TDCs).
#' A \code{formula} must be provided where the right hand side (RHS) contains
#' the structure of the TDCs
#'
#' @inheritParams tidyr::nest
#' @inheritParams split_data
#' @param data A suitable data structure (e.g. unnested data frame with
#' concurrent TDCs or a list wehre each element is a data frame, potentially
#' containing TDCs as specified in the RHS of \code{formula}).
#' Only TDCs present in \code{formula} will be returned.
#' @param formula A two sided formula with a two part RHS, where the second
#' part indicates the structure of the TDC structure.
#' @param ... Further arguments passed to methods.
#' @import checkmate dplyr
#' @importFrom tidyr nest
#' @importFrom purrr map map_int reduce
#' @export
#' @keywords internal
nest_tdc <- function(data, formula,...) {
  UseMethod("nest_tdc", data)
}

#' @inherit nest_tdc
#' @rdname nest_tdc
#' @param vars A character vector of TDCs that will be nested.
#' @param id A character giving the name of the ID column.
#' @export
nest_tdc.default <- function(data, formula, ...) {

  dots <- list(...)
  id <- dots$id

  tdc_vars     <- names_tdc(formula)
  outcome_vars <- names_lhs(formula)
  time_var     <- outcome_vars[1]
  tdc_vars     <- setdiff(tdc_vars, outcome_vars)

  if (!any(colnames(data) %in% tdc_vars)) {
    return(data)
  } else {
    nested_df <- map(tdc_vars,
        ~nest(data=data[, c(id, .)], -one_of(id), .key = !!.)) %>%
      reduce(left_join)# Would be better to have numeric vecotrs within each list element
      class(nested_df) <- c("nested_fdf", class(nested_df))
  }

  nested_df

}



#' @inherit nest_tdc
#' @rdname nest_tdc
#' @export
nest_tdc.list <- function(data, formula, ...) {

  dots         <- list(...)
  tdc_vars     <- names_tdc(formula)
  outcome_vars <- names_lhs(formula)
  time_var     <- outcome_vars[1]
  tdc_vars     <- setdiff(tdc_vars, outcome_vars)

  nested_df <- map(data, ~nest_tdc(., formula = formula, id=dots$id,...)) %>%
    reduce(left_join)

  ## add atrributes
  if(!is.null(dots$cut)) {
    cut <- dots$cut
  } else {
    cut <- nested_df %>% pull(time_var) %>% unique() %>% sort()
  }
  id_n <- nested_df %>% pull(time_var) %>%
    pmin(max(cut)) %>%
    map_int(findInterval, vec=cut, left.open=TRUE, rightmost.closed=TRUE)

  attr_old <- attributes(nested_df)
  attributes(nested_df) <- c(attr_old, list(
    id_var     = dots$id,
    time_var   = time_var,
    status_var = outcome_vars[2],
    tdc_vars   = tdc_vars,
    breaks     = cut,
    id_n       = id_n,
    id_tseq    = id_n %>% map(seq_len) %>% unlist(),
    id_teseq   = rep(seq_along(nested_df[[dots$id]]), times = id_n)))

  class(nested_df) <- c("nested_fdf", class(nested_df))

  nested_df

}


names_tdc <- function(formula, specials = "func") {

  f2      <- formula(Formula(formula), lhs=FALSE, rhs = 2)
  terms_f <- terms(f2, specials = specials)
  all.vars(terms_f)

}

names_lhs <- function(formula, specials="Surv") {
  all.vars(formula(Formula(formula), lsh=TRUE, rhs=FALSE))
}
