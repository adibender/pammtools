#' Transform data to Piece-wise Exponential Data (PED)
#'
#' @inherit split_data
#' @export
as_ped <- function(data, formula, ...) {
  UseMethod("as_ped", data)
}

#' @inherit as_ped
#' @export
as_ped.data.frame <- function(data, formula, ...) {

  dots         <- list(...)
  dots$data    <- data
  dots$formula <- formula
  do.call(split_data, dots)

}

#' @inherit as_ped
#' @export
as_ped.nested_fdf <- function(data, formula, ...) {

  Form <- Formula(formula)
  form_elra <- formula(Form, lhs=FALSE, rhs = 2)
  vars_elra <- all.vars(form_elra)
  vars_elra <- vars_elra[!(vars_elra == attr(data, "time_var"))]
  func_components <- get_func(data, form_elra)
  te_vars <- func_components[["te_vars"]] %>% unlist()
  te      <- func_components[["te"]]
  ll_funs <- func_components[["ll_funs"]]
  names(te_vars) <- names(ll_funs) <- te_vars

  form_ped <- formula(Form, lhs=1, rhs = 1)
  ped <- select(data, -one_of(vars_elra)) %>%
    as_ped.data.frame(formula = form_ped, ...)

  func_components <- func_components$func_mats
  for(i in seq_along(func_components)) {
    ped[[names(func_components)[i]]] <- func_components[[i]]
  }
  attr(ped, "ll_funs")  <- ll_funs
  attr(ped, "te")      <- te
  attr(ped, "te_vars") <- te_vars

  ped

}

#' @inherit as_ped
#' @export
as_ped.list <- function(data, formula, ...) {

  assert_class(data, "list")
  assert_class(formula, "formula")

  nl    <- length(data)
  form  <- Formula(formula)
  n_rhs <- length(form)[2]

  if (nl == 1 & n_rhs == 1) {
    ped <- data[[1]] %>% as_ped(formula=form, ...)
  } else {
    if (nl == 2 & n_rhs == 1) {
    stop("Two data sets provided in 'data' but no specification of
      time-dependent covariate effects in 'formula'")
    } else {

      nested_fdf <- nest_tdc(data, form, ...)
      ped <- as_ped(nested_fdf, formula, ...)
    }
  }

  ped

}
#' @rdname as_ped
#' @param x any R object.
#' @export
is.ped <- function(x) inherits(x, "ped")
