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

  status_error(data, formula)

  dots         <- list(...)
  dots$data    <- data
  dots$formula <- formula(Formula(formula), lhs=1, rhs=1)
  ped <- do.call(split_data, dots)
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  ped

}

#' @inherit as_ped
#' @export
as_ped.nested_fdf <- function(data, formula, ...) {

  status_error(data, formula)

  dots <- list(...)
  # update interval break points (if neccessary)
  cut <- dots$cut
  if(is.null(cut)) {
    cut <- attr(data, "breaks")
  }
  ccr_breaks <- attr(data, "ccr_breaks")
  cut <- union(cut, ccr_breaks) %>% sort()

  ped <- data %>%
    select_if(is.atomic) %>%
    as_ped.data.frame(
      formula  = formula,
      id       = dots$id,
      cut      = cut,
      max_time = dots$max_time)

  # replace updated attributes
  attr(data, "breaks") <- attr(ped, "breaks")
  attr(data, "id_n") <- ped %>% group_by(!!sym(attr(data, "id_var"))) %>%
    summarize(id_n = n()) %>% pull("id_n") %>% as_vector()
  attr(data, "id_tseq") <- ped %>% group_by(!!sym(attr(data, "id_var"))) %>%
    transmute(id_tseq = row_number()) %>% pull("id_tseq") %>% as_vector()
  attr(data, "id_teseq") <- rep(seq_len(nrow(data)), times=attr(data, "id_n"))

  if(has_special(formula, "concurrent")) {
    ped <- ped %>% add_concurrent(data=data, id_var=dots$id)
  }

  if(has_special(formula, "cumulative")) {
    ped <- add_cumulative(ped, data=data, formula=formula)
    class(ped) <- c("fped", class(ped))
  }
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  ped

}

#' @inherit as_ped
#' @export
as_ped.list <- function(data, formula, ...) {

  assert_class(data, "list")
  assert_class(formula, "formula")

  status_error(data[[1]], formula)

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
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  ped

}
#' @rdname as_ped
#' @param x any R object.
#' @export
is.ped <- function(x) inherits(x, "ped")
