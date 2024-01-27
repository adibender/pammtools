#' Create nested data frame from data with time-dependent covariates
#'
#' Provides methods to nest data with time-dependent covariates (TDCs).
#' A \code{formula} must be provided where the right hand side (RHS) contains
#' the structure of the TDCs
#'

#' @param data A suitable data structure (e.g. unnested data frame with
#' concurrent TDCs or a list where each element is a data frame, potentially
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
nest_tdc <- function(data, formula, ...) {
  UseMethod("nest_tdc", data)
}

#' @rdname nest_tdc
#' @export
nest_tdc.default <- function(data, formula, ...) {

  dots <- list(...)
  id <- dots$id

  tdc_vars     <- dots$tdc_vars
  outcome_vars <- dots$outcome_vars
  if (is.null(tdc_vars)) {
    tdc_vars <- get_tdc_vars(formula)
  }
  if (is.null(outcome_vars)) {
    outcome_vars <- get_lhs_vars(formula)
  }

  tdc_vars <- setdiff(tdc_vars, outcome_vars)

  if (!any(colnames(data) %in% tdc_vars) | !has_tdc(data, id)) {
    vars_to_exclude <- intersect(colnames(data), tdc_vars)
    return(data %>% select(-one_of(vars_to_exclude)))
  } else {
    df_list <- map(
      tdc_vars,
      ~ tidyr::nest(.data = data[, c(id, .x)], {{.x}} := one_of(.x)))

    suppressMessages(nested_df <- df_list %>% reduce(left_join)) # better: numeric vectors in each list element
    class(nested_df) <- c("nested_fdf", class(nested_df))
  }

  nested_df %>% as_tibble()

}


#' @rdname nest_tdc
#' @export
nest_tdc.list <- function(data, formula, ...) {

  dots <- list(...)
  cut  <- dots[["cut"]]

  data_dummy <- suppressMessages(
    map(data, ~.x[1,]) %>% do.call(what = left_join))

  # preprocess information on time-dependent covariates
  lgl_concurrent <- has_special(formula, special = "concurrent")
  lgl_cumulative <- has_special(formula)
  tdc_vars <- character(0)
  form_tdc <- get_tdc_form(formula, data = data_dummy)

  if (lgl_concurrent) {
    # obtain information on concurrent effects
    ccr <- prep_concurrent(data, form_tdc)
    # update cut points
    ccr_time <- ccr[["ccr_time"]]
    # update vector of TDCs
    ccr_tz_vars <- map_chr(ccr[["ccr_list"]], ~.x[["tz_var"]]) %>% unique()
    ccr_vars <- map(ccr[["ccr_list"]], ~.x[["col_vars"]]) %>% unlist()
    tdc_vars <- c(tdc_vars, ccr_tz_vars, ccr_vars)

  } else {
    ccr      <- NULL
    ccr_time <- NULL
  }

  if (lgl_cumulative) {
    func_list    <- eval_special(form_tdc, data = data[[2]])
    func_vars    <- map(func_list, ~.x[["col_vars"]]) %>% unlist()
    func_tz_vars <- map_chr(func_list, ~.x[["tz_var"]]) %>% unique()
    tdc_vars     <- c(tdc_vars, func_vars, func_tz_vars) %>% unique()
  } else {
    func_list <- NULL
  }
  # remove outcome vars from TDCs vector
  outcome_vars <- get_lhs_vars(formula)
  time_var     <- outcome_vars[1]
  tdc_vars     <- setdiff(tdc_vars, outcome_vars)
  suppressMessages(
  nested_df <- map(data, ~nest_tdc(.x, formula = formula, id = dots$id,
      tdc_vars = tdc_vars, outcome_vars = outcome_vars)) %>%
    reduce(left_join) %>% as_tibble()
    )

  ## add atrributes
  cut <- get_cut(nested_df, formula, cut = dots$cut, max_time = dots$max_time)

  id_n <- nested_df %>%
    pull(time_var) %>%
    pmin(max(cut)) %>%
    map_int(findInterval, vec = cut, left.open = TRUE, rightmost.closed = TRUE)

  attr_old <- attributes(nested_df)
  attributes(nested_df) <- c(attr_old, list(
    id_var     = dots$id,
    time_var   = time_var,
    status_var = outcome_vars[2],
    tdc_vars   = tdc_vars,
    breaks     = cut,
    ccr        = ccr,
    ccr_breaks = ccr_time,
    func_list  = func_list,
    id_n       = id_n,
    id_tseq    = id_n %>% map(seq_len) %>% unlist(),
    id_tz_seq   = rep(seq_along(nested_df[[dots$id]]), times = id_n)))

  class(nested_df) <- c("nested_fdf", class(nested_df))

  nested_df

}
