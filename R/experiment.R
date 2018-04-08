#' Create new data for cumulative effects models
#'
#' Experimental new implementation of \code{make_newdata} for \code{fped}
#' objects.
#' @param x An object of class \code{fped}.
#' @param ... Specification of variables and their values. Functions can be used,
#' e.g., \code{x1 = seq_range(x1)}.
#'
#' @importFrom tidyr crossing
#' @export
#' @keywords internal
mk_ndf <- function(x, ...) {

  expressions    <- quos(...)
  expr_evaluated <- map(expressions, lazyeval::f_eval, data=x)
  lgl_atomic     <- map_lgl(expr_evaluated, is_atomic)

  partI  <- expr_evaluated[lgl_atomic] %>% cross_df()
  partII <- do.call(combine_df, expr_evaluated[!lgl_atomic])
  ndf    <- combine_df(partI, partII)

  cumulative_vars <- attr(x, "func_mat_names") %>%
    setdiff(names(ndf))
  time_part <- x %>%
    group_by(tend) %>%
    slice(1) %>%
    select(one_of(cumulative_vars)) %>%
    map(simplify_mat) %>%
    bind_cols()

  rest   <- x %>% select(-one_of(c(colnames(ndf), colnames(time_part))))
  si     <- sample_info.fped(rest)
  out_df <- combine_df(si, time_part, ndf)
  int_df <- int_info(attr(x, "breaks"))
  right_join(int_df, out_df) %>%
    select(intersect(colnames(x), names(.)))

}

simplify_mat <- function(x) {
  if(is.matrix(x)) {
    return(x[, 1])
  } else {
    return(x)
  }
}
