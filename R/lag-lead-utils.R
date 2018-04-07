#' Construct data that represents a lag-lead window
#'
#' Constructs lag-lead window data set from raw inputs or from data objects
#' with suitable information stored in attributes.
#'
#' @param x Either a numeric vector of follow-up cut points or a suitable object.
#' @param ... Further arguments passed to methods.
#' @export
get_laglead <- function(x, ...) {
  UseMethod("get_laglead", x)
}

#' @rdname get_laglead
#' @inherit get_laglead
#' @param te A vector of exposure times
#' @param ll_fun Function with that indicates how the lag-lead matrix
#' should be contructed. First argument is the follow up time
#' second argument is the time of exposure.s
#' @importFrom dplyr mutate
#' @importFrom tidyr crossing
#' @export
get_laglead.default <- function(x, te, ll_fun, ...) {

  LL_df <- crossing(t=x, te=te) %>% mutate(LL = ll_fun(t, te)*1L)
  class(LL_df) <- c("LL_df", class(LL_df))

  LL_df

}

#' @rdname get_laglead
#' @inherit get_laglead
#' @importFrom purrr map2_dfr
#' @export
get_laglead.data.frame <- function(x, ...) {

  t       <- attr(x, "breaks")
  te      <- attr(x, "te")
  ll_funs <- attr(x, "ll_funs")

  LL_df <- map2_dfr(te, ll_funs,
      ~get_laglead.default(t, .x, ll_fun=.y), .id="te_var")
  if(!inherits(LL_df, "LL_df")) {
    class(LL_df) <- c("LL_df", class(LL_df))
  }

  LL_df

}

#' Plot lag lead window from appropriate DF
#'
#' Plot Lag-Lead windows
#'
#' Given a matrix defining a lag lead window, returns respective plot as a
#' \code{ggplot2} object.
#'
#' @inheritParams get_laglead
#' @param high_col Color used to highlight exposure times within the lag-lead window.
#' @param low_col Color of exposure times outside the lag-lead window.
#' @param grid_col Color of grid lines.
#' @import checkmate ggplot2
#' @examples
#' ## Example 1: supply t, te, ll_fun directly
#' gg_laglead(1:10, 1:5, ll_fun = function(t,te) {t >= te})
#'
#' ## Example 2: extract information on t, te, ll_from data with respective attributes
#' data("simdf_elra", package = "pammtools")
#' gg_laglead(simdf_elra)
#' @export
gg_laglead <- function(x, ...) {
  UseMethod("gg_laglead", x)
}

#' @rdname gg_laglead
#' @inherit gg_laglead
#' @export
gg_laglead.default <- function(x, te, ll_fun, ...) {

  LL_df <- get_laglead(x, te, ll_fun = ll_fun)
  gg_laglead(LL_df, ...)

}

#' @inherit gg_laglead
#' @rdname gg_laglead
#' @export
gg_laglead.LL_df <- function(
  x,
  high_col   = "grey20",
  low_col    = "smokewhite",
  grid_col   = "lightgrey",
  ...) {

  x <- left_join(x, int_info(unique(x$t)), by = c("t" = "tend"))
  x <- x %>% filter(!is.na(interval)) %>%
    mutate(
      te = as.factor(te),
      interval = factor(interval, levels = rev(levels(interval))) )
  gg_ll <- ggplot(x, aes_string(y = "interval", x = "te")) +
    geom_tile(aes_string(fill = "LL"), colour = grid_col) +
    scale_fill_gradient(low = low_col, high = high_col) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), position = c("left"))
    xlab(expression(t[e])) + ylab(expression(t)) +
    theme(legend.position = "none")

  if(!is.null(x$te_var)) {
    gg_ll <- gg_ll + facet_wrap(~te_var, scales="free_x")
  }

  gg_ll

}

#' @rdname gg_laglead
#' @inherit gg_laglead
#' @export
gg_laglead.nested_fdf <- function(x, ...) {

  LL_df <- get_laglead(x)

  gg_laglead(LL_df, ...)

}
