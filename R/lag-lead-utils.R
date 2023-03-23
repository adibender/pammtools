#' Construct or extract data that represents a lag-lead window
#'
#' Constructs lag-lead window data set from raw inputs or from data objects
#' with suitable information stored in attributes, e.g., objects created
#' by \code{\link{as_ped}}.
#'
#' @param x Either a numeric vector of follow-up cut points or a suitable object.
#' @param ... Further arguments passed to methods.
#' @examples
#' get_laglead(0:10, tz=-5:5, ll_fun=function(t, tz) { t >= tz + 2 & t <= tz + 2 + 3})
#' gg_laglead(0:10, tz=-5:5, ll_fun=function(t, tz) { t >= tz + 2 & t <= tz + 2 + 3})
#' @export
get_laglead <- function(x, ...) {
  UseMethod("get_laglead", x)
}

#' @rdname get_laglead
#' @param tz A vector of exposure times
#' @param ll_fun Function that specifies how the lag-lead matrix
#' should be constructed. First argument is the follow up time
#' second argument is the time of exposure.
#' @importFrom dplyr mutate
#' @importFrom tidyr crossing
#' @export
get_laglead.default <- function(x, tz, ll_fun, ...) {

  LL_df <- crossing(t = x, tz = tz) %>%
    mutate(LL = ll_fun(.data$t, .data$tz) * 1L) %>%
    group_by(tz) %>%
    mutate(LL = lag(.data$LL, default = 0)) %>%
    ungroup()
  class(LL_df) <- c("LL_df", class(LL_df))

  LL_df

}

#' @rdname get_laglead
#' @importFrom purrr map2_dfr
#' @export
get_laglead.data.frame <- function(x, ...) {

  time    <- attr(x, "breaks")
  tz      <- attr(x, "tz")
  ll_funs <- attr(x, "ll_funs")

  LL_df <- map2_dfr(tz, ll_funs,
      ~get_laglead.default(time, .x, ll_fun = .y), .id = "tz_var")
  if (!inherits(LL_df, "LL_df")) {
    class(LL_df) <- c("LL_df", class(LL_df))
  }

  LL_df

}


#' Plot Lag-Lead windows
#'
#' Given data defining a Lag-lead window, returns respective plot as a
#' \code{ggplot2} object.
#'
#' @inheritParams get_laglead
#' @param high_col Color used to highlight exposure times within the lag-lead window.
#' @param low_col Color of exposure times outside the lag-lead window.
#' @param grid_col Color of grid lines.
#' @import checkmate ggplot2
#' @examples
#' ## Example 1: supply t, tz, ll_fun directly
#'  gg_laglead(1:10, tz=-5:5,
#'   ll_fun=function(t, tz) { t >= tz + 2 & t <= tz + 2 + 3})
#'
#' ## Example 2: extract information on t, tz, ll_from data with respective attributes
#' data("simdf_elra", package = "pammtools")
#' gg_laglead(simdf_elra)
#' @export
#' @seealso get_laglead
gg_laglead <- function(x, ...) {
  UseMethod("gg_laglead", x)
}

#' @rdname gg_laglead
#' @export
gg_laglead.default <- function(x, tz, ll_fun, ...) {

  LL_df <- get_laglead(x, tz, ll_fun = ll_fun)
  gg_laglead(LL_df, ...)

}

#' @rdname gg_laglead
#' @export
gg_laglead.LL_df <- function(
  x,
  high_col   = "grey20",
  low_col    = "whitesmoke",
  grid_col   = "lightgrey",
  ...) {

  x <- left_join(x, int_info(unique(x$t)), by = c("t" = "tend"))
  x <- x %>% filter(!is.na(.data$interval)) %>%
    mutate(
      tz = factor(.data$tz, levels = sort(unique(.data$tz),
        decreasing = FALSE)),
      interval = factor(.data$interval, levels = levels(.data$interval)) )
  gg_ll <- ggplot(x, aes(x = .data[["interval"]], y = .data[["tz"]])) +
    geom_tile(aes(fill = .data[["LL"]]), colour = grid_col) +
    scale_fill_gradient(low = low_col, high = high_col) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    ylab(expression(t[z])) + xlab(expression(t))

  if (!is.null(x[["tz_var"]])) {
    gg_ll <- gg_ll + facet_wrap(~tz_var, scales = "free_y")
  }

  gg_ll + theme(legend.position = "none")

}

#' @rdname gg_laglead
#' @inherit gg_laglead
#' @export
gg_laglead.nested_fdf <- function(x, ...) {

  LL_df <- get_laglead(x)

  gg_laglead(LL_df, ...)

}
