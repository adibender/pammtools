# Stolen from the \code{RmcdrPlugin.KMggplot2} (slightly modified)

#' Step ribbon plots.
#'
#' \code{geom_stepribbon} is an extension of the \code{geom_ribbon}, and
#' is optimized for Kaplan-Meier plots with pointwise confidence intervals
#' or a confidence band. The default \code{direction}-argument \code{"hv"} is 
#' appropriate for right-continuous step functions like the hazard rates etc 
#' returned by \code{pammtools}.
#'
#' @seealso
#'   \code{\link[ggplot2]{geom_ribbon}} \code{geom_stepribbon}
#'   inherits from \code{geom_ribbon}.
#' @inheritParams ggplot2:::geom_ribbon
#' @inheritParams ggplot2:::geom_step
#' @examples
#' library(ggplot2)
#' huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
#' h <- ggplot(huron, aes(year))
#' h + geom_stepribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#'     geom_step(aes(y = level))
#' h + geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#'     geom_line(aes(y = level))
#' @rdname geom_stepribbon
#' @importFrom ggplot2 layer GeomRibbon
#' @export
geom_stepribbon <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
  direction   = "hv", 
  na.rm       = FALSE,
  show.legend = NA,
  inherit.aes = TRUE, ...) {

  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepribbon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, direction = direction, ... )
  )

}

#' @rdname geom_stepribbon
#' @importFrom ggplot2 ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomStepribbon <- ggproto(
  "GeomStepribbon", GeomRibbon,

  extra_params = c("na.rm"),

  draw_group = function(data, panel_scales, coord, na.rm = FALSE, direction = "hv") {

    if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
    data   <- rbind(data, data)
    data   <- data[order(data$x), ]
    data   <- ggplot2_stairstep(data[complete.cases(data["x"]), ], 
                               direction = direction)
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = na.rm)
  }

)


# code adapted from 
# https://github.com/tidyverse/ggplot2/blob/9741da5050f81b7b5c012c56d02f45fc93d68f89/R/geom-path.r#L320
ggplot2_stairstep <- function(data, direction =  c("hv", "vh", "mid")) {
  direction <- match.arg(direction)
  data <- as.data.frame(data)[order(data$x), ]
  n <- nrow(data)
  if (n <= 1) {
    return(data[0, , drop = FALSE])
  }
  if (direction == "vh") {
    xs <- rep(1:n, each = 2)[-2 * n]
    ys <- c(1, rep(2:n, each = 2))
  }
  if (direction == "hv") {
    xs <- c(1, rep(2:n, each = 2))
    ys <- rep(1:n, each = 2)[-2 * n]
  }
  if (direction == "mid") {
    xs <- rep(1:(n - 1), each = 2)
    ys <- rep(1:n, each = 2)
  }
  
  ymin <- c(data$ymin[ys])
  ymax <- c(data$ymax[ys])
  if (direction == "mid") {
    gaps <- data$x[-1] - data$x[-n]
    mid_x <- data$x[-n] + gaps/2
    x <- c(data$x[1], mid_x[xs], data$x[n])
    data_attr <- data[c(1, xs, n), 
                      setdiff(names(data), c("x", "ymin", "ymax"))]
  } else {
    x <- data$x[xs]
    ymin <- data$ymin[ys]
    ymax <- data$ymax[ys]
    data_attr <- data[xs, setdiff(names(data), c("x", "ymin", "ymax"))]
  }
  cbind(data.frame(x = x, ymin = ymin, ymax = ymax), data_attr)
}
