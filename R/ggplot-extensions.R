# Stolen from the \code{RmcdrPlugin.KMggplot2} (slightly modified)

#' Step ribbon plots.
#'
#' \code{geom_stepribbon} is an extension of the \code{geom_ribbon}, and
#' is optimized for Kaplan-Meier plots with pointwise confidence intervals
#' or a confidence band.
#'
#' @section Aesthetics:
#' \Sexpr[results=rd,stage=build]{ggplot2:::rd_aesthetics("geom", "ribbon")}
#'
#' @seealso
#'   \code{\link[ggplot2]{geom_ribbon}} \code{geom_stepribbon}
#'   inherits from \code{geom_ribbon}.
#' @inheritParams ggplot2:::geom_ribbon
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
#' @author Triad Sou <triadsou@@gmail.com>
#' @export
geom_stepribbon <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
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
    params      = list(na.rm = na.rm, ... )
  )

}

#' @rdname geom_stepribbon
#' @format NULL
#' @usage NULL
#' @export

GeomStepribbon <- ggproto(
  "GeomStepribbon", GeomRibbon,

  extra_params = c("na.rm"),

  draw_group = function(data, panel_scales, coord, na.rm = FALSE) {

    if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
    data   <- rbind(data, data)
    data   <- data[order(data$x), ]
    data$x <- c(data$x[2:nrow(data)], NA)
    data   <- data[complete.cases(data["x"]), ]
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)

  }

)
