# Stolen from the \code{RmcdrPlugin.KMggplot2} (slightly modified)

#' PAMM (Cumulative) (Step-)Hazard Plots.
#'
#' \code{geom_hazard} is an extension of the \code{geom_line}, and
#' is optimized for (cumulative) hazard plots. Essentially, it add a (0,0)
#' row to the data, if not already the case.
#'
#' @section Aesthetics:
#' \Sexpr[results=rd,stage=build]{ggplot2:::rd_aesthetics("geom", "line")}
#'
#' @seealso
#'   \code{\link[ggplot2]{geom_line}}.
#' @inheritParams ggplot2::geom_line
#' @rdname geom_hazard
#' @importFrom ggplot2 layer GeomStep
#' @export
geom_hazard <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
  na.rm       = FALSE,
  show.legend = NA,
  inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHazard,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_hazard
#' @format NULL
#' @usage NULL
#' @export

GeomHazard <- ggproto(
  "GeomHazard", GeomLine,
  setup_data = function(data, params) {
    row1   <- data %>% group_by(group) %>% slice(1)
    row1$x <- 0
    row1$y <- 0
    data   <- bind_rows(row1, data)
    data[order(data$group, data$x), ]
  }
)

#' @inheritParams ggplot2::geom_step
#' @rdname geom_hazard
#' @importFrom ggplot2 layer GeomStep
#' @export
geom_stephazard <- function(mapping = NULL, data = NULL, stat = "identity",
                      position = "identity", direction = "vh",
                      na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepHazard,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      direction = direction,
      na.rm     = na.rm,
      ...
    )
  )
}

#' @rdname geom_hazard
#' @format NULL
#' @usage NULL
#' @export
GeomStepHazard <- ggproto(
  "GeomStepHazard",
  GeomStep,
  draw_panel = function(data, panel_params, coord, direction = "vh") {
    data <- plyr::ddply(data, "group", stairstep, direction = direction)
    GeomPath$draw_panel(data, panel_params, coord)
  },
  setup_data = function(data, params) {
    row1   <- data %>% group_by(group) %>% slice(1)
    row1$x <- 0
    row1$y <- 0
    data   <- bind_rows(row1, data)
    data[order(data$PANEL, data$group, data$x), ]
  }
)

# copied from https://github.com/tidyverse/ggplot2/blob/master/R/geom-path.r
# @keyword internal
stairstep <- function(data, direction="hv") {
  direction <- match.arg(direction, c("hv", "vh"))
  data <- as.data.frame(data)[order(data$x), ]
  n <- nrow(data)

  if (n <= 1) {
    # Need at least one observation
    return(data[0, , drop = FALSE])
  }

  if (direction == "vh") {
    xs <- rep(1:n, each = 2)[-2*n]
    ys <- c(1, rep(2:n, each = 2))
  } else {
    ys <- rep(1:n, each = 2)[-2*n]
    xs <- c(1, rep(2:n, each = 2))
  }

  data.frame(
    x = data$x[xs],
    y = data$y[ys],
    data[xs, setdiff(names(data), c("x", "y"))]
  )
}
