#' (Cumulative) (Step-) Hazard Plots.
#'
#' \code{geom_hazard} is an extension of the \code{geom_line}, and
#' is optimized for (cumulative) hazard plots. Essentially, it adds a (0,0)
#' row to the data, if not already the case. Stolen from the
#' \code{RmcdrPlugin.KMggplot2} (slightly modified).
#'
#' @seealso
#'   \code{\link[ggplot2]{geom_line}},
#'   \code{\link[ggplot2]{geom_step}}.
#' @inheritParams ggplot2::geom_line
#' @rdname geom_hazard
#' @importFrom ggplot2 layer GeomLine
#' @examples
#' library(ggplot2)
#' library(pammtools)
#' ped <- tumor[10:50,] %>% as_ped(Surv(days, status)~1)
#' pam <- mgcv::gam(ped_status ~ s(tend), data=ped, family = poisson(), offset = offset)
#' ndf <- make_newdata(ped, tend = unique(tend)) %>% add_hazard(pam)
#' # piece-wise constant hazards
#' ggplot(ndf, aes(x = tend, y = hazard)) +
#'  geom_vline(xintercept = c(0, ndf$tend[c(1, (nrow(ndf)-2):nrow(ndf))]), lty = 3) +
#'  geom_hline(yintercept = c(ndf$hazard[1:3], ndf$hazard[nrow(ndf)]), lty = 3) +
#'  geom_stephazard() +
#'  geom_step(col=2) +
#'  geom_step(col=2, lty = 2, direction="vh")
#'
#' # comulative hazard
#' ndf <- ndf %>% add_cumu_hazard(pam)
#' ggplot(ndf, aes(x = tend, y = cumu_hazard)) +
#'  geom_hazard() +
#'  geom_line(col=2) # doesn't start at (0, 0)
#'
#' # survival probability
#' ndf <- ndf %>% add_surv_prob(pam)
#' ggplot(ndf, aes(x = tend, y = surv_prob)) +
#'  geom_surv() +
#'  geom_line(col=2) # doesn't start at c(0,1)
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
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomHazard,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, ... )
  )
}


#' @rdname geom_hazard
#'
#' @format NULL
#' @usage NULL
#' @import ggplot2
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
geom_stephazard <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
  direction   = "vh",
  na.rm       = FALSE,
  show.legend = NA,
  inherit.aes = TRUE, ...) {
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
  setup_data = function(data, params) {
    row1   <- data %>% group_by(group) %>% slice(1)
    row1$x <- 0
    row1$y <- row1$y
    data   <- bind_rows(row1, data)
    data[order(data$PANEL, data$group, data$x), ]
  }
)



#' @inheritParams ggplot2::geom_line
#' @rdname geom_hazard
#' @importFrom ggplot2 layer GeomLine
#' @export
geom_surv <- function(
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
    geom = GeomSurv,
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
GeomSurv <- ggproto(
  "GeomSurv", GeomLine,
  setup_data = function(data, params) {
    row1   <- data %>% group_by(group) %>% slice(1)
    row1$x <- 0
    row1$y <- 1
    data   <- bind_rows(row1, data)
    data[order(data$group, data$x), ]
  }
)
