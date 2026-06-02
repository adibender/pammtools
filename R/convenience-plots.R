#' Plot smooth 1d terms of gam objects
#'
#' Given a gam model this convenience function returns a plot of its univariate
#' smooth terms. If \code{terms} is not specified, all univariate smooths are
#' plotted; otherwise only the requested ones (see \code{\link{get_terms}} for how
#' terms are matched). Different smooths are faceted. Smooths that are indexed by a
#' factor -- a factor \code{by}-variable or a factor-smooth interaction
#' (\code{bs = "fs"}/\code{"sz"}) -- are drawn in a single facet with one
#' coloured/filled curve per factor level.
#'
#' @param x A data frame or object of class \code{ped}.
#' @param ... Further arguments passed to \code{\link{get_terms}} (e.g.
#' \code{terms}).
#' @import ggplot2
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples
#' g1 <- mgcv::gam(Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data=iris)
#' gg_smooth(iris, g1, terms=c("Sepal.Width", "Petal.Length"))
#' # all univariate smooths (terms omitted)
#' gg_smooth(iris, g1)
#' # factor-by smooth: one coloured curve per Species
#' g2 <- mgcv::gam(Sepal.Length ~ s(Sepal.Width, by = Species), data = iris)
#' gg_smooth(iris, g2, terms = "Sepal.Width")
#' @export
#' @seealso get_terms
gg_smooth <- function(x, ...) {
  UseMethod("gg_smooth", x)
}

#' @rdname gg_smooth
#' @param fit A model object.
#' @export
gg_smooth.default <- function(x, fit, ...) {
  sobj <- get_terms(data = x, fit = fit, ...)

  if (nrow(sobj) == 0) {
    stop("No univariate smooth terms to plot.", call. = FALSE)
  }

  # smooths without a factor index (level == NA) are drawn as plain (black) curves,
  # factor-indexed smooths (factor `by`-variable, bs = "fs"/"sz") get one
  # coloured/filled curve per level. These are split into separate layers so that
  # the colour/fill scale and legend only contain the actual factor levels.
  plain_df <- sobj[is.na(sobj[["level"]]), , drop = FALSE]
  level_df <- sobj[!is.na(sobj[["level"]]), , drop = FALSE]

  ggsmooth <- ggplot(mapping = aes(x = .data[["x"]], y = .data[["eff"]])) +
    geom_hline(yintercept = 0, lty = 3)

  if (nrow(plain_df) > 0) {
    ggsmooth <- ggsmooth +
      geom_ribbon(
        data = plain_df,
        aes(
          ymin  = .data[["ci_lower"]],
          ymax  = .data[["ci_upper"]],
          group = .data[["term"]]
        ),
        alpha = 0.2
      ) +
      geom_line(data = plain_df, aes(group = .data[["term"]]))
  }

  if (nrow(level_df) > 0) {
    ggsmooth <- ggsmooth +
      geom_ribbon(
        data = level_df,
        aes(
          ymin  = .data[["ci_lower"]],
          ymax  = .data[["ci_upper"]],
          fill  = .data[["level"]],
          group = interaction(.data[["term"]], .data[["level"]])
        ),
        alpha = 0.2
      ) +
      geom_line(
        data = level_df,
        aes(
          colour = .data[["level"]],
          group  = interaction(.data[["term"]], .data[["level"]])
        )
      ) +
      scale_colour_discrete(name = "level") +
      scale_fill_discrete(name = "level")
  }

  ggsmooth +
    facet_wrap(~term, scales = "free_x") +
    ylab(expression(f[p](x[p]))) +
    xlab(expression(x[p]))
}

#' Plot tensor product effects
#'
#' Given a gam model this convenience function returns a \code{ggplot2} object
#' depicting 2d smooth terms specified in the model as heat/contour plots. If
#' more than one 2d smooth term is present individual terms are faceted.
#' @inheritParams tidy_smooth2d
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @examples
#' g <- mgcv::gam(Sepal.Length ~ te(Sepal.Width, Petal.Length), data=iris)
#' gg_tensor(g)
#' gg_tensor(g, ci=TRUE)
#' gg_tensor(update(g, .~. + te(Petal.Width, Petal.Length)))
#' @seealso \code{\link{tidy_smooth2d}}
#' @export
gg_tensor <- function(x, ci = FALSE, ...) {
  df2d <- tidy_smooth2d(x, ci = ci, se = ci, ...)
  if (ci) {
    df2d <- df2d %>%
      gather("type", "fit", .data$fit, .data$ci_lower, .data$ci_upper) %>%
      mutate(
        type = factor(.data$type, levels = c("ci_lower", "fit", "ci_upper"))
      )
  }

  gg2d <- ggplot(
    df2d,
    aes(x = .data[["x"]], y = .data[["y"]], z = .data[["fit"]])
  ) +
    geom_raster(aes(fill = .data[["fit"]])) +
    scale_fill_gradient2(
      name = expression(f(list(x, y))),
      low = "steelblue",
      high = "firebrick2"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_contour(col = "grey30")
  if (ci) {
    gg2d + facet_grid(main ~ type, scales = "free")
  } else {
    gg2d + facet_wrap(~main, scales = "free")
  }
}


#' Plot Normal QQ plots for random effects
#'
#' @inherit tidy_re
#' @import ggplot2
#' @examples
#' library(pammtools)
#' data("patient")
#' ped <- patient %>%
#'  dplyr::slice(1:100) %>%
#'  as_ped(Surv(Survdays, PatientDied)~ ApacheIIScore + CombinedicuID, id="CombinedID")
#' pam <- mgcv::gam(ped_status ~ s(tend) + ApacheIIScore + s(CombinedicuID, bs="re"),
#'  data=ped, family=poisson(), offset=offset)
#' gg_re(pam)
#' plot(pam, select = 2)
#' @seealso \code{\link{tidy_re}}
#' @export
gg_re <- function(x, ...) {
  re <- tidy_re(x, ...)
  ggplot(re, aes(sample = .data[["fit"]])) +
    geom_abline(aes(
      intercept = .data[["qqintercept"]],
      slope = .data[["qqslope"]]
    )) +
    geom_qq(distribution = stats::qnorm) +
    facet_wrap(~main)
}


#' Forrest plot of fixed coefficients
#'
#' @inherit tidy_fixed
#' @param intercept Logical, indicating whether intercept term should be included.
#' Defaults to \code{FALSE}.
#' @import ggplot2
#' @seealso \code{\link{tidy_fixed}}
#' @examples
#' g <- mgcv::gam(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
#'  data=iris)
#' gg_fixed(g, intercept=TRUE)
#' gg_fixed(g)
#' @export
gg_fixed <- function(x, intercept = FALSE, ...) {
  fixed_df <- tidy_fixed(x, intercept = intercept, ...)

  ggplot(
    fixed_df,
    aes(
      x = .data[["variable"]],
      y = .data[["coef"]],
      ymin = .data[["ci_lower"]],
      ymax = .data[["ci_upper"]]
    )
  ) +
    geom_hline(yintercept = 0, lty = 3) +
    geom_pointrange() +
    coord_flip() +
    ylab(expression(hat(beta) %+-% 1.96 %.% SE)) +
    xlab("")
}
