#' Plot smooth 1d terms of gam objects
#'
#' Given a gam model this convenience function returns a plot of all
#' smooth terms contained in the model. If more than one smooth is present, the
#' different smooth are faceted.
#'
#' @inheritParams get_term
#' @param x A data frame or object of class \code{pamm}.
#' @param ... Further arguments passed to \code{\link{get_terms}}
#' @import ggplot2
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples
#' g1 <- mgcv::gam(Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data=iris)
#' gg_smooth(iris, g1, terms=c("Sepal.Width", "Petal.Length"))
#' @export
#' @seealso get_terms
gg_smooth <- function(x, ...) {
  UseMethod("gg_smooth", x)
}

#' @rdname gg_smooth
#' @inherit gg_smooth
#' @param fit A model object.
#' @export
gg_smooth.default <- function(x, fit, ...) {

  sobj <- get_terms(data = x, fit = fit, ...)

  ggsmooth <- ggplot(sobj, aes_string(x = "x", y = "eff", group = "term")) +
    geom_hline(yintercept = 0, lty = 3) +
    geom_line() +
    geom_ribbon(aes_string(ymin = "ci_lower", ymax = "ci_upper"), alpha = 0.2) +
    facet_wrap(~term, scales = "free_x") +
    ylab(expression(f[p](x[p]))) + xlab(expression(x[p]))

  return(ggsmooth)

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
      mutate(type = factor(.data$type,
        levels = c("ci_lower", "fit", "ci_upper")))
  }

  gg2d <- ggplot(df2d, aes_string(x = "x", y = "y", z = "fit")) +
    geom_raster(aes_string(fill = "fit")) +
    scale_fill_gradient2(
      name = expression(f(list(x, y))),
      low  = "steelblue", high = "firebrick2") +
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
#' data("lung", package="survival")
#' lung$inst <- as.factor(lung$inst) # for mgcv
#' ped <- lung %>% as_ped(Surv(time, status)~ph.ecog + inst, id="id")
#' pam <- mgcv::gam(ped_status ~ s(tend) + ph.ecog + s(inst, bs="re"),
#'  data=ped, family=poisson(), offset=offset)
#' gg_re(pam)
#' @seealso \code{\link{tidy_re}}
#' @export
gg_re <- function(x, ...) {

  re <- tidy_re(x, ...)
  ggplot(re, aes_string(sample = "fit")) +
    geom_abline(aes_string(intercept = "qqintercept", slope = "qqslope")) +
    geom_qq(distribution = "qnorm") +
    facet_wrap(~main) +
    theme_set(theme_bw())

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
gg_fixed <- function(x, intercept=FALSE, ...) {

  fixed_df <- tidy_fixed(x, intercept = intercept, ...)

  ggplot(fixed_df, aes_string(x = "variable", y = "coef", ymin = "ci_lower",
      ymax = "ci_upper")) +
    geom_hline(yintercept = 0, lty = 3) +
    geom_pointrange() +
    coord_flip() +
    ylab(expression(hat(beta) %+-% 1.96 %.% SE)) +
    xlab("")

}
