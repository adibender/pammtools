##' Plot smooth 1d terms of gam objects
#' 
#' Given a gam model this convenience function returns a plot of all 
#' smooth terms contained in the model. If more than one smooth is present, the 
#' different smooth are faceted. 
#' 
#' @inheritParams get_term
#' @param ... Further arguments passed to \code{\link{get_terms}}
#' @import ggplot2
#' @return A \code{\link[ggplot2]{ggplot2}} object.
#' @examples
#' g1 <- mgcv::gam(Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data=iris)
#' gg_smooth(iris, g1, terms=c("Sepal.Width", "Petal.Length"))
#' @export 
gg_smooth <- function(data, fit, ...) {

	sobj <- get_terms(data=data, fit=fit, ...)

	ggsmooth <- ggplot(sobj, aes(x=x, y=eff, group=term)) + 
		geom_hline(yintercept = 0, lty=3) +
		geom_line() + 
		geom_ribbon(aes(ymin=ci.lower, ymax=ci.upper), alpha=0.2) + 
		facet_wrap(~term, scales="free_x") + 
		ylab(expression(f[j](x[j]))) + xlab(expression(x[j]))

	return(ggsmooth)

}

#' Plot tensor product effects 
#' 
#' Given a gam model this convenience function returns a \code{ggplot2} object 
#' depicting 2d smooth terms specified in the model as heat/contour plots. If 
#' more than one 2d smooth term is present individual terms are faceted.
#' @inheritParams get_term
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom magrittr "%<>%"
#' @examples
#' library(mgcv)
#' g <- gam(Sepal.Length ~ te(Sepal.Width, Petal.Length), data=iris)
#' gg_tensor(g)
#' gg_tensor(g, ci=TRUE)
#' gg_tensor(update(g, .~. + te(Petal.Width, Petal.Length)))
#' @seealso \code{\link{tidy_smooth2d}}
#' @export
gg_tensor <- function(fit, ci=FALSE, ...) {

	df2d <- tidy_smooth2d(fit, ci=ci, se=ci, ...)
	if (ci) {
		df2d %<>% gather(type, fit, fit, low, high) %>% 
			mutate(
				type = factor(
					type, 
					levels = c("low", "fit", "high"), 
					labels = c("lower", "fit", "upper")))
	}

	gg2d <- ggplot(df2d, aes_string(x="x", y="y", z="fit")) + 
		geom_raster(aes(fill=fit)) + 
		scale_fill_gradient2(
			name = expression(f(list(x,y))),
			low  = "steelblue", high = "firebrick2") +
		geom_contour(col="grey30")
		if(ci) {
			gg2d + facet_grid(main ~ type, scales="free")
		} else {
			gg2d + facet_wrap(~main, scales="free")
		}

}

