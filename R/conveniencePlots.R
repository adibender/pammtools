##' Plot smooth terms of gam objects
#' 
#' Given a gam model this convenience functions returns a plot of all 
#' smooth terms contained in the model. If more than one smooth is present, the 
#' different smooth are faceted. 
#' 
#' @inheritParams get_terms
#' @import ggplot2
#' @seealso \code{\link[mgcvtools]{tidy_smooth}}
#' @return A \code{\link[ggplot2]{ggplot2}} object.
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