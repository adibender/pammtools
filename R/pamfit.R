#' Fit a piece-wise exponential additive model
#' 
#' Basically a wrapper around \code{\link[mgcv]{gam}}. However, we set 
#' \code{family=poisson()}, \code{offset=data$offset} and \code{method="REML"} 
#' by default. The first two can not be overriden. The \code{method} argument 
#' can be specified as usually, but defaults to \code{GCV.cp} in \code{\link[mgcv]{gam}}.
#' 
#' @import mgcv
#' @inherit mgcv::gam
#' @param ... Further arguments passed to \code{\link[mgcv]{gam}}.
#' @importFrom stats poisson
#' @rdname pam
#' @export
pam <- function(formula, data=list(), method="REML", ...) {

	dots <- list(...)
	dots$formula = formula
	dots$family  = poisson()
	dots$data    = data
	dots$offset  = data$offset

	pamfit         <- do.call(gam, dots)
	class(pamfit)  <- c("pam", class(pamfit))
	pamfit$cut     <- attr(data, "cut")
	pamfit$intvars <- attr(data, "intvars")

	pamfit

}

#' @rdname pam
#' @param x any R object.
#' @export
is.pam <- function(x) inherits(x, "pam")
