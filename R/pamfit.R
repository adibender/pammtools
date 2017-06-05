unpam <- function(pam) {
  class(pam) <- class(pam)[-1]
  pam
}
repam <- function(x) {
  class(x) <- c("ped", class(x))
  x
}


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

#' Check if object is of class pam
#' @export
#' @param x Any R object.
is.pam <- function(x) inherits(x, "pam")


#' @param x An object of class \code{pam} as returned by \code{\link{pam}}.
#' @rdname pam
#' @export 
print.pam <- function(x, ...) {

	print(unpam(x), ...)

}

#' @rdname pam
#' @param object An object of class \code{pam} as returned by \code{\link{pam}}.
#' @export
summary.pam <- function(object, ...) {

	summary(unpam(object), ...)

}

#' @rdname pam
#' @export
plot.pam <- function(x, ...) {

	plot(unpam(x), ...)

}