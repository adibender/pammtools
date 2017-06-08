unpam <- function(pamm) {
  class(pamm) <- class(pamm)[-1]
  pamm
}
repam <- function(x) {
  class(x) <- c("pamm", class(x))
  x
}

append_ped_attr <- function(pamm, ped) {

	attr_ped <- ped_attr(ped)
	pamm[names(attr_ped)] <- attr_ped

	pamm

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
#' @rdname pamm
#' @export
pamm <- function(formula, data=list(), method="REML", ...) {

	dots <- list(...)
	dots$formula = formula
	dots$family  = poisson()
	dots$data    = data
	dots$offset  = data$offset

	pamm_fit        <- do.call(gam, dots)
	class(pamm_fit) <- c("pamm", class(pamm_fit))
	pamm_fit        <- append_ped_attr(pamm_fit, data)

	pamm_fit

}


#' Check if object is of class pamm
#' 
#' @export
#' @param x Any R object.
is.pamm <- function(x) inherits(x, "pamm")


#' @param x An object of class \code{pamm} as returned by \code{\link{pamm}}.
#' @rdname pamm
#' @export 
print.pamm <- function(x, ...) {

	print(unpam(x), ...)

}

#' @rdname pamm
#' @param object An object of class \code{pamm} as returned by \code{\link{pamm}}.
#' @export
summary.pamm <- function(object, ...) {

	summary(unpam(object), ...)

}

#' @rdname pamm
#' @export
plot.pamm <- function(x, ...) {

	plot(unpam(x), ...)

}