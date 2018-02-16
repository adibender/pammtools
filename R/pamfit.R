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
#' @inheritParams mgcv::gam
#' @param ... Further arguments passed to \code{\link[mgcv]{gam}}.
#' @import mgcv
#' @importFrom stats poisson
#' @rdname pamm
#' @seealso \code{\link[mgcv]{gam}}
#' @export
pamm <- function(formula, data=list(), method="REML", ..., trafo.args=NULL) {

	dots <- list(...)
	dots$formula = formula
	dots$family  = poisson()
  if (!is.null(trafo.args)) {
    trafo.args$data <- data
    data <- do.call(split_data, trafo.args)
  }
	dots$data    = data
	dots$offset  = data$offset

	pamm_fit        <- do.call(gam, dots)
	class(pamm_fit) <- c("pamm", class(pamm_fit))
	pamm_fit        <- append_ped_attr(pamm_fit, data)

	pamm_fit

}


#' Check if object is of class pamm
#'
#' @param x Any R object.
#' @rdname pamm
#' @keywords internal
#' @export
is.pamm <- function(x) inherits(x, "pamm")


#' @rdname pamm
#' @keywords internal
#' @export
print.pamm <- function(x, ...) {

	print(unpam(x), ...)

}

#' @rdname pamm
#' @param object An object of class \code{pamm} as returned by \code{\link{pamm}}.
#' @keywords internal
#' @export
summary.pamm <- function(object, ...) {

	summary(unpam(object), ...)

}

#' @rdname pamm
#' @keywords internal
#' @export
plot.pamm <- function(x, ...) {

	plot(unpam(x), ...)

}
