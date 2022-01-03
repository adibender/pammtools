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
  pamm[["attr_ped"]] <- attr_ped

  pamm

}


#' Fit a piece-wise exponential additive model
#'
#' A thin wrapper around \code{\link[mgcv]{gam}}, however, some arguments are
#' prespecified:
#' \code{family=poisson()} and \code{offset=data$offset}.
#' These two can not be overwritten. In many cases it will also be advisable to
#' set \code{method="REML"}.
#'
#' @inheritParams mgcv::gam
#' @param ... Further arguments passed to \code{engine}.
#' @param trafo_args A named list. If data is not in PED format, \code{as_ped}
#' will be called internally with arguments provided in \code{trafo_args}.
#' @param engine Character name of the function that will be called to fit the
#' model. The intended entries are either \code{"gam"} or \code{"bam"}
#' (both from package \code{mgcv}).
#' @import mgcv
#' @importFrom stats poisson
#' @rdname pamm
#' @seealso \code{\link[mgcv]{gam}}
#' @examples
#' ped <- tumor[1:100, ] %>%
#'  as_ped(Surv(days, status) ~ complications, cut = seq(0, 3000, by = 50))
#' pam <- pamm(ped_status ~ s(tend) + complications, data = ped)
#' summary(pam)
#' ## Alternatively
#' pamm(
#'  ped_status ~ s(tend) + complications,
#'  data = tumor[1:100, ],
#' trafo_args = list(formula = Surv(days, status)~complications))
#' @export
pamm <- function(
  formula,
  data       = list(),
  ...,
  trafo_args = NULL,
  engine     = "gam") {

  dots <- list(...)
  dots$formula <- formula
  dots$family  <- poisson()
  if (!is.null(trafo_args)) {
    trafo_args$data <- data
    data <- do.call(split_data, trafo_args)
  }
  dots$data   <- data
  dots$offset <- data$offset

  pamm_fit        <- do.call(engine, dots)
  class(pamm_fit) <- c("pamm", class(pamm_fit))
  # pamm_fit        <- append_ped_attr(pamm_fit, data)
  pamm_fit[["trafo_args"]] <- attr(data, "trafo_args")
  ind_attr_keep <- !(names(attributes(data)) %in%
    c("names", "row.names", "trafo_args", "class"))
  pamm_fit[["attr_ped"]] <- attributes(data)[ind_attr_keep]

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
