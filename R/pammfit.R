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
#' @param trafo_args \strong{Deprecated.} A named list passed to \code{as_ped}
#' for inline data transformation. Convert your data with \code{as_ped()} before
#' calling \code{pamm()} instead.
#' @param engine Character name of the function that will be called to fit the
#' model. The intended entries are \code{"gam"} or \code{"bam"}
#' (both from package \code{mgcv}) or \code{"scam"} (from package \code{scam},
#' for shape-constrained PAMMs, e.g. monotone baseline hazards).
#' @import mgcv
#' @importFrom stats poisson
#' @rdname pamm
#' @seealso \code{\link[mgcv]{gam}}
#' @examples
#' ped <- tumor[1:100, ] %>%
#'  as_ped(Surv(days, status) ~ complications, cut = seq(0, 3000, by = 50))
#' pam <- pamm(ped_status ~ s(tend) + complications, data = ped)
#' summary(pam)
#' ## Deprecated: trafo_args inline transformation (use as_ped() instead)
#' # ped2 <- as_ped(tumor[1:100, ], Surv(days, status) ~ complications)
#' # pamm(ped_status ~ s(tend) + complications, data = ped2)
#' @export
pamm <- function(
  formula,
  data = list(),
  ...,
  trafo_args = NULL,
  engine = "gam"
) {
  dots <- list(...)
  dots$formula <- formula
  dots$family <- poisson()
  if (!is.null(trafo_args)) {
    .Deprecated(
      msg = paste0(
        "The 'trafo_args' argument of pamm() is deprecated and will be removed ",
        "in a future version.\n",
        "Please convert your data first:\n",
        "  ped <- as_ped(data, formula = ...)\n",
        "  pamm(formula, data = ped)"
      )
    )
    trafo_args$data <- data
    data <- do.call(split_data, trafo_args)
  }

  dots$data <- data
  dots$offset <- data$offset

  if (is.null(data$offset)) {
    warning(paste0(
      deparse(substitute(data)),
      " does not contain an offset. PAMM assumes a risk time of 1 for all subjects"
    ))
  }

  engine_fun <- switch(
    engine,
    gam = mgcv::gam,
    bam = mgcv::bam,
    scam = scam::scam,
    match.fun(engine)
  )
  pamm_fit <- do.call(engine_fun, dots)
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
