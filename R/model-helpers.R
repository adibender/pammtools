## Internal generics that abstract over the supported model classes
## (glm, mgcv::gam/bam, scam::scam). All downstream calculations are based on
## the triplet (X, coefs, V) with X = make_X(object, newdata), coefs =
## get_coefs(object) and V = get_Vp(object), such that
## X %*% coefs is the linear predictor (minus offset) and
## sqrt(rowSums((X %*% V) * X)) its standard error.

#' Extract model coefficients on the scale of the design matrix
#'
#' Returns the coefficient vector \code{coefs} such that
#' \code{make_X(object, newdata) \%*\% coefs} yields the linear predictor.
#' For most models this is simply \code{coef(object)}. For \code{scam} models,
#' however, \code{coef()} returns the coefficients on the underlying
#' unconstrained scale, while the linear predictor is calculated from the
#' re-parametrized (partially exponentiated) coefficients
#' \code{object$coefficients.t}.
#'
#' @param object A fitted model object.
#' @param ... Further arguments passed to methods.
#' @keywords internal
get_coefs <- function(object, ...) {
  UseMethod("get_coefs", object)
}

#' @rdname get_coefs
#' @importFrom stats coef
#' @keywords internal
get_coefs.default <- function(object, ...) {
  coef(object)
}

#' @rdname get_coefs
#' @keywords internal
get_coefs.scam <- function(object, ...) {
  object$coefficients.t
}

#' Extract the (Bayesian) covariance matrix of the model coefficients
#'
#' Returns the covariance matrix that matches the coefficients returned by
#' \code{\link{get_coefs}}. For \code{mgcv} models this is the Bayesian
#' posterior covariance matrix \code{object$Vp}, for \code{scam} models the
#' covariance matrix of the re-parametrized coefficients \code{object$Vp.t}
#' and \code{vcov(object)} otherwise.
#'
#' @inheritParams get_coefs
#' @keywords internal
get_Vp <- function(object, ...) {
  UseMethod("get_Vp", object)
}

#' @rdname get_Vp
#' @importFrom stats vcov
#' @keywords internal
get_Vp.default <- function(object, ...) {
  vcov(object)
}

#' @rdname get_Vp
#' @keywords internal
get_Vp.gam <- function(object, ...) {
  object$Vp
}

#' @rdname get_Vp
#' @keywords internal
get_Vp.scam <- function(object, ...) {
  object$Vp.t
}

#' Draw coefficients from their approximate posterior distribution
#'
#' Simulation based confidence intervals are calculated by drawing coefficient
#' vectors from their asymptotic (posterior) distribution, a multivariate
#' normal with mean \code{\link{get_coefs}} and covariance \code{\link{get_Vp}}.
#' For \code{scam} models this means that draws are obtained on the scale of
#' the re-parametrized (partially exponentiated) coefficients, i.e., based on
#' the same normal approximation that underlies the reported standard errors of
#' the model (the exact posterior of the constrained coefficients is not
#' Gaussian, so individual draws may violate the shape constraints slightly).
#'
#' @inheritParams get_coefs
#' @param nsim Number of draws.
#' @return A matrix with \code{nsim} rows, one coefficient vector per row, on
#' the scale of the design matrix returned by \code{\link{make_X}}.
#' @importFrom mvtnorm rmvnorm
#' @keywords internal
sample_coefs <- function(object, nsim, ...) {
  UseMethod("sample_coefs", object)
}

#' @rdname sample_coefs
#' @keywords internal
sample_coefs.default <- function(object, nsim, ...) {
  mvtnorm::rmvnorm(nsim, mean = get_coefs(object), sigma = get_Vp(object))
}

#' Draw hazard trajectories from a model's sampling distribution
#'
#' Internal seam used by the simulation-based confidence interval helpers
#' (\code{\link{get_sim_ci}}, \code{get_sim_ci_cumu}, \code{get_sim_ci_surv}).
#' It returns a matrix of \code{nsim} draws of the (response-scale) hazard, one
#' column per draw and one row per row of \code{newdata}. The default method
#' draws coefficient vectors via \code{\link{sample_coefs}} and evaluates the
#' linear predictor \code{make_X(object, newdata) \%*\% z}; other backends (e.g.
#' a bootstrap ensemble that has no coefficient covariance) can provide their own
#' method to obtain simulation-based intervals from the same machinery.
#'
#' @inheritParams get_coefs
#' @param newdata A data frame for which hazards are predicted.
#' @param nsim Number of draws.
#' @param sim_coef_mat Optional pre-drawn coefficient matrix (as returned by
#' \code{\link{sample_coefs}}); used to share one set of draws across groups.
#' @return A numeric matrix with \code{nrow(newdata)} rows and \code{nsim}
#' columns of hazard draws on the response scale.
#' @keywords internal
sim_hazard <- function(object, newdata, nsim = 100L, ...) {
  UseMethod("sim_hazard", object)
}

#' @rdname sim_hazard
#' @keywords internal
sim_hazard.default <- function(object, newdata, nsim = 100L, sim_coef_mat = NULL, ...) {
  X <- make_X(object, newdata, ...)
  if (is.null(sim_coef_mat)) {
    sim_coef_mat <- sample_coefs(object, nsim)
  }
  matrix(
    apply(sim_coef_mat, 1, function(z) exp(drop(X %*% z))),
    nrow = nrow(newdata)
  )
}
