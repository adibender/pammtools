% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-functions.R
\name{get_cumu_hazard}
\alias{get_cumu_hazard}
\title{Calculate cumulative hazard}
\usage{
get_cumu_hazard(
  newdata,
  object,
  ci = TRUE,
  ci_type = c("default", "delta", "sim"),
  time_var = NULL,
  se_mult = 2,
  interval_length = "intlen",
  nsim = 100L,
  ...
)
}
\arguments{
\item{newdata}{ A data frame or list containing the values of the model covariates at which predictions
  are required. If this is not provided then predictions corresponding to the
                original data are returned. If \code{newdata} is provided then
                it should contain all the variables needed for prediction: a
                warning is generated if not. See details for use with \code{link{linear.functional.terms}}. }

\item{object}{ a fitted \code{gam} object as produced by \code{gam()}.
                }

\item{ci}{\code{logical}. Indicates if confidence intervals should be
calculated. Defaults to \code{TRUE}.}

\item{ci_type}{The method by which standard errors/confidence intervals
will be calculated. Default transforms the linear predictor at
respective intervals. \code{"delta"} calculates CIs based on the standard
error calculated by the Delta method. \code{"sim"} draws the
property of interest from its posterior based on the normal distribution of
the estimated coefficients. See \href{https://adibender.github.io/simpamm/confidence-intervals.html}{here}
for details and empirical evaluation.}

\item{time_var}{Name of the variable used for the baseline hazard. If
not given, defaults to \code{"tend"} for \code{\link[mgcv]{gam}} fits, else
\code{"interval"}. The latter is assumed to be a factor, the former
numeric.}

\item{se_mult}{Factor by which standard errors are multiplied for calculating
the confidence intervals.}

\item{interval_length}{The variable in newdata containing the interval lengths.
Can be either bare unquoted variable name or character. Defaults to \code{"intlen"}.}

\item{...}{Further arguments passed to \code{\link[mgcv]{predict.gam}} and
\code{\link{get_hazard}}}
}
\description{
Calculate cumulative hazard
}
\keyword{internal}
