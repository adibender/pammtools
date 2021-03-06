% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as-ped.R
\name{get_event_types}
\alias{get_event_types}
\title{Exctract event types}
\usage{
get_event_types(data, formula, censor_code)
}
\arguments{
\item{data}{Either an object inheriting from data frame or in case of
time-dependent covariates a list of data frames (of length 2), where the first data frame
contains the time-to-event information and static covariates while the second
(and potentially further data frames) contain information on time-dependent
covariates and the times at which they have been observed.}

\item{formula}{A two sided formula with a \code{\link[survival]{Surv}} object
on the left-hand-side and covariate specification on the right-hand-side (RHS).
The RHS can be an extended formula, which specifies how TDCs should be transformed
using specials \code{concurrent} and \code{cumulative}. The left hand-side can
be in start-stop-notation. This, however, is only used to create left-truncated
data and does not support the full functionality.}

\item{censor_code}{Specifies the value of the status variable that indicates censoring.
Often this will be \code{0}, which is the default.}
}
\description{
Given a formula that specifies the status variable of the outcome, this function
extracts the different event types (except for censoring, specified by
\code{censor_code}).
}
\keyword{internal}
