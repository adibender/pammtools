#' Time until extubation
#'
#' This is a preprocessed
#' subset of the data discussed in Heyard, et. al 2018 (and provided in a slightly
#' different format as \code{VAP_data} in the package \code{TBFmultinomial}).
#' In this package, the data is split in two parts, \code{extub_event} contains
#' time-to-event data and time-constant covariates and \code{extub_tdc} contains
#' the information on the time-dependent covariate SOFA score.
#' The data contains the following variables:
#' \describe{
#'  \item{ID}{Unique patient ID}
#'  \item{gender}{The patients' gender}
#'  \item{type}{Type of admission, either \code{Medical} or \code{Surgical}}
#'  \item{SAPSadmission}{SAPS score at admission}
#'  \item{time}{Time (days) until extubation}
#'  \item{extubation}{0 = no extubation/censoring, 1 = extubation}
#'  \item{day}{Exposure time, i.e., time at which the SOFA score was observed}
#' \item{SOFA}{The SOFA score at respective \code{day}s}
#'}
#' @aliases extub_tdc
"extub_event"

#' @rdname extub_event
"extub_tdc"
