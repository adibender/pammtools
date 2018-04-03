#' Survival data of critically ill ICU patients
#'
#' A data set containing the survival time (or hospital release time) among
#' other covariates. This is a subset of the data discussed in
#' \href{https://doi.org/10.1093/biostatistics/kxy003}{Bender et. al., 2018}.
#' The full data is available \href{https://github.com/adibender/elra-biostats}{here}.
#' The following variables are provided:
#' \describe{
#'  \item{Year}{The year of ICU Admission}
#'  \item{CombinedicuID}{Intensive Care Unit (ICU) ID}
#'  \item{CombinedID}{Patient identificator}
#'  \item{Survdays}{Survival time of patients. Here it is assumed that patients
#'    survive until t=30 if released from hospital.}
#'  \item{PatientDied}{Status indicator; 1=death, 0=censoring}
#'  \item{survhosp}{Survival time in hospital. Here it is assumed that patients
#'    are censored at time of hospital release (potentially informative)}
#'  \item{Gender}{Male or female}
#'  \item{Age}{The patients age at Admission}
#'  \item{AdmCatID}{Admission category: medical, surgical elective or surgical emergency}
#'  \item{ApacheIIScore}{The patient's Apache II Score at Admission}
#'  \item{BMI}{Patient's Body Mass Index}
#'  \item{DiagID2}{Diagnosis at admission in 9 categories}  }
"patient"


#' Time-depedent covariates of the \code{\link{patient}} data set.
#'
#' This data set containt the time-dependent covariates (TDCs) for the \code{\link{patient}}
#' data set. Note that nutrition was protocoled for at most 12 days after
#' ICU admission. The data set includes:
#' \describe{
#'  \item{CombinedID}{Unique patient identifier. Can be used to merge with
#'  \code{\link{patient}} data}
#'  \item{Study_Day}{The calendar (!) day at which calories (or proteins) were
#'  administered}
#'  \item{caloriesPercentage}{The percentage of target calories supplied to the
#'    patient by the ICU staff}
#'  \item{proteinGproKG}{The amount of protein supplied to the patient by the
#'    ICU staff}}
"daily"


#' Simulated data with cumulative effects
#'
"simdf_elra"
