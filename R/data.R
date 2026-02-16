#' Survival data of critically ill ICU patients
#'
#' A data set containing the survival time (or hospital release time) among
#' other covariates.
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


#' Time-dependent covariates of the \code{\link{patient}} data set.
#'
#' This data set contains the time-dependent covariates (TDCs) for the \code{\link{patient}}
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
#' This is data simulated using the \code{\link[pammtools]{sim_pexp}} function.
#' It contains two time-constant and two time-dependent covariates (observed
#' on different exposure time grids). The code used for simulation is
#' contained in the examples of \code{?sim_pexp}.
#'
"simdf_elra"


#' Stomach area tumor data
#'
#' Information on patients treated for a cancer disease
#' located in the stomach area.
#' The data set includes:
#' \describe{
#'  \item{days}{Time from operation until death in days.}
#'  \item{status}{Event indicator (0 = censored, 1 = death).}
#'  \item{age}{The subject's age.}
#'  \item{sex}{The subject's sex (male/female).}
#'  \item{charlson_score}{Charlson comorbidity score, 1-6.}
#'  \item{transfusion}{Has subject received transfusions (no/yes).}
#'  \item{complications}{Did major complications occur during operation (no/yes).}
#'  \item{metastases}{Did the tumor develop metastases? (no/yes).}
#'  \item{resection}{Was the operation accompanied by a major resection (no/yes).}
#' }
#'
"tumor"


#' Time until staphylococcus aureaus infection in children, with possible recurrence
#'
#' This dataset originates from the Drakenstein child health study.
#' The data contains the following variables:
#' \describe{
#'  \item{id}{Randomly generated unique child ID}
#'  \item{t.start}{The time at which the child enters the risk set for the $k$-th event}
#'  \item{t.stop}{Time of $k$-th infection or censoring}.
#'  \item{enum}{Event number. Maximum of 6.}
#'  \item{hiv}{}
#' }
"staph"
