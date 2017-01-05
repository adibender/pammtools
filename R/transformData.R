#' Function to transform data of simple structure into piece-wise exponential data
#' format
#'
#' @inheritParams survival::survSplit
#' @param ... Further arguments passed to \code{\link[survival]{survSplit}}
#' @param max.end logical. Should the last interval span until the last 
#' observed censoring or event time (if larger than the largest specified 
#' cut point).
#' @import survival checkmate
#' @importFrom dplyr left_join
#' @return A data frame class \code{ped} in piece-wise exponential data format.
#' @examples
#' data("leuk2", package="bpcp")
#' head(leuk2)
#' leuk.ped <- split_data(Surv(time, status)~treatment, data=leuk2, 
#'    cut=c(0:5, 10, 40))
#' head(leuk.ped)
#' nrow(leuk.ped)
#' class(leuk.ped) # class ped (piece-wise exponential data)
#' @seealso \code{\link[survival]{survSplit}}
#' @export

split_data <- function(formula, data, cut=NULL, ..., max.end=FALSE) {

  ## assert that inputs have correct formats
  assert_class(formula, "formula")
  assert_data_frame(data, min.rows=2, min.cols=2)
  assert_numeric(cut, lower=0, finite=TRUE, any.missing=FALSE, min.len=1, 
    null.ok=TRUE)
  assert_flag(max.end)

  ## extract names for event time and status variables
  tvars     <- all.vars(update(formula, .~0))
  if(length(tvars)!=2) {
    stop(
      "The left hand side of the formula has only one variable.\n 
      At the moment a formula of the form Surv(time, event)~., is required.\n
      See ?Surv for more details.")
  }
  ## standardize event time and status names 
  data <- rename_(data, 
    .dots=setNames(c(tvars), c("time", "status")))
  formula <- as.formula(paste0("Surv(time, status)", paste0(formula[-2], collapse = "")))

  if(is.null(cut)) {
    cut <- unique(data[["time"]][data[["status"]]==1])
  }
  max.fail <- max(data[["time"]][data[["status"]]==1])
  max.time <- max(max(data[["time"]]), max(cut))


  # sort breaks in case they are not (so that interval factor variables will 
  # be in correct ordering)
  cut <- sort(cut)
  split.data <- survSplit(formula=formula, data=data, cut=cut)
  rm(data)

  ## Add variables for piece-wise exponential (additive) model
  # add last observation to cut if necessary
  if(max.end & (max.time > max(cut))) {
    cut <- c(cut, max.time)
  } else {
    split.data  <- split.data %>%
      mutate(
        status = ifelse(status == 1 & time > max(cut), 0, status),
        time   = pmin(time, max(cut)),
        offset = log(time - tstart)) %>%
      filter(!(tstart==time))
  }

  ## combine data with general interval info
  split.data <- left_join(split.data, int_info(brks=cut), by=c("tstart"="tstart"))

  ## set class and return
  class(split.data) <- c("ped", class(split.data))

  return(split.data)

}
