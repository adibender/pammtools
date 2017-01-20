#' Function to transform data without time-dependent covariates into piece-wise 
#' exponential data format
#'
#' @inheritParams survival::survSplit
#' @param ... Further arguments passed to \code{\link[survival]{survSplit}}
#' @param max.end logical. Should the last interval span until the last 
#' observed censoring or event time (if larger than the largest specified 
#' cut point).
#' @import survival checkmate dplyr
#' @importFrom stats as.formula setNames update
#' @return A data frame class \code{ped} in piece-wise exponential data format.
#' @examples
#' data("leuk2", package="bpcp")
#' head(leuk2)
#' leuk.ped <- split_data(Surv(time, status)~treatment, data=leuk2, 
#'    cut=c(0:5, 10, 40), id="id")
#' head(leuk.ped)
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
  vars <- if("." %in% all.vars(formula)) {
      colnames(data)
    } else {
      all.vars(formula)
    } 
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


  # sort interval cut points in case they are not (so that interval factor 
  # variables will be in correct ordering)
  cut <- sort(cut)
  # add last observation to cut if necessary
  if(max.end & (max.time > max(cut))) {
    cut <- c(cut, max.time)
  }

  ## crate argument list to be passed to survSplit 
  dots <- list(...)
  dots$data    <- data
  dots$formula <- formula
  dots$cut     <- cut
  rm(data) 

  # if id allready in the data set, remove id variable from dots but keep 
  # id variable for later rearrangment 
  if(!is.null(dots$id)) {
    id.var <- dots$id
    if(id.var %in% names(dots$data) & id.var %in% vars) {
      dots$id <- NULL
    }
  }

  # crate data in ped format
  split.data <- do.call(survSplit, args=dots)

  ## Add variables for piece-wise exponential (additive) model
  split.data  <- split.data %>%
  mutate(
    status = ifelse(status == 1 & time > max(cut), 0, status),
    time   = pmin(time, max(cut)),
    offset = log(time - tstart)) %>%
  filter(!(tstart==time))

  ## combine data with general interval info
  split.data <- left_join(split.data, int_info(brks=cut), by=c("tstart"="tstart"))

  ## rearrange columns 
  move <- c("tstart", "tend", "interval", "intmid", "intlen", "offset",
    "time", "status")
  if(exists("id.var")) move <- c(id.var, move)
  split.data <- select(split.data, one_of(move), everything())

  ## set class and and attributes
  class(split.data) <- c("ped", class(split.data))
  attr(split.data, "cut") <- cut 
  attr(split.data, "intvars") <- c("id", "tstart", "tend", "intlen", "intmid", 
    "interval", "offset", "time", "status")

  return(split.data)

}
