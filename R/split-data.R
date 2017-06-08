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
#' @importFrom magrittr "%<>%"
#' @return A data frame class \code{ped} in piece-wise exponential data format.
#' @examples
#' data("veteran", package="survival")
#' head(veteran)
#' ped <- split_data(Surv(time, status)~trt + age, data=veteran, 
#'    cut=c(0, 100, 400), id="id")
#' head(ped)
#' class(ped) # class ped (piece-wise exponential data)
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
  tvars <- all.vars(update(formula, .~0))
  vars <- if("." %in% all.vars(formula)) {
      names(data)
    } else {
      all.vars(formula)
    } 
  uvars <- union(tvars, vars)
  if(!all(uvars %in% vars)) {
    stop(paste("Variables provided in formula not in data set:", 
      paste0(setdiff(uvars, vars), collapse=", ")))
  }

  
  if(length(tvars)!=2) {
    stop(
      "Currently a formula of the form Surv(time, event)~., is required.\n
      See ?Surv for more details.")
  }
  ## standardize event time and status names 
  proposed.names <- c("ped_time", "ped_status")
  if(any(proposed.names %in% names(data))) {
    stop(paste0("Error in attempt to rename provided time/status variables: Variables 
      ", intersect(proposed.names, names(data)), " allready in data set."))
  }
  data <- rename_(data, 
    .dots=setNames(c(tvars), c(proposed.names)))
  formula <- as.formula(paste0("Surv(ped_time, ped_status)", paste0(formula[-2], collapse = "")))

  if(is.null(cut)) {
    cut <- unique(data[["ped_time"]][data[["ped_status"]]==1])
  }
  max.fail <- max(data[["ped_time"]][data[["ped_status"]]==1])
  max.time <- max(max(data[["ped_time"]]), max(cut))


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
    id_var <- dots$id
  } else {
    id_var  <- "id"
    dots$id <- id_var
  }

  if (id_var %in% names(dots$data)) {
    if (length(unique(dots$data[[id_var]])) != nrow(dots$data)) {
      stop(paste0("Specified ID variable (", id_var, ") must have same number of 
        unique values as number of rows in 'data'."))
    }
    if (id_var %in% vars) {
      dots$id <- NULL 
    }
  }

  # create data in ped format
  split_df <- do.call(survSplit, args=dots)

  ## Add variables for piece-wise exponential (additive) model
  split_df  <- split_df %>%
    mutate(
      ped_status = ifelse(ped_status == 1 & ped_time > max(cut), 0, ped_status),
      ped_time   = pmin(ped_time, max(cut)),
      offset     = log(ped_time - tstart)) %>%
    filter(!(tstart==ped_time))


  ## combine data with general interval info
  split_df <- left_join(split_df, int_info(cut), by=c("tstart"="tstart"))

  ## rearrange columns 
  move <- c(id_var, "tstart", "tend", "interval", "intmid", "intlen", "offset",
    "ped_time", "ped_status")
  split_df %<>% select(one_of(move), everything(), -intmid, -intlen, -ped_time)
  

  ## set class and and attributes
  class(split_df) <- c("ped", class(split_df))
  attr(split_df, "cut") <- cut 
  attr(split_df, "id_var") <- id_var
  attr(split_df, "intvars") <- c(id_var, "tstart", "tend", "interval", "offset", 
    "ped_status")

  return(split_df)

}
