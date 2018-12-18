#' Function to transform data without time-dependent covariates into piece-wise
#' exponential data format
#'
#' @inheritParams as_ped
#' @import survival checkmate dplyr
#' @importFrom stats as.formula update
#' @importFrom purrr set_names
#' @examples
#' data("veteran", package="survival")
#' head(veteran)
#' ped <- split_data(Surv(time, status)~trt + age, data=veteran,
#'    cut=c(0, 100, 400), id="id")
#' head(ped)
#' class(ped) # class ped (piece-wise exponential data)
#' @seealso \code{\link[survival]{survSplit}}
#' @export
#' @keywords internal
split_data <- function(
  formula,
  data,
  cut      = NULL,
  max_time = NULL,
  ...) {

  ## assert that inputs have correct formats
  assert_class(formula, "formula")
  assert_data_frame(data, min.rows = 1, min.cols = 2)
  assert_numeric(cut, lower = 0, finite = TRUE, any.missing = FALSE,
    min.len = 1, null.ok = TRUE)
  assert_number(max_time, lower = 0, finite = TRUE, null.ok = TRUE)


  ## extract names for event time and status variables
  surv_vars <- all.vars(update(formula, .~0))
  vars <- if ("." %in% all.vars(formula)) {
      names(data)
    } else {
      all.vars(formula)
    }
  uvars <- union(surv_vars, vars)
  if (!all(uvars %in% vars)) {
    stop(paste("Variables provided in formula not in data set:",
      paste0(setdiff(uvars, vars), collapse = ", ")))
  }


  if (length(surv_vars) != 2) {
    stop(
      "Currently a formula of the form Surv(time, event)~., is required.\n
      See ?Surv for more details.")
  }
  ## standardize event time and status names
  proposed.names <- c("ped_time", "ped_status")
  if (any(proposed.names %in% names(data))) {
    stop(paste0("Error in attempt to rename provided time/status variables:
      Variables ",
      intersect(proposed.names, names(data)), " allready in data set."))
  }
  data    <- rename(data, !!!set_names(surv_vars, as.list(proposed.names)))
  formula <- as.formula(
    paste0("Surv(ped_time, ped_status)",
      paste0(formula[-2], collapse = "")))

  # obtain interval breaks points
  cut <- get_cut(data, formula, cut = cut, max_time = max_time)

  ## crate argument list to be passed to survSplit
  dots         <- list(...)
  dots$data    <- data
  dots$formula <- formula
  dots$cut     <- cut
  rm(data)

  # if id allready in the data set, remove id variable from dots but keep
  # id variable for later rearrangment
  if (!is.null(dots$id)) {
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
    } else {
      dots$id <- NULL
      dots$formula <- update(dots$formula, paste0("~ . + ", id_var))
    }
  }

  # create data in ped format
  split_df <- do.call(survSplit, args = dots)

  # Add variables for piece-wise exponential (additive) model
  split_df  <- split_df %>%
    mutate(
      ped_status = ifelse(.data$ped_status == 1 & .data$ped_time > max(cut),
          0L, .data$ped_status),
      ped_time   = pmin(.data$ped_time, max(cut)),
      offset     = log(.data$ped_time - .data$tstart)) %>%
    filter(!(.data$tstart == .data$ped_time))


  ## combine data with general interval info
  split_df <- left_join(split_df, int_info(cut), by = c("tstart" = "tstart"))

  ## rearrange columns
  move <- c(id_var, "tstart", "tend", "interval", "intmid", "intlen", "offset",
    "ped_time", "ped_status")
  split_df <- split_df %>%
    select(one_of(move), everything(),
      -one_of(c("intmid", "intlen", "ped_time")))

  ## set class and and attributes
  class(split_df) <- c("ped", class(split_df))
  attr(split_df, "breaks") <- cut
  attr(split_df, "id_var") <- id_var
  attr(split_df, "intvars") <- c(id_var, "tstart", "tend", "interval", "offset",
    "ped_status")

  split_df

}
