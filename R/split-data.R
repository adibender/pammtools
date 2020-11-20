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
#' \dontrun{
#'  data("cgd", package = "frailtyHL")
#' cgd2 <- cgd %>%
#'  filter(id %in% c(1:5)) %>%
#'  select(id, tstart, tstop, enum, status, age)
#' cgd2 <- filter(cgd2, enum == 2)
#' ped_re <- split_data(
#'   formula = Surv(tstart, tstop, status) ~ age + enum,
#'   data = cgd2)
#' }
#' @seealso \code{\link[survival]{survSplit}}
#' @export
#' @keywords internal
split_data <- function(
  formula,
  data,
  cut      = NULL,
  max_time = NULL,
  ...) {

  dots_in         <- list(...)
  dots_in$formula <- formula

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


  ## standardize event time and status names
  proposed_names <- c("ped_start", "ped_time", "ped_status")
  ind <- ifelse(length(surv_vars) == 2, 2, 1):3
  proposed_names <- proposed_names[ind]
  if (any(proposed_names %in% names(data))) {
    stop(paste0("Error in attempt to rename provided time/status variables:
      Variables ",
      intersect(proposed_names, names(data)), " allready in data set."))
  }
  data <- rename(
    data,
    !!!set_names(
      surv_vars,
      as.list(proposed_names)))
  formula_cut <- update_formula(formula, proposed_names)

  # obtain interval breaks points
  cut <- get_cut(data, formula_cut, cut = cut, max_time = max_time)

  ## crate argument list to be passed to survSplit
  dots         <- list(...)
  dots$data    <- data
  dots$formula <- update_formula(formula, proposed_names)
  dots$cut     <- dots_in$cut <- cut
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
  if("ped_start" %in% colnames(split_df)) {
    split_df <- rename(split_df, !!!set_names("ped_start", "tstart"))
  }


  # Add variables for piece-wise exponential (additive) model
  split_df  <- split_df %>%
    mutate(
      ped_status = ifelse(.data$ped_status == 1 & .data$ped_time > max(cut),
          0L, .data$ped_status),
      ped_time   = pmin(.data$ped_time, max(cut)),
      offset     = log(.data$ped_time - .data$tstart)) %>%
    filter(!(.data$tstart == .data$ped_time))


  ## combine data with general interval info
  if(length(surv_vars) == 3) {
    info_cut <- split_df %>%
      select(one_of(c("tstart", "ped_time"))) %>% unique()
  } else {
    info_cut <- cut
  }
  int_info <- int_info(info_cut)
  split_df <- left_join(split_df, int_info, by = c("tstart" = "tstart"))

  ## rearrange columns
  move <- c(id_var, "tstart", "tend", "interval", "intmid", "intlen", "offset",
    "ped_time", "ped_status")
  split_df <- split_df %>%
    select(one_of(move), everything(),
      -one_of(c("intmid", "intlen", "ped_time")))

  ## set class and and attributes
  class(split_df) <- c("ped", class(split_df))
  attr(split_df, "breaks") <- cut
  attr(split_df, "id_var") <- dots_in$id <- id_var
  attr(split_df, "intvars") <- c(id_var, "tstart", "tend", "interval", "offset",
    "ped_status")
  attr(split_df, "trafo_args") <- dots_in

  split_df

}
