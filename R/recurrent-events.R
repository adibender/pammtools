#' Split data to obtain recurrent event data in PED format
#'
#' Currently, the input data must be in start-stop notation for each spell and
#' contain a colum that indicates the spell (event number).
#' @inherit split_data
#' @inheritParams get_cut
#' @param episode_var A character indicating the column in data that indicates the
#' event/episode number for recurrent events.
#' @param event The value that encodes the occurrence of an event in the data set.
#' @param timescale Defines the timescale for the recurrent event data transformation.
#' Defaults to \code{"gaptime"}.
#' @param min_events Minimum number of events for each event number.
#' @examples
#' \dontrun{
#' data("cgd", package = "frailtyHL")
#' cgd2 <- filter(cgd, enum %in% c(1:2)) %>%
#'  select(id, age, tstart, tstop, status, enum)
#' cgd2 <- filter(cgd2, id %in% c(1:5))
#' # gaptime scale
#' ped_gt <- split_data_recurrent(
#'  formula = Surv(tstart, tstop, status) ~ age + enum,
#'  data = cgd2,
#'  episode_var = "enum", max_time = 388)
#' ped_ct <- split_data_recurrent(Surv(tstart, tstop, status)~enum + age, data = cgd2,
#'    episode_var = "enum", timescale = "calendar", max_time = 388)
#' }
#' @export
split_data_recurrent <- function(
  formula,
  data,
  episode_var    = character(),
  cut        = NULL,
  max_time   = NULL,
  event      = 1L,
  min_events = 1L,
  timescale = c("gap", "calendar"),
  ...) {

  assert_character(episode_var, min.chars = 1L, min.len = 1L, any.missing = FALSE,
    len = 1L)
  assert_integer(min_events, lower = 1L, len = 1L)
  assert_character(timescale)
  timescale <- match.arg(timescale)

  dots_in <- list(...)
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

  ## obtain interval breaks points for each spell
  if(timescale == "gap") {
    data <- mutate(data, !!!list(.time = quo(!!as.name(surv_vars[2]) - !!as.name(surv_vars[1]))))
    formula <- update(formula, Surv(.time, status)~.)
    formula <- update_formula(formula, proposed_names = c(".time", surv_vars[3]))
  }
  # split data for each spell
  data_list <- split(data, data[[episode_var]])
  rm(data)
  # only keep spells with minimum number of events
  data_list <- data_list[map_dbl(data_list, nrow) > min_events]
  cuts <- get_cut(data_list, formula, cut = cut, max_time = max_time,
    event = event, min_events = min_events)

  ## create argument list to be passed to survSplit
  dots <- list(...)

  # if id allready in the data set, remove id variable from dots but keep
  # id variable for later rearrangment
  if (!is.null(dots$id)) {
    id_var <- dots$id
  } else {
    id_var  <- "id"
    dots$id <- id_var
  }

  split_df_list <- map2(
    .x = data_list,
    .y = ifelse(is.list(cuts), cuts, list(cuts)),
    .f = ~ {
      dots$data    <- .x
      dots$formula <- formula
      dots$cut     <- .y
      split_df     <- do.call(split_data, dots)
    }
  )

  split_df <- bind_rows(split_df_list)

  ## set class and and attributes
  # class(split_df) <- c("ped", class(split_df))
  # attr(split_df, "breaks") <- cut
  # attr(split_df, "id_var") <- dots_in$id <- id_var
  # attr(split_df, "intvars") <- c(id_var, "tstart", "tend", "interval", "offset",
  #   "ped_status")
  # attr(split_df, "trafo_args") <- dots_in

  split_df

}


#' Recurrent events trafo
#'
#' @examples
#' \dontrun{
#' data("cgd", package = "frailtyHL")
#' cgd <- cgd %>%
#'  select(id, tstart, tstop, enum, status, age)
#' ped_re <- split_data(
#'   formula = Surv(tstart, tstop, status) ~ age + enum,
#'   data = cgd2)
#' }
#' @rdname as_ped
#' @export
#' @keywords internal
as_ped_recurrent <- function(
  data,
  formula,
  cut          = NULL,
  max_time     = NULL,
  tdc_specials = c("concurrent", "cumulative"),
  censor_code  = 0L,
  episode_var  = character(),
  timescale    = c("gap", "calendar"),
  min_events   = 1L,
  ...
) {

  assert_character(episode_var, min.chars = 1L, min.len = 1L, any.missing = FALSE,
    len = 1L)
  assert_integer(min_events, lower = 1L, len = 1L)

  status_error(data, formula)
  assert_subset(tdc_specials, c("concurrent", "cumulative"))


  dots             <- list(...)
  dots$data        <- data
  dots$formula     <- get_ped_form(formula, data = data, tdc_specials = tdc_specials)
  dots$cut         <- cut
  dots$max_time    <- max_time
  dots$episode_var <- episode_var
  dots$min_events  <- min_events
  dots$timescale   <- timescale

  ped <- do.call(split_data_recurrent, dots)
  attr(ped, "time_var")   <- get_lhs_vars(dots$formula)[1]
  attr(ped, "status_var") <- get_lhs_vars(dots$formula)[2]

  ped

}
