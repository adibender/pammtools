#' Transform data to Piece-wise Exponential Data (PED)
#'
#' This is the general data transformation function provided by the
#' \code{pammtools} package. Two main applications must be distinguished:
#' \enumerate{
#'  \item Transformation of standard time-to-event data.
#'  \item Transformation of left-truncated time-to-event data.
#'  \item Transformation of time-to-event data with time-dependent covariates (TDC).
#' }
#' For the latter, the type of effect one wants to estimate is also
#' important for the data transformation step.
#' In any case, the data transformation is specified by a two sided formula.
#' In case of TDCs, the right-hand-side of the formula can contain formula specials
#' \code{concurrent} and \code{cumulative}.
#' See the \href{https://adibender.github.io/pammtools//articles/data-transformation.html}{data-transformation}
#' vignette for details.
#'
#'
#' @rdname as_ped
#' @param data Either an object inheriting from data frame or in case of
#' time-dependent covariates a list of data frames (of length 2), where the first data frame
#' contains the time-to-event information and static covariates while the second
#' (and potentially further data frames) contain information on time-dependent
#' covariates and the times at which they have been observed.
#' @param formula A two sided formula with a \code{\link[survival]{Surv}} object
#' on the left-hand-side and covariate specification on the right-hand-side (RHS).
#' The RHS can be an extended formula, which specifies how TDCs should be transformed
#' using specials \code{concurrent} and \code{cumulative}. The left hand-side can
#' be in start-stop-notation. This, however, is only used to create left-truncated
#' data and does not support the full functionality.
#' @param cut Split points, used to partition the follow up into intervals.
#' If unspecified, all unique event times will be used.
#' @param max_time If \code{cut} is unspecified, this will be the last
#' possible event time. All event times after \code{max_time}
#' will be administratively censored at \code{max_time}.
#' @param tdc_specials A character vector. Names of potential specials in
#' \code{formula} for concurrent and or cumulative effects.
#' @param censor_code Specifies the value of the status variable that indicates censoring.
#' Often this will be \code{0}, which is the default.
#' @param ... Further arguments passed to the \code{data.frame} method and
#' eventually to \code{\link[survival]{survSplit}}
#' @importFrom Formula Formula
#' @examples
#' tumor[1:3, ]
#' tumor[1:3, ] %>% as_ped(Surv(days, status)~ age + sex, cut = c(0, 500, 1000))
#' tumor[1:3, ] %>% as_ped(Surv(days, status)~ age + sex)
#' @return A data frame class \code{ped} in piece-wise exponential data format.
#' @export
as_ped <- function(data, ...) {
  UseMethod("as_ped", data)
}

#' @rdname as_ped
#' @export
as_ped.data.frame <- function(
  data,
  formula,
  cut          = NULL,
  max_time     = NULL,
  tdc_specials = c("concurrent", "cumulative"),
  censor_code  = 0L,
  transition   = character(),
  timescale    = c("gap", "calendar"),
  min_events   = 1L,
  ...) {

  status_error(data, formula, censor_code)
  assert_subset(tdc_specials, c("concurrent", "cumulative"))

  if (test_character(transition, min.chars = 1L, min.len = 1L)) {
    ped <- as_ped_multistate(data = data, formula = formula, cut = cut,
      max_time = max_time, tdc_specials = tdc_specials, censor_code = censor_code,
      transition = transition, timescale = timescale, min_events = min_events, ... )
    return(ped)
  }

  event_types <- get_event_types(data, formula, censor_code)
  if (length(event_types) > 1) {

    ped <- as_ped_cr(data = data, formula = formula, cut = cut, max_time = max_time,
      tdc_specials = tdc_specials, censor_code = censor_code, ...)

  } else {

    dots          <- list(...)
    dots$data     <- data
    dots$formula  <- get_ped_form(formula, data = data, tdc_specials = tdc_specials)
    dots$cut      <- cut
    dots$max_time <- max_time

    ped <- do.call(split_data, dots)
    attr(ped, "time_var") <- get_lhs_vars(dots$formula)[1]
    attr(ped, "status_var") <- get_lhs_vars(dots$formula)[2]

  }

  ped

}

#' @rdname as_ped
#' @export
as_ped.nested_fdf <- function(
  data,
  formula,
  ...) {

  dots <- list(...)
  # update interval break points (if necessary)
  cut <- dots$cut
  if (is.null(cut)) {
    cut <- attr(data, "breaks")
  }
  ccr_breaks <- attr(data, "ccr_breaks")
  cut <- union(cut, ccr_breaks[ccr_breaks <= max(cut)]) %>% sort()

  # ped <- data %>%
  #   select_if(is.atomic) %>%
  #   as.data.frame() %>%
  #   as_ped(
  #     formula  = formula,
  #     id       = dots$id,
  #     cut      = cut,
  #     max_time = dots$max_time,
  #     ...)
  dots$formula <- formula
  dots$data    <- as.data.frame(select_if(data, is.atomic))
  dots$cut     <- cut
  ped          <- do.call(as_ped, dots)

  # replace updated attributes
  attr(data, "breaks") <- attr(ped, "breaks")
  attr(data, "id_n") <- ped %>% group_by(!!sym(attr(data, "id_var"))) %>%
    summarize(id_n = n()) %>% pull("id_n") %>% as_vector()
  attr(data, "id_tseq") <- ped %>% group_by(!!sym(attr(data, "id_var"))) %>%
    transmute(id_tseq = row_number()) %>% pull("id_tseq") %>% as_vector()
  attr(data, "id_tz_seq") <- rep(seq_len(nrow(data)),
    times = attr(data, "id_n"))

  if (has_special(formula, "concurrent")) {
    ped <- ped %>% add_concurrent(data = data, id_var = dots$id)
  }

  if (has_special(formula, "cumulative")) {
    ped <- add_cumulative(ped, data = data, formula = formula)
    attr(ped, "ll_weights") <- imap(attr(ped, "tz"),
      ~bind_cols(!!.y := .x, ll_weight = c(mean(abs(diff(.x))), abs(diff(.x)))))
    class(ped) <- c("fped", class(ped))
  }
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  attr(ped, "func_mat_names") <- make_mat_names(
  attr(ped, "func"),
  attr(ped, "time_var"))

  ped

}

#' @rdname as_ped
#' @export
as_ped.list <- function(
  data,
  formula,
  tdc_specials = c("concurrent", "cumulative"),
  censor_code = 0L,
  ...) {

  assert_class(data, "list")
  assert_class(formula, "formula")

  status_error(data[[1]], formula, censor_code)

  nl    <- length(data)
  # form  <- Formula(formula)
  has_tdc <- has_tdc_form(formula, tdc_specials = tdc_specials)

  if (nl == 1 & !has_tdc) {
    ped <- data[[1]] %>% as_ped(formula = formula, tdc_specials = tdc_specials, ...)
  } else {
    if (nl == 2 & !has_tdc) {
    stop("Two data sets provided in 'data' but no specification of
      time-dependent covariate effects in 'formula'")
    } else {

      nested_fdf <- nest_tdc(data, formula, ...)
      ped <- as_ped(nested_fdf, formula, ...)

    }
  }
  lhs_vars <- get_lhs_vars(formula)
  attr(ped, "time_var") <- lhs_vars[1]
  attr(ped, "trafo_args")$formula <- formula

  ped

}

#' @rdname as_ped
#' @param x any R object.
#' @export
is.ped <- function(x) inherits(x, "ped")


#' @rdname as_ped
#' @param newdata A new data set (\code{data.frame}) that contains the same
#' variables that were used to create the PED object (\code{data}).
#' @export
as_ped.ped <- function(data, newdata, ...) {

  if (is.ped(newdata)) {
    stop("newdata already in ped format.")
  }

  trafo_args <- attr(data, "trafo_args")
  trafo_args[["data"]] <- newdata
  do.call(as_ped,  trafo_args)

}



#' @rdname as_ped
#' @export
as_ped.pamm <- function(data, newdata, ...) {

  if (is.ped(newdata)) {
    stop("newdata already in ped format.")
  }
  trafo_args      <- data[["trafo_args"]]
  trafo_args$data <- newdata
  do.call(split_data, trafo_args)

}

## Competing risks

#' Competing risks trafo
#'
#' @inherit as_ped
#' @importFrom rlang .env
#'
#' @keywords internal
as_ped_cr <- function(
  data,
  formula,
  cut          = NULL,
  max_time     = NULL,
  tdc_specials = c("concurrent", "cumulative"),
  censor_code  = 0L,
  combine      = TRUE,
  ...) {

  lhs_vars <- get_lhs_vars(formula)
  n_lhs <- length(lhs_vars)
  event_types <- get_event_types(data, formula, censor_code)
  n_events <- sum(event_types != censor_code)

  cut <- map2(
    event_types,
    if(is.list(cut)) cut else list(cut),
    function(.event, .cut) {
      get_cut(data, formula = formula, cut = .cut, max_time = NULL, event = .event)
    }
  )
  if(length(cut) > 1 & combine) {
    cut <- list(reduce(cut, union))
  }

  ped <- map2(
    event_types,
    cut,
    function(.event, .cut) {
      ped_i <- data %>%
        mutate(!!lhs_vars[n_lhs] := 1L * (.data[[lhs_vars[n_lhs]]] == .env[[".event"]])) %>%
        as_ped(
          formula      = formula,
          cut          = .cut,
          max_time     = max_time,
          tdc_specials = tdc_specials,
          ...)
      ped_i$cause <- .event
      ped_i
    })

  if (combine) {
    ped <- do.call(rbind, ped)
    class(ped) <- c("ped_cr_union", "ped_cr", class(ped))
    attr(ped, "intvars") <- c(attr(ped, "intvars"), "cause")
    attr(ped, "breaks") <- if (length(cut) ==1) unlist(cut) else cut
  } else {
    class(ped) <- c("ped_cr_list", "ped_cr", "ped", class(ped))
    names(ped) <- paste0("cause = ", event_types)
    attributes(ped)$trafo_args$id <- attributes(ped[[1]])$trafo_args$id
    attributes(ped)$trafo_args$formula <- formula
  }

  attr(ped, "trafo_args")[["cut"]] <- if (length(cut) ==1) unlist(cut) else cut
  attr(ped, "trafo_args")[["combine"]] <- combine
  attr(ped, "trafo_args")[["censor_code"]] <- censor_code
  attr(ped, "risks") <- event_types

  ped

}

#' Exctract event types
#'
#' Given a formula that specifies the status variable of the outcome, this function
#' extracts the different event types (except for censoring, specified by
#' \code{censor_code}).
#'
#' @inheritParams as_ped
#'
#' @keywords internal
get_event_types <- function(data, formula, censor_code) {

  lhs_vars <- get_lhs_vars(formula)
  status_values <- unique(data[[lhs_vars[length(lhs_vars)]]]) %>% sort()
  status_values[status_values != censor_code]

}



#' Recurrent events trafo
#'
#' @examples
#' \dontrun{
#' data("cgd", package = "frailtyHL")
#' cgd2 <- cgd %>%
#'  select(id, tstart, tstop, enum, status, age) %>%
#'  filter(enum %in% c(1:2))
#' ped_re <- as_ped_multistate(
#'   formula = Surv(tstart, tstop, status) ~ age + enum,
#'   data = cgd2,
#'  transition = "enum",
#'  timescale = "calendar")
#' }
#' @rdname as_ped
#' @export
#' @keywords internal
as_ped_multistate <- function(
  data,
  formula,
  cut          = NULL,
  max_time     = NULL,
  tdc_specials = c("concurrent", "cumulative"),
  censor_code  = 0L,
  transition  = character(),
  timescale    = c("gap", "calendar"),
  min_events   = 1L,
  ...
) {

  assert_character(transition, min.chars = 1L, min.len = 1L, any.missing = FALSE,
    len = 1L)
  assert_integer(min_events, lower = 1L, len = 1L)

  status_error(data, formula, censor_code)
  assert_subset(tdc_specials, c("concurrent", "cumulative"))

  rhs_vars <- get_rhs_vars(formula)
  if (!(transition %in% rhs_vars)) {
    formula <- add_to_rhs(formula, transition)
  }

  dots            <- list(...)
  dots$data       <- data
  dots$formula    <- get_ped_form(formula, data = data, tdc_specials = tdc_specials)
  dots$cut        <- sort(unique(cut))
  dots$max_time   <- max_time
  dots$transition <- transition
  dots$min_events <- min_events
  dots$timescale  <- timescale
  
  ped <- do.call(split_data_multistate, dots)
  attr(ped, "time_var")   <- get_lhs_vars(dots$formula)[1]

  return(ped)

}
