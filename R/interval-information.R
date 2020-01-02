#' Create start/end times and interval information
#'
#' Given interval breaks points, returns data frame with information on
#' interval start time, interval end time, interval length and a factor
#' variable indicating the interval (left open intervals). If an object of class
#' \code{ped} is provided, extracts unique interval information from object.
#'
#' @param x A numeric vector of cut points in which the follow-up should be
#' partitioned in or object of class \code{ped}.
#' @param ... Currently ignored.
#' @rdname int_info
#' @return A data frame containing the start and end times of the
#' intervals specified by the \code{x} argument. Additionally, the interval
#' length, interval mid-point and a factor variable indicating the intervals.
#' @export
int_info <- function(x, ...) {
  UseMethod("int_info",  x)
}


#' @param min_time Only intervals that have lower borders larger than
#' this value will be included in the resulting data frame.
#' @import checkmate dplyr
#' @examples
#' ## create interval information from cut points
#' int_info(c(1, 2.3, 5))
#'
#' @rdname int_info
#' @export
int_info.default <- function(
  x,
  min_time = 0L, ...) {

  # check inputs
  assert_numeric(x, lower = 0, any.missing = FALSE)
  assert_numeric(min_time, lower  = 0L)

  # sort x and add origin if necessary
  if (is.unsorted(x)) {
    x <- sort(x)
  }
  if (min(x != 0)) {
    x <- c(0, x)
  }

  intlen <- diff(x)
  tstart <- x[-length(x)]
  tend   <- tstart + intlen

  tdf <- data.frame(
    tstart = tstart,
    tend   = tend,
    intlen = intlen) %>%
    mutate(
      intmid = tstart + intlen / 2,
      interval = paste0("(", tstart, ",", tend, "]"),
      interval = factor(.data$interval, levels = .data$interval))

  filter(tdf, tstart >= min_time)

}

#' @import dplyr
#' @rdname int_info
#' @examples
#' ## extract interval information used to create ped object
#' tdf <- data.frame(time=c(1, 2.3, 5), status=c(0, 1, 0))
#' ped <- tdf %>% as_ped(Surv(time, status)~., id="id")
#' int_info(ped)
#'
#' @seealso as_ped ped_info
#' @export
int_info.ped <- function(x, ...) {

  int_info(attr(x, "breaks"), ...)

}

#' @rdname int_info
#' @export
#' @keywords internal
int_info.pamm <- function(x, ...) {

  int_info(x[["trafo_args"]][["cut"]])

}


#' Information on intervals in which times fall
#'
#' @inheritParams int_info
#' @param x An object from which interval information can be obtained,
#' see \code{\link{int_info}}.
#' @param times A vector of times for which corresponding interval information
#' should be returned.
#' @param ... Further arguments passed to \code{\link[base]{findInterval}}.
#' @import dplyr
#' @return A \code{data.frame} containing information on intervals in which
#' values of \code{times} fall.
#' @seealso \code{\link[base]{findInterval}} \code{\link{int_info}}
#' @rdname get_intervals
#' @export
#' @examples
#' set.seed(111018)
#' brks <- c(0, 4.5, 5, 10, 30)
#' int_info(brks)
#' x <- runif (3, 0, 30)
#' x
#' get_intervals(brks, x)
get_intervals <- function(x, times, ...) {
  UseMethod("get_intervals", x)
}

#' @rdname get_intervals
#' @inheritParams base::findInterval
#' @export
get_intervals.default <- function(
  x,
  times,
  left.open        = TRUE,
  rightmost.closed = TRUE,
  ...) {

  # check inputs
  assert_numeric(times, lower = 0, finite = TRUE, all.missing = FALSE)

  int_df <- int_info(x)
  int    <- findInterval(
    x                = times,
    vec              = sort(union(int_df$tstart, int_df$tend)),
    left.open        = left.open,
    rightmost.closed = rightmost.closed)

  int_df %>%
    slice(int) %>%
    mutate(times = times) %>%
    # arrange(times) %>%
    select(times, everything())

}


#' Extract interval information and median/modus values for covariates
#'
#' Given an object of class \code{ped}, returns data frame with one row for each
#' interval containing interval information, mean values for numerical
#' variables and modus for non-numeric variables in the data set.
#'
#' @param ped An object of class \code{ped} as returned by
#' \code{\link[pammtools]{as_ped}}.
#' @import checkmate dplyr
#' @examples
#' ped <- tumor[1:4,] %>% as_ped(Surv(days, status)~ sex + age)
#' ped_info(ped)
#' @export
#' @return A data frame with one row for each unique interval in \code{ped}.
#' @seealso \code{\link[pammtools]{int_info}}, \code{\link[pammtools]{sample_info}}
ped_info <- function(ped) {

  assert_class(ped, classes = "ped")

  int_df <- int_info(ped)
  sdf    <- sample_info(ped)
  if (is.null(sdf)) {
    return(int_df)
  } else {
    bind_cols(
      int_df %>% slice(rep(seq_len(nrow(int_df)), times = nrow(sdf))),
      sdf %>% slice(rep(seq_len(nrow(sdf)), each = nrow(int_df)))) %>%
      grouped_df(vars = group_vars(sdf))
  }

}

#' Extract risk set information for each interval.
#'
#' The columns \code{ped_riskset, ped_events, ped_censored} provide the
#' size of the risk set at the beginning of each interval as well as the number
#' of events and censorings that occurred in the interval, respectively.
#'
#' @inheritParams ped_info
#' @import checkmate dplyr
#' @examples
#' ped <- tumor[1:4,] %>% as_ped(Surv(days, status)~ .)
#' riskset_info(ped)
#' @keywords internal
#' @export
#' @return A data frame with one row for each interval in \code{ped}.
#' @seealso \code{\link[pammtools]{int_info}}, \code{\link[pammtools]{sample_info}}
riskset_info <- function(ped) {
  assert_class(ped, classes = "ped")

  # how often is interval the last row for a given id when status == 0?
  censored <- ped %>% group_by(id, add = TRUE) %>%
    arrange(.data$tend) %>% slice(n()) %>%
    filter(.data$ped_status == 0) %>% ungroup(id) %>%
    grouped_df(vars = c(group_vars(ped), "interval")) %>%
    summarize(ped_censored = n())

  join_vars <- c(group_vars(ped), "interval")
  ped %>% group_by(.data$interval, add = TRUE) %>%
    summarize(
      ped_riskset = n(),
      ped_events = sum(.data$ped_status)) %>%
    left_join(censored, by = join_vars) %>%
    mutate(ped_censored = ifelse(is.na(.data$ped_censored),
        0, .data$ped_censored))
}
