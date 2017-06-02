#' Create start/end times and interval information
#'
#' Given interval breaks points, returns data frame with information on
#' interval start time, interval end time, interval length and a interval factor
#' variable (left open intervals). If object of class ped is provided, extracts
#' unique interval information from object.
#'
#' @param x A numeric vector of cut points in which the follow-up should be
#' partitioned in or object of class \code{ped}.
#' @param ... Currently ignored.
#' @rdname int_info
#' @return data.frame. A data frame containing the start and end times of the
#' intervals specified by the \code{x} argument. Additionally the interval
#' length, interval mid-point and a factor variable of the intervals themselves.
#' @export
int_info <- function(x, ...) {
  UseMethod("int_info",  x)
}


#' @inheritParams int_info
#' @param min.time Only intervals that have lower borders larger than
#' this value will be included in the resulting data frame.
#' @import checkmate dplyr
#' @export
#' @examples
#' int_info(c(1, 2.3, 5))
#' @rdname int_info
int_info.default <- function(
  x,
  min.time = 0L, ...) {

  # check inputs
  assert_numeric(x, lower = 0, any.missing = FALSE)
  assert_numeric(min.time, lower  = 0L)

  # sort x and add origin if necessary
  if(is.unsorted(x)) {
    x <- sort(x)
  }
  if(min(x!=0)) {
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
      intmid = tstart + intlen/2,
      interval = paste0("(", tstart, ",", tend, "]"),
      interval = factor(interval, levels=interval))

  filter(tdf, tstart >= min.time)

}

#' @inheritParams int_info
#' @import dplyr
#' @rdname int_info
#' @examples
#' tdf <- data.frame(time=c(1, 2.3, 5), status=c(0, 1, 0))
#' ped <- split_data(Surv(time, status)~., data=tdf, id="id", max.end=TRUE)
#' int_info(ped)
#' @export
#' @seealso split_data
int_info.ped <- function(x, ...) {

  int_info(attr(x, "cut"), ...)

}


#' Given breaks, return intervals in which times vector falls
#'
#' @inheritParams int_info
#' @param brks Vector of values for which interval information should be returned.
#' @param ... Further arguments passed to \code{\link[base]{findInterval}}.
#' @import dplyr
#' @return A \code{data.frame} containing information on intervals in which
#' values of x fall.
#' @examples
#' set.seed(111018)
#' brks <- c(0, 4.5, 5, 10, 30)
#' int_info(brks)
#' x <- runif(3, 0, 30)
#' get_intervals(brks, x, left.open=TRUE)
#' @export
#' @seealso findInterval int_info

get_intervals <- function(brks, x, ...) {

  # check inputs
  assert_numeric(brks, lower = 0, any.missing = FALSE)
  assert_numeric(x, finite = TRUE, all.missing = FALSE)

  int_df <- int_info(brks)
  int    <- findInterval(x, union(int_df$tstart, int_df$tend), ...)

  int_df %>%
    slice(int)      %>%
    mutate(x = x)   %>%
    arrange(tstart) %>%
    select(x, everything())

}


#' Extract interval information and median/modus values for covariates
#'
#' Given an object of class \code{ped}, returns data frame with one row for each
#' interval containing interval information, median values for numerical
#' variables and modi for non-numeric variables in the data set.
#'
#' @param ped An object of class \code{ped} as returned by \code{\link[pam]{split_data}}.
#' @import checkmate dplyr
#' @examples
#' data("veteran", package="survival")
#' ped <- split_data(Surv(time, status)~ trt + age, data=veteran, id="id")
#' ped_info(ped) # note that trt is coded 1/2, should be fixed beforehand
#' @export
#' @return A data frame with one row for each interval in \code{ped}.
#' @seealso \code{\link[pam]{int_info}}, \code{\link[pam]{sample_info}}
ped_info <- function(ped) {

  assert_class(ped, classes="ped")

  int_df <- int_info(ped)
  sdf    <- sample_info(ped)
  bind_cols(
    int_df %>% slice(rep(seq_len(nrow(int_df)), times = nrow(sdf))),
    sdf %>% slice(rep(seq_len(nrow(sdf)), each = nrow(int_df)))) %>%
    grouped_df(vars = groups(sdf))
}

#' Extract risk set information for each interval.
#'
#' The columns \code{ped_riskset, ped_events, ped_censored} provide the
#' size of the riskset at the beginning of each interval as well as the number
#' of events and censorings that occured in the interval, respectively.
#' @param ped An object of class \code{ped} as returned by \code{\link[pam]{split_data}}.
#' @import checkmate dplyr
#' @examples
#' data("veteran", package="survival")
#' ped <- split_data(Surv(time, status)~ ., data = veteran, id = "id",
#'   cut = seq(0,400, by = 100))
#' riskset_info(ped)
#' (riskset_celltype <- riskset_info(group_by(ped, celltype)))
#' ## add descriptive statistics for riskset at beginning of each interval:
#' # left_join(riskset_celltype,
#' #           group_by(ped, celltype, interval) %>% sample_info())
#' @export
#' @return A data frame with one row for each interval in \code{ped}.
#' @seealso \code{\link[pam]{int_info}}, \code{\link[pam]{sample_info}}
riskset_info <- function(ped) {
  assert_class(ped, classes="ped")

  # how often is interval the last row for a given id when status == 0?
  censored <- ped %>% group_by(id, add = TRUE) %>%
    arrange(tend) %>% slice(n()) %>%
    filter(ped_status == 0) %>% ungroup(id) %>%
    grouped_df(vars = c(groups(ped), as.symbol("interval"))) %>%
    summarize(ped_censored = n())

  join_vars <- unlist(c(sapply(groups(ped), deparse), "interval"))
  ped %>% group_by(interval, add = TRUE) %>%
    summarize(
      ped_riskset = n(),
      ped_events = sum(ped_status)) %>%
    left_join(censored, by = join_vars) %>%
    mutate(ped_censored = ifelse(is.na(ped_censored), 0, ped_censored))
}

#' Extract information for plotting step functions
#'
#'
#' @param pinfo A data frame as returned by \code{\link[pam]{ped_info}} and
#' potentially additional information from predictions, etc.
#' @examples
#' data("veteran", package="survival")
#' leuk.ped <- split_data(Surv(time, status)~., data=veteran, id="id")
#' pem <- glm(ped_status ~ interval, data = leuk.ped, family=poisson(), offset=offset)
#' pinfo <- ped_info(leuk.ped)
#' pinfo$basehaz <- predict(pem, newdata=pinfo, type="response")
#' plot.inf <- plot_df(pinfo)
#' @import dplyr
#' @export
plot_df <- function(pinfo) {

  bind_rows(pinfo, pinfo[nrow(pinfo), ]) %>%
    mutate(
      tend   = lag(tend, default   = min(tstart)),
      intlen = lag(intlen, default = intlen[1])) %>%
    select(-one_of("interval", "tstart")) %>%
    rename(time=tend)

}
