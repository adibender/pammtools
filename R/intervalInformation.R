# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2016-12-12 16:11:27
# @Last Modified by:   andreas.bender@stat.uni-muenchen.de
# @Last Modified time: 2017-01-19 21:03:23

#' Create start/end times and interval information
#'
#' Given interval breaks points, retuns data frame with information on 
#' interval start time, interval end time, interval length and a interval factor 
#' variable (left open intervals). If object of class ped is provided, extracts 
#' unique interval information from object. 
#' 
#' @param brks A numeric vector of cut points in which the follow-up should be
#' partitioned in or object of class \code{ped}. 
#' @param ... Currently ignored.
#' @rdname int_info
#' @return data.frame. A data frame containing the start and end times of the
#' intervals specified by the \code{brks} argument. Additionally the interval
#' length, interval mid-point and a factor variable of the intervals themselfs.
#' @export
int_info <- function(brks, ...) {
  UseMethod("int_info",  brks)
}


#' @inheritParams int_info
#' @param min.time Only intervals that have lower borders larger than
#' this value will be included in the resulting data frame.
#' @import checkmate dplyr
#' @rdname int_info
#' @export
int_info.numeric <- function(
  brks,
  min.time = 0L, ...) {

  # check inputs
  assert_numeric(brks, lower = 0, any.missing = FALSE)
  assert_numeric(min.time, lower  = 0L)

  # sort brks and add origin if necessary
  if(is.unsorted(brks)) {
    brks <- sort(brks)
  }
  if(min(brks!=0)) {
    brks <- c(0, brks)
  }

  intlen <- diff(brks)
  tstart <- brks[-length(brks)]
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
#' @seealso split_data
int_info.ped <- function(brks, ...) {

  brks %>% select(one_of(setdiff(
      attr(brks, "intvars"), 
      c("id", "offset", "time", "status")))) %>%
    unique()

}



#' Given breaks, return intervals in which times vector falls
#' 
#' @inheritParams int_info
#' @param x Vector of values for which interval information should be returned.
#' @param ... Further arguments passed to \code{\link[base]{findInterval}}.
#' @import dplyr
#' @return A \code{data.frame} containing information on intervals in which 
#' values of x fall
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

  int.df <- int_info(brks)
  int <- findInterval(x, union(int.df$tstart, int.df$tend), ...)

  int.df %>% 
    slice(int)      %>%
    mutate(x = x)   %>%
    arrange(tstart) %>%
    select(x, everything())

}


#' Extract inetrval information and median/modus values vor covariates
#' 
#' Given an object of class \code{ped}, returns data frame with interval information, 
#' median values for numerical variables and modi for non-numeric variables 
#' in the data set. 
#' 
#' @param ped An object of class \code{ped} as returned by \code{\link[pam]{split_data}}.
#' @import checkmate dplyr
#' @examples 
#' data("leuk2", package="bpcp")
#' leuk.ped <- split_data(Surv(time, status)~., data=leuk2, id="id")
#' ped_info(leuk.ped)
#' @export 
#' @seealso \code{\link[pam]{int_info}}, \code{\link[pam]{sample_info}}
ped_info <- function(ped) {

  assert_class(ped, classes="ped")

  int.df <- int_info(ped)
  sdf <- sample_info(ped)

  bind_cols(int.df, sdf[rep(1, nrow(int.df)), ])

}

#' Extract information for plotting step functions
#' 
#' 
#' @param pinf A data frame as returned by \code{\link[pam]{ped_info}} and 
#' potentially additional information from predictions, etc. 
#' @examples
#' data("leuk2", package="bpcp")
#' leuk.ped <- split_data(Surv(time, status)~., data=leuk2, id="id")
#' pem <- glm(status ~ interval, data = leuk.ped, family=poisson(), offset=offset)
#' pinf <- ped_info(leuk.ped)
#' pinf$basehaz <- predict(pem, newdata=pinf, type="response")
#' plot.inf <- plot_info(pinf)
#' 
#' 

plot_df <- function(pinf) {

  pinf <- bind_rows(pinf, pinf[nrow(pinf), ]) %>% 
    mutate(
      tend = lag(tend, default = min(tstart)), 
      intlen = lag(intlen, default = intlen[1])) %>% 
    select(-one_of("interval", "tstart")) %>% 
    rename(time=tend)

}

