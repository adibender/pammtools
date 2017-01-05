#' Given breaks, create start/end times and interval information
#'
#' @param brks numeric. A vector of cut point in which the follow-up should be
#' partitioned in.
#' @param min.time numeric. Only intervals that have lower boarders larger than
#' this value will be included in the resulting data frame.
#' @return data.frame. A data frame containing the start and end times of the
#' intervals specified by the \code{brks} argument. Additionally the interval
#' length, interval mid-point and a factor variable of the intervals themselfs.
#' @import checkmate dplyr
#' @export

int_info <- function(
  brks,
  min.time = 0L) {

  # check inputs
  assert_numeric(brks, lower = 0, any.missing = FALSE)
  assert_numeric(min.time, lower  = 0L)

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
