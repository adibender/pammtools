t_lead <- function(te, te.add.lead=0, lead.const=4, lead.factor=2) {
  lead.const + (te + te.add.lead)*lead.factor
}

get_end <- function(start, lead, max.end) {

  end <- start + lead
  end <- ifelse(end > max.end, max.end, end)
}

find_interval <- function(
  x, 
  cut = c(0L:12L, seq(15, 55, by = 5), 61),
  labels  = FALSE, ...) {

  
  ind <- cut(x, cut, labels=FALSE, ...)
  # if(any(is.na(ind))) ind[is.na(ind)] <- max(ind, na.rm=TRUE)
  if(labels) {
    ind <- int_info2(min.time=0, x=cut)[ind, "interval"]   
  }

  ind

}


##' Given interval break-points, provides data frame with information regarding
#' interval start-/end-points, lengths, midpoints,etc.
#'
#' @inheritParams int_info
#' @importFrom dplyr mutate
int_info2 <- function(
  x    = c(0:12, seq(15, 55, by=5), 61),
  min.time = 4) {

  intlen <- diff(x)
  tstart <- c(0, cumsum(intlen)[-length(intlen)])
  tend   <- tstart + intlen

  tdf <- data.frame(
    tstart = tstart,
    tend   = tend,
    intlen = intlen)
  tdf          <- mutate(tdf, intmid = tstart + intlen/2)
  tdf$interval <- paste0("(", tdf$tstart, ",", tdf$tend, "]")
  tdf$interval <- factor(tdf$interval, levels=tdf$interval, labels=tdf$interval)

  ind.keep <- which(tstart >= min.time)
  subset(tdf, tstart >= min.time)

}

lag_lead_df <- function(
  te,
  cut,
  te.add.lag   = 0,
  t.lag        = 1,
  te.add.lead  = 0,
  lead.const   = 0,
  lead.factor  = 0,
  labels       = TRUE, 
  ...) {

  if(is.null(lead.const)) {
    lead.const <- max(cut)
  }
  lead    = t_lead(te, te.add.lead = te.add.lead, lead.const = lead.const, 
    lead.factor = lead.factor)
  w.begin = (te + te.add.lag) + t.lag
  w.end   = get_end(w.begin, lead, max.end=max(cut))

  data.frame(
    te        = te,
    lag       = t.lag,
    lead      = lead,
    w.begin   = w.begin,
    w.end     = w.end,
    int.begin = find_interval(w.begin, cut=cut, labels = labels, ...),
    int.end   = find_interval(w.end, cut=cut, labels = labels, ...)
    )
}
#' creates one instance of Lag/Lead mat
#' @param te Numeric/Integer vector specifying the times at which exposure occurred.
#' @param te.add.lag A numeric constant added to te before application of lag time
#' @param t.lag A numeric constant, specifying the time (from \code{te}) before 
#' \code{te} can affect hazard.
#' @param te.add.lead A numeric constant, added to te before application of lead time.
#' @param lead.const A numeric constant, specifying the constant amount of time 
#' after \code{te + t.lag}, in which \code{te} can still affect hazard.
#' @param lead.factor If the lead time is dynamic, this factor can be set different 
#' to zero, such that \code{t.lead=lead.const + lead.factor*te}.
#' @param cut The break points dividing the follow up into intervals.
#' @param t.min If some intervals are not of interest only intervals for t > t.min are 
#' returned.
#' @return A data frame with intervals as first column and \code{length(te)} 
#' columns specifying the lag/lead for each \code{te}.
#' @import checkmate dplyr 
create_Lmat <- function(
  te,
  cut,
  t.lag        = 1,
  lead.const   = NULL,
  te.add.lag   = 0,
  te.add.lead  = 0,
  lead.factor  = 0,
  t.min        = 0) {
  
  assert_integer(te,         lower = 1, any.missing = FALSE, unique    = TRUE)
  assert_numeric(cut,        lower = 0, any.missing = FALSE, min.len   = 2)
  assert_number(te.add.lag,  lower = 0, upper = max(cut), finite = TRUE)
  assert_number(t.lag,       lower = 0, upper = max(cut), finite = TRUE)
  assert_number(lead.const,  lower = 0, upper = max(cut), finite = TRUE, null.ok = TRUE)
  assert_number(lead.factor, lower = 0, upper = max(cut), finite = TRUE)
  assert_number(t.min,       lower = 0, upper = max(cut), finite = TRUE)

  # create lag-lead information matrix
  ldf <- lag_lead_df(te=te, te.add.lag=te.add.lag, te.add.lead=te.add.lead, 
    t.lag=t.lag, lead.const=lead.const, lead.factor=lead.factor, 
    cut=cut, left.open=TRUE, rightmost.closed=TRUE)

  ind.begin <-find_interval(ldf$w.begin, cut=cut, left.open=TRUE, 
    rightmost.closed = TRUE)
  ind.end   <-find_interval(ldf$w.end, cut=cut, left.open=TRUE, 
    rightmost.closed = TRUE)

  int.info  <- int_info2(x=cut, min.time=0)
  int.keep  <- int.info$interval[which(int.info$tstart >= t.min)]

  ints <- apply(cbind(ind.begin, ind.end), 1, function(z) {
    z.i <- int.info$interval[z[1]:z[2]]
    int.info$interval %in% z.i
  }) * 1

  ints <- data.frame(intsL=int.info$interval, Lcols=ints)
  filter(ints, intsL %in% int.keep)


}

#' Create one instance of the Lag-Lead matrix
#' @inheritParams as_fped
#' @param te_vec The vector of exposure times (times at which Change in TDC occured)
#' @param t_vec The vector of event/split-times on the scale of the follow-up.
Lsimp <- function(te_vec, t_vec, ll_fun=function(te, t) { te <= t}) {
  t(outer(te_vec, t_vec, ll_fun))*1
}

#' Extend instance of Lag-Lead matrix to whole data set 
#' @param LL An instance of the Lag-Lead matrix (for all time points of the follow-up)
#' @param id The id vector of the data set in PED format. 
Lmat <- function(LL, id) {

  rle_id <- rle(id)
  LL_ind <- unlist(sapply(rle_id$lengths, seq_len))

  LL[LL_ind, ]


}