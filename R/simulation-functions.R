#' Construct full data set with time-dependent covariate (TDC)
#'
#' @param n Number of subjects to simulate.
#' @param m the number of follow up intervals (\code{t_end = 1:m})
#' @param te the vector of timepoints at which the TDC was observed
#' @param h0 A constant component of the baseline hazard.
#' @param f0 A function of time \code{t}, representing the smooth part of the
#'   baseline hazard.
#' @param fz A function of time \code{t}, exposure time \code{te} and the TDC
#'   \code{z}, representing the partial effect of \code{z} at \code{te} on the
#'   hazard at \code{t}.
#' @param fwindow  A function of time \code{t} and exposure time \code{te}
#'   defining the "window of effectiveness" during which exposures affect the
#'   hazard at time t. Should return \code{TRUE/FALSE}.
#' @param rng_z RNG for TDC $z(t_e)$, a function of \code{te}. Must return a
#'   vector of same length as te.
#' @import checkmate
#' @importFrom stats runif arima.sim
#' @keywords internal
#' @export
make_X <- function(
	n   = 500L,
	m   = 30,
  te  = -29:30,
	h0  = -2,
	f0  = function(t) {0 * t},
	fz = function(t, te, z) {(-(te - min(te))/10 + 0.005 * (te - min(te))^2) * z},
  fwindow = function(t, te) {(te <= t) & (te >= t - 30)},
  rng_z = function(te) {arima.sim(n = length(te), list(ar = c(.8, -.6)))}) {

	## check inputs
	assert_int(n, lower = 2)
	assert_int(m, lower = 1)

	t <- seq_len(m)

	assert_numeric(t, min.len = 1, any.missing = FALSE)
	assert_true(all(diff(t) == 1))
	assert_numeric(te, min.len = 1, any.missing = FALSE)
	assert_true(all(diff(te) == 1))
	assert_number(h0, lower=-5, upper=5)
	assert_function(f0)
	assert_function(fz)

	me <- length(te)


	# create time-dependent covariate
	Z <- t(replicate(n, rng_z(te)))

	## create data frame (start stop format)
	id <- rep(seq_len(n), each = m)
	Xdf <- data.frame(id = id)
	Xdf$tstart <- rep(t - 1, times = n)
	Xdf$tend   <- rep(t, times = n)
	Xdf$intmid <- Xdf$tstart + (Xdf$tend -Xdf$tstart)/2
	Xdf$Z      <- Z[id, ]
	Xdf$Fz <- NA * Xdf$Z
	Xdf$LL <- NA * Xdf$Z
  for (i in seq_len(nrow(Xdf))) {
    Xdf$Fz[i,] <- fz(t = Xdf$intmid[i], te = te, z = Xdf$Z[i,])
    Xdf$LL[i,] <- c(mean(diff(te)), diff(te)) * fwindow(t = Xdf$intmid[i],  te = te)
  }
	colnames(Xdf$Fz) <- paste0("V", seq_len(me))
	# Xdf$z_vec  <- as.vector(t(Z))

	## add baseline + cumulative effects
	Xdf$eta_base <- h0 + f0(Xdf$intmid)
	Xdf$eta_wce  <- rowSums(Xdf$Fz * Xdf$LL)
	Xdf$eta      <- Xdf$eta_base + Xdf$eta_wce

	Xdf$te_df   <- matrix(te, nrow = nrow(Xdf), ncol = me, byrow = TRUE)
	Xdf$time_df <- matrix(Xdf$tend, nrow = nrow(Xdf), ncol = me)

	Xdf
}


#' Simulate survival data using the piece-wise exponential distribution
#'
#' @rdname sim_wce
#' @param Xdf A data frame containing full (time-dependent) covariate information
#' for each subject
#' @param tmax The maximal time of follow-up. Survival times after \code{tmax}
#' will be administratively censored.
#' @import checkmate
#' @importFrom msm rpexp
#' @keywords internal
#' @export
sim_wce <- function(Xdf, tmax=30) {

	assert_data_frame(Xdf, all.missing=FALSE, min.rows=2, min.cols=3)
	assert_number(tmax, lower=0, finite=TRUE)

	tdf  <- data.frame(rate = exp(Xdf$eta), t = Xdf$tstart)
	rt.l <- split(tdf, f=Xdf$id)
	rt.l <- lapply(rt.l, as.list)

	## simulate new survival
	new.times                   <- sapply(rt.l, do.call, what = rpexp)
	status                      <- (new.times <= tmax)*1
	new.times[new.times > tmax] <- tmax # administrative censoring at t = tmax
	Xdf$time                    <- new.times[Xdf$id]
	Xdf$status                  <- 0

	Xdf_sub      <- subset(Xdf, tstart < time)
	Xdf_sub$time <- round(Xdf_sub$time, 2) + 0.01# round to have less splits,
	# +0.01 to avoid 0.00 times
	sub_ind                 <- cumsum(rle(Xdf_sub$id)$lengths)
	Xdf_sub$status[sub_ind] <- status
	Xdf_sub$offset          <- 0
	Xdf_sub$offset[sub_ind] <- log(new.times - Xdf_sub$tstart[sub_ind])

	Xdf_sub

}

#' @rdname sim_wce
#' @inherit sim_wce
#' @importFrom dplyr mutate left_join
#' @keywords internal
#' @export
sim_wce2 <- function(Xdf, tmax=30) {

	assert_data_frame(Xdf, all.missing=FALSE, min.rows=2, min.cols=3)
	assert_number(tmax, lower=0, finite=TRUE)

	tdf  <- data.frame(rate = exp(Xdf$eta), t = Xdf$tstart)
	rt.l <- split(tdf, f=Xdf$id)
	rt.l <- lapply(rt.l, as.list)

	## simulate new survival
	new_df <- data.frame(
			id   = unique(Xdf$id),
			time = round(sapply(rt.l, do.call, what = rpexp), 2) + 0.01)	%>% # +0.01 to avoid 0.00 times
		mutate(
			status = (time <= tmax)*1,
			time   = ifelse(time > tmax, tmax, time))
	first        <- !(duplicated(Xdf$id))
	new_df$Z     <- Xdf$Z[first, ]
	new_df$te_df <- Xdf$te_df[first, ]

	new_df

}
