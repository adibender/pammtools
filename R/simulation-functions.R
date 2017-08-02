#' Construct full data set with time-dependent covariate (TDC)
#' 
#' @param n Number of subjects to simulate.
#' @param m = the number of time-points at which TDC is observed.
#' @param h0 A constant component of the baseline hazard.
#' @param f0 A function of t, representing the smooth part of the baseline hazard.
#' @param fte A function of the exposure time, representing the functional 
#' effect of the timing of the TDC.
#' @import checkmate
#' @importFrom stats runif
#' @keywords internal
#' @export 
make_X <- function(
	n   = 500L,
	m   = 30L,
	h0  = -3,
	f0  = function(t) (0*t),
	fte = function(te) (-2* te + 2.5*te^2)) {

	## check inputs
	assert_int(n, lower=2)
	assert_int(m, lower=1)
	assert_number(h0, lower=-5, upper=5)
	assert_function(f0)
	assert_function(fte)

	# create time-dependent covariate
	Z <- lapply(seq_len(n), function(z) {
		runif(m, 0, 1)
	})
	Z <- do.call(rbind, Z)
	# weight by weight function
	fZ <- lapply(seq_len(n), function(z) {
		fte(seq_len(m))
	})
	Fz <- do.call(rbind, fZ)
	Fz <- Z * Fz
	colnames(Fz) <- paste0("V", seq_len(ncol(Fz)))
	## craete data frame (start stop format)
	id <- rep(1:nrow(Fz), each=m)

	## output data set 
	Xdf        <- data.frame(id=id)
	Xdf$tstart <- rep(0:(m-1), times = n)
	Xdf$tend   <- rep(1:m, times = n)
	Xdf$intmid <- Xdf$tstart + 0.5
	Xdf$Z      <- Z[id, ]
	Xdf$Fz     <- Fz[id, ]
	Xdf$z_vec  <- as.vector(t(Z))

	## Lag Lead Matrix
	Lmat <- matrix(1, nrow=m, ncol = m)
	Lmat[upper.tri(Lmat)] <- 0
	Xdf$LL <- Lmat[rep(1:m, times=n), ]

	## add baseline + cumulative effects
	Xdf$eta_base <- h0 + f0(Xdf$intmid, tmax=max(Xdf$intmid))
	Xdf$eta_wce  <- rowSums(Xdf$Fz*Xdf$LL)
	Xdf$eta      <- Xdf$eta_base + Xdf$eta_wce

	Xdf$te_df   <- matrix(1:m, nrow=nrow(Xdf), ncol=m, byrow=TRUE)
	Xdf$time_df <- matrix(1:m, nrow=nrow(Xdf), ncol=m)

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