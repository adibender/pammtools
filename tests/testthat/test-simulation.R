context("Simulation functions")
## number of subjects 
n <- 10
## number of time-points (m) at which time-dependent covariate was observed 
# (here equal to maximal follow-up time (tmax))
tmax <- 30
## constant baseline hazard contribution 
h0 <- -3
## smooth baseline hazard contribution
f0 <- function(t, ...) {
	exp(-.1*(t))
}
## weight function lattency
# te = time of exposure
fz_dlnm <- function(t, te, z, tscale=1, zscale=1,....) {

	lattency <- (t-te)/length(te)
  (tscale*(-.6 + 2.5*lattency - 2*lattency^2 + 0.5*lattency^3))*
  	(zscale*(0.5*z -0.5*z^2))

}
te <- -30:30

test_that("Simulation functions produce correct dimensions", {

	cnames_X <- c("id" ,"tstart", "tend", "intmid", "Z", "Fz", "LL", "eta_base", 
		"eta_wce", "eta", "te_df", "time_df")

	X <- make_X(
		n       = n,
		m       = tmax,
		te      = te,
		h0      = h0,
		f0      = f0,
		fz      = fz_dlnm,
		fwindow = function(te, t) { t-te <= tmax & t-te >=0},
		rng_z   = function(te) {runif(length(te), 0, 1)})

	assert_data_frame(X, any.missing=FALSE, nrows=tmax*n, ncols=12)
	assert_subset(colnames(X), cnames_X)
	for(i in cnames_X[c(5:7, 11, 12)]) {
		assert_matrix(X[[i]], nrows=tmax*n, ncols=length(te))
	}
# ggplot(Xdf_full, aes(x=tstart, y=eta_base, group=id)) + 
#   geom_line(alpha = 0.1)

	instance <- as_fped(Surv(time, status)~., 
		data=sim_wce2(X), 
		data_orig=X, 
		cut=0:30, max.end=T,
		ll_fun = function(te, t) { (t-te <= 30) & (t-te >=0)})

	for(i in cnames_X[c(5:7, 11, 12)]) {
		assert_matrix(instance[[i]], ncols=length(te))
	}

})