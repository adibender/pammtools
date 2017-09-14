## Functions to extract effect information from fitted models



#' Calculate confidence intervals 
#' 
#' Given 2 column matrix or data frame, returns 3 column data.frame 
#' with coefficient estimate plus lower and upper borders of the 
#' 95% confidence intervals.
#' 
#' @param ftab A table with two columns, containing coefficients in the first 
#' column and standard-errors in the second column. 
#' @importFrom tibble as_tibble
calc_ci <- function(ftab) {

	colnames(ftab) <- c("coef", "se")
	rnames         <- rownames(ftab)
	ftab           <- as_tibble(ftab)
	ftab$variable  <- rnames
	ftab$lower     <- ftab$coef - 2*ftab$se
	ftab$upper     <- ftab$coef + 2*ftab$se
	ftab$se        <- NULL

	ftab[, c("variable", "coef", "lower", "upper")]

}

#' Extract fixed coefficient table from model object
#' 
#' Given a model object, returns data frame with columns \code{variable},
#' \code{coef} (coefficient), \code{lower} (lower 95\% CI) and 
#' \code{upper} (upper 95% CI). 
#' 
#' @param x A model object.
#' @param ... Currently not used.
#' @export 
tidy_fixed <- function(x, ...) {
	UseMethod("tidy_fixed", x)
}

#' @inheritParams tidy_fixed
#' @param intercept Should intercept also be returned?
#' @rdname tidy_fixed
#' @export 
tidy_fixed.gam <- function(x, intercept=FALSE, ...) {

  ftab <- summary(x)[["p.table"]][, 1:2]
  if(!intercept) {
  	ftab <- ftab[!grepl("Intercept", rownames(ftab)), ]
  }
  calc_ci(ftab)

}

#' @inheritParams tidy_fixed
#' @importFrom tibble as_tibble
#' @rdname tidy_fixed
#' @export
tidy_fixed.coxph <- function(x, ...) {

  ftab <- summary(x)[["coefficients"]][, c(1, 3)]
  calc_ci(ftab)

}


#' Extract 1d smooth objects in tidy data format.
#'
#' @inheritParams get_plotinfo
#' @param keep A vector of variables to keep.
#' @param ci A logical value indicating whether confidence intervals should be
#' calculated and returned. Defaults to \code{TRUE}.
#' @importFrom dplyr bind_rows
#' @export
tidy_smooth <- function(
	x,
	keep = c("x", "fit", "se", "xlab", "ylab"),
	ci   = TRUE,
	...) {

	po <- get_plotinfo(x, ...)
	# index of list elements that are 1d smooths and not random effects
	ind1d <- vapply(
		X         = po,
		FUN       = function(z) !is.null(z[["x"]]) & is.null(z[["main"]]),
		FUN.VALUE = logical(1))
	# keep only variables of interest
	po <- lapply(po[ind1d], "[", i=keep, drop=TRUE)

	# transform to data.frame
	po <- lapply(po, function(z) {
		z[["fit"]] <- as.vector(z[["fit"]])
		temp <- as_tibble(z)
		if(ci) {
			temp %<>% mutate(
				low  = fit - se,
				high = fit + se)
		}
		temp
	})

	return(bind_rows(po))

}


#' Extract 2d smooth objects in tidy format.
#'
#' @inheritParams tidy_smooth
#' @importFrom purrr cross_df
#' @importFrom tibble as_tibble
#' @import dplyr
#' @export
tidy_smooth2d <- function(
	x,
	keep = c("x", "y", "fit","se", "xlab", "ylab", "main"),
	ci   = FALSE,
	...) {

	po <- get_plotinfo(x, ...)

	ind2d <- vapply(
		X         = po,
		FUN       = function(z) !is.null(z[["x"]]) & !is.null(z[["y"]]),
		FUN.VALUE = logical(1))

	# keep only variables of interes
	po <- lapply(po[ind2d], "[", i=keep, drop=TRUE)

	# transform to data.frame
	po <- lapply(po, function(z) {
		z[["fit"]] <- as.vector(z[["fit"]])
		p1 <- as_tibble(z[setdiff(keep, c("x", "y"))])
		xy <- cross_df(z[c("x", "y")])
		xy <- bind_cols(xy, p1)
		if(ci) {
			xy %<>% mutate(
				low  = fit - se,
				high = fit + se)
		}
		xy
	})

	return(bind_rows(po))

}


#' Extract random effects objects in tidy data format.
#'
#' @inheritParams tidy_smooth
#' @importFrom dplyr bind_rows
#' @importFrom stats ppoints qnorm quantile
#' @rdname tidy_smooth
#' @seealso \code{\link[stats]{qqline}}
#' @export
tidy_re <- function(x, keep=c("fit", "main", "xlab", "ylab"), ...) {

	po <- get_plotinfo(x, ...)
	ind.re <- vapply(
		X         = po,
		FUN       = function(z) !is.null(z[["main"]]) & z[["xlab"]] == "Gaussian quantiles",
		FUN.VALUE = logical(1))

	po <- lapply(po[ind.re], "[", i=keep, drop=TRUE)
	po <- lapply(po, function(z) {
		re.df = do.call(cbind.data.frame, c(z, stringsAsFactors=FALSE))
		re.df$x = qnorm(ppoints(length(re.df$fit))[order(order(re.df$fit))])
		# code to calculate qqslope and qqintercept from ?stats::qqline
		yl <- quantile(re.df$fit, probs=c(0.1, 0.9), type=7, names=FALSE)
		xl <- qnorm(c(0.1, 0.9))
		re.df$qqslope <- diff(yl)/diff(xl)
		re.df$qqintercept <- yl[1L] - re.df$qqslope*xl[1L]

		re.df

	})

	return(bind_rows(po))

}