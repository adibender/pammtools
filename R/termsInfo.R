#' Add info about term effects to data set 
#' 
#' @inheritParams mgcv::predict.gam
#' @param term A string or regular expression indicating for which term 
#' information should be extracted and added to data set. 
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}}
#' @import checkmate dplyr mgcv 
#' @importFrom magrittr %<>%
#' @importFrom stats predict
#' @examples
#' library(mgcv)
#' data("leuk2", package="bpcp")
#' leuk.ped <- split_data(Surv(time, status)~., data=leuk2, id="id")
#' pam <- gam(status ~ ti(tend) + ti(tend, by=as.ordered(treatment)), 
#' 	data = leuk.ped, family=poisson(), offset=offset)
#' pinfo <- ped_info(leuk.ped)
#' add_term(pinfo, pam, term="treatment")
#' @export 
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pam]{add_hazard}}
add_term <- function(newdata, object, term, se.fit=TRUE, type="terms", ...) {

	assert_data_frame(newdata, all.missing=FALSE) 
	assert_class(object, classes = c("gam", "glm", "lm"))
	assert_string(term)
	pred <- predict(object=object, newdata=newdata, se.fit=se.fit, type=type, ...)
	ind.term <- grep(term, colnames(pred$fit))
	stopifnot(nrow(pred$fit) == nrow(newdata))

	newdata %>% mutate(
		tmp.fit    = pred$fit[, ind.term],
		tmp.se     = pred$se.fit[, ind.term],
		term       = tmp.fit,
		term.lower = tmp.fit - 2*tmp.se,
		term.upper = tmp.fit + 2*tmp.se) %>%
	select(-tmp.fit, -tmp.se)

}



#' Extract plot information for all special model terms
#' 
#' Given a \code{mgcv} \code{\link[mgcv]{gamObject}}, returns the information 
#' used for the default plots produced by \code{\link[mgcv]{plot.gam}}.
#' 
#' @inheritParams mgcv::plot.gam
#' @param ... Further arguments passed to \code{\link[mgcv]{plot.gam}}
#' @import mgcv
#' @importFrom checkmate assert_class
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot
#' @export 
get_plotinfo <- function(x, ...) {

	assert_class(x, c("gam", "glm", "lm"))

	tmp <- paste0(tempfile(), ".png")
	png(tmp)
	po <- plot(x, page=1, ...)
	dev.off()
	if(file.exists(tmp)) file.remove(tmp)

	class(po) <- c("mgcv.plotlist", class(po))

	return(po)

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
	ci = TRUE, 
	...) {

	po <- get_plotinfo(x, ...)
	# index of list elements that are 1d smooths and not random effects 
	ind.1d <- vapply(
		X         = po,
		FUN       = function(z) !is.null(z$x) & is.null(z$main),
		FUN.VALUE = logical(1))
	# keep only variables of interes
	po <- lapply(po[ind.1d], "[", i=keep, drop=TRUE)
	# use cbind.data.frame here, b/c as_data_frame does not work here 
	po <- lapply(po, function(z) do.call(cbind.data.frame, c(z, stringsAsFactors=FALSE)))
	if(ci) {
		po <- lapply(po, function(z) {
			z$low  = z$fit - z$se
			z$high = z$fit + z$se
			z
		})
	}

	return(bind_rows(po))

}


#' Extract random effects objects in tidy data format.
#' 
#' @inheritParams tidy_smooth
#' @importFrom dplyr bind_rows
#' @importFrom stats ppoints qnorm quantile
#' @rdname tidy_smooth
#' @export 
tidy_re <- function(x, keep=c("fit", "main", "xlab", "ylab"), ...) {

	po <- get_plotinfo(x, ...)
	ind.re <- vapply(
		X         = po,
		FUN       = function(z) !is.null(z$main) & z$xlab == "Gaussian quantiles",
		FUN.VALUE = logical(1))

	po <- lapply(po[ind.re], "[", i=keep, drop=TRUE)
	po <- lapply(po, function(z) {
		re.df = do.call(cbind.data.frame, c(z, stringsAsFactors=FALSE))
		re.df$x = qnorm(ppoints(length(re.df$fit))[order(order(re.df$fit))])
		yl <- quantile(re.df$fit, probs=c(0.1, 0.9), type=7, names=FALSE)
		xl <- qnorm(c(0.1, 0.9))
		re.df$qqslope <- diff(yl)/diff(xl)
		re.df$qqintercept <- yl[1L] - re.df$qqslope*xl[1L]

		re.df

	})

	return(bind_rows(po))

}
