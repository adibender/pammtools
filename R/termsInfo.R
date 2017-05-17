#' Add info about term effects to data set
#'
#' Adds the contribution (plus confidence intervals) of a specific term to the
#' linear predictor to the provided data.
#' Largely a wrapper to \code{\link[mgcv]{predict.gam}}, with \code{type="terms"}.
#' Thus most arguments and their documentation below is from \code{predict.gam}.
#'
#'
#' @inheritParams mgcv::predict.gam
#' @param term A character (vector) or regular expression indicating for
#' which term(s) information should be extracted and added to data set.
#' @param se.mult The factor by which standard errors are multiplied to form
#' confidence intervals.
#' @param relative If \code{TRUE}, calculates relative risk contribution,
#' that is \eqn{(X-\bar{X})'\beta} and respective confidence intervals
#' if \code{se.fit = TRUE}. Defaults to \code{FALSE}.
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}}
#' @import checkmate dplyr mgcv
#' @importFrom magrittr %<>%
#' @importFrom stats predict
#' @examples
#' library(mgcv)
#' data("veteran", package="survival")
#' ped <- split_data(Surv(time, status)~ age, data=veteran, id="id")
#' pam <- gam(ped_status ~ s(tend),
#' 	data = ped, family=poisson(), offset=offset)
#' pinf <- ped_info(ped)
#' add_term(pinf, pam, term="tend")
#' @export
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pam]{add_hazard}}
#' @importFrom stats model.matrix vcov
add_term <- function(
	newdata,
	object,
	term,
	se.fit   = TRUE,
	type     = "terms",
	se.mult  = 2,
	relative = FALSE,
	...) {

	assert_data_frame(newdata, all.missing=FALSE)
	assert_character(term, min.chars=1, any.missing=FALSE, min.len=1)

	col_ind <- lapply(term, grep, x=names(object$coefficients)) %>%
		unlist() %>% unique() %>% sort()
  is_pam <- inherits(object, "gam")

  X <- if (is_pam) {
    predict(object, newdata = newdata, type = "lpmatrix", ...)[,col_ind, drop=FALSE]
  } else  {
    model.matrix(object$formula[-2], data = newdata)[,col_ind, drop=FALSE]
  }
  if (relative) {
    if (is.null(object$model)) {
      stop("Relative risk can only be calculated when original data is present in 'object'!")
    }
    data_bar <- sample_info(object$model)[rep(1, nrow(newdata)), ]
    X_bar <- if (is_pam) {
      predict(
        object,
        newdata = data_bar,
        type = "lpmatrix")[, col_ind, drop = FALSE]
  	} else {
  	  model.matrix(object$formula[-2], data = data_bar)[, col_ind, drop = FALSE]
  	}
  	X <- X - X_bar
  }

  newdata[["fit"]] <- drop(X %*% object$coefficients[col_ind])
  if (se.fit) {
  	cov.coefs <- if(is_pam) {
  	  object$Vp[col_ind, col_ind]
  	} else {
  	  vcov(object)[col_ind, col_ind]
  	}
  	se <- sqrt(drop(diag(X %*% cov.coefs %*% t(X))))
		newdata %<>% mutate(
			low = fit - se.mult * se,
			high = fit + se.mult * se)
  }

  return(newdata)

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
	ind1d <- vapply(
		X         = po,
		FUN       = function(z) !is.null(z[["x"]]) & is.null(z[["main"]]),
		FUN.VALUE = logical(1))
	# keep only variables of interest
	po <- lapply(po[ind1d], "[", i=keep, drop=TRUE)

	if(ci) {
		po <- lapply(po, function(z) {
			z$low  = z$fit - z$se
			z$high = z$fit + z$se
			z
		})
	}

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
	ci = FALSE,
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
