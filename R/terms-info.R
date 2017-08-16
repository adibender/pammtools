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
		unlist() %>% 
		unique() %>% 
		sort()
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
        type    = "lpmatrix")[, col_ind, drop = FALSE]
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
			low  = fit - se.mult * se,
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