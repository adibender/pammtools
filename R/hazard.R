#' Add predicted hazard to data set
#'
#' @inheritParams mgcv::predict.gam
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}}
#' @param ci Logical indicating whether to iclude confidence intervals. Defaults
#' to \code{TRUE}
#' @param se.mult Factor by which standard errors are multiplied for calculating
#' the confidence intervals.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to \code{FALSE}. If \code{TRUE}, columns with names
#' \code{c("hazard", "se", "lower", "upper")} will be overwritten.
#' @import checkmate dplyr mgcv
#' @importFrom magrittr %<>%
#' @importFrom stats predict
#' @examples
#' library(mgcv)
#' data("leuk2", package="bpcp")
#' leuk.ped <- split_data(Surv(time, status)~., data=leuk2, id="id")
#' pam <- gam(ped_status ~ s(tend), data = leuk.ped, family=poisson(), offset=offset)
#' pinfo <- ped_info(leuk.ped)
#' add_hazard(pinfo, pam)
#' @export
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pam]{add_cumhazard}}
#' @rdname add_hazard
add_hazard <- function(
	newdata,
	object,
	type      = c("response", "link"),
	ci        = TRUE,
	se.mult   = 2,
	overwrite = FALSE,
	...)  {

	if(!overwrite) {
		if("hazard" %in% names(newdata)) {
			stop("Data set already contains 'hazard' column. Set `overwrite=TRUE` to overwrite")
		}
	} else {
			rm.vars <- intersect(c("hazard", "se", "lower", "upper"), names(newdata))
			newdata %<>% select(-one_of(rm.vars))
	}

	pred <- get_hazard(newdata, object, ci=ci, type=type, se.mult=se.mult, ...)
	stopifnot(nrow(pred) == nrow(newdata))

	newdata %<>% bind_cols(rm_grpvars(pred))

	return(newdata)

}

#' Calculate predicted hazard
#'
#' @inheritParams add_hazard
#' @rdname add_hazard
get_hazard <- function(
	newdata,
	object,
	ci      = TRUE,
	type    = c("response", "link"),
	se.mult = 2,
	...)  {

	assert_data_frame(newdata, all.missing=FALSE)
	assert_class(object, classes = "glm")
	type <- match.arg(type)

	pred <-
	  predict(object=object, newdata = newdata, se.fit = TRUE, type=type,
	    ...)[c("fit", "se.fit")] %>%
		bind_cols() %>%
		rename(hazard = fit, se = se.fit) %>%
		mutate(
			hazard = as.numeric(hazard),
			se     = as.numeric(se))
	stopifnot(nrow(pred) == nrow(newdata))

	if(ci) {
		pred %<>%
			mutate(
				lower = hazard - se.mult*se,
				upper = hazard + se.mult*se)
	}

	# it is necessary to include the grouping variables here, otherwise
	# functions calculating the cumulative hazard will cumulate over all rows
	# instead of group wise
	if(is.grouped_df(newdata)) {
		group.df <- select_(newdata, .dots=unlist(groups(newdata)))
		pred     <- bind_cols(group.df, pred)
	}

	return(pred)

}

#' Add cumulative hazard estimate to data set
#'
#' @inheritParams add_hazard
#' @export
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pam]{add_hazard}}
#' @rdname add_hazard
add_cumhazard <- function(
	newdata,
	object,
	type      = c("response"),
	ci        = TRUE,
	se.mult   = 2,
	overwrite = FALSE,
	...)  {

	if(!overwrite) {
		if("cumhazard" %in% names(newdata)) {
			stop("Data set already contains 'hazard' column. Set `overwrite=TRUE` to overwrite")
		}
	} else {
			rm.vars <- intersect(c("cumhazard", "cumlower", "cumupper"), names(newdata))
			newdata %<>% select(-one_of(rm.vars))
	}

	pred <- get_cumhazard(newdata, object, ci=ci, type=type, se.mult=se.mult, ...)

	newdata %<>% bind_cols(rm_grpvars(pred))

	return(newdata)

}

#' Calculate cumulative hazard
#'
#' @inheritParams add_cumhazard
#' @rdname add_hazard
get_cumhazard <- function(
	newdata,
	object,
	ci   = TRUE,
	type = "response",
	...) {

	hazard.df <- get_hazard(newdata, object, ci=TRUE, type=type)
	hazard.df %>%
		bind_cols(select_(rm_grpvars(newdata), .dots="intlen"))%>%
		mutate(
			cumhazard = cumsum(hazard * intlen),
			cumlower  = cumsum(lower  * intlen),
			cumupper  = cumsum(upper  * intlen)) %>%
		select(-one_of(c("hazard", "se", "lower", "upper", "intlen")))

}


