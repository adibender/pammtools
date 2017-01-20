#' Add predicted hazard to data set
#' 
#' @inheritParams mgcv::predict.gam
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}}
#' @import checkmate dplyr mgcv 
#' @importFrom magrittr %<>%
#' @importFrom stats predict
#' @examples
#' library(mgcv)
#' data("leuk2", package="bpcp")
#' leuk.ped <- split_data(Surv(time, status)~., data=leuk2, id="id")
#' pam <- gam(status ~ s(tend), data = leuk.ped, family=poisson(), offset=offset)
#' pinfo <- ped_info(leuk.ped)
#' add_hazard(pinfo, pam)
#' @export 
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pam]{add_cumhazard}}
#' @rdname add_hazard
add_hazard <- function(newdata, object, se.fit=TRUE, type=c("response", "link"), ...)  {

	assert_data_frame(newdata, all.missing=FALSE) 
	assert_class(object, classes = c("gam", "glm", "lm"))
	type <- match.arg(type)

	pred <- predict(object=object, newdata=newdata, se.fit=se.fit, type=type, ...)
	stopifnot(length(pred$fit) == nrow(newdata))

	fit <- as.numeric(pred$fit)
	newdata %<>% mutate(hazard=fit)
	if(se.fit) {
		se <- as.numeric(pred$se.fit) 
		newdata %<>% mutate(
			lower  = hazard - 2*se, 
			upper  = hazard + 2*se)
	}

	return(newdata)
	
}


#' Add cumulative hazard estimate to data set 
#' 
#' @inheritParams add_hazard
#' @examples
#' library(mgcv)
#' data("leuk2", package="bpcp")
#' leuk.ped <- split_data(Surv(time, status)~., data=leuk2, id="id")
#' pam <- gam(status ~ s(tend), data = leuk.ped, family=poisson(), offset=offset)
#' pinfo <- ped_info(leuk.ped)
#' add_cumhazard(pinfo, pam)
#' @export 
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pam]{add_hazard}}
#' @rdname add_hazard
add_cumhazard <- function(newdata, object, se.fit=TRUE, type="response", ...) {

	assert_data_frame(newdata, all.missing=FALSE) 
	assert_class(object, classes = c("gam", "glm", "lm"))
	pred <- predict(object=object, newdata=newdata, se.fit=se.fit, type=type, ...)
	stopifnot(length(pred$fit) == nrow(newdata))

	newdata %>% mutate(
		cumhazard = cumsum(pred$fit*intlen),
		cumlower  = cumsum(pred$fit*intlen - 2*pred$se.fit),
		cumupper  = cumsum(pred$fit*intlen + 2*pred$se.fit))

}



