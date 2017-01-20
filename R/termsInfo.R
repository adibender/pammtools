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

