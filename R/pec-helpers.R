#' Predict survival probabilities at specified time points
#' 
#' This function is needed to work with the \code{pec} package for estimation 
#' of the time-dependent C-Index and Brier-Score (prediction error)
#' 
#' @param object A PAM model object. 
#' @param newdata Data containing covariate information only, for which to 
#' predict survival probabilities. 
#' @param times Time points for which survival probabilities should be predicted.
#' @param ... Arguments passed to \code{\link[pam]{int_info}}.
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom tidyr spread_
predictSurvProb.pamm <- function(
	object, 
	newdata, 
	times, ...) {

	assert_class(object, "pamm")
	## to calculate survival probabilities
	int_df <- int_info(object) %>% 
		filter(tstart <= max(times))
	newdata %<>% mutate(.id = row_number())

	## calculate survival probabilities
	surv_df <- combine_df(int_df, newdata) %>% 
		group_by_(".id") %>% 
		add_survprob(object) %>% 
		select(one_of(".id", "tend", "survprob"))

	out_df <- get_intervals(object, times) %>% 
		select(one_of("times", "tend"))

	left_join(out_df, surv_df) %>% 
		select(one_of(".id", "times", "survprob")) %>% 
			spread_("times", "survprob") %>% 
			select(-one_of(".id")) %>% 
			as.matrix()

}


#' k-fold cross validated prediction error curve (pec)
#' 
#' @inheritParams split_data
#' @param formula_ped The formula used to transform data to PED (see 
#' \code{\link{split_data}}).
#' @param formula_pam The formula used to fit the PAM (see \code{\link{pam}}).
#' @param data The data in standard time-to-event data format which can be tansformed 
#' in PED (see \code{\link{split_data}}).
#' @param method The function to which \code{formula_pam} will be passed to fit the 
#' model in the k-th fold. 
#' @param k The number of folds to be used in cross-validation. 
#' @param formula_ipcw The formula used to calculate inverse probability censoring 
#' weights (IPCW). See \code{\link[pec]{ipcw}} for details.
#' @param ... Further arguments passed to \code{\link{split_data}}.
#' @import dplyr 
#' @importFrom modelr crossv_kfold
#' @importFrom magrittr "%<>%"
#' @importFrom purrr map map2
#' @importFrom tidyr unnest
#' @seealso \code{\link[pec]{pec}}, \code{\link[pec]{ipcw}}
#' @export
pec_cv <- function(
	data, 
	formula_pam, 
	formula_ped  = Surv(time, status)~.,
	k            = 10L,
	formula_ipcw = Surv(time, status)~1,
	cut          = NULL,
	method       = pamm,
	...)  {

	assert_data_frame(data)
	assert_int(k, lower=2)

	cv_df <- crossv_kfold(data, k=k)

	cv_df %>% 
		mutate(
			idx_train = map(train, as.integer),
			idx_test  = map(test, as.integer),
			pec       = map2(
				idx_train,
				idx_test,
				.f           = pec_pamm,
				formula_ped  = formula_ped,
				formula_pam  = formula_pam,
				data         = data,
				cut          = cut,
				formula_ipcw = formula_ipcw,
				... )) %>% 
		select(one_of(".id", "pec")) %>% unnest()
	
}


#' Fit PAM to train data and calculate PEC on test data 
#' 
#' @inherit pec_cv
#' @param idx_train Index of rows that indicate training set. 
#' @param idx_test Index of rows that indicate test set.
#' @param times The time points at which prediction error will be evaluated
#' @param n If \code{times} argument not specified, will create a sequence 
#' from minimal to maximal time of this length. 
#' @param ... further arguments passed to \code{\link{split_data}}.
#' @importFrom pec pec
#' @importFrom modelr seq_range
#' @importFrom prodlim Hist
pec_pamm <- function(
	data,
	idx_train, 
	idx_test, 
	formula_ped, 
	formula_pam, 
	formula_ipcw = Surv(time, status)~1,
	times = NULL,
	cut   = NULL,
	n     = 50L,
	...) {

	train <- data %>% slice(idx_train)
	test  <- data %>% slice(idx_test)

	if(is.null(times)) {
		times = c(0, round(modelr::seq_range(train_ped$tend, n=n), 2))
	}

	# Necessary, b/c pec package internally performs the below operation and 
	# then would throw error b/c prediction matrix (nrow(newdata)x tiems) and 
	# requested newdata (nrow(newdata) x new times)
	# would have different dimensions and pec would throw an error. 
	times <- times[times <= max(test$time)]

	train_ped <- split_data(formula_ped, train, cut=cut, ...)
	## either provided cut or default cut created by split_data
	## (makes sure same cut is used for test data)
	train_pam <- pamm(formula_pam, data=train_ped)
	
	pred <- predictSurvProb.pamm(train_pam, newdata=test, times=times)

	
	pec_df <- pec(
			pred, 
			data    = test,
			times   = times,
			exact   = FALSE,
			formula = formula_ipcw) %>%
		tidy_pec()

}

#' Tranform data from pec object to tidy format
#' 
#' @param pec An object of class pec (see \code{\link[pec]{pec}}).
#' @importFrom checkmate assert_class
#' @importFrom dplyr mutate
#' @importFrom  tidyr gather
tidy_pec <- function(pec) {

	assert_class(pec, "pec")

	as.data.frame(pec[c( "time", "AppErr")]) %>% 
		gather(method, brier, -time) %>% 
		mutate(method = sub("AppErr.", "", method, fixed=TRUE))

}