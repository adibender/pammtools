#' Transform data with one row per subject to functional PED format
#' 
#' @inherit split_data 
#' @param ll_fun Function with two arguments that indicates how lag-lead matrix 
#' should be contructed. First argument is the time on the scale of the exposure, 
#' second argument is the time on the scale of the follow-up. Internally the 
#' lag lead matrix will be constructed by calling \code{t(outer(te_vec, t_vec, ll_fun))}.
#' @importFrom purrr compose
#' @keywords internal 
#' @export
as_fped <- function(
	formula, 
	data, 
  data_orig,
	cut     = NULL,
	..., 
	max.end = FALSE,
	ll_fun  = function(te, t) {te <= t}, 
  keep =  c("id", "tstart", "eta_base", "eta_wce", "eta")) {

  ## assert that inputs have correct formats
  assert_class(formula, "formula")
  assert_data_frame(data, min.rows=2, min.cols=2)
  assert_data_frame(data_orig, min.rows=nrow(data), min.cols=ncol(data))
  assert_numeric(cut, lower=0, finite=TRUE, any.missing=FALSE, min.len=1, 
    null.ok=TRUE)
  assert_flag(max.end)

  raw_df  <- data %>% select_if(compose("!", is.matrix))
  surv_df <- split_data(formula, raw_df, cut=cut, max.end=max.end)
  surv_df <- left_join(surv_df, data_orig[,keep]) %>% arrange(id, tstart)
  ## create functional exporsure variables
  surv_df$Z     <- data$Z[surv_df$id, ]
  surv_df$Fz    <- data$Z[surv_df$id, ]
  surv_df$te_df <- data$te_df[surv_df$id, ]

  ## crate Lag-Lead matrix
  times      <- attr(surv_df, "cut")
  te         <- surv_df$te_df[1,]
  LL         <- outer(attr(surv_df, "cut"), surv_df$te_df[1,], ll_fun)
  surv_df$LL <- Lmat(LL, surv_df$id)

 	## create time_df 
 	surv_df$time_df <- matrix(surv_df$tend, ncol=1)[, rep(1, ncol(surv_df$Z))]

  ## create intmid var
  surv_df$intmid <- surv_df$tstart + (surv_df$tend - surv_df$tstart)/2

 	surv_df

}


#' Exctract id with most observations from fped object
#' 
#' @param fped A functional PED. 
#' @keywords internal
#' @export
extract_one <- function(fped) {
  rle_id <- rle(fped$id)
  max_id <- rle_id$values[[which.max(rle_id$lengths)]]

  subset(fped, id==max_id)

}