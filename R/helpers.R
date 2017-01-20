#' Calculate the modus
#' 
#' @param var A atomic vector
#' @importFrom checkmate assert_atomic_vector
modus <- function(var) {

	# input checks 
	assert_atomic_vector(var, all.missing=FALSE, min.len=1)

	# calculate modus 
  freqs <- table(var)
  mod <- names(freqs)[which.max(freqs)]
  
  # factors should be returned as factors with all factor levels
  if(is.factor(var)) {
  	mod <- factor(mod, levels=levels(var))
  }

  return(mod)

}


#' Extract information of the sample contained in a data set 
#' 
#' Given a data set and grouping variables, this function returns median values 
#' for numeric variables and modus for characters and factors.
#' 
#' @param x A data frame (or object that inherits from \code{data.frame}).
#' @param ... Further arguments passed to specialized methods.
#' @importFrom stats median
#' @export
#' @rdname sample_info
sample_info <- function(x, ...) {
	UseMethod("sample_info", x)
}

#' @inheritParams sample_info
#' @import checkmate dplyr
#' @importFrom magrittr %<>%
#' @rdname sample_info
sample_info.data.frame <- function(x, ...) {

	assert_data_frame(x, all.missing=FALSE, min.rows=1, min.cols=1)

	bind_cols(
		summarize_if(x, .predicate=function(column) is.numeric(column), funs(median(., na.rm=TRUE))),
		summarize_if(x, .predicate=function(column) !is.numeric(column), modus))

}


#' @inheritParams sample_info
#' @import checkmate dplyr
#' @importFrom magrittr %<>%
#' @rdname sample_info
#' @seealso \code{\link[pam]{split_data}}
sample_info.ped <- function(x, ...) {

	# remove "noise" information on interval variables 
	x %<>% select(-one_of(attr(x, "intvars")))
	sample_info.data.frame(x, ...)

}