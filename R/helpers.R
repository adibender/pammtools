# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2017-01-19 17:53:32
# @Last Modified by:   andreas.bender@stat.uni-muenchen.de
# @Last Modified time: 2017-01-19 19:02:32

#' Returns the modus of a variable
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
#' for numeric variables and modus for characters and factors. One row 
#' 
#' @param x A data frame (or object that inherits from \code{data.frame}).
#' @param ... Further arguments passed to specialized methods.
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
		summarize_if(x, .predicate=function(x) is.numeric(x), funs(median(., na.rm=TRUE))),
		summarize_if(x, .predicate=function(x) !is.numeric(x), modus))

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