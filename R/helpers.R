#' Calculate the modus
#' 
#' @param var A atomic vector
#' @importFrom checkmate assert_atomic_vector
modus <- function(var) {

	# input checks 
	assert_atomic_vector(var, all.missing=FALSE, min.len=1)

	# calculate modus 
  freqs <- table(var)
  mod   <- names(freqs)[which.max(freqs)]
  
  # factors should be returned as factors with all factor levels
  if(is.factor(var)) {
  	mod <- factor(mod, levels=levels(var))
  }

  return(mod)

}

#' Creates sequence from minimum to maximum
#' 
#' @param x A numeric or integer vector.
#' @inheritParams base::seq
#' @import checkmate
seq_range <- function(x, length.out=100L) {

  assert_numeric(x, finite=TRUE, all.missing=FALSE, min.len=2)
  assert_integer(length.out, lower=2)

  range.x <- range(x)
  seq(range.x[1], range.x[2], length.out=length.out)

}


#' Combines multiple data frames
#' 
#' @importFrom dplyr slice bind_cols combine
#' @importFrom purrr map2 transpose cross
#' @importFrom checkmate test_data_frame
#' @param ... Data frames that should be combined to one data frame. 
#' Elements of first df vary fastest, elements of last df vary slowest
#' @examples 
#' combine_df(
#'   data.frame(x=1:3, y=3:1), 
#'   data.frame(x1=c("a", "b"), x2=c("c", "d")),
#'   data.frame(z=c(0, 1)))
#' @export
combine_df <- function(...) {

  dots <- list(...)
  if(!all(sapply(dots, test_data_frame))) {
    stop("All elements in ... must inherit from data.frame!")
  }
  seq.list <- lapply(dots, function(z) seq_len(nrow(z)))
  ind.list <- cross(seq.list) %>% transpose() %>% lapply(combine)
  
  map2(dots, ind.list, function(.x, .y) slice(.x, .y)) %>% bind_cols()

}


#' Construct a data frame suitable for prediction
#' 
#' Given a data set of class \code{ped}, returns a data frame that can be used 
#' as \code{newdata} argument in a call to \code{predict} and similar functions. 
#' When \code{expand==NULL} falls back to \code{\link{ped_info}}.
#' 
#' @inheritParams ped_info
#' @param expand A charachter vector of column names in \code{ped}.
#' @importFrom checkmate assert_data_frame
construct_newdata <- function(ped, expand=NULL) {

  assert_data_frame(ped, all.missing=FALSE, min.rows=2, min.cols=1)
  assert_character(min.chars=1, any.missing=FALSE)

  if(is.null(expand)) {
    return(ped_info(ped))
  }

  if(!all(expand %in% names(ped)) | any(is.null(expand))) {
    stop("All arguments provided in ... must be named and have names equal to 
      column names of ped object")
  }

}


#' extract character vector of grouping variables 
#' 
#' @importFrom dplyr groups
#' @param data A data frame (potentially `grouped_df`) from which to extract 
#' names of grouping variables. 
#' @return A character vector (of length 0 if no grouping variables present). 
get_grpvars <- function(data) {
  vapply(groups(data), as.character, character(1))
}

#' return ungrouped data frame without grouping variables 
#' 
#' @param data A data frame (potentially `grouped_df`) from which grouping 
#' variables will be removed. 
#' @import dplyr
#' @return Returns \code{data} without grouping variables (and group properties). 
rm_grpvars <- function(data) {

  if(is.grouped_df(data)) {
    grp.vars <- get_grpvars(data)
    data %<>% ungroup() %>% select(-one_of(grp.vars))
  } 
  return(data)

}