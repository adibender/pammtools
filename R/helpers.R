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


#' Extract character vector of grouping variables 
#' 
#' @importFrom dplyr groups
#' @param data A data frame (potentially `grouped_df`) from which to extract 
#' names of grouping variables. 
#' @return A character vector (of length 0 if no grouping variables present). 
get_grpvars <- function(data) {
  vapply(groups(data), as.character, character(1))
}

#' Return ungrouped data frame without grouping variables 
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