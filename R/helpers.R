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