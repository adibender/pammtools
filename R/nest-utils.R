# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2017-06-29 13:01:13
# @Last Modified by:   andreas.bender@stat.uni-muenchen.de
# @Last Modified time: 2017-06-29 13:45:29

#' Transform a nested data frame to data frame with matrix columns
#' 
#' @param nested_df A nested data frame with list columns that contain 
#' a matrix per row. 
#' @importFrom dplyr select_if slice
#' @importFrom purrr compose
#' @importFrom magrittr "%<>%" "%>%"
as_matdf <- function(nested_df) {

	ndf  <- select_if(nested_df, compose("!", is.list))
	ldf  <- select_if(nested_df, is.list)
	nrep <- sapply(ldf[[1]], function(z) nrow(z))
	rind <- lapply(seq_len(nrow(ndf)), function(z) rep(z, nrep[z])) %>%
		unlist()

	ldf <- lapply(ldf, function(z) do.call(rbind, z))
	ndf %<>% slice(rind) %>% as.data.frame()
	for(i in seq_along(ldf)) {
		ndf[[names(ldf)[i]]] <- ldf[[i]]
	}

	return(ndf)

}