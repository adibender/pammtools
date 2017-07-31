#' Function to check if an object is of class \code{ped} or 
#' transform data into piece-wise exponential data (PED) format if possible.
#'
#' @inheritParams survival::survSplit
#' @param tdc_df A data set containing information on time-dependent covariates (TDC).
#' @param ... Further arguments passed to methods. 
#' @return A data frame class \code{ped} in piece-wise exponential data (PED) format.
as_ped <- function(x, data, tdc_df=NULL, ...) {
	UseMethod("as_ped", x)
}

as_ped.formula <- function(x, data, tdc_df=NULL, ...) {

	if (is.null(tdc_df)) {
		split_data(x, data, ...)
	}

}

as_ped.resample <- function(x, data, ...) {

	left_join(as_tibble(x), data)

}

#' @rdname as_ped
#' @param x any R object.
is.ped <- function(x) inherits(x, "ped")