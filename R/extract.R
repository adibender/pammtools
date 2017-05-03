#' Calculate confidence intervals 
#' 
#' Given 2 column matrix or data frame, returns 3 column data.frame 
#' with coefficient estimate plus lower and upper borders of the 
#' 95% confidence intervals.
#' 
#' @param ftab A table with two columns, containing coefficients in the first 
#' column and standard-errors in the second column. 
#' @importFrom tibble as_tibble
calc_ci <- function(ftab) {

	colnames(ftab) <- c("coef", "se")
	rnames         <- rownames(ftab)
	ftab           <- as_tibble(ftab)
	ftab$variable <- rnames
	ftab$lower     <- ftab$coef - 2*ftab$se
	ftab$upper     <- ftab$coef + 2*ftab$se
	ftab$se        <- NULL

	ftab[, c("variable", "coef", "lower", "upper")]

}

#' Extract fixed coefficient table from model object
#' 
#' Given a model object, returns data frame with columns \code{variable},
#' \code{coef} (coefficient), \code{lower} (lower 95\% CI) and 
#' \code{upper} (upper 95% CI). 
#' 
#' @param x A model object.
#' @param ... Currently not used.
#' @export 
extract_fixed <- function(x, ...) {
	UseMethod("extract_fixed", x)
}

#' @inheritParams extract_fixed
#' @rdname extract_fixed
#' @export 
extract_fixed.gam <- function(x, ...) {

  ftab <- summary(x)[["p.table"]][, 1:2]
  calc_ci(ftab)

}

#' @inheritParams extract_fixed
#' @importFrom tibble as_tibble
#' @rdname extract_fixed
#' @export
extract_fixed.coxph <- function(x, ...) {

  ftab <- summary(x)[["coefficients"]][, c(1, 3)]
  calc_ci(ftab)

}
