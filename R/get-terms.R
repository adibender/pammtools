#' Extract partial effects for specified model terms
#'
#' @param data Any data frame containing variables used to fit the model. Only
#' first row will be used.
#' @param fit A fitted object of class \code{\link[mgcv]{gam}}.
#' @param term The (non-linear) model term of interest.
#' @param ... Currently ignored.
#' @import magrittr dplyr
#' @importFrom stats predict setNames
get_term <- function(data, fit, term, ...) {

	range.term <- range(data[[term]], na.rm=TRUE)
	seq.term   <- seq(range.term[1], range.term[2], length.out = 100)

	newdf <- data[1, ]
	rm(data)
	gc()

	newdf         <- newdf[rep(1, length(seq.term)), ]
	newdf[[term]] <- seq.term
	pred.term     <- predict(fit, newdata=newdf, type="terms", se.fit=TRUE)
	ind.term      <- grep(term, colnames(pred.term$fit), value=TRUE)

	newdf %>% mutate(
		term     = term,
		eff      = as.numeric(pred.term$fit[, ind.term]),
		se       = as.numeric(pred.term$se.fit[, ind.term]),
		ci.lower = eff - 2*se,
		ci.upper = eff + 2*se) %>%
	select_(.dots=c("term", term, "eff", "se", "ci.lower", "ci.upper")) %>%
	rename_(.dots=setNames(term, "x"))

}

#' Extract the partial effects of non-linear model terms
#'
#' @inheritParams get_term
#' @param terms A character vector (can be length one). Specifies the terms
#' for which partial effects will be returned
#' @import checkmate
#' @return A data frame with 5 columns.
#' @seealso \code{\link[survival]{coxph}}
#' @export
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ pspline(karno, df=4), data=veteran)
#' term.karno <- get_terms(veteran, fit, terms="karno")
get_terms <- function(data, fit, terms, ...) {

  # check inputs
  assert_class(data, "data.frame")
  assert_character(terms, min.len=1, unique=TRUE)

  cols.term <- sapply(terms, grep, x=colnames(data), value=TRUE)

  # apply get_term to each element of terms
	term.list <- lapply(terms, function(z) {
		get_term(fit=fit, data=data, term=z, ...)
	})

	bind_rows(term.list)

}