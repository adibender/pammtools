# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2017-03-20 17:57:02

# @Last Modified by:   fabian.scheipl@stat.uni-muenchen.de
# @Last Modified time: 2017-04-19


unped <- function(ped) {
  class(ped) <- class(ped)[-1]
  ped
}
reped <- function(.data) {
  class(.data) <- c("ped", class(.data))
  .data
}


#' @name dplyr_verbs
#' @title \code{dplyr} Verbs for \code{ped}-Objects
#' @param .data an  object of class \code{ped}, see \code{\link{split_data}}.
#' @param ... see \code{dplyr} documentation
#' @param .dots see \code{dplyr} documentation
#' @description See \code{dplyr} documentation of the respective functions for
#'   description and examples.
#' @return a modified \code{ped} object (except for \code{do})
#' @import dplyr
#FIXME: replace deprecated "underscore" verbs
NULL


#' @export
#' @export do
#' @rdname dplyr_verbs
do.ped <- function(.data, ...) {
  do(unped(.data), ...)
}

#' @export
#' @export filter
#' @rdname dplyr_verbs
filter.ped <- function(.data, ...) {
  reped(filter(unped(.data), ...))
}

#' @export
#' @export slice
#' @rdname dplyr_verbs
slice.ped <- function(.data, ...) {
  reped(slice(unped(.data), ...))
}

#' @export
#' @export arrange
#' @rdname dplyr_verbs
arrange.ped <- function(.data, ...) {
  reped(arrange(unped(.data), ...))
}

#' @export
#' @export select
#' @rdname dplyr_verbs
select.ped <- function(.data, ...) {
  reped(select(unped(.data), ...))
}

#' @export
#' @export select_
#' @rdname dplyr_verbs
select_.ped <- function(.data, ..., .dots = list()) {
  reped(select_(unped(.data), ..., .dots = .dots))
}

#' @export
#' @export rename
#' @rdname dplyr_verbs
rename.ped <- function(.data, ...) {
  reped(rename(unped(.data), ...))
}

#' @export
#' @export summarise
#' @rdname dplyr_verbs
summarise.ped <- function(.data, ...) {
  reped(summarise(unped(.data), ...))
}
#' @export
#' @rdname dplyr_verbs
summarize.ped <- summarise.ped

#' @export
#' @export rename_
#' @rdname dplyr_verbs
rename_.ped <- function(.data, ..., .dots = list()) {
  reped(rename_(unped(.data), ..., .dots = .dots))
}


#' @export
#' @export distinct_
#' @rdname dplyr_verbs
distinct_.ped <- function(.data, ..., .dots = list()) {
  reped(distinct_(unped(.data), ..., .dots = .dots))
}


#' @param keep.attributes conserve attributes? defaults to \code{TRUE}
#' @export
#' @export mutate
#' @rdname dplyr_verbs
mutate.ped <- function(.data, ..., keep.attributes=TRUE) {
	if (keep.attributes) {
		.data.attr   <- attributes(.data)
		attr.names <- setdiff(names(.data.attr), c("class", "row.names", "names"))
	}
  .data <- reped(mutate(unped(.data), ...))
	if (keep.attributes) {
		attributes(.data) <- c(attributes(.data), .data.attr[attr.names])
	}
	return(.data)
}

#' @export
#' @export transmute
#' @rdname dplyr_verbs
transmute.ped <- function(.data, ...) {
  reped(transmute(unped(.data), ...))
}


#' @export
#' @export sample_n
#' @inheritParams dplyr::sample_n
#' @rdname dplyr_verbs
sample_n.ped <- function(tbl, size, replace = FALSE, weight = NULL,
  .env = NULL) {
  reped(sample_n(unped(tbl), size, replace, weight, .env))
}


#' @export
#' @export sample_frac
#' @inheritParams dplyr::sample_frac
#' @rdname dplyr_verbs
sample_frac.ped <- function(tbl, size = 1, replace = FALSE, weight = NULL,
  .env = NULL) {
  reped(sample_frac(unped(tbl), size, replace, weight, .env))
}


#' @inheritParams dplyr::left_join
#' @export
#' @export left_join
#' @rdname dplyr_verbs
left_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  reped(left_join(unped(tbl), y, by, copy, suffix, ...))
}



#' @inheritParams dplyr::group_by_
#' @export
#' @export group_by_
#' @rdname dplyr_verbs
group_by_.ped <- function(.data, ..., .dots = list(), add = FALSE) {
  reped(group_by_(unped(.data), ..., .dots = .dots, add = add))
}


#' @export
#' @export ungroup
#' @rdname dplyr_verbs
ungroup.ped <- function(x, ...) {
  reped(ungroup(unped(x), ...))
}
