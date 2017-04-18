# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2017-03-20 17:57:02
# @Last Modified by:   fabian.scheipl@stat.uni-muenchen.de
# @Last Modified time: 2017-04-18

use_dplyr_verb <- function(dplyr_verb, .data, ...) {
  .data.class  <- class(.data)[1]
  class(.data) <- class(.data)[-1]
  args <- list(.data)
  args <- c(args, ...)
  .data        <- do.call(dplyr_verb, args)
  class(.data) <- c(.data.class, class(.data))
  .data
}

#' @name dplyr_verbs
#' @title \code{dplyr} Verbs for \code{ped}-Objects
#' @param .data an  object of class \code{ped}, see \code{\link{split_data}}.
#' @param ... see \code{dplyr} documentation
#' @description See \code{dplyr} documentation of the respective function for
#'   details and examples.
#' @return a modified \code{ped} object
NULL

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
filter.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::filter,.data, ...)
}


#' @import dplyr
#' @export
#' @rdname dplyr_verbs
filter_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::filter_,.data, ...)
}


#' @import dplyr
#' @export
#' @rdname dplyr_verbs
slice.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::slice,.data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
slice_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::slice_,.data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
arrange.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::arrange,.data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
arrange_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::arrange_,.data, ...)
}


#' @import dplyr
#' @export
#' @rdname dplyr_verbs
select_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::select_,.data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
rename_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::rename_,.data, ...)
}


#' @import dplyr
#' @export
#' @rdname dplyr_verbs
distinct_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::distinct_,.data, ...)
}


#' @import dplyr
#' @param keep.attributes conserve attributes? defaults to \code{TRUE}
#' @export
#' @rdname dplyr_verbs
mutate.ped <- function(.data, keep.attributes=TRUE, ...) {
	if (keep.attributes) {
		.data.attr   <- attributes(.data)
		attr.names <- setdiff(names(.data.attr), c("class", "row.names", "names"))
	}
  .data <- use_dplyr_verb(dplyr::mutate,.data, ...)
	if (keep.attributes) {
		attributes(.data) <- c(attributes(.data), .data.attr[attr.names])
	}
	return(.data)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
mutate_.ped <- function(.data, keep.attributes=TRUE, ...) {

	if (keep.attributes) {
	  .data.attr   <- attributes(.data)
		attr.names <- setdiff(names(.data.attr), c("class", "row.names", "names"))
	}
  .data <- use_dplyr_verb(dplyr::mutate_,.data, ...)
	if (keep.attributes) {
		attributes(.data) <- c(attributes(.data), .data.attr[attr.names])
	}
	return(.data)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
transmute.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::transmute,.data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
transmute_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::transmute_,.data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
summarise_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::summarise_,.data, ...)
}

#' @import dplyr
#' @inherit dplyr::sample_n params
#' @export
#' @rdname dplyr_verbs
sample_n.ped <- function(.data, ...) {
  .data.class  <- class(.data)[1]
  class(.data) <- class(.data)[-1]
  .data        <- sample_n(.data, ...)
  class(.data) <- c(.data.class, class(.data))
  return(.data)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
sample_frac.ped <- function(.data, ...) {
  .data.class  <- class(.data)[1]
  class(.data) <- class(.data)[-1]
  .data        <- sample_frac(.data, ...)
  class(.data) <- c(.data.class, class(.data))
  return(.data)
}

#' @import dplyr
#' @inherit dplyr::left_join params
#' @export
#' @rdname dplyr_verbs
left_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  use_dplyr_verb(dplyr::left_join, x, y, by, copy, suffix, ...)
}


#' @import dplyr
#' @export
#' @rdname dplyr_verbs
group_by_.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::group_by_, .data, ...)
}

#' @import dplyr
#' @export
#' @rdname dplyr_verbs
ungroup.ped <- function(.data, ...) {
  use_dplyr_verb(dplyr::ungroup, .data, ...)
}
