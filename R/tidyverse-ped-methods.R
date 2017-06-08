
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
#' @param tbl an  object of class \code{ped}, see \code{\link{split_data}}.
#' @param x an  object of class \code{ped}, see \code{\link{split_data}}.
#' @param funs see \code{\link[dplyr]{summarize_all}}
#' @param ... see \code{dplyr} documentation
#' @param .dots see \code{dplyr} documentation
#' @description See \code{dplyr} documentation of the respective functions for
#'   description and examples.
#' @return a modified \code{ped} object (except for \code{do})
#' @import dplyr
#' @aliases arrange distinct_ filter full_join group_by group_by_ inner_join left_join mutate mutate_each rename rename_ right_join sample_frac sample_n select select_ slice summarise summarise_each transmute ungroup
# FIXME: replace deprecated "underscore" verbs, [summarise|mutate]_each
NULL

#-------------------------------------------------------------------------------
# single table: grouping/sorting

#' @export
#' @export arrange
#' @rdname dplyr_verbs
arrange.ped <- function(.data, ...) {
  reped(arrange(unped(.data), ...))
}

#' @inheritParams dplyr::group_by
#' @export
#' @export group_by
#' @rdname dplyr_verbs
group_by.ped <- function(.data, ..., add = FALSE) {
  reped(group_by(unped(.data), ..., add = add))
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

#-------------------------------------------------------------------------------
# single table: row ops

#' @export
#' @export distinct_
#' @rdname dplyr_verbs
distinct_.ped <- function(.data, ..., .dots = list()) {
  reped(distinct_(unped(.data), ..., .dots = .dots))
}

#' @export
#' @export filter
#' @rdname dplyr_verbs
filter.ped <- function(.data, ...) {
  reped(filter(unped(.data), ...))
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

#' @export
#' @export slice
#' @rdname dplyr_verbs
slice.ped <- function(.data, ...) {
  reped(slice(unped(.data), ...))
}

#-------------------------------------------------------------------------------
# single table: column ops

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

#' @inheritParams dplyr::mutate_all
#' @export
#' @export mutate_each
#' @rdname dplyr_verbs
mutate_each.ped <- function(tbl, funs, ..., keep.attributes=TRUE) {
  if (keep.attributes) {
    .data.attr   <- attributes(tbl)
    attr.names <- setdiff(names(.data.attr), c("class", "row.names", "names"))
  }
  tbl <- reped(mutate_each(unped(tbl), funs, ...))
  if (keep.attributes) {
    attributes(tbl) <- c(attributes(tbl), .data.attr[attr.names])
  }
  return(tbl)
}

#' @export
#' @export rename
#' @rdname dplyr_verbs
rename.ped <- function(.data, ...) {
  reped(rename(unped(.data), ...))
}

#' @export
#' @export rename_
#' @rdname dplyr_verbs
rename_.ped <- function(.data, ..., .dots = list()) {
  reped(rename_(unped(.data), ..., .dots = .dots))
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
#' @export summarise_each
#' @rdname dplyr_verbs
summarise_each.ped <- function(tbl, funs, ...) {
  reped(summarise_each(unped(tbl), funs,...))
}

#' @export
#' @export transmute
#' @rdname dplyr_verbs
transmute.ped <- function(.data, ...) {
  reped(transmute(unped(.data), ...))
}

#-------------------------------------------------------------------------------
# joins

#' @inheritParams dplyr::inner_join
#' @export
#' @export inner_join
#' @rdname dplyr_verbs
inner_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  reped(inner_join(unped(x), y, by, copy, suffix, ...))
}

#' @inheritParams dplyr::full_join
#' @export
#' @export full_join
#' @rdname dplyr_verbs
full_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  reped(full_join(unped(x), y, by, copy, suffix, ...))
}

#' @inheritParams dplyr::left_join
#' @export
#' @export left_join
#' @rdname dplyr_verbs
left_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  reped(left_join(unped(x), y, by, copy, suffix, ...))
}

#' @inheritParams dplyr::right_join
#' @export
#' @export right_join
#' @rdname dplyr_verbs
right_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  reped(right_join(unped(x), y, by, copy, suffix, ...))
}



#' @name tidyr_verbs
#' @title \code{tidyr} Verbs for \code{ped}-Objects
#' @param data an  object of class \code{ped}, see \code{\link{split_data}}.
#' @return a modified \code{ped} object.
#' @import tidyr
NULL

#' @rdname tidyr_verbs
#' @inheritParams tidyr::fill
#' @importFrom tidyr fill
#' @export fill
#' @export
fill.ped <- function(data, ..., .direction=c("down", "up")) {

  reped(fill(unped(data), ..., .direction))

}

#' @rdname tidyr_verbs
#' @inheritParams tidyr::fill_
#' @importFrom tidyr fill_
#' @export fill_
#' @export
fill_.ped <- function(data, fill_cols, .direction=c("down", "up")) {

  reped(fill_(unped(data), fill_cols, .direction))

}


