
unped <- function(ped) {
  class(ped) <- class(ped)[-1]
  ped
}
reped <- function(.data) {
  class(.data) <- c("ped", class(.data))
  .data
}

ped_attr <- function(ped) {
  attributes(ped)[c("cut", "id_var", "intvars")]
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

#' @param keep_attributes conserve attributes? defaults to \code{TRUE}
#' @export
#' @export mutate
#' @rdname dplyr_verbs
mutate.ped <- function(.data, ..., keep_attributes=TRUE) {
  if (keep_attributes) {
    data_attr   <- ped_attr(.data)
  }
  .data <- reped(mutate(unped(.data), ...))
  if (keep_attributes) {
    attributes(.data) <- c(attributes(.data), data_attr)
  }
  return(.data)
}

#' @inheritParams dplyr::mutate_all
#' @export
#' @export mutate_each
#' @rdname dplyr_verbs
mutate_each.ped <- function(tbl, funs, ..., keep_attributes=TRUE) {
  if (keep_attributes) {
    data_attr   <- ped_attr(tbl)
  }
  tbl <- reped(mutate_each(unped(tbl), funs, ...))
  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
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
  ..., keep_attributes=TRUE) {

  if (keep_attributes) {
    data_attr   <- ped_attr(x)
  }
  #FIXME?
  tbl <- reped(left_join(unped(x), y, by, copy, suffix, ...))

  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }

  return(tbl)

}

#' @inheritParams dplyr::right_join
#' @export
#' @export right_join
#' @rdname dplyr_verbs
right_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ..., keep_attributes=TRUE) {

  if (keep_attributes) {
    data_attr   <- ped_attr(y)
  }
  #FIXME?
  tbl <- reped(right_join(unped(x), y, by, copy, suffix, ...))

  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }

  return(tbl)
}



#' @name tidyr_verbs
#' @title \code{tidyr} Verbs for \code{ped}-Objects
#' @param data an  object of class \code{ped}, see \code{\link{split_data}}.
#' @return A modified \code{ped} object.
#' @importFrom tidyr fill fill_
#' @description See \code{tidyr} documentation of the respective functions for
#'   description and examples.
#' @aliases fill fill_
NULL

#' @inheritParams tidyr::fill
#' @param keep_attributes conserve attributes? defaults to \code{TRUE}.
#' @export
#' @export fill
#' @rdname tidyr_verbs
fill.ped <- function(data, ..., .direction=c("down", "up"), keep_attributes=TRUE) {
  if (keep_attributes) {
    data_attr   <- ped_attr(data)
  }
  tbl <- reped(fill(unped(data), ..., .direction=.direction))
  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }

  return(tbl)

}

# #' @inheritParams tidyr::fill_
# #' @export fill_
# #' @export
# #' @rdname tidyr_verbs
# fill_.ped <- function(data, fill_cols, .direction=c("down", "up")) {

#   data_attr   <- ped_attr(data)
#   tbl <- reped(fill_(unped(data), fill_cols, .direction)) 
#   attributes(tbl) <- c(attributes(tbl), data_attr)

#   return(tbl)

# }