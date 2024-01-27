ped_classes <- function(ped) {

  ind_ped <- class(ped) %in% c("ped", "ped_cr", "ped_cr_union", "fped", "nested_fdf")
  class(ped)[ind_ped]

}

re_attribute <- function(.data, attr2) {

  attr1 <- attributes(.data)
  attributes(.data) <- c(attr1,
    attr2[setdiff(names(attr1), names(attr2))])

  .data

}




#' @importFrom purrr discard
unped <- function(ped, classes_ped = "ped") {

  class(ped) <- setdiff(class(ped), classes_ped)

  ped

}

reped <- function(.data, ped_classes = "ped") {

  class(.data) <- c(ped_classes, class(.data))
  .data

}

ped_attr <- function(
  ped,
  ped_attributes = c("breaks", "id_var", "intvars", "combine", "censor_code", "risks")
) {

  attr_ped <- attributes(ped)
  ped_attr_avail <- intersect(names(attr_ped), ped_attributes)

  attr_ped[ped_attr_avail]

}

unfped <- function(fped) {
  class(fped) <- class(fped) %>% discard(~.=="fped")
  fped
}

refped <- function(.data) {
  class(.data) <- c("fped", class(.data))
  .data
}

fped_attr <- function(fped) {
  attributes(fped)[c("breaks", "id_var", "intvars")]
}


#' @name dplyr_verbs
#' @title \code{dplyr} Verbs for \code{ped}-Objects
#' @param .data an  object of class \code{ped}, see \code{\link{as_ped}}.
#' @param tbl an  object of class \code{ped}, see \code{\link{as_ped}}.
#' @param x an  object of class \code{ped}, see \code{\link{as_ped}}.
#' @param ... see \code{dplyr} documentation
#' @description See \code{dplyr} documentation of the respective functions for
#'   description and examples.
#' @return a modified \code{ped} object (except for \code{do})
#' @import dplyr
#' @aliases arrange filter distinct full_join group_by inner_join left_join mutate rename right_join sample_frac sample_n select slice summarise transmute ungroup
#' @keywords internal
NULL

#-------------------------------------------------------------------------------
# single table: grouping/sorting

#' @export
#' @export arrange
#' @rdname dplyr_verbs
arrange.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- arrange(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}


#' @export
#' @export group_by
#' @rdname dplyr_verbs
group_by.ped <- function(.data, ..., .add = FALSE) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- dplyr::group_by(unped(.data, classes_ped), ..., .add = .add)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}



#' @export
#' @export ungroup
#' @rdname dplyr_verbs
ungroup.ped <- function(x, ...) {

  classes_ped <- ped_classes(x)
  attr_ped    <- attributes(x)
  x           <- ungroup(unped(x, classes_ped), ...)
  x           <- reped(x, classes_ped)

  re_attribute(x, attr_ped)

}


#' @export
#' @export distinct
#' @rdname dplyr_verbs
distinct.ped <- function(.data, ..., .keep_all = FALSE) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data           <- distinct(unped(.data, classes_ped), ..., .keep_all = .keep_all)
  .data           <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}

#-------------------------------------------------------------------------------
# single table: row ops

#' @export
#' @export filter
#' @rdname dplyr_verbs
filter.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- filter(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}

#' @export
#' @export sample_n
#' @inheritParams dplyr::sample_n
#' @rdname dplyr_verbs
sample_n.ped <- function(tbl, size, replace = FALSE, weight = NULL, .env = NULL, ...) {

  classes_ped <- ped_classes(tbl)
  attr_ped    <- attributes(tbl)
  tbl         <- sample_n(unped(tbl, classes_ped), size, replace, weight, .env, ...)
  tbl         <- reped(tbl, classes_ped)

  re_attribute(tbl, attr_ped)

}

#' @export
#' @export sample_frac
#' @inheritParams dplyr::sample_frac
#' @rdname dplyr_verbs
sample_frac.ped <- function(tbl, size = 1, replace = FALSE, weight = NULL, .env = NULL, ...) {

  classes_ped <- ped_classes(tbl)
  attr_ped    <- attributes(tbl)
  tbl         <- sample_n(unped(tbl, classes_ped), size, replace, weight, .env, ...)
  tbl         <- reped(tbl, classes_ped)

  re_attribute(tbl, attr_ped)

}

#' @export
#' @export slice
#' @rdname dplyr_verbs
slice.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- slice(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}

#-------------------------------------------------------------------------------
# single table: column ops

#' @export
#' @export select
#' @rdname dplyr_verbs
select.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- select(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}

#' @export
#' @export mutate
#' @rdname dplyr_verbs
mutate.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- mutate(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}


#' @export
#' @export rename
#' @rdname dplyr_verbs
rename.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- rename(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}

#' @export
#' @export summarise
#' @rdname dplyr_verbs
summarise.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- summarise(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}
#' @export
#' @rdname dplyr_verbs
summarize.ped <- summarise.ped

#' @export
#' @export transmute
#' @rdname dplyr_verbs
transmute.ped <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- transmute(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}

#-------------------------------------------------------------------------------
# joins

#' @inheritParams dplyr::inner_join
#' @export
#' @export inner_join
#' @rdname dplyr_verbs
inner_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {

  classes_ped_x <- ped_classes(x)
  classes_ped_y <- ped_classes(y)
  attr_ped_x    <- attributes(x)
  .data         <- inner_join(unped(x, classes_ped_x), unped(y, classes_ped_y), by, copy, suffix, ...)
  .data         <- reped(.data, classes_ped_x)

  re_attribute(.data, attr_ped_x)

}

#' @inheritParams dplyr::full_join
#' @export
#' @export full_join
#' @rdname dplyr_verbs
full_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {

  classes_ped_x <- ped_classes(x)
  classes_ped_y <- ped_classes(y)
  attr_ped_x    <- attributes(x)
  .data         <- full_join(unped(x, classes_ped_x), unped(y, classes_ped_y), by, copy, suffix, ...)
  .data         <- reped(.data, classes_ped_x)

  re_attribute(.data, attr_ped_x)

}

#' @inheritParams dplyr::left_join
#' @export
#' @export left_join
#' @rdname dplyr_verbs
left_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {

  classes_ped_x <- ped_classes(x)
  classes_ped_y <- ped_classes(y)
  attr_ped_x    <- attributes(x)
  .data         <- left_join(unped(x, classes_ped_x), unped(y, classes_ped_y), by, copy, suffix, ...)
  .data         <- reped(.data, classes_ped_x)

  re_attribute(.data, attr_ped_x)

}

#' @inheritParams dplyr::right_join
#' @export
#' @export right_join
#' @rdname dplyr_verbs
right_join.ped <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ...) {

  classes_ped_x <- ped_classes(x)
  classes_ped_y <- ped_classes(y)
  attr_ped_x    <- attributes(x)
  .data         <- inner_join(unped(x, classes_ped_x), unped(y, classes_ped_y), by, copy, suffix, ...)
  .data         <- reped(.data, classes_ped_x)

  re_attribute(.data, attr_ped_x)

}
