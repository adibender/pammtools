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

ped_attr <- function(ped) {
  attributes(ped)[c("breaks", "id_var", "intvars", "combine", "censor_code", "risks")]
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
#' @param funs see \code{\link[dplyr]{summarize_all}}
#' @param ... see \code{dplyr} documentation
#' @param .dots see \code{dplyr} documentation
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
  attr_ped <- attributes(x)
  x <- ungroup(unped(x, classes_ped), ...)
  x <- reped(x, classes_ped)

  re_attribute(x, attr_ped)

}

#-------------------------------------------------------------------------------
# single table: row ops

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
  .env = NULL, ...) {
  reped(sample_n(unped(tbl), size, replace, weight, .env, ...))
}

#' @export
#' @export sample_frac
#' @inheritParams dplyr::sample_frac
#' @rdname dplyr_verbs
sample_frac.ped <- function(tbl, size = 1, replace = FALSE, weight = NULL,
  .env = NULL, ...) {
  reped(sample_frac(unped(tbl), size, replace, weight, .env, ...))
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
  out <- NextMethod()

  # I don't know if this is the correct thing to do here; you'll
  # need to think about what attributes should be preserved
  attr <- attributes(.data)
  attr$names <- names(out)
  attributes(out) <- attr

  out
}

#' @param keep_attributes conserve attributes? defaults to \code{TRUE}
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

## ' @param keep_attributes conserve attributes? defaults to \code{TRUE}
## ' @export
## ' @export mutate
## ' @rdname dplyr_verbs
## mutate.ped_cr <- function(.data, ..., keep_attributes = TRUE) {
#   if (keep_attributes) {
#     data_attr   <- ped_attr(.data)
#     data_class <- class(.data)
#   }
#   .data <- mutate(unped(.data), ...)
#   if (keep_attributes) {
#     attributes(.data) <- c(attributes(.data), data_attr)
#   }
#   class(.data) <- data_class

#   return(.data)

# }
##' @param keep_attributes conserve attributes? defaults to \code{TRUE}
##' @export
##' @export mutate
##' @rdname dplyr_verbs
# mutate.ped_cr_union <- function(.data, ..., keep_attributes = TRUE) {
#   if (keep_attributes) {
#     data_attr   <- ped_attr(.data)
#     data_class <- class(.data)
#   }
#   .data <- mutate(unped(.data), ...)
#   if (keep_attributes) {
#     attributes(.data) <- c(attributes(.data), data_attr)
#   }
#   class(.data) <- data_class

#   return(.data)

# }

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



#' @name tidyr_verbs
#' @title \code{tidyr} Verbs for \code{ped}-Objects
#' @param data an  object of class \code{ped}, see \code{\link{as_ped}}.
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
fill.ped <- function(
  data,
  ...,
  .direction = c("down", "up", "downup", "updown"),
  keep_attributes = TRUE) {

  data_class <- class(data)

  if (keep_attributes) {
    data_attr   <- ped_attr(data)
  }
  .direction = match.arg(.direction)
  tbl <- fill(unped(data), ..., .direction = .direction)
  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }
  class(tbl) <- data_class

  return(tbl)

}

#' @importFrom purrr discard
un_nested_df <- function(nested_fdf) {
  class(nested_fdf) <- class(nested_fdf) %>% discard(~.=="nested_fdf")
  nested_fdf
}
re_nested_df <- function(.data) {
  class(.data) <- c("nested_fdf", class(.data))
  .data
}

nested_fdf_attr <- function(nested_fdf) {
  attributes(nested_fdf)[c("breaks", "id_var", "intvars")]
}



#' @export
#' @export arrange
#' @rdname dplyr_verbs
#' @keywords internal
arrange.nested_fdf <- function(.data, ...) {

  classes_ped <- ped_classes(.data)
  attr_ped    <- attributes(.data)
  .data       <- arrange(unped(.data, classes_ped), ...)
  .data       <- reped(.data, classes_ped)

  re_attribute(.data, attr_ped)

}



#' @export
#' @export group_by
#' @rdname dplyr_verbs
group_by.nested_fdf <- ifelse(
  utils::packageVersion("dplyr") > "0.8.5",
  function(.data, ..., .add = FALSE) {
    re_nested_df(group_by(un_nested_df(.data), ..., .add = .add))
  },
  function(.data, ..., add = FALSE) {
    re_nested_df(group_by(un_nested_df(.data), ..., add = add))
  }
)

# #' @export
# #' @export ungroup
# #' @rdname dplyr_verbs
# ungroup.nested_fdf <- function(x, ...) {
#   re_nested_df(ungroup(un_nested_df(x), ...))
# }

#-------------------------------------------------------------------------------
# single table: row ops

#' @export
#' @export filter
#' @rdname dplyr_verbs
filter.nested_fdf <- function(.data, ...) {
  re_nested_df(filter(un_nested_df(.data), ...))
}

#' @export
#' @export sample_n
#' @inheritParams dplyr::sample_n
#' @rdname dplyr_verbs
sample_n.nested_fdf <- function(tbl, size, replace = FALSE, weight = NULL,
  .env = NULL, ...) {
  re_nested_df(sample_n(un_nested_df(tbl), size, replace, weight, .env, ...))
}

#' @export
#' @export sample_frac
#' @inheritParams dplyr::sample_frac
#' @rdname dplyr_verbs
sample_frac.nested_fdf <- function(tbl, size = 1, replace = FALSE, weight = NULL,
  .env = NULL, ...) {
  re_nested_df(sample_frac(un_nested_df(tbl), size, replace, weight, .env, ...))
}

#' @export
#' @export slice
#' @rdname dplyr_verbs
slice.nested_fdf <- function(.data, ...) {
  re_nested_df(slice(un_nested_df(.data), ...))
}

#-------------------------------------------------------------------------------
# single table: column ops

#' @export
#' @export select
#' @rdname dplyr_verbs
select.nested_fdf <- function(.data, ...) {
  re_nested_df(select(un_nested_df(.data), ...))
}


#' @param keep_attributes conserve attributes? defaults to \code{TRUE}
#' @export
#' @export mutate
#' @rdname dplyr_verbs
mutate.nested_fdf <- function(.data, ..., keep_attributes=TRUE) {
  if (keep_attributes) {
    data_attr   <- nested_fdf_attr(.data)
  }
  .data <- re_nested_df(mutate(un_nested_df(.data), ...))
  if (keep_attributes) {
    attributes(.data) <- c(attributes(.data), data_attr)
  }
  return(.data)
}

#' @export
#' @export rename
#' @rdname dplyr_verbs
rename.nested_fdf <- function(.data, ...) {
  re_nested_df(rename(un_nested_df(.data), ...))
}


#' @export
#' @export summarise
#' @rdname dplyr_verbs
summarise.nested_fdf <- function(.data, ...) {
  re_nested_df(summarise(un_nested_df(.data), ...))
}
#' @export
#' @rdname dplyr_verbs
summarize.nested_fdf <- summarise.nested_fdf

#' @export
#' @export transmute
#' @rdname dplyr_verbs
transmute.nested_fdf <- function(.data, ...) {
  re_nested_df(transmute(un_nested_df(.data), ...))
}

#-------------------------------------------------------------------------------
# joins

#' @inheritParams dplyr::inner_join
#' @export
#' @export inner_join
#' @rdname dplyr_verbs
inner_join.nested_fdf <- function(x, y, by = NULL, copy = FALSE,
  suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  re_nested_df(inner_join(un_nested_df(x), y, by, copy, suffix, ...))
}

#' @inheritParams dplyr::full_join
#' @export
#' @export full_join
#' @rdname dplyr_verbs
full_join.nested_fdf <- function(x, y, by = NULL, copy = FALSE,
  suffix = c(".x", ".y"),
  ...) {
  #FIXME?
  re_nested_df(full_join(un_nested_df(x), y, by, copy, suffix, ...))
}

#' @inheritParams dplyr::left_join
#' @export
#' @export left_join
#' @rdname dplyr_verbs
left_join.nested_fdf <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ..., keep_attributes=TRUE) {

  if (keep_attributes) {
    data_attr   <- nested_fdf_attr(x)
  }
  #FIXME?
  tbl <- re_nested_df(left_join(un_nested_df(x), y, by, copy, suffix, ...))

  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }

  return(tbl)

}

#' @inheritParams dplyr::right_join
#' @export
#' @export right_join
#' @rdname dplyr_verbs
right_join.nested_fdf <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),
  ..., keep_attributes=TRUE) {

  if (keep_attributes) {
    data_attr   <- nested_fdf_attr(y)
  }
  #FIXME?
  tbl <- re_nested_df(right_join(un_nested_df(x), y, by, copy, suffix, ...))

  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }

  return(tbl)
}


#' @inheritParams tidyr::fill
#' @inheritParams dplyr::filter
#' @export
#' @export fill
#' @rdname tidyr_verbs
#' @keywords internal
fill.nested_fdf <- function(
  data,
  ...,
  .direction = c("down", "up"),
  keep_attributes = TRUE) {
  if (keep_attributes) {
    data_attr   <- nested_fdf_attr(data)
  }
  tbl <- re_nested_df(fill(un_nested_df(data), ..., .direction = .direction))
  if (keep_attributes) {
    attributes(tbl) <- c(attributes(tbl), data_attr)
  }

  return(tbl)

}


# vectrs methods ----------------------------------------------------------

#' @importFrom vctrs vec_cast vec_ptype2
NULL

ped_ptype2 <- function(x, y, ...) {
  reped(vctrs::df_ptype2(x, y, ...))
}
ped_cast <- function(x, y, ...) {
  reped(vctrs::df_cast(x, y, ...))
}

#' @export
vec_ptype2.ped.ped <- function(x, y, ...)        ped_ptype2(x, y, ...)
#' @export
vec_ptype2.ped.data.frame <- function(x, y, ...) ped_ptype2(x, y, ...)
#' @export
vec_ptype2.data.frame.ped <- function(x, y, ...) ped_ptype2(x, y, ...)

#' @export
vec_cast.ped.ped <- function(x, y, ...)          ped_cast(x, y, ...)
#' @export
vec_cast.ped.data.frame <- function(x, y, ...)   ped_cast(x, y, ...)
#' @export
vec_cast.data.frame.ped <- function(x, y, ...)   vctrs::df_cast(x, y, ...)
