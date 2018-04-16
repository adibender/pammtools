#' Extract information of the sample contained in a data set
#'
#' Given a data set and grouping variables, this function returns mean values
#' for numeric variables and modus for characters and factors.
#'
#' @rdname sample_info
#' @param x A data frame (or object that inherits from \code{data.frame}).
#' @importFrom stats median
#' @return A data frame containing sample information (for each group).
#' If applied to an object of class \code{ped}, the sample means of the
#' original data is returned.
#' Note: When applied to a \code{ped} object, that doesn't contain covariates
#' (only interval information), returns data frame with 0 columns.
#' @export
sample_info <- function(x, ...) {
  UseMethod("sample_info", x)
}

#' @inheritParams sample_info
#' @import checkmate dplyr
#' @importFrom purrr compose
#' @export
#' @rdname sample_info
sample_info.data.frame <- function(x, ...) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 1, min.cols = 1)

  cn  <- colnames(x)
  num <- summarize_if(x, .predicate = is.numeric, ~mean(., na.rm = TRUE))
  fac <- summarize_if(x, .predicate = compose("!", is.numeric), modus)

  nnames <- intersect(names(num), names(fac))

  if(length(nnames) != 0) {
    x <- left_join(num, fac) %>% grouped_df(vars = lapply(nnames, as.name))
  } else {
    x <- bind_cols(num, fac)
  }

  return(select(x, one_of(cn)))

}

#' @rdname sample_info
#' @inheritParams sample_info
#' @import checkmate dplyr
#' @importFrom rlang sym
#' @export
sample_info.ped <- function(x) {
  # is.grouped_df
  # remove "noise" information on interval variables
  grps <- group_vars(x)
  iv <- attr(x, "intvars")
  id_var <- attr(x, "id_var")
  x <- x %>%
    group_by(!!sym(id_var)) %>%
    slice(1) %>%
    ungroup() %>%
    grouped_df(grps) %>%
    select(-one_of(iv))
  if (test_data_frame(x, min.rows = 1, min.cols = 1)) {
    sample_info.data.frame(x)
  } else {
    NULL
  }

}

#' @rdname sample_info
#' @inherit sample_info
#' @export
sample_info.fped <- function(x) {
  grps   <- group_vars(x)
  iv     <- attr(x, "intvars")
  id_var <- attr(x, "id_var")

  x %>% select_if(~!is.matrix(.x)) %>% sample_info.ped()

}


#' Combines multiple data frames
#'
#' Works like \code{\link[base]{expand.grid}} but for data frames.
#'
#' @importFrom dplyr slice bind_cols combine
#' @importFrom purrr map map_lgl map2 transpose cross
#' @importFrom checkmate test_data_frame
#' @param ... Data frames that should be combined to one data frame.
#' Elements of first df vary fastest, elements of last df vary slowest.
#' @examples
#' combine_df(
#'   data.frame(x=1:3, y=3:1),
#'   data.frame(x1=c("a", "b"), x2=c("c", "d")),
#'   data.frame(z=c(0, 1)))
#' @export
combine_df <- function(...) {

  dots <- list(...)
  if(!all(sapply(dots, test_data_frame))) {
    stop("All elements in ... must inherit from data.frame!")
  }
  ind_seq   <- map(dots, ~ seq_len(nrow(.x)))
  not_empty <- map_lgl(ind_seq, ~ length(.x) > 0)
  ind_list  <- ind_seq[not_empty] %>% cross() %>% transpose() %>% map(combine)

  map2(dots[not_empty], ind_list, function(.x, .y) slice(.x, .y)) %>% bind_cols()

}


#' Construct a data frame suitable for prediction
#'
#' Given a data set, returns a data frame type object that can be used
#' as \code{newdata} argument in a call to \code{predict} and similar functions.
#'
#' @rdname newdata
#' @aliases make_newdata
#' @inheritParams sample_info
#' @param ... Covariate specifications (expressions) that will be evaulated
#' by looking for variables in \code{x}. Must be of form \code{z = f(z)}
#' where \code{z} is a variable in the data set \code{x} and \code{f} a known
#' function that can be usefully applied to $\code{age}$. See examples below.
#' @import dplyr
#' @importFrom checkmate assert_data_frame assert_character
#' @importFrom purrr map cross_df
#' @details Depening on the class of \code{x}, mean or modus values will be
#' used for variables not specified in ellipsis. If x an object inherits class
#' \code{ped}, useful data set completion will be performed depending on variables
#' specified in ellipsis.
#' @examples
#' library(dplyr)
#' tumor %>% make_newdata()
#' tumor %>% make_newdata(age=c(50))
#' tumor %>% make_newdata(days=seq_range(days, 3), age=c(50, 55))
#' tumor %>% make_newdata(days=seq_range(days, 3), status=unique(status), age=c(50, 55))
#' # mean/modus values of unspecified variables are calculated over whole data
#' tumor %>% make_newdata(sex=unique(sex))
#' tumor %>% group_by(sex) %>% make_newdata()
#' # You can also pass complete data sets to make_newdata
#' purrr::cross_df(list(days = c(0, 500, 1000), sex = c("male", "female"))) %>%
#'   make_newdata(x=tumor)
#'
#' # Examples for PED data
#' ped <- tumor %>% slice(1:3) %>% as_ped(Surv(days, status)~., cut = c(0, 500, 1000))
#' ped %>% make_newdata(age=c(50, 55))
#' # if time information is specified, other time variables will be specified
#' # accordingly and offset calculated correctly
#' ped %>% make_newdata(tend = c(1000), age = c(50, 55))
#' ped %>% make_newdata(tend = unique(tend))
#' ped %>% group_by(sex) %>% make_newdata(tend = unique(tend))
#' @export
make_newdata <- function(x, ...) {
  UseMethod("make_newdata", x)
}


#' @inherit make_newdata
#' @rdname newdata
#' @param expand A character vector of column names in \code{ped}.
#' @param n If \code{expand} specified, respective variables will be expanded
#' in \code{n} values from minimum to maximum.
#' @export
make_newdata.default <- function(x, ...) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 2, min.cols = 1)

  orig_names <- names(x)

  expressions    <- quos(...)
  expr_evaluated <- map(expressions, lazyeval::f_eval, data=x)

  # construct data parts depending on input type
  lgl_atomic     <- map_lgl(expr_evaluated, is_atomic)
  partI  <- expr_evaluated[lgl_atomic] %>% cross_df()
  partII <- do.call(combine_df, expr_evaluated[!lgl_atomic])

  ndf    <- combine_df(partI, partII)
  rest   <- x %>% select(-one_of(c(colnames(ndf))))
  si     <- sample_info(rest) %>% ungroup()
  out_df <- combine_df(si, ndf)

  out_df %>% select(one_of(orig_names))

}

#' @rdname newdata
#' @inherit make_newdata.default
#' @export
make_newdata.ped <- function(x, ...) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 2, min.cols = 1)

  # prediction time points have to be interval end points so that piece-wise
  # constancy of predicted hazards is respected. If user overrides this, warn.

  orig_vars <- names(x)
  int_df <- int_info(x)

  expressions <- quos(...)
  dot_names   <- names(expressions)
  int_names   <- names(int_df)
  x <- select(x, -one_of(setdiff(int_names, c(dot_names, "intlen", "intmid"))))

  ndf <- make_newdata(unped(x), ...)

  if (any(names(int_df) %in% names(ndf))) {
    ndf <- right_join(int_df, ndf)
  } else {
    ndf <- combine_df(int_df[1,], ndf)
  }

  ndf %>% select(one_of(orig_vars)) %>%
    mutate(offset = log(tend - tstart))

}
