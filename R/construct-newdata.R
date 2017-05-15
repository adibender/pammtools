#' Extract information of the sample contained in a data set
#'
#' Given a data set and grouping variables, this function returns median values
#' for numeric variables and modus for characters and factors.
#'
#' @param x A data frame (or object that inherits from \code{data.frame}).
#' @param ... Further arguments passed to specialized methods.
#' @importFrom stats median
#' @export
#' @rdname sample_info
sample_info <- function(x, ...) {
  UseMethod("sample_info", x)
}

#' @inheritParams sample_info
#' @import checkmate dplyr
#' @importFrom magrittr %<>%
#' @importFrom purrr compose
#' @export
#' @rdname sample_info
sample_info.data.frame <- function(x, ...) {

  assert_data_frame(x, all.missing=FALSE, min.rows=1, min.cols=1)

  cn <- colnames(x)
  num <- summarize_if(x, .predicate=is.numeric, funs(mean(., na.rm=TRUE)))
  fac <- summarize_if(x, .predicate=compose("!", is.numeric), modus)

  nnames <- intersect(names(num), names(fac))

  if(length(nnames) != 0) {
    x <- left_join(num, fac) %>% grouped_df(vars=lapply(nnames, as.name))
  } else {
    x <- bind_cols(num, fac)
  }

  return(select(x, one_of(cn)))

}


#' @inheritParams sample_info
#' @import checkmate dplyr
#' @importFrom magrittr %<>%
#' @export
#' @rdname sample_info
#' @seealso \code{\link[pam]{split_data}}
sample_info.ped <- function(x, ...) {
  # is.grouped_df
  # remove "noise" information on interval variables
  iv <- attr(x, "intvars")
  x %<>% select(-one_of(iv))
  sample_info.data.frame(x, ...)

}


#' Creates sequence from minimum to maximum
#'
#' @param x A numeric or integer vector.
#' @inheritParams base::seq
#' @import checkmate
#' @export
seq_range <- function(x, length.out=100L) {

  assert_numeric(x, finite=TRUE, all.missing=FALSE, min.len=2)
  assert_number(length.out, lower=2, finite=TRUE)

  range.x <- range(x)
  seq(range.x[1], range.x[2], length.out=length.out)

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
#' @inheritParams sample_info
#' @inheritParams seq_range
#' @param ... Further specifications of variables that should be set
#' to a specific value.
#' @param expand A character vector of column names in \code{ped}.
#' @import dplyr
#' @importFrom checkmate assert_data_frame assert_character
#' @importFrom purrr map cross_df
#' @details Extracts information from \code{ped}, using
#'   \code{\link{sample_info}}. If variables are specified with specific values
#'   in \code{...}, the values in from \code{sample_info} will be overwritten.
#'   If variables are provided in \code{expand}, these will be expanded from \code{min} to
#'   \code{max} using in \code{length.out} equidistant steps.
#' @examples
#' \dontrun{
#' library(dplyr)
#' iris %>% make_newdata()
#' iris %>% make_newdata(Sepal.Length=5)
#' iris %>% make_newdata(Sepal.Length=c(5, 10), Sepal.Width=c(5, 5.1, 5.2))
#' iris %>% make_newdata(Sepal.Length=c(5, 10), expand="Sepal.Width", length.out=5)
#' iris %>% group_by(Species) %>%
#'   make_newdata(Sepal.Length=c(5, 10), expand="Sepal.Width", length.out=5) %>%
#'   print(n=30)
#' }
#' @export
make_newdata <- function(x, ..., expand=NULL, length.out=50L) {
  UseMethod("make_newdata", x)
}


#' @inherit make_newdata
#' @importFrom magrittr "%<>%"
#' @export
make_newdata.default <- function(
  x,
  ...,
  expand     = NULL,
  length.out = 50L) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 2, min.cols = 1)
  assert_character(expand, min.chars = 1, any.missing = FALSE, null.ok = TRUE)

  orig_names <- names(x)

  si      <- sample_info(x) %>% ungroup()
  dots_df <- cross_df(list(...))
  assert_subset(names(dots_df), orig_names)
  # return(dots_df)
  if (!is.null(expand)) {
    if (!all(expand %in% names(x))) {
      stop("All arguments provided in 'expand' must be quoted variable names
      that are  equal to column names of x object")
    } else {
      expanded_df <- x %>%
        ungroup() %>% # ungroup here to obtain sequence from min to max for all data
        select(one_of(expand)) %>% as.list() %>%
        map(seq_range, length.out=length.out) %>%
        cross_df()
    }
  } else {
    expanded_df <- data.frame()
  }

  si_names    <- intersect(names(si), c(names(dots_df), names(expanded_df)))
  si %<>% select(-one_of(si_names))

  combine_df(si, dots_df, expanded_df) %>%
    select(one_of(orig_names))

}

#' @inherit make_newdata
#' @importFrom magrittr "%<>%"
#' @export
make_newdata.ped <- function(
  x,
  ...,
  expand     = NULL,
  length.out = 50L) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 2, min.cols = 1)
  assert_character(expand, min.chars = 1, any.missing = FALSE, null.ok = TRUE)

  # prediction time points have to be interval end points so that piece-wise
  # constancy of predicted hazards is respected. If user overrides this, warn.
  user_override_dots <- any(c("tstart", "tend") %in% names(list(...)))
  user_override_expand <- any(c("tstart", "tend") %in% expand)
  if (user_override_dots | user_override_expand) {
    warning("Setting interval start or end points to values not used for",
      " original 'ped'-data can invalidate PEM assumption and yield incorrect",
      " predictions. Proceed with caution!")
  }

  g_vars <- group_vars(x)
  x %<>% ungroup() %>% group_by(id) %>% slice(1) %>% ungroup(id) %>% unped() %>%
    group_by_(.dots = g_vars)

  make_newdata(x, ..., expand = expand, length.out = length.out)
}
