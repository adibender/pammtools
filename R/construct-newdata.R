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
  num <- summarize_if(x, .predicate=is.numeric, funs(median(., na.rm=TRUE)))
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
#' Given a data set of class \code{ped}, returns a data frame that can be used 
#' as \code{newdata} argument in a call to \code{predict} and similar functions. 
#' When \code{expand==NULL} falls back to \code{\link{ped_info}}.
#' 
#' @inheritParams ped_info
#' @inheritParams seq_range
#' @param ... Further specifications of variables that should be set 
#' to a specific value. 
#' @param expand A charachter vector of column names in \code{ped}.
#' @import dplyr 
#' @importFrom checkmate assert_data_frame assert_character
#' @importFrom purrr map cross_df
#' @details Details 
#' Extracts information from \code{ped}, using \code{\link{sample_info}}. 
#' If variables are specified with specific values in \code{...}, the values
#' in from \code{sample_info} will be overwritten. If variables provided 
#' in expand, these will be exanded from \code{min} to \code{max} using 
#' in \code{length.out} equidistant steps.
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
make_newdata <- function(
  ped, 
  ..., 
  expand     = NULL,
  length.out = 25L) {

  assert_data_frame(ped, all.missing = FALSE, min.rows = 2, min.cols = 1)
  assert_character(expand, min.chars = 1, any.missing = FALSE, null.ok = TRUE)

  orig_names <- names(ped)

  si      <- sample_info(ped) %>% ungroup()
  dots_df <- cross_df(list(...))
  assert_subset(names(dots_df), orig_names)
  # return(dots_df)
  if (!is.null(expand)) {
    if (!all(expand %in% names(ped))) {
      stop("All arguments provided in 'expand' must be quoted variable names 
      that are  equal to column names of ped object")
    } else {
      expanded_df <- ped %>% 
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