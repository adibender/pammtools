#' Extract information of the sample contained in a data set
#'
#' Given a data set and grouping variables, this function returns mean values
#' for numeric variables and modus for characters and factors. Usually
#' this function should not be called directly but will rather be called
#' as part of a call to \code{make_newdata}.
#'
#' @rdname sample_info
#' @param x A data frame (or object that inherits from \code{data.frame}).
#' @importFrom stats median
#' @return A data frame containing sample information (for each group).
#' If applied to an object of class \code{ped}, the sample means of the
#' original data is returned.
#' Note: When applied to a \code{ped} object, that doesn't contain covariates
#' (only interval information), returns data frame with 0 columns.
#'
#' @export
#' @keywords internal
sample_info <- function(x) {
  UseMethod("sample_info", x)
}

#' @import checkmate dplyr
#' @importFrom purrr compose
#' @export
#' @rdname sample_info
sample_info.data.frame <- function(x) {

  cn  <- colnames(x)
  num <- summarize_if (x, .predicate = is.numeric, ~mean(., na.rm = TRUE))
  fac <- summarize_if (x, .predicate = compose("!", is.numeric), modus)

  nnames <- intersect(names(num), names(fac))

  if (length(nnames) != 0) {
    suppressMessages(
      x <- left_join(num, fac) %>% group_by(!!!lapply(nnames, as.name))
    )
  } else {
    x <- bind_cols(num, fac)
  }

  return(select(x, one_of(cn)))

}

#' @rdname sample_info
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
#' @export
sample_info.fped <- function(x) {

  x %>% select_if(~!is.matrix(.x)) %>% sample_info.ped()

}


#' Create a data frame from all combinations of data frames
#'
#' Works like \code{\link[base]{expand.grid}} but for data frames.
#'
#' @importFrom dplyr slice bind_cols
#' @importFrom vctrs vec_c
#' @importFrom purrr map map_lgl
#' @importFrom checkmate test_data_frame
#' @param ... Data frames that should be combined to one data frame.
#' Elements of first df vary fastest, elements of last df vary slowest.
#' @examples
#' combine_df(
#'   data.frame(x=1:3, y=3:1),
#'   data.frame(x1=c("a", "b"), x2=c("c", "d")),
#'   data.frame(z=c(0, 1)))
#' @export
#' @keywords internal
combine_df <- function(...) {

  dots <- list(...)
  if (!all(sapply(dots, test_data_frame))) {
    stop("All elements in ... must inherit from data.frame!")
  }
  ind_seq   <- map(dots, ~ seq_len(nrow(.x)))
  not_empty <- map_lgl(ind_seq, ~ length(.x) > 0)

  ord <- lapply(dots[not_empty], function(z) colnames(z)) |> unlist()
  out <- do.call(expand_grid, rev(dots[not_empty]))
  out <- out[, ord]

}


#' Construct a data frame suitable for prediction
#'
#' This functions provides a flexible interface to create a data set that
#' can be plugged in as \code{newdata} argument to a suitable  \code{predict}
#' function (or similar).
#' The function is particularly useful in combination with one of the
#' \code{add_*} functions, e.g., \code{\link[pammtools]{add_term}},
#' \code{\link[pammtools]{add_hazard}}, etc.
#'
#' @rdname newdata
#' @aliases make_newdata
#' @importFrom tidyr expand_grid
#' @inheritParams sample_info
#' @param ... Covariate specifications (expressions) that will be evaluated
#' by looking for variables in \code{x}. Must be of the form \code{z = f(z)}
#' where \code{z} is a variable in the data set and \code{f} a known
#' function that can be usefully applied to \code{z}. Note that this is also
#' necessary for single value specifications (e.g. \code{age = c(50)}).
#' For data in PED (piece-wise exponential data) format, one can also specify
#' the time argument, but see "Details" an "Examples" below.
#' @import dplyr
#' @importFrom checkmate assert_data_frame assert_character
#' @importFrom purrr map
#' @importFrom tidyr expand_grid
#' @details Depending on the type of variables in \code{x}, mean or modus values
#' will be used for variables not specified in ellipsis
#' (see also \code{\link[pammtools]{sample_info}}). If \code{x} is an object
#' that inherits from class \code{ped}, useful data set completion will be
#' attempted depending on variables specified in ellipsis. This is especially
#' useful, when creating a data set with different time points, e.g. to
#' calculate survival probabilities over time (\code{\link[pammtools]{add_surv_prob}})
#' or to calculate a time-varying covariate effects (\code{\link[pammtools]{add_term}}).
#' To do so, the time variable has to be specified in \code{...}, e.g.,
#' \code{tend = seq_range(tend, 20)}. The problem with this specification is that
#' not all values produced by \code{seq_range(tend, 20)} will be actual values
#' of \code{tend} used at the stage of estimation (and in general, it will
#' often be tedious to specify exact \code{tend} values). \code{make_newdata}
#' therefore finds the correct interval and sets \code{tend} to the respective
#' interval endpoint. For example, if the intervals of the PED object are
#' \eqn{(0,1], (1,2]} then \code{tend = 1.5} will be set to \code{2} and the
#' remaining time-varying information (e.g. offset) completed accordingly.
#' See examples below.
#' @examples
#' # General functionality
#' tumor %>% make_newdata()
#' tumor %>% make_newdata(age=c(50))
#' tumor %>% make_newdata(days=seq_range(days, 3), age=c(50, 55))
#' tumor %>% make_newdata(days=seq_range(days, 3), status=unique(status), age=c(50, 55))
#' # mean/modus values of unspecified variables are calculated over whole data
#' tumor %>% make_newdata(sex=unique(sex))
#' tumor %>% group_by(sex) %>% make_newdata()
#'
#' # Examples for PED data
#' ped <- tumor %>% slice(1:3) %>% as_ped(Surv(days, status)~., cut = c(0, 500, 1000))
#' ped %>% make_newdata(age=c(50, 55))
#'
#' # if time information is specified, other time variables will be specified
#' # accordingly and offset calculated correctly
#' ped %>% make_newdata(tend = c(1000), age = c(50, 55))
#' ped %>% make_newdata(tend = unique(tend))
#' ped %>% group_by(sex) %>% make_newdata(tend = unique(tend))
#'
#' # tend is set to the end point of respective interval:
#' ped <- tumor %>% as_ped(Surv(days, status)~.)
#' seq_range(ped$tend, 3)
#' make_newdata(ped, tend = seq_range(tend, 3))
#' @export
make_newdata <- function(x, ...) {
  UseMethod("make_newdata", x)
}


#' @inherit make_newdata
#' @rdname newdata
#' @export
make_newdata.default <- function(x, ...) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 2, min.cols = 1)

  orig_names <- names(x)

  expressions    <- quos(...)
  expr_evaluated <- map(expressions, lazyeval::f_eval, data = x) |>
    map(c)

  # construct data parts depending on input type
  lgl_atomic <- map_lgl(expr_evaluated, is_atomic)
  # part1 <- expr_evaluated[lgl_atomic] |> cross_df()
  part1 <- do.call(tidyr::expand_grid, rev(expr_evaluated[lgl_atomic]))
  part2 <- do.call(combine_df, expr_evaluated[!lgl_atomic])

  ndf  <- combine_df(part1, part2)
  rest <- x %>% select(-one_of(c(colnames(ndf))))
  if (ncol(rest) > 0) {
    si  <- sample_info(rest) %>% ungroup()
    ndf <- combine_df(si, ndf)

  }

  ndf %>% select(one_of(orig_names))

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
  # x <- select(x, -one_of(setdiff(int_names, c(dot_names, "intlen", "intmid"))))
  ndf <- x %>%
    select(-one_of(setdiff(int_names, c(dot_names, "intlen", "intmid")))) %>%
    unped() %>%
    make_newdata(...)

  if (any(names(int_df) %in% names(ndf))) {
    int_tend <- get_intervals(x, ndf$tend)$tend
    if (!all(ndf$tend == int_tend)) {
      message("Some values of 'tend' have been set to the respective interval end-points")
    }
    ndf$tend <- int_tend
    suppressMessages(
      ndf <- ndf %>% left_join(int_df)
      )
  } else {
    ndf <- combine_df(int_df[1, ], ndf)
  }

  int_names <- intersect(int_names, c("intlen", orig_vars))
  ndf %>% select(one_of(c(int_names, setdiff(orig_vars, int_names)))) %>%
    mutate(
      intlen = .data$tend - .data$tstart,
      offset = log(.data$tend - .data$tstart),
      ped_status = 0)

}

#' @rdname newdata
#' @inherit make_newdata.ped
#' @importFrom rlang quos
#' @export
make_newdata.fped <- function(x, ...) {

  assert_data_frame(x, all.missing = FALSE, min.rows = 2, min.cols = 1)

  # prediction time points have to be interval end points so that piece-wise
  # constancy of predicted hazards is respected. If user overrides this, warn.
  expressions <- quos(...)
  dot_names   <- names(expressions)
  orig_vars   <- names(x)
  cumu_vars   <- setdiff(unlist(attr(x, "func_mat_names")), dot_names)
  cumu_smry   <- smry_cumu_vars(x, attr(x, "time_var")) %>%
    select(one_of(cumu_vars))

  int_names   <- attr(x, "intvars")
  ndf <- x %>%
    select(one_of(setdiff(names(x), cumu_vars))) %>%
    unfped() %>%
    make_newdata(...)

  out_df <- do.call(combine_df, compact(list(cumu_smry, ndf)))
  int_df <- int_info(attr(x, "breaks"))
  suppressMessages(
    out_df <- out_df %>% left_join(int_df) %>%
      select(-one_of(c("intmid"))) %>% as_tibble()
      )

  # adjust lag-lead indicator
  out_df <- adjust_ll(out_df, x)

  out_df

}


smry_cumu_vars <- function(data, time_var) {

  cumu_vars <- unlist(attr(data, "func_mat_names"))
  func_list <- attr(data, "func")
  z_vars    <- map(func_list, ~get_zvars(.x, time_var, length(func_list))) %>%
    unlist()
  smry_z <- select(data, one_of(z_vars)) %>%
    map(~ .x[1, ]) %>% map(~mean(unlist(.x))) %>% bind_cols()
  smry_time <- select(data, setdiff(cumu_vars, z_vars)) %>% map(~.x[1, 1])

  bind_cols(smry_z, smry_time)

}

get_zvars <- function(func, time_var, n_func) {

  col_vars <- func$col_vars
  all_vars <- make_mat_names(c(col_vars, "LL"), func$latency_var, func$tz_var,
    func$suffix, n_func)
  time_vars <- make_mat_names(c(time_var, func$tz_var, "LL"),
    func$latency_var, func$tz_var, func$suffix, n_func)

  setdiff(all_vars, time_vars)

}


## apply ll_fun to newly created data
adjust_ll <- function(out_df, data) {

  func_list <- attr(data, "func")
  n_func    <- length(func_list)
  LL_names <- grep("LL", unlist(attr(data, "func_mat_names")), value = TRUE)

  for (i in LL_names) {
    ind_ll <- map_lgl(names(attr(data, "ll_funs")), ~grepl(.x, i))
    if (any(ind_ll)) {
      ind_ll <- which(ind_ll)
    } else {
      ind_ll <- 1
    }

    func   <- func_list[[ind_ll]]
    ll_i   <- attr(data, "ll_funs")[[ind_ll]]
    tz_var <- attr(data, "tz_vars")[[ind_ll]]
    tz_var <- make_mat_names(tz_var, func$latency_var, func$tz_var, func$suffix,
      n_func)
    if (func$latency_var == "") {
      out_df[[i]] <- ll_i(out_df[["tend"]], out_df[[tz_var]]) * 1L
    } else {
      out_df[[i]] <- ll_i(out_df[["tend"]], out_df[["tend"]] -
        out_df[[tz_var]]) * 1L
    }
  }

  out_df

}

# All variables that represent follow-up time should have the same values
# adjust_time_vars <- function(out_df, data, dot_names) {

#   time_vars <- c("tend",
#     grep(attr(data, "time_var"), unlist(attr(data, "func_mat_names")), value=TRUE))
#   time_vars_dots <- c(grep("tend", dot_names, value=TRUE),
#     grep(attr(data, "time_var"), dot_names, value=TRUE))
#   if (length(time_vars_dots) == 0) {
#     time_vars_dots <- "tend"
#   } else {
#     if (length(time_vars_dots) > 1) {
#       warning(paste0("Only one of ", paste0(time_vars_dots, collapse=", "),
#         "must be specified. Only the first one will be used!"))
#       time_vars_dots <- time_vars_dots[1]
#     }
#   }
#   for (i in setdiff(time_vars, time_vars_dots)) {
#     out_df[[i]] <- out_df[[time_vars_dots]]
#   }

#   out_df

# }
