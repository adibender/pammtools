#' Transform data to Piece-wise Exponential Data (PED)
#'
#' This is the general data transformation function provided by the
#' \code{pammtools} package. Two main applications must be distinguished:
#' \enumerate{
#'  \item Transformation of standard time-to-event data.
#'  \item Transformation of time-to-event data with time-dependent covariates (TDC).
#' }
#' For the latter, the type of effect one wants to estimate is also
#' important for the data transformation step.
#' In any case, the data transformation is specified by a two sided formula.
#' In case of TDCs, the right-hand-side of the formula can contain formula specials
#' \code{concurrent} and \code{cumulative}.
#' See the \href{https://adibender.github.io/pammtools//articles/data-transformation.html}{data-transformation}
#' vignette for details.
#'
#'
#' @rdname as_ped
#' @param data Either an object inheriting from data frame or in case of
#' time-dependent covariates a list of data frames, where the first data frame
#' contains the time-to-event information and static covariates while the second
#' (and potentially further data frames) contain information on time-dependent
#' covariates and the times at which they have been observed.
#' @param formula A two sided formula with a \code{\link[survival]{Surv}} object
#' on the left-hand-side and covariate specification on the right-hand-side (RHS).
#' The RHS can be an extended formula, which specifies how TDCs should be transformed
#' using specials \code{concurrent} and \code{cumulative}.
#' @inheritParams survival::survSplit
#' @param cut Break points, used to partition the follow up into intervals.
#' If unspecified, all unique event times will be used.
#' @param max_time If \code{cut} is unspecified, this will be the last
#' possible event time. All event times after \code{max_time}
#' will be administratively censored at \code{max_time}.
#' @param ... Further arguments passed to the \code{data.frame} method and
#' eventually to \code{\link[survival]{survSplit}}
#' @importFrom Formula Formula
#' @examples
#' tumor[1:3, ]
#' tumor[1:3, ] %>% as_ped(Surv(days, status)~ age + sex, cut = c(0, 500, 1000))
#' tumor[1:3, ] %>% as_ped(Surv(days, status)~ age + sex)
#' @return A data frame class \code{ped} in piece-wise exponential data format.
#' @export
as_ped <- function(data, formula, ...) {
  UseMethod("as_ped", data)
}

#' @rdname as_ped
#' @inherit as_ped
#' @export
as_ped.data.frame <- function(
  data,
  formula,
  cut      = NULL,
  max_time = NULL,
  ...) {

  status_error(data, formula)

  dots          <- list(...)
  dots$data     <- data
  dots$formula  <- formula(Formula(formula), lhs = 1, rhs = 1)
  dots$cut      <- cut
  dots$max_time <- max_time
  ped <- do.call(split_data, dots)
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  ped

}

#' @rdname as_ped
#' @inherit as_ped
#' @export
as_ped.nested_fdf <- function(data, formula, ...) {

  status_error(data, formula)

  dots <- list(...)
  # update interval break points (if necessary)
  cut <- dots$cut
  if (is.null(cut)) {
    cut <- attr(data, "breaks")
  }
  ccr_breaks <- attr(data, "ccr_breaks")
  cut <- union(cut, ccr_breaks[ccr_breaks <= max(cut)]) %>% sort()

  ped <- data %>%
    select_if (is.atomic) %>%
    as_ped.data.frame(
      formula  = formula,
      id       = dots$id,
      cut      = cut,
      max_time = dots$max_time)

  # replace updated attributes
  attr(data, "breaks") <- attr(ped, "breaks")
  attr(data, "id_n") <- ped %>% group_by(!!sym(attr(data, "id_var"))) %>%
    summarize(id_n = n()) %>% pull("id_n") %>% as_vector()
  attr(data, "id_tseq") <- ped %>% group_by(!!sym(attr(data, "id_var"))) %>%
    transmute(id_tseq = row_number()) %>% pull("id_tseq") %>% as_vector()
  attr(data, "id_tz_seq") <- rep(seq_len(nrow(data)),
    times = attr(data, "id_n"))

  if (has_special(formula, "concurrent")) {
    ped <- ped %>% add_concurrent(data = data, id_var = dots$id)
  }

  if (has_special(formula, "cumulative")) {
    ped <- add_cumulative(ped, data = data, formula = formula)
    attr(ped, "ll_weights") <- imap(attr(ped, "tz"),
      ~bind_cols(!!.y := .x, ll_weight = c(mean(abs(diff(.x))), abs(diff(.x)))))
    class(ped) <- c("fped", class(ped))
  }
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  attr(ped, "func_mat_names") <- make_mat_names(
    attr(ped, "func"),
    attr(ped, "time_var"))

  ped

}

#' @rdname as_ped
#' @inherit as_ped
#' @export
as_ped.list <- function(data, formula, ...) {

  assert_class(data, "list")
  assert_class(formula, "formula")

  status_error(data[[1]], formula)

  nl    <- length(data)
  form  <- Formula(formula)
  n_rhs <- length(form)[2]

  if (nl == 1 & n_rhs == 1) {
    ped <- data[[1]] %>% as_ped(formula = form, ...)
  } else {
    if (nl == 2 & n_rhs == 1) {
    stop("Two data sets provided in 'data' but no specification of
      time-dependent covariate effects in 'formula'")
    } else {

      nested_fdf <- nest_tdc(data, form, ...)
      ped <- as_ped(nested_fdf, formula, ...)
    }
  }
  attr(ped, "time_var") <- get_lhs_vars(formula)[1]
  ped

}

#' @rdname as_ped
#' @param x any R object.
#' @export
is.ped <- function(x) inherits(x, "ped")
