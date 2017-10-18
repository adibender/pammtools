#' Add predicted hazard to data set
#'
#' @inheritParams mgcv::predict.gam
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}} and
#'   \code{\link{get_hazard}}
#' @param ci Logical indicating whether to iclude confidence intervals. Defaults
#' to \code{TRUE}
#' @param se.mult Factor by which standard errors are multiplied for calculating
#' the confidence intervals.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to \code{FALSE}. If \code{TRUE}, columns with names
#' \code{c("hazard", "se", "lower", "upper")} will be overwritten.
#' @param time_variable Name of the variable used for the baseline hazard. If
#'   not given, defaults to \code{"tend"} for \code{\link[mgcv]{gam}} fits, else
#'   \code{"interval"}. The latter is assumed to be a factor, the former
#'   numeric.
#' @import checkmate dplyr mgcv
#' @importFrom stats predict
#' @examples
#' \dontrun{
#' library(mgcv)
#' data("veteran", package="survival")
#' ped <- split_data(Surv(time, status)~. cut=seq(0, 500, by=100), data=veteran,
#'  id="id")
#' pam <- gam(ped_status ~ s(tend, k=5), data = ped, family=poisson(), offset=offset)
#' pinfo <- ped_info(ped)
#' add_hazard(pinfo, pam)
#' }
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pammtools]{add_cumu_hazard}}
#' @rdname add_hazard
#' @export
add_hazard <- function(
  newdata,
  object,
  type          = c("response", "link"),
  ci            = TRUE,
  se.mult       = 2,
  overwrite     = FALSE,
  time_variable = NULL,
  ...)  {

  if (!overwrite) {
    if ("hazard" %in% names(newdata)) {
      stop("Data set already contains 'hazard' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
      rm.vars <- intersect(
        c("hazard", "se", "ci_lower", "ci_upper"),
        names(newdata))
      newdata <- newdata %>% select(-one_of(rm.vars))
  }

  pred <- get_hazard(newdata, object, ci = ci, type = type, se.mult = se.mult,
    time_variable = time_variable, ...)
  stopifnot(nrow(pred) == nrow(newdata))

  newdata <- newdata %>% bind_cols(rm_grpvars(pred))

  return(newdata)

}

#' Calculate predicted hazard
#'
#' @inheritParams add_hazard
#' @importFrom stats model.frame
#' @keywords internal
get_hazard <- function(
  newdata,
  object,
  ci            = TRUE,
  type          = c("response", "link"),
  se.mult       = 2,
  time_variable = NULL,
  ...)  {

  assert_data_frame(newdata, all.missing = FALSE)
  assert_class(object, classes = "glm")
  type <- match.arg(type)

  is_pam <- inherits(object, "gam")
  if (is.null(time_variable)) {
    time_variable <- ifelse(is_pam, "tend", "interval")
  } else {
    assert_string(time_variable)
    assert_choice(time_variable, colnames(newdata))
  }

  original_intervals <- if (is_pam) {
    unique(model.frame(object)[[time_variable]])
  } else levels(model.frame(object)[[time_variable]])
  prediction_intervals <- if (is_pam) {
    unique(newdata[[time_variable]])
  } else levels(factor(newdata[[time_variable]]))
  new_ints <- which(!(prediction_intervals %in% original_intervals))
  if (length(new_ints)) {
   message <- paste0("Intervals in <newdata> contain values (",
     paste(prediction_intervals[new_ints], collapse = ","),
     ") not used in original fit.",
     " Setting intervals to values not used for original fit in <object>",
     "can invalidate the PEM assumption and yield incorrect predictions.")
   if (is_pam) warning(message) else stop(message)
  }

  pred <- predict(
    object  = object,
    newdata = newdata,
    se.fit  = TRUE,
    type    = "link", ...)
  pred <- pred[c("fit", "se.fit")] %>%
    bind_cols() %>%
    rename(
      hazard = "fit",
      se     = "se.fit") %>%
    mutate(
      hazard = as.numeric(hazard),
      se     = as.numeric(se))
  stopifnot(nrow(pred) == nrow(newdata))

  mutate_vars <- c("hazard")
  if (ci) {
    pred <- pred %>%
      mutate(
        ci_lower = hazard - se.mult * se,
        ci_upper = hazard + se.mult * se)
    mutate_vars <- c(mutate_vars, "ci_lower", "ci_upper")
  }

  if (type == "response") {
    pred <- pred %>%
      mutate_at(mutate_vars, funs(exp(.)))
  }

  # it is necessary to include the grouping variables here, otherwise
  # functions calculating the cumulative hazard will cumulate over all rows
  # instead of group wise
  if (is.grouped_df(newdata)) {
    group.df <- newdata %>% select(one_of(unlist(groups(newdata))))
    pred     <- bind_cols(group.df, pred)
  }

  return(pred)

}

#' Add cumulative hazard estimate to data set
#'
#' @rdname add_hazard
#' @inheritParams add_hazard
#' @param interval_length \code{quosure} providing the name of the variable in
#'  newdata containing the interval lengths. Defaults to \code{intlen}.
#' @importFrom dplyr bind_cols
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pammtools]{add_hazard}}
#' @export
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pammtools]{add_hazard}}
add_cumu_hazard <- function(
  newdata,
  object,
  ci              = TRUE,
  se.mult         = 2,
  overwrite       = FALSE,
  time_variable   = NULL,
  interval_length = quo(intlen),
  ...)  {

  if (!overwrite) {
    if ("cumu_hazard" %in% names(newdata)) {
      stop(
        "Data set already contains 'hazard' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
      rm.vars <- intersect(c("cumu_hazard", "cumu_lower", "cumu_upper"),
        names(newdata))
      newdata <- newdata %>% select(-one_of(rm.vars))
  }

  pred <- get_cumu_hazard(newdata, object, ci = ci, se.mult = se.mult,
    time_variable = time_variable, interval_length = interval_length, ...)

  newdata <- newdata %>% bind_cols(rm_grpvars(pred))

  return(newdata)

}

#' Calculate cumulative hazard
#'
#' @inheritParams add_cumu_hazard
#' @import checkmate dplyr
#' @importFrom rlang quo UQ UQS
#' @keywords internal
get_cumu_hazard <- function(
  newdata,
  object,
  ci              = TRUE,
  time_variable   = NULL,
  interval_length = quo(intlen),
  ...) {

  assert_class(interval_length, "quosure")
  assert_choice(as.character(interval_length)[2], colnames(newdata))

  lengths <- select(rm_grpvars(newdata), !!interval_length)

  mutate_args  <- list(cumu_hazard = quo(cumsum(hazard * (!!interval_length))))
  vars_exclude <- c("hazard", "se")
  if (ci) {
    mutate_args <- mutate_args %>%
      append(
        list(
          cumu_lower = quo(cumsum(ci_lower * (!!interval_length))),
          cumu_upper = quo(cumsum(ci_upper * (!!interval_length)))))
    vars_exclude <- c(vars_exclude, "ci_lower", "ci_upper")
  }
  get_hazard(newdata, object, ci = ci, type = "response",
    time_variable = time_variable, ...) %>%
    bind_cols(lengths) %>%
    mutate(!!!mutate_args) %>%
    select(-one_of(vars_exclude)) %>%
    select(-!!interval_length)
}


#' Add survival probability estimates to data set
#'
#' @inherit add_cumu_hazard
#' @export
add_surv_prob <- function(
  newdata,
  object,
  ci              = TRUE,
  se.mult         = 2,
  overwrite       = FALSE,
  time_variable   = NULL,
  interval_length = quo(intlen),
  ...)  {

  if (!overwrite) {
    if ("surv_prob" %in% names(newdata)) {
      stop("Data set already contains 'surv_prob' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
      rm.vars <- intersect(
        c("surv_prob", "surv_lower", "surv_upper"),
        names(newdata))
      newdata <- newdata %>% select(-one_of(rm.vars))
  }

  pred <- get_surv_prob(newdata, object, ci = ci, se.mult = se.mult,
    time_variable = time_variable, interval_length = interval_length, ...)

  newdata %>% bind_cols(rm_grpvars(pred))

}


#' Calculate survival probabilities
#'
#' @inheritParams add_surv_prob
#' @keywords internal
get_surv_prob <- function(
  newdata,
  object,
  ci              = TRUE,
  time_variable   = NULL,
  interval_length = "intlen",
  ...) {

  assert_class(interval_length, "quosure")
  assert_choice(as.character(interval_length)[2], colnames(newdata))

  lengths      <- select(rm_grpvars(newdata), !!interval_length)
  exclude_vars <- c("cumu_hazard")
  mutate_args  <- list(surv_prob = quo(exp(-cumu_hazard)))
  if (ci) {
    mutate_args <- mutate_args %>%
      append(
        list(
          surv_lower = quo(exp(-cumu_lower)),
          surv_upper = quo(exp(-cumu_upper))))
    exclude_vars <- c(exclude_vars, "cumu_lower", "cumu_upper")
  }

  newdata %>%
    get_cumu_hazard(object, ci = ci, time_variable = time_variable, ...) %>%
    mutate(!!!mutate_args) %>%
    select(-one_of(exclude_vars))

}
