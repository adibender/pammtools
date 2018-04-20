#' Add info about term effects to data set
#'
#' Adds the contribution (plus confidence intervals) of a specific term to the
#' linear predictor to the provided data.
#' Largely a wrapper to \code{\link[mgcv]{predict.gam}}, with \code{type="terms"}.
#' Thus most arguments and their documentation below is from \code{predict.gam}.
#'
#'
#' @inheritParams mgcv::predict.gam
#' @param term A character (vector) or regular expression indicating for
#' which term(s) information should be extracted and added to data set.
#' @param se.mult The factor by which standard errors are multiplied to form
#' confidence intervals.
#' @param relative If \code{TRUE}, calculates relative risk contribution,
#' that is \eqn{(X-\bar{X})'\beta} and respective confidence intervals
#' if \code{se.fit = TRUE}. Defaults to \code{FALSE}.
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}}
#' @import checkmate dplyr mgcv
#' @importFrom stats predict
#' @importFrom purrr map
#' @examples
#' ped <- tumor[1:50,] %>% as_ped(Surv(days, status)~ age)
#' pam <- mgcv::gam(ped_status ~ s(tend)+age, data = ped, family=poisson(), offset=offset)
#' ped_info(ped) %>% add_term(pam, term="tend")
#' @export
#' @importFrom stats model.matrix vcov
add_term <- function(
  newdata,
  object,
  term,
  se.fit   = TRUE,
  type     = "terms",
  se.mult  = 2,
  relative = FALSE,
  ...) {

  assert_data_frame(newdata, all.missing=FALSE)
  assert_character(term, min.chars=1, any.missing=FALSE, min.len=1)

  col_ind <- map(term, grep, x=names(object$coefficients)) %>%
    unlist() %>% unique() %>% sort()
  is_pam <- inherits(object, "gam")

  X <- if (is_pam) {
    predict(object, newdata = newdata, type = "lpmatrix", ...)[, col_ind, drop=FALSE]
  } else  {
    model.matrix(object$formula[-2], data = newdata)[,col_ind, drop=FALSE]
  }
  if (relative) {
    if (is.null(object$model)) {
      stop("Relative risk can only be calculated when original data is present in 'object'!")
    }
    data_bar <- sample_info(object$model)[rep(1, nrow(newdata)), ]
    X_bar <- if (is_pam) {
      predict(
        object,
        newdata = data_bar,
        type    = "lpmatrix")[, col_ind, drop = FALSE]
    } else {
      model.matrix(object$formula[-2], data = data_bar)[, col_ind, drop = FALSE]
    }
    X <- X - X_bar
  }

  newdata[["fit"]] <- drop(X %*% object$coefficients[col_ind])
  if (se.fit) {
    cov.coefs <- if(is_pam) {
      object$Vp[col_ind, col_ind]
    } else {
      vcov(object)[col_ind, col_ind]
    }
    se <- sqrt(drop(diag(X %*% cov.coefs %*% t(X))))
    newdata <- newdata %>%
      mutate(
        ci_lower = .data$fit - se.mult * se,
        ci_upper = .data$fit + se.mult * se)
  }

  return(newdata)

}

add_term2 <- function(
  newdata,
  object,
  term,
  se.fit   = TRUE,
  type     = "terms",
  se.mult  = 2,
  reference = NULL,
  ...) {

  assert_data_frame(newdata, all.missing=FALSE)
  assert_character(term, min.chars=1, any.missing=FALSE, min.len=1)

  col_ind <- map(term, grep, x=names(object$coefficients)) %>%
    unlist() %>% unique() %>% sort()
  is_pam <- inherits(object, "gam")

  X <- if (is_pam) {
    predict(object, newdata = newdata, type = "lpmatrix", ...)[, col_ind, drop=FALSE]
  } else  {
    model.matrix(object$formula[-2], data = newdata)[,col_ind, drop=FALSE]
  }
  if (!is.null(reference)) {
    reference <- newdata %>% mutate(!!!reference)
    X_ref <- if (is_pam) {
      predict(
        object,
        newdata = reference,
        type    = "lpmatrix")[, col_ind, drop = FALSE]
    } else {
      model.matrix(object$formula[-2], data = reference)[, col_ind, drop = FALSE]
    }
    X <- X - X_ref
  }

  newdata[["fit"]] <- drop(X %*% object$coefficients[col_ind])
  if (se.fit) {
    cov.coefs <- if(is_pam) {
      object$Vp[col_ind, col_ind]
    } else {
      vcov(object)[col_ind, col_ind]
    }
    se <- sqrt(drop(diag(X %*% cov.coefs %*% t(X))))
    newdata <- newdata %>%
      mutate(
        ci_lower = .data$fit - se.mult * se,
        ci_upper = .data$fit + se.mult * se)
  }

  return(newdata)

}


#' Add predicted (cumulative) hazard to data set
#'
#' Add (cumulative) hazard based on the provided data set and model.
#' If \code{ci=TRUE} confidence intervals are also provided. Their width can
#' be controlled via the \code{se.mult} argument.
#'
#' @rdname add_hazard
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
#' ped <- tumor[1:50,] %>% as_ped(Surv(days, status)~ age)
#' pam <- mgcv::gam(ped_status ~ s(tend)+age, data = ped, family=poisson(), offset=offset)
#' ped_info(ped) %>% add_hazard(pam, type="link")
#' ped_info(ped) %>% add_hazard(pam, type="response")
#' ped_info(ped) %>% add_cumu_hazard(pam)
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
      hazard = as.numeric(.data$hazard),
      se     = as.numeric(.data$se))
  stopifnot(nrow(pred) == nrow(newdata))

  mutate_vars <- c("hazard")
  if (ci) {
    pred <- pred %>%
      mutate(
        ci_lower = .data$hazard - se.mult * .data$se,
        ci_upper = .data$hazard + se.mult * .data$se)
    mutate_vars <- c(mutate_vars, "ci_lower", "ci_upper")
  }

  if (type == "response") {
    pred <- pred %>% mutate_at(mutate_vars, exp)
  }

  # it is necessary to include the grouping variables here, otherwise
  # functions calculating the cumulative hazard will cumulate over all rows
  # instead of group wise
  if (is.grouped_df(newdata)) {
    group.df <- newdata %>% select(one_of(group_vars(newdata)))
    pred     <- bind_cols(group.df, pred)
  }

  return(pred)

}

#' Add cumulative hazard estimate to data set
#'
#' @rdname add_hazard
#' @inheritParams add_hazard
#' @param interval_length The variable in newdata containing the interval lengths.
#' Can be either bare unquoted variable name or character. Defaults to \code{"intlen"}.
#' @importFrom dplyr bind_cols
#' @seealso \code{\link[mgcv]{predict.gam}}, \code{\link[pammtools]{add_hazard}}
#' @export
add_cumu_hazard <- function(
  newdata,
  object,
  ci              = TRUE,
  se.mult         = 2,
  overwrite       = FALSE,
  time_variable   = NULL,
  interval_length = "intlen",
  ...)  {

  interval_length <- quo_name(enquo(interval_length))

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
    time_variable = time_variable, interval_length = interval_length,
    ...)

  newdata <- newdata %>% bind_cols(rm_grpvars(pred))

  return(newdata)

}

#' Calculate cumulative hazard
#'
#' @inheritParams add_cumu_hazard
#' @import checkmate dplyr
#' @importFrom rlang UQ sym quo_name .data
#' @keywords internal
get_cumu_hazard <- function(
  newdata,
  object,
  ci              = TRUE,
  time_variable   = NULL,
  interval_length = "intlen",
  ...) {

  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))

  interval_length <- sym(interval_length)

  lengths <- select(rm_grpvars(newdata), !!interval_length)

  mutate_args  <- list(cumu_hazard = quo(cumsum(.data$hazard * (!!interval_length))))
  vars_exclude <- c("hazard", "se")
  if (ci) {
    mutate_args <- mutate_args %>%
      append(
        list(
          cumu_lower = quo(cumsum(.data$ci_lower * (!!interval_length))),
          cumu_upper = quo(cumsum(.data$ci_upper * (!!interval_length)))))
    vars_exclude <- c(vars_exclude, "ci_lower", "ci_upper")
  }
  get_hazard(newdata, object, ci = ci, type = "response",
    time_variable = time_variable, ...) %>%
    bind_cols(lengths) %>%
    mutate(!!!mutate_args) %>%
    select(-one_of(vars_exclude)) %>%
    select(-!!interval_length)
}


#' Add survival probability estimates
#'
#' Given suitable data (i.e. data with all columns used for estimation of the model),
#' this functions adds a column \code{surv_prob} containing survival probabilities
#' for the specified covariate and follow-up information (and CIs
#' \code{surv_lower}, \code{surv_upper} if \code{ci=TRUE}).
#'
#' @inherit add_cumu_hazard
#' @examples
#' ped <- tumor[1:50,] %>% as_ped(Surv(days, status)~ age)
#' pam <- mgcv::gam(ped_status ~ s(tend)+age, data=ped, family=poisson(), offset=offset)
#' ped_info(ped) %>% add_surv_prob(pam, ci=TRUE)
#' @export
add_surv_prob <- function(
  newdata,
  object,
  ci              = TRUE,
  se.mult         = 2,
  overwrite       = FALSE,
  time_variable   = NULL,
  interval_length = "intlen",
  ...)  {

  interval_length <- quo_name(enquo(interval_length))

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
    time_variable = time_variable, interval_length = interval_length,
    ...)

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

  assert_character(interval_length)
  assert_choice(interval_length, colnames(newdata))

  exclude_vars <- c("cumu_hazard")
  mutate_args  <- list(surv_prob = quo(exp(-.data$cumu_hazard)))
  if (ci) {
    mutate_args <- mutate_args %>%
      append(
        list(
          surv_lower = quo(exp(-.data$cumu_lower)),
          surv_upper = quo(exp(-.data$cumu_upper))))
    exclude_vars <- c(exclude_vars, "cumu_lower", "cumu_upper")
  }

  newdata %>%
    get_cumu_hazard(object, ci = ci, time_variable = time_variable,
      interval_length = interval_length, ...) %>%
    mutate(!!!mutate_args) %>%
    select(-one_of(exclude_vars))

}
