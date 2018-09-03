#' Add info about term effects to data set
#'
#' Adds the contribution of a specific term to the
#' linear predictor to the data specified by \code{newdata}.
#' Essentially a wrapper to \code{\link[mgcv]{predict.gam}}, with \code{type="terms"}.
#' Thus most arguments and their documentation below is from \code{\link[mgcv]{predict.gam}}.
#'
#'
#' @inheritParams mgcv::predict.gam
#' @param term A character (vector) or regular expression indicating for
#' which term(s) information should be extracted and added to data set.
#' @param se_mult The factor by which standard errors are multiplied to form
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
  se_mult  = 2,
  relative = FALSE,
  ...) {

  assert_data_frame(newdata, all.missing = FALSE)
  assert_character(term, min.chars = 1, any.missing = FALSE, min.len = 1)

  col_ind <- map(term, grep, x = names(object$coefficients)) %>%
    unlist() %>% unique() %>% sort()
  is_gam <- inherits(object, "gam")

  X <- if (is_gam) {
    predict.gam(object, newdata = newdata,
      type = "lpmatrix", ...)[, col_ind, drop = FALSE]
  } else  {
    model.matrix(object$formula[-2], data = newdata)[, col_ind, drop = FALSE]
  }
  if (relative) {
    if (is.null(object$model)) {
      stop("Relative risk can only be calculated when original data is present in 'object'!")
    }
    data_bar <- sample_info(object$model)[rep(1, nrow(newdata)), ]
    X_bar <- if (is_gam) {
      predict.gam(
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
    cov.coefs <- if (is_gam) {
      object$Vp[col_ind, col_ind]
    } else {
      vcov(object)[col_ind, col_ind]
    }
    se <- sqrt(rowSums( (X %*% cov.coefs) * X ))
    newdata <- newdata %>%
      mutate(
        ci_lower = .data$fit - se_mult * se,
        ci_upper = .data$fit + se_mult * se)
  }

  return(newdata)

}

add_term2 <- function(
  newdata,
  object,
  term,
  se.fit    = TRUE,
  type      = "terms",
  se_mult   = 2,
  reference = NULL,
  ...) {

  assert_data_frame(newdata, all.missing = FALSE)
  assert_character(term, min.chars = 1, any.missing = FALSE, min.len = 1)

  col_ind <- map(term, grep, x = names(object$coefficients)) %>%
    unlist() %>% unique() %>% sort()
  is_gam <- inherits(object, "gam")

  X <- if (is_gam) {
    predict.gam(object, newdata = newdata,
      type = "lpmatrix", ...)[, col_ind, drop = FALSE]
  } else  {
    model.matrix(object$formula[-2], data = newdata)[, col_ind, drop = FALSE]
  }
  if (!is.null(reference)) {
    reference <- newdata %>% mutate(!!!reference)
    X_ref <- if (is_gam) {
      predict.gam(
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
    cov.coefs <- if (is_gam) {
      object$Vp[col_ind, col_ind]
    } else {
      vcov(object)[col_ind, col_ind]
    }
    se <- sqrt(rowSums( (X %*% cov.coefs) * X ))
    newdata <- newdata %>%
      mutate(
        ci_lower = .data$fit - se_mult * se,
        ci_upper = .data$fit + se_mult * se)
  }

  return(newdata)

}


#' Add predicted (cumulative) hazard to data set
#'
#' Add (cumulative) hazard based on the provided data set and model.
#' If \code{ci=TRUE} confidence intervals are also returned. Their width can
#' be controlled via the \code{se_mult} argument.
#'
#' @rdname add_hazard
#' @inheritParams mgcv::predict.gam
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}} and
#'   \code{\link{get_hazard}}
#' @param ci Logical indicating whether to include confidence intervals. Defaults
#' to \code{TRUE}
#' @param ci_type The method by which standard errors/confidence intervals
#' will be calculated. Default transforms the linear predictor at
#' respective intervals. \code{"delta"} calculates CIs based on the standard
#' error calculated by the Delta method. \code{"sim"} draws the
#' property of interest from its posterior based on the normal distribution of
#' the estimated coefficients. CIs are given by respective quantiles.
#' @param se_mult Factor by which standard errors are multiplied for calculating
#' the confidence intervals.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to \code{FALSE}. If \code{TRUE}, columns with names
#' \code{c("hazard", "se", "lower", "upper")} will be overwritten.
#' @param time_var Name of the variable used for the baseline hazard. If
#'   not given, defaults to \code{"tend"} for \code{\link[mgcv]{gam}} fits, else
#'   \code{"interval"}. The latter is assumed to be a factor, the former
#'   numeric.
#' @import checkmate dplyr mgcv
#' @importFrom stats predict
#' @examples
#' ped <- tumor[1:50,] %>% as_ped(Surv(days, status)~ age)
#' pam <- mgcv::gam(ped_status ~ s(tend)+age, data = ped, family=poisson(), offset=offset)
#' ped_info(ped) %>% add_hazard(pam, type="link")
#' ped_info(ped) %>% add_hazard(pam, type = "response")
#' ped_info(ped) %>% add_cumu_hazard(pam)
#' @export
add_hazard <- function(
  newdata,
  object,
  type          = c("response", "link"),
  ci            = TRUE,
  se_mult       = 2,
  ci_type       = c("default", "delta", "sim"),
  overwrite     = FALSE,
  time_var = NULL,
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

  get_hazard(newdata, object, ci = ci, type = type, se_mult = se_mult,
    ci_type = ci_type, time_var = time_var, ...)

}

#' Calculate predicted hazard
#'
#' @inheritParams add_hazard
#' @importFrom stats model.frame
#' @importFrom mgcv predict.gam predict.bam
#' @keywords internal
get_hazard <- function(
  newdata,
  object,
  ci            = TRUE,
  type          = c("response", "link"),
  ci_type       = c("default", "delta", "sim"),
  time_var = NULL,
  se_mult        = 2,
  ...)  {

  assert_data_frame(newdata, all.missing = FALSE)
  assert_class(object, classes = "glm")
  type    <- match.arg(type)
  ci_type <- match.arg(ci_type)

  is_gam <- inherits(object, "gam")
  if (is.null(time_var)) {
    time_var <- ifelse(is_gam, "tend", "interval")
  } else {
    assert_string(time_var)
    assert_choice(time_var, colnames(newdata))
  }

  # throw warning or error if evaluation time points/intervals do not correspond
  # to evaluation time-points/intervals do not correspond to the ones used for
  # estimation
  warn_about_new_time_points(newdata, object, time_var)

  if (is_gam) {
    X <- predict.gam(object, newdata, type = "lpmatrix")
  } else {
    X <- model.matrix(object$formula[-2], data = newdata)
  }
  coefs <- coef(object)
  newdata$hazard <- unname(drop(X %*% coefs))
  if (ci) {
    newdata <- newdata %>% add_ci(object, X, type = type, ci_type = ci_type,
      se_mult = se_mult, ...)
  }
  if (type == "response") {
    newdata <- newdata %>% mutate(hazard = exp(.data$hazard))
  }

  newdata

}


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
  se_mult         = 2,
  overwrite       = FALSE,
  time_var   = NULL,
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

  get_cumu_hazard(newdata, object, ci = ci, se_mult = se_mult,
    time_var = time_var, interval_length = interval_length, ...)

}

#' Calculate cumulative hazard
#'
#' @inheritParams add_cumu_hazard
#' @import checkmate dplyr
#' @importFrom rlang UQ sym quo_name .data
#' @importFrom purrr map_lgl
#' @keywords internal
get_cumu_hazard <- function(
  newdata,
  object,
  ci              = TRUE,
  ci_type         = c("default", "delta", "sim"),
  time_var   = NULL,
  se_mult         = 2,
  interval_length = "intlen",
  nsim            = 100L, ...)  {

  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)
  assert_class(object, classes = "glm")

  ci_type <- match.arg(ci_type)

  interval_length <- sym(interval_length)

  mutate_args  <- list(cumu_hazard = quo(cumsum(.data$hazard *
    (!!interval_length))))
  haz_vars_in_data <- map(c("hazard", "se", "ci_lower", "ci_upper"),
    ~ grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
  vars_exclude <- c("hazard")

  if (ci) {
    if (ci_type == "default" | ci_type == "delta") {
      vars_exclude <- c(vars_exclude, "se", "ci_lower", "ci_upper")
      newdata <- get_hazard(newdata, object, type = "response", ci = ci,
        ci_type = ci_type, time_var = time_var, se_mult = se_mult, ...)
      if (ci_type == "default") {
        mutate_args <- mutate_args %>%
          append(list(
            cumu_lower = quo(cumsum(.data$ci_lower * (!!interval_length))),
            cumu_upper = quo(cumsum(.data$ci_upper * (!!interval_length)))))
      } else {
        # ci delta rule
        newdata <- split(newdata, group_indices(newdata)) %>%
            map_dfr(add_delta_ci_cumu, object = object, se_mult = se_mult)
      }
    } else {
      if (ci_type == "sim") {
        newdata <- get_hazard(newdata, object, type = "response", ci = FALSE,
          time_var = time_var, ...)
        ci_df <- split(newdata, group_indices(newdata)) %>%
          map_dfr(get_sim_ci_cumu, object = object, nsim = nsim)
        newdata <- newdata %>% bind_cols(ci_df)
      }
    }
  } else {
    newdata <- newdata %>%
      get_hazard(object, type = "response", ci = ci,
        ci_type = ci_type, time_var = time_var, se_mult = se_mult, ...)
  }
  newdata <- newdata %>%
    mutate(!!!mutate_args)

  vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
  if (length(vars_exclude) != 0 ) {
    newdata <- newdata %>% select(-one_of(vars_exclude))
  }

  newdata

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
  se_mult         = 2,
  overwrite       = FALSE,
  time_var   = NULL,
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

  get_surv_prob(newdata, object, ci = ci, se_mult = se_mult,
    time_var = time_var, interval_length = interval_length, ...)

}


#' Calculate survival probabilities
#'
#' @inheritParams add_surv_prob
#' @keywords internal
get_surv_prob <- function(
  newdata,
  object,
  ci              = TRUE,
  ci_type         = c("default", "delta", "sim"),
  se_mult         = 2L,
  time_var   = NULL,
  interval_length = "intlen",
  nsim            = 100L,
  ...) {

  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)
  assert_class(object, classes = "glm")

  ci_type <- match.arg(ci_type)

  interval_length <- sym(interval_length)

  mutate_args  <- list(surv_prob = quo(exp(-cumsum(.data$hazard *
    (!!interval_length)))))
  haz_vars_in_data <- map(c("hazard", "se", "ci_lower", "ci_upper"),
    ~grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
  vars_exclude <- c("hazard")

  if (ci) {
    if (ci_type == "default" | ci_type == "delta") {
      vars_exclude <- c(vars_exclude, "se", "ci_lower", "ci_upper")
      newdata <- get_hazard(newdata, object, type = "response", ci = ci,
        ci_type = ci_type, time_var = time_var,  se_mult = se_mult, ...)
      if (ci_type == "default") {
        mutate_args <- mutate_args %>%
          append(list(
            surv_upper = quo(exp(-cumsum(.data$ci_lower * (!!interval_length)))),
            surv_lower = quo(exp(-cumsum(.data$ci_upper * (!!interval_length))))))
      } else {
        # ci delta rule
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(add_delta_ci_surv, object = object, se_mult = se_mult)
      }
    } else {
      if (ci_type == "sim") {
        newdata <- get_hazard(newdata, object, type = "response", ci = FALSE,
          time_var = time_var)
        ci_df <- split(newdata, group_indices(newdata)) %>%
          map_dfr(get_sim_ci_surv, object = object, nsim = nsim)
        newdata <- newdata %>% bind_cols(ci_df)
      }
    }
  } else {
    newdata <- newdata %>%
      get_hazard(object = object, type = "response", ci = FALSE,
        time_var = time_var)
  }

  newdata <- newdata %>%
    mutate(!!!mutate_args)

  vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
  if (length(vars_exclude) != 0 ) {
    newdata <- newdata %>% select(-one_of(vars_exclude))
  }

  newdata

}

add_ci <- function(
  newdata,
  object,
  X,
  type    = c("response", "link"),
  se_mult = 2,
  ci_type = c("default", "delta", "sim"),
  nsim = 100) {

  ci_type <- match.arg(ci_type)

  is_gam <- inherits(object, "gam")
  if (is_gam) {
    V <- object$Vp
  } else {
    V <- vcov(object)
  }
  se <- unname(sqrt(rowSums( (X %*% V) * X) ))
  newdata$se <- se
  if (type == "link") {
    newdata <- newdata %>%
      mutate(
        ci_lower = .data$hazard - se_mult * .data$se,
        ci_upper = .data$hazard + se_mult * .data$se)
  }

  if (type != "link") {
    if (ci_type == "default") {
      newdata <- newdata %>%
        mutate(
          ci_lower = exp(.data$hazard - se_mult * .data$se),
          ci_upper = exp(.data$hazard + se_mult * .data$se))
    } else {
      if (ci_type == "delta") {
        newdata <- split(newdata, group_indices(newdata)) %>%
            map_dfr(add_delta_ci, object = object, se_mult = se_mult)
      } else {
        if (ci_type == "sim") {
          ci_df <- split(newdata, group_indices(newdata)) %>%
            map_dfr(get_sim_ci, object = object)
          newdata <- newdata %>% bind_cols(ci_df)
        }
      }
    }
  }
  newdata
}

add_delta_ci <- function(newdata, object, se_mult = 2) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  V     <- object$Vp

  Jacobi <- diag(exp(newdata$hazard)) %*% X
  newdata %>%
    mutate(
      se       = sqrt(rowSums( (Jacobi %*% V) * Jacobi )),
      ci_lower = exp(.data$hazard) - .data$se * se_mult,
      ci_upper = exp(.data$hazard) + .data$se * se_mult)

}

add_delta_ci_cumu <- function(newdata, object, se_mult = 2) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  V     <- object$Vp

  Delta  <- lower.tri(diag(nrow(X)), diag = TRUE) %*% diag(newdata$intlen)
  Jacobi <- diag(newdata$hazard) %*% X
  LHS <- Delta %*% Jacobi
  newdata %>%
    mutate(
      se       = sqrt(rowSums( (LHS %*% V) * LHS )),
      cumu_lower = cumsum(.data$intlen * .data$hazard) - .data$se * se_mult,
      cumu_upper = cumsum(.data$intlen * .data$hazard) + .data$se * se_mult)

}

add_delta_ci_surv <- function(newdata, object, se_mult = 2) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  V     <- object$Vp

  Delta  <- lower.tri(diag(nrow(X)), diag = TRUE) %*% diag(newdata$intlen)
  Jacobi <- diag(newdata$hazard) %*% X
  LHS <- -diag(exp(-rowSums(Delta %*% diag(newdata$hazard)))) %*%
    (Delta %*% Jacobi)
  newdata %>%
    mutate(
      se       = sqrt(rowSums( (LHS %*% V) * LHS)),
      surv_lower = exp(-cumsum(.data$hazard * .data$intlen)) - .data$se * se_mult,
      surv_upper = exp(-cumsum(.data$hazard * .data$intlen)) + .data$se * se_mult)

}

#' Calculate simulation based confidence intervals
#' @keywords internal
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef
get_sim_ci <- function(newdata, object, alpha = 0.05, nsim = 100L) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  V     <- object$Vp
  coefs <- coef(object)

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))

  ci_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  ci_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)

  tibble(ci_lower = ci_lower, ci_upper = ci_upper)

}


get_sim_ci_cumu <- function(newdata, object, alpha = 0.05, nsim = 100L) {

  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  V     <- object$Vp
  coefs <- coef(object)

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigm = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    cumsum(newdata$intlen * exp(X %*% z)))

  cumu_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  cumu_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)

  tibble(cumu_lower = cumu_lower, cumu_upper = cumu_upper)

}

get_sim_ci_surv <- function(newdata, object, alpha = 0.05, nsim = 100L) {

  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  V     <- object$Vp
  coefs <- coef(object)

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    exp(-cumsum(newdata$intlen * exp(X %*% z))))

  surv_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  surv_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)

  tibble(surv_lower = surv_lower, surv_upper = surv_upper)

}
