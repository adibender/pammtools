#' Embeds the data set with the specified (relative) term contribution
#'
#' Adds the contribution of a specific term to the
#' linear predictor to the data specified by \code{newdata}.
#' Essentially a wrapper to \code{\link[mgcv]{predict.gam}}, with \code{type="terms"}.
#' Thus most arguments and their documentation below is from \code{\link[mgcv]{predict.gam}}.
#'
#' @inheritParams mgcv::predict.gam
#' @param term A character (vector) or regular expression indicating for
#' which term(s) information should be extracted and added to data set.
#' @param ci \code{logical}. Indicates if confidence intervals should be
#' calculated. Defaults to \code{TRUE}.
#' @param se_mult The factor by which standard errors are multiplied to form
#' confidence intervals.
#' @param reference A data frame with number of rows equal to \code{nrow(newdata)} or
#' one, or a named list with (partial) covariate specifications. See examples.
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}}
#' @import checkmate dplyr mgcv
#' @importFrom stats predict
#' @importFrom purrr map
#' @importFrom stats model.matrix vcov
#' @examples
#' library(ggplot2)
#' ped <- as_ped(tumor, Surv(days, status)~ age, cut = seq(0, 2000, by = 100))
#' pam <- mgcv::gam(ped_status ~ s(tend) + s(age), family = poisson(),
#'   offset = offset, data = ped)
#' #term contribution for sequence of ages
#' s_age <- ped %>% make_newdata(age = seq_range(age, 50)) %>%
#'   add_term(pam, term = "age")
#' ggplot(s_age, aes(x = age, y = fit)) + geom_line() +
#'   geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .3)
#' # term contribution relative to mean age
#' s_age2 <- ped %>% make_newdata(age = seq_range(age, 50)) %>%
#'   add_term(pam, term = "age", reference = list(age = mean(.$age)))
#' ggplot(s_age2, aes(x = age, y = fit)) + geom_line() +
#'   geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .3)
#' @export
add_term <- function(
  newdata,
  object,
  term,
  reference = NULL,
  ci        = TRUE,
  se_mult   = 2,
  ...) {

  assert_data_frame(newdata, all.missing = FALSE)
  assert_character(term, min.chars = 1, any.missing = FALSE, min.len = 1)

  col_ind <- map(term, grep, x = names(object$coefficients)) %>%
    unlist() %>% unique() %>% sort()
  is_gam <- (inherits(object, "gam") | inherits(object, "scam"))

  X <- prep_X(object, newdata, reference, ...)[, col_ind, drop = FALSE]

  newdata[["fit"]] <- unname(drop(X %*% object$coefficients[col_ind]))
  if (ci) {
    cov.coefs <- if (is_gam) {
      object$Vp[col_ind, col_ind]
    } else {
      vcov(object)[col_ind, col_ind]
    }
    se <- unname(sqrt(rowSums( (X %*% cov.coefs) * X )))
    newdata <- newdata %>%
      mutate(
        ci_lower = .data[["fit"]] - se_mult * se,
        ci_upper = .data[["fit"]] + se_mult * se)
  }

  return(newdata)

}


#' Create design matrix from a suitable object
#'
#' @keywords internal
#' @param object A suitable object from which a design matrix can be generated.
#' Often a model object.
make_X <- function(object, ...) {

  UseMethod("make_X", object)

}

#' @inherit make_X
#' @keywords internal
#' @rdname make_X
#' @inherit make_X
#' @param newdata A data frame from which design matrix will be constructed
make_X.default <- function(object, newdata, ...) {

  model.matrix(object$formula[-2], data = newdata, ...)

}

#' @inherit make_X
#' @inherit make_X.default
#' @rdname make_X
#' @keywords internal
make_X.gam <- function(object, newdata, ...) {

  predict.gam(object, newdata = newdata, type = "lpmatrix", ...)

}

#' @inherit make_X
#' @importFrom scam predict.scam
#' @keywords internal
make_X.scam <- function(object, newdata, ...) {

  X <- predict.scam(object, newdata = newdata, type = "lpmatrix", ...)

}

prep_X <- function(object, newdata, reference = NULL, ...) {

  X <- make_X(object, newdata, ...)

  if (!is.null(reference)) {
    reference <- preproc_reference(reference, colnames(newdata), nrow(newdata))
    reference <- newdata %>% mutate(!!!reference)
    X_ref <- make_X(object, reference, ...)
    X <- X - X_ref
  }

  X
}


preproc_reference <- function(reference, cnames, n_rows) {

  # check that provided variables contained in newdata
  names_ref <- names(reference)
  if (!check_subset(names_ref, cnames)) {
    stop(paste0("Columns in 'reference' but not in 'newdata':",
      paste0(setdiff(names_ref, cnames), collapse = ",")))
  }
  # transform to list if inherits from data frame, so it can be processed
  # in mutate via !!!
  if (inherits(reference, "data.frame")) {
    if (!(nrow(reference) == n_rows || nrow(reference) == 1)) {
      stop("If reference is provided as data frame, number of rows must be
        either 1 or the number of rows in newdata.")
    }
    reference <- as.list(reference)
  }

  reference

}

#' Add predicted (cumulative) hazard to data set
#'
#' Add (cumulative) hazard based on the provided data set and model.
#' If \code{ci=TRUE} confidence intervals (CI) are also added. Their width can
#' be controlled via the \code{se_mult} argument. The method by which the
#' CI are calculated can be specified by \code{ci_type}.
#' This is a wrapper around
#' \code{\link[mgcv]{predict.gam}}. When \code{reference} is specified, the
#' (log-)hazard ratio is calculated.
#'
#' @rdname add_hazard
#' @inheritParams mgcv::predict.gam
#' @inheritParams add_term
#' @param type Either \code{"response"} or \code{"link"}. The former calculates
#' hazard, the latter the log-hazard.
#' @param ... Further arguments passed to \code{\link[mgcv]{predict.gam}} and
#'   \code{\link{get_hazard}}
#' @param ci_type The method by which standard errors/confidence intervals
#' will be calculated. Default transforms the linear predictor at
#' respective intervals. \code{"delta"} calculates CIs based on the standard
#' error calculated by the Delta method. \code{"sim"} draws the
#' property of interest from its posterior based on the normal distribution of
#' the estimated coefficients. See \href{https://adibender.github.io/simpamm/confidence-intervals.html}{here}
#' for details and empirical evaluation.
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
add_hazard <- function(newdata, object, ...) {
  UseMethod("add_hazard", object)
}

#' @rdname add_hazard
#' @export
add_hazard.default <- function(
  newdata,
  object,
  reference = NULL,
  type      = c("response", "link"),
  ci        = TRUE,
  se_mult   = 2,
  ci_type   = c("default", "delta", "sim"),
  overwrite = FALSE,
  time_var  = NULL,
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

  get_hazard(object, newdata, reference = reference,
    ci = ci, type = type, se_mult = se_mult, ci_type = ci_type,
    time_var = time_var, ...)

}

#' Calculate predicted hazard
#'
#' @inheritParams add_hazard
#' @importFrom stats model.frame
#' @importFrom mgcv predict.gam predict.bam
#' @keywords internal
get_hazard <- function(object, newdata, ...) {
  UseMethod("get_hazard", object)
}

#' @rdname get_hazard
get_hazard.default <- function(
  object,
  newdata,
  reference = NULL,
  ci        = TRUE,
  type      = c("response", "link"),
  ci_type   = c("default", "delta", "sim"),
  time_var  = NULL,
  se_mult   = 2,
  ...)  {

  assert_data_frame(newdata, all.missing = FALSE)
  assert_class(object, classes = "glm")
  type    <- match.arg(type)
  ci_type <- match.arg(ci_type)

  is_gam <- (inherits(object, "gam") | inherits(object, "scam"))
  if (is.null(time_var)) {
    time_var <- ifelse(is_gam, "tend", "interval")
  } else {
    assert_string(time_var)
    assert_choice(time_var, colnames(newdata))
  }

  # throw warning or error if evaluation time points/intervals do not correspond
  # to evaluation time-points/intervals do not correspond to the ones used for
  # estimation
  #warn_about_new_time_points(object, newdata, time_var)

  X <- prep_X(object, newdata, reference, ...)
  coefs <- coef(object)
  newdata$hazard <- unname(drop(X %*% coefs))
  if (ci) {
    newdata <- newdata %>%
      add_ci(object, X, type = type, ci_type = ci_type, se_mult = se_mult, ...)
  }
  if (type == "response") {
    newdata <- newdata %>% mutate(hazard = exp(.data[["hazard"]]))
  }

  newdata %>% arrange(.data[[time_var]], .by_group = TRUE)

}


#' @rdname add_hazard
#' @inheritParams add_hazard
#' @param interval_length The variable in newdata containing the interval lengths.
#' Can be either bare unquoted variable name or character. Defaults to \code{"intlen"}.
#' @importFrom dplyr bind_cols
#' @seealso \code{\link[mgcv]{predict.gam}},
#' \code{\link[pammtools]{add_surv_prob}}
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

  is_gam <- (inherits(object, "gam") | inherits(object, "scam"))
  if (is.null(time_var)) {
    time_var <- ifelse(is_gam, "tend", "interval")
  } else {
    assert_string(time_var)
    assert_choice(time_var, colnames(newdata))
  }


  trafo_args <- attr(newdata, "trafo_args")
  intvars    <- attr(newdata, "intvars")

  times <- setdiff(sort(unique(newdata[[time_var]])), c(0))
  brks <- setdiff(trafo_args[["cut"]][trafo_args[["cut"]]<= max(times)], c(0))

  # if selected time points contain all times already, do not extend newdata
  if (all(brks %in% times)) {
    joindata <- newdata
  } else {
    if (length(groups(newdata))!=0) {
      old_groups <- dplyr::groups(newdata)
      joindata <- group_split(newdata) |>
        map(newdata, .f = ~ expand_df(.x, object, trafo_args, intvars, time_var))|> #expand uses distinct, hence need to regroup
        map(newdata, .f = ~ group_by(.x, !!!old_groups)) |>
        bind_rows()
    } else {
      joindata <- newdata %>% expand_df(object, trafo_args, intvars)
    }
  }

  joindata <- get_cumu_hazard(joindata, object, ci = ci, se_mult = se_mult,
                              time_var = time_var, interval_length = interval_length, ...)

  suppressMessages(
    newdata %>% left_join(joindata)
  )

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

  mutate_args  <- list(cumu_hazard = quo(cumsum(.data[["hazard"]] *
    (!!interval_length))))
  haz_vars_in_data <- map(c("hazard", "se", "ci_lower", "ci_upper"),
    ~ grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
  vars_exclude <- c("hazard")

  if (ci) {
    if (ci_type == "default" | ci_type == "delta") {
      vars_exclude <- c(vars_exclude, "se", "ci_lower", "ci_upper")
      newdata <- get_hazard(object, newdata, type = "response", ci = ci,
        ci_type = ci_type, time_var = time_var, se_mult = se_mult, ...)
      if (ci_type == "default") {
        mutate_args <- mutate_args %>%
          append(list(
            cumu_lower = quo(cumsum(.data[["ci_lower"]] * (!!interval_length))),
            cumu_upper = quo(cumsum(.data[["ci_upper"]] * (!!interval_length)))))
      } else {
        # ci delta rule
        newdata <- split(newdata, group_indices(newdata)) %>%
            map_dfr(add_delta_ci_cumu, object = object, se_mult = se_mult, ...)
      }
    } else {
      if (ci_type == "sim") {
        newdata <- get_hazard(object, newdata, type = "response", ci = FALSE,
          time_var = time_var, ...)
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(get_sim_ci_cumu, object = object, nsim = nsim, ...)
      }
    }
  } else {
    newdata <-
      get_hazard(object, newdata, type = "response", ci = ci,
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

  mutate_args  <- list(surv_prob = quo(exp(-cumsum(.data[["hazard"]] *
    (!!interval_length)))))
  haz_vars_in_data <- map(c("hazard", "se", "ci_lower", "ci_upper"),
    ~grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
  vars_exclude <- c("hazard")

  if (ci) {
    if (ci_type == "default" | ci_type == "delta") {
      vars_exclude <- c(vars_exclude, "se", "ci_lower", "ci_upper")
      newdata <- get_hazard(object, newdata, type = "response", ci = ci,
        ci_type = ci_type, time_var = time_var,  se_mult = se_mult, ...)
      if (ci_type == "default") {
        mutate_args <- mutate_args %>%
          append(list(
            surv_upper = quo(exp(-cumsum(.data[["ci_lower"]] * (!!interval_length)))),
            surv_lower = quo(exp(-cumsum(.data[["ci_upper"]] * (!!interval_length))))))
      } else {
        # ci delta rule
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(add_delta_ci_surv, object = object, se_mult = se_mult, ...)
      }
    } else {
      if (ci_type == "sim") {
        newdata <- get_hazard(object, newdata, type = "response", ci = FALSE,
          time_var = time_var, ...)
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(get_sim_ci_surv, object = object, nsim = nsim, ...)
      }
    }
  } else {
    newdata <-
      get_hazard(object = object, newdata, type = "response", ci = FALSE,
        time_var = time_var, ...)
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
  nsim = 100, ...) {

  ci_type <- match.arg(ci_type)

  is_gam <- (inherits(object, "gam") | inherits(object, "scam"))
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
        ci_lower = .data[["hazard"]] - se_mult * .data[["se"]],
        ci_upper = .data[["hazard"]] + se_mult * .data[["se"]])
  }

  if (type != "link") {
    if (ci_type == "default") {
      newdata <- newdata %>%
        mutate(
          ci_lower = exp(.data[["hazard"]] - se_mult * .data[["se"]]),
          ci_upper = exp(.data[["hazard"]] + se_mult * .data[["se"]]))
    } else {
      if (ci_type == "delta") {
        newdata <- split(newdata, group_indices(newdata)) %>%
            map_dfr(add_delta_ci, object = object, se_mult = se_mult, ...)
      } else {
        if (ci_type == "sim") {
          newdata <- split(newdata, group_indices(newdata)) %>%
            map_dfr(get_sim_ci, object = object, nsim = nsim, ...)
        }
      }
    }
  }
  newdata
}

add_delta_ci <- function(newdata, object, se_mult = 2, ...) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp

  Jacobi <- diag(exp(newdata$hazard)) %*% X
  newdata %>%
    mutate(
      se       = sqrt(rowSums( (Jacobi %*% V) * Jacobi )),
      ci_lower = exp(.data[["hazard"]]) - .data[["se"]] * se_mult,
      ci_upper = exp(.data[["hazard"]]) + .data[["se"]] * se_mult)

}

add_delta_ci_cumu <- function(newdata, object, se_mult = 2, ...) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp

  Delta  <- lower.tri(diag(nrow(X)), diag = TRUE) %*% diag(newdata$intlen)
  Jacobi <- diag(newdata$hazard) %*% X
  LHS <- Delta %*% Jacobi
  newdata %>%
    mutate(
      se       = sqrt(rowSums( (LHS %*% V) * LHS )),
      cumu_lower = cumsum(.data[["intlen"]] * .data[["hazard"]]) - .data[["se"]] * se_mult,
      cumu_upper = cumsum(.data[["intlen"]] * .data[["hazard"]]) + .data[["se"]] * se_mult)

}

add_delta_ci_surv <- function(newdata, object, se_mult = 2, ...) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp

  Delta  <- lower.tri(diag(nrow(X)), diag = TRUE) %*% diag(newdata$intlen)
  Jacobi <- diag(newdata$hazard) %*% X
  LHS <- -diag(exp(-rowSums(Delta %*% diag(newdata$hazard)))) %*%
    (Delta %*% Jacobi)
  newdata %>%
    mutate(
      se       = sqrt(rowSums( (LHS %*% V) * LHS)),
      surv_lower = exp(-cumsum(.data[["hazard"]] * .data[["intlen"]])) - .data[["se"]] * se_mult,
      surv_upper = exp(-cumsum(.data[["hazard"]] * .data[["intlen"]])) + .data[["se"]] * se_mult)

}

#' Calculate simulation based confidence intervals
#'
#' @keywords internal
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef
get_sim_ci <- function(newdata, object, alpha = 0.05, nsim = 100L, ...) {
  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp
  coefs <- coef(object)

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))

  newdata$ci_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  newdata$ci_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)

  newdata

}


get_sim_ci_cumu <- function(newdata, object, alpha = 0.05, nsim = 100L, ...) {

  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp
  coefs <- coef(object)

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    cumsum(newdata$intlen * exp(X %*% z)))

  newdata$cumu_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  newdata$cumu_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)

  newdata

}

get_sim_ci_surv <- function(newdata, object, alpha = 0.05, nsim = 100L, ...) {

  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix", ...)
  V     <- object$Vp
  coefs <- coef(object)

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    exp(-cumsum(newdata$intlen * exp(X %*% z))))

  newdata$surv_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2)
  newdata$surv_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2)

  newdata

}


## Cumulative Incidence Function (CIF) for competing risks data


#' Add cumulative incidence function to data
#'
#' @inheritParams add_hazard
#' @param alpha The alpha level for confidence/credible intervals.
#' @param nsim Number of simulations (draws from posterior of estimated coefficients)
#' on which estimation of CIFs and their confidence/credible intervals will be
#' based on.
#' @param cause_var Character. Column name of the 'cause' variable.
#'
#' @export
add_cif <- function(
  newdata,
  object,
  ...) {

  UseMethod("add_cif", object)

}


#' @rdname add_cif
#' @export
add_cif.default <- function(
  newdata,
  object,
  ci        = TRUE,
  overwrite = FALSE,
  alpha     = 0.05,
  nsim      = 500L,
  cause_var = "cause",
  time_var  = NULL,
  ...) {

  coefs        <- coef(object)
  V            <- object$Vp
  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)

  map_dfr(
    split(newdata, group_indices(newdata)),
    ~get_cif(
      newdata = .x, object = object, ci = ci, alpha = alpha, nsim = nsim,
      cause_var = cause_var, coefs = coefs, V = V, sim_coef_mat = sim_coef_mat,
      time_var = time_var, ...)
  )

}

#' Calculate CIF for one cause
#'
#' @keywords internal
get_cif <- function(newdata, object, ...) {

  UseMethod("get_cif", object)

}

#' @rdname get_cif
#' @keywords internal
get_cif.default <- function(
  newdata,
  object,
  ci,
  time_var,
  alpha,
  nsim,
  cause_var,
  coefs,
  V,
  sim_coef_mat,
  ...) {

  is_gam <- (inherits(object, "gam") | inherits(object, "scam"))
  if (is.null(time_var)) {
    time_var <- ifelse(is_gam, "tend", "interval")
  } else {
    assert_string(time_var)
    assert_choice(time_var, colnames(newdata))
  }


  causes_model <- as.factor(object$attr_ped$risks)
  cause_data   <- unique(newdata[[cause_var]])

  if(length(cause_data) > 1) {
    stop("Did you forget to group by cause?")
  }

  hazards <- map(
    causes_model,
    ~ {
        .df <- mutate(newdata, cause = .x) %>%
          arrange(.data[[time_var]], .by_group = TRUE)
        X <- predict(object, .df, type = "lpmatrix")
        apply(sim_coef_mat, 1, function(z) exp(X %*% z))
      }
    )
  overall_survivals <- apply(
    Reduce("+", hazards),
    2,
    function(z) exp(-cumsum(z * newdata[["intlen"]])))
  names(hazards) <- causes_model
  # calculate cif
  hazard           <- hazards[[cause_data]]
  # Value of survival just prior to time-point
  survival         <- overall_survivals - 1e-20
  hps              <- hazard * survival
  cifs             <- apply(hps, 2, function(z) cumsum(z * newdata[["intlen"]]))
  newdata[["cif"]] <- rowMeans(cifs)
 if(ci) {
    newdata[["cif_lower"]] <- apply(cifs, 1, quantile, alpha/2)
    newdata[["cif_upper"]] <- apply(cifs, 1, quantile, 1-alpha/2)
  }

  newdata

}

## Transition Probability Matrix for multi-state data
#' @keywords internal
get_trans_prob <- function(
    newdata,
    # object,
    time_var        = NULL,
    interval_length = "intlen",
    transition      = "transition",
    tend            = "tend",
    cumu_hazard     = "cumu_hazard",
    ...) {


  # interval_length
  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  # transition
  assert_character(transition)
  assert_subset(transition, colnames(newdata))
  # time
  assert_character(tend)
  assert_subset(tend, colnames(newdata))
  # cumu_hazard
  assert_character(cumu_hazard)
  assert_subset(cumu_hazard, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)

  interval_length <- sym(interval_length)
  transition <- sym(transition)
  # FIXME: create time_var variable (="tend"), if time_var == NULL
  # use time_var rather then "tend" below

  # include from and to, to obtain transition probability in multidim array
  newdata <- newdata %>%
    mutate(
      from = as.numeric(gsub("->.*", "", !!transition)),
      to = as.numeric(gsub(".*->", "", !!transition))
    )

  # get unique transitions to build transition matrix
  unique_transition <- newdata %>%
    select(!!transition, "from", "to") %>%
    unique() %>%
    data.frame()
  # get unique time points
  unique_tend <- newdata %>%
    ungroup(!!transition) %>%
    select(!!tend) %>%
    unique() %>%
    data.frame()


  # FIXME could think about writing a separate function that performs the below
  # calculations (or use a function e.g. from etm that does that)
  # transition matrix
  m <- max(newdata$to) + 1 #transition starts at 0, integer of matrix at 1
  M <- array(0, dim=c(max(m), max(m), nrow(unique_transition)))

  # create transition matrices to be used at every time point,
  # multiply matrices with "scalar" alpha_ij_k which is the delta cumu hazard at time t_k for transition i->j

  for (iter in seq_len(nrow(unique_transition))) {
    M[unique_transition$from[iter] + 1, unique_transition$to[iter] + 1,iter] <- 1
    M[unique_transition$from[iter] + 1, unique_transition$from[iter] + 1,iter] <- -1
  }

  # add cumu hazards to dataset
  newdata <- newdata %>%
    group_by(!!transition) %>%
    mutate(delta_cumu_hazard = cumu_hazard - ifelse(is.na(lag(cumu_hazard)), 0, lag(cumu_hazard)))

  # create dA array, to calculate transition probabilities
  alpha <- array(
    rep(0, nrow(unique_tend)*nrow(unique_transition)),
    dim=c(nrow(unique_tend), nrow(unique_transition)))
  I     <- array(
    rep(diag(max(m)), nrow(unique_tend)),
    dim=c( max(m), max(m), nrow(unique_tend)))
  A     <- array(0, dim=c(max(m), max(m), nrow(unique_tend)))
  cum_A <- array(0, dim=c(max(m), max(m), nrow(unique_tend)))

  # calculate differences in hazards
  alpha <- sapply(seq_len(nrow(unique_transition)),
    function(iter) {
      val <- newdata %>%
        ungroup() %>%
        filter(transition == unique_transition[iter,1]) %>%
        arrange(tend)
      val$delta_cumu_hazard
  })

  for (t in seq_len(nrow(unique_tend))) {
    for (trans in seq_len(nrow(unique_transition))) {
      A[,,t] <- A[,,t] + M[,,trans] * alpha[t, trans]
    }
  }

  # # for debugging
  # print(A[,,nrow(unique_tend)])

  # prepare transition probabilities
  A <- I + A

  for (iter in seq_len(nrow(unique_tend))) {
    if (iter == 1) {
      cum_A[,,iter] = A[,,iter]
    } else {
      cum_A[,,iter] = cum_A[,,iter-1] %*% A[,,iter] #use matrix multiplikation
    }
  }

  # transform array so that transition probability can be joined via tend and transition
  tmp <- cbind(
    unique_tend,
    sapply(
      seq_len(nrow(unique_transition)),
      function(row) {
                cum_A[unique_transition$from[row] + 1, unique_transition$to[row] + 1, ]
              }
    )
  )
  # FIXME: replace "tend" with time_var
  colnames(tmp) <- c("tend", as.character(unique_transition$transition))
  trans_prob_df <- tmp %>%
    pivot_longer(
      cols      = c(as.character(unique_transition$transition)),
      names_to  = "transition",
      values_to = "trans_prob") %>%
    mutate(trans_prob = pmin(pmax(.data$trans_prob, 0), 1))

  # join probabilities and return matrix
  newdata <- newdata %>%
    left_join(trans_prob_df, by=c("tend", "transition")) %>%
    select(-one_of(c("delta_cumu_hazard", "from", "to")))

  return(newdata)

}

#' Add transition probabilities
#'
#' @inherit add_hazard
#' @inherit add_cif
#' @export
add_trans_prob <- function(
    newdata
    , object
    , overwrite       = FALSE
    , ci              = FALSE
    , alpha           = 0.05
    , nsim            = 100L
    , time_var        = NULL
    , interval_length = "intlen",
    ...
) {


  if (!overwrite) {
    if ("trans_prob" %in% names(newdata)) {
      stop("Data set already contains 'trans_prob' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
    rm.vars <- intersect(
      c("trans_prob"
        , "trans_lower"
        , "trans_upper"
      ),
      names(newdata))
    newdata <- newdata %>% select(-one_of(rm.vars))
  }

  # add confidence intervals if wanted
  has_cumu = "cumu_hazard" %in% colnames(newdata)
  if (!has_cumu) {
    newdata <- newdata |>
    add_cumu_hazard(object, ci = FALSE)
  }
  if (ci) {
    newdata <- newdata |>
      add_trans_ci(object)
  }

  old_groups <- group_vars(newdata)
  out_df <- newdata |>
    ungroup("transition") |>
    group_split() |>
      map(.f = ~ group_by(.x, transition)) |>
      map(.f = ~ get_trans_prob(.x)) |>
      bind_rows() |>
      group_by(across(all_of(old_groups)))

  if (!has_cumu) {
    out_df[["cumu_hazard"]] <- NULL
  }

  out_df

}

#' helper function for add_trans_ci
#' @keywords internal
get_sim_cumu <- function(newdata, ...) {

  newdata$cumu_hazard <- cumsum(newdata$intlen * newdata$hazard)

  newdata

}

#' Add transition probabilities confidence intervals
#' @keywords internal
add_trans_ci <- function(newdata, object, nsim=100L, alpha=0.05, ...) {

  X     <- predict.gam(object, newdata = newdata, type = "lpmatrix")
  coefs <- coef(object)
  V     <- object$Vp

  # define groups: 1. all grouping variables -> cumu hazards, 2. all but transition -> trans_prob
  groups_array <- group_indices(newdata)
  groups_trans <- newdata %>% ungroup("transition") %>% group_indices()

  sim_coef_mat <- mvtnorm::rmvnorm(nsim, mean = coefs, sigma = V)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z)
    exp(X %*% z))

  # create list with replicated newdata
  nlst <- as.list(replicate(nsim, newdata[,c("tend", "transition", "intlen")], simplify=F))

  # add cumu-hazard in each element and calculate trans_prob with perturbed hazards
  nlst <- lapply(1:nsim, function(i) {
    nlst[[i]] <- cbind(nlst[[i]], hazard = sim_fit_mat[, i]) # add hazard
    # split by group and calculate cumu hazard
    nlst[[i]] <- split(nlst[[i]], groups_array) %>% map_dfr(get_sim_cumu)
    nlst[[i]] <- split(nlst[[i]], groups_trans) %>% map_dfr(get_trans_prob)

    nlst[[i]]
  })

  sim_trans_probs <- do.call(cbind, lapply(nlst, function(df) df$trans_prob))
  newdata$trans_lower <- apply(sim_trans_probs, 1, quantile, probs = alpha / 2)
  newdata$trans_upper <- apply(sim_trans_probs, 1, quantile, probs = 1 - alpha / 2)

  newdata
}


