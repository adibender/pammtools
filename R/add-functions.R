#' Embeds the data set with the specified (relative) term contribution
#'
#' Adds the contribution of a specific term to the
#' linear predictor to the data specified by \code{newdata}.
#' Essentially a wrapper to \code{\link[mgcv]{predict.gam}}, with \code{type="terms"}.
#' Thus most arguments and their documentation below is from \code{\link[mgcv]{predict.gam}}.
#' Shape-constrained additive models fit with \code{\link[scam]{scam}} are
#' supported as well.
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
  ci = TRUE,
  se_mult = 2,
  ...
) {
  assert_data_frame(newdata, all.missing = FALSE)
  assert_character(term, min.chars = 1, any.missing = FALSE, min.len = 1)

  col_ind <- map(term, grep, x = names(object$coefficients)) %>%
    unlist() %>%
    unique() %>%
    sort()
  if (length(col_ind) == 0) {
    stop(
      paste0(
        "No model coefficients matched `term`: ",
        paste(term, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  X <- prep_X(object, newdata, reference, ...)[, col_ind, drop = FALSE]

  newdata[["fit"]] <- unname(drop(X %*% get_coefs(object)[col_ind]))
  if (ci) {
    cov.coefs <- get_Vp(object)[col_ind, col_ind]
    se <- unname(sqrt(rowSums((X %*% cov.coefs) * X)))
    newdata <- newdata %>%
      mutate(
        ci_lower = .data[["fit"]] - se_mult * se,
        ci_upper = .data[["fit"]] + se_mult * se
      )
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
  predict.scam(object, newdata = newdata, type = "lpmatrix", ...)
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
    stop(paste0(
      "Columns in 'reference' but not in 'newdata':",
      paste0(setdiff(names_ref, cnames), collapse = ",")
    ))
  }
  # transform to list if inherits from data frame, so it can be processed
  # in mutate via !!!
  if (inherits(reference, "data.frame")) {
    if (!(nrow(reference) == n_rows || nrow(reference) == 1)) {
      stop(
        "If reference is provided as data frame, number of rows must be
        either 1 or the number of rows in newdata."
      )
    }
    reference <- as.list(reference)
  }

  reference
}

resolve_time_var <- function(time_var, object, newdata) {
  if (is.null(time_var)) {
    is_gam <- inherits(object, "gam") || inherits(object, "scam")
    time_var <- if (is_gam) "tend" else "interval"
  } else {
    assert_string(time_var)
  }
  assert_choice(time_var, colnames(newdata))

  time_var
}

#' Add predicted (cumulative) hazard to data set
#'
#' Add (cumulative) hazard based on the provided data set and model.
#' If \code{ci=TRUE} confidence intervals (CI) are also added. Their width can
#' be controlled via the \code{se_mult} argument. The method by which the
#' CI are calculated can be specified by \code{ci_type}.
#' This is a wrapper around
#' \code{\link[mgcv]{predict.gam}}. When \code{reference} is specified, the
#' (log-)hazard ratio is calculated. In addition to models fit with
#' \code{\link[mgcv]{gam}}/\code{\link[mgcv]{bam}} or \code{\link[stats]{glm}},
#' shape-constrained additive models fit with \code{\link[scam]{scam}} are
#' supported (e.g., for monotone baseline hazards). For \code{scam} models all
#' calculations (including delta-method and simulation based confidence
#' intervals) are based on the re-parametrized coefficients and their
#' covariance matrix, i.e., on the same normal approximation that underlies
#' the standard errors reported by \code{scam} itself.
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
#' for details and empirical evaluation. For \code{ci_type = "sim"}, interval
#' bounds are empirical quantiles (type 6, see \code{\link[stats]{quantile}})
#' of \code{nsim} posterior draws (default \code{nsim = 100L}, passed via
#' \code{...}). Type-6 quantiles avoid the systematic inward bias that the
#' \code{\link[stats]{quantile}} default (type 7) exhibits for small
#' \code{nsim}, but at the default \code{nsim = 100} the bounds are estimated
#' from few tail draws and thus noisy; increase \code{nsim} (e.g., to 500 or
#' more) for more stable interval bounds. Very small \code{nsim}
#' (\code{nsim < 2 / alpha - 1}, i.e., below 39 for \code{alpha = 0.05})
#' cannot achieve the nominal level at all.
#' @param se_mult Factor by which standard errors are multiplied for calculating
#' the confidence intervals.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to \code{FALSE}. If \code{TRUE}, columns with names
#' \code{c("hazard", "se", "lower", "upper")} will be overwritten.
#' @param time_var Name of the variable used for the baseline hazard. Defaults
#'   to \code{"tend"}.
#'
#' @details
#' When computing cumulative hazards or survival probabilities across groups,
#' the input data must be grouped via \code{group_by()} prior to calling
#' \code{add_cumu_hazard()} or \code{add_surv_prob()}. Omitting
#' \code{group_by()} will not produce an error or warning but will return
#' silently incorrect results, as the cumulative hazard will be accumulated
#' over the entire dataset rather than within each group.
#' See the \href{https://adibender.github.io/pammtools/articles/convenience.html#cumulative-hazard}{workflow vignette}
#' for a worked example.
#'
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
  type = c("response", "link"),
  ci = TRUE,
  se_mult = 2,
  ci_type = c("default", "delta", "sim"),
  overwrite = FALSE,
  time_var = NULL,
  ...
) {
  if (!overwrite) {
    if ("hazard" %in% names(newdata)) {
      stop(
        "Data set already contains 'hazard' column.
        Set `overwrite=TRUE` to overwrite"
      )
    }
  } else {
    rm.vars <- intersect(
      c("hazard", "se", "ci_lower", "ci_upper"),
      names(newdata)
    )
    newdata <- newdata %>% select(-one_of(rm.vars))
  }

  get_hazard(
    object,
    newdata,
    reference = reference,
    ci = ci,
    type = type,
    se_mult = se_mult,
    ci_type = ci_type,
    time_var = time_var,
    ...
  )
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
  ci = TRUE,
  type = c("response", "link"),
  ci_type = c("default", "delta", "sim"),
  time_var = NULL,
  se_mult = 2,
  ...
) {
  assert_data_frame(newdata, all.missing = FALSE)
  assert_multi_class(object, classes = c("glm", "scam"))
  type <- match.arg(type)
  ci_type <- match.arg(ci_type)

  time_var <- resolve_time_var(time_var, object, newdata)

  # throw warning or error if evaluation time points/intervals do not correspond
  # to evaluation time-points/intervals do not correspond to the ones used for
  # estimation
  #warn_about_new_time_points(object, newdata, time_var)

  X <- prep_X(object, newdata, reference, ...)
  coefs <- get_coefs(object)
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
#' @export
add_cumu_hazard <- function(newdata, object, ...) {
  UseMethod("add_cumu_hazard", object)
}

#' @rdname add_hazard
#' @inheritParams add_hazard
#' @param interval_length The variable in newdata containing the interval lengths.
#' Can be either bare unquoted variable name or character. Defaults to \code{"intlen"}.
#' @param boundary Logical. If \code{TRUE} (default), a boundary row at
#'   \code{time = 0} with cumulative hazard \code{0} is prepended (per group),
#'   so that cumulative hazards start at the natural origin (consistent with
#'   \code{\link{add_surv_prob}}, \code{\link{add_cif}} and
#'   \code{\link{add_trans_prob}}).
#' @importFrom dplyr bind_cols
#' @seealso \code{\link[mgcv]{predict.gam}},
#' \code{\link[pammtools]{add_surv_prob}}
#' @export
add_cumu_hazard.default <- function(
  newdata,
  object,
  ci = TRUE,
  se_mult = 2,
  overwrite = FALSE,
  time_var = NULL,
  interval_length = "intlen",
  boundary = TRUE,
  ...
) {
  interval_length <- quo_name(enquo(interval_length))

  if (!overwrite) {
    if ("cumu_hazard" %in% names(newdata)) {
      stop(
        "Data set already contains 'hazard' column.
        Set `overwrite=TRUE` to overwrite"
      )
    }
  } else {
    rm.vars <- intersect(
      c("cumu_hazard", "cumu_lower", "cumu_upper"),
      names(newdata)
    )
    newdata <- newdata %>% select(-one_of(rm.vars))
  }

  time_var <- resolve_time_var(time_var, object, newdata)
  # The boundary is a continuous-time row at time == 0 on the resolved time axis
  # (`time_var`, i.e. "tend" or a renamed continuous time variable). It is only
  # well defined for models predicted on that axis (gam/scam/pamm). Interval-
  # factor models (glm/PEM) have no time == 0 level, so no boundary is added --
  # consistent with add_surv_prob/add_cif, whose boundary also targets the
  # continuous time axis.
  boundary <- boundary &&
    (inherits(object, "gam") || inherits(object, "scam"))

  if (boundary) {
    newdata <- drop_cumulative_boundary(newdata, time_var)
  }

  joindata <- reconstruct_cutpoints(newdata, object, time_var, interval_length)
  joindata <- get_cumu_hazard(
    joindata,
    object,
    ci = ci,
    se_mult = se_mult,
    time_var = time_var,
    interval_length = interval_length,
    ...
  )

  out <- suppressMessages(
    newdata %>% left_join(joindata)
  )
  out <- restore_prediction_attrs(out, newdata)

  if (boundary) {
    out <- add_cumulative_boundary(
      out,
      time_var = time_var,
      values = c(cumu_hazard = 0, cumu_lower = 0, cumu_upper = 0),
      interval_length = interval_length
    )
  }

  out
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
  ci = TRUE,
  ci_type = c("default", "delta", "sim"),
  time_var = NULL,
  se_mult = 2,
  interval_length = "intlen",
  nsim = 100L,
  ...
) {
  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)
  assert_multi_class(object, classes = c("glm", "scam"))

  ci_type <- match.arg(ci_type)

  interval_length_name <- interval_length
  interval_length <- sym(interval_length)

  mutate_args <- list(
    cumu_hazard = quo(cumsum(
      .data[["hazard"]] *
        (!!interval_length)
    ))
  )
  haz_vars_in_data <- map(
    c("hazard", "se", "ci_lower", "ci_upper"),
    ~ grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)
  ) %>%
    flatten_chr()
  vars_exclude <- c("hazard")

  if (ci) {
    if (ci_type == "default" | ci_type == "delta") {
      vars_exclude <- c(vars_exclude, "se", "ci_lower", "ci_upper")
      newdata <- get_hazard(
        object,
        newdata,
        type = "response",
        ci = ci,
        ci_type = ci_type,
        time_var = time_var,
        se_mult = se_mult,
        ...
      )
      if (ci_type == "default") {
        mutate_args <- mutate_args %>%
          append(list(
            cumu_lower = quo(cumsum(.data[["ci_lower"]] * (!!interval_length))),
            cumu_upper = quo(cumsum(.data[["ci_upper"]] * (!!interval_length)))
          ))
      } else {
        # ci delta rule
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(
            add_delta_ci_cumu,
            object = object,
            se_mult = se_mult,
            interval_length = interval_length_name,
            ...
          )
      }
    } else {
      if (ci_type == "sim") {
        newdata <- get_hazard(
          object,
          newdata,
          type = "response",
          ci = FALSE,
          time_var = time_var,
          ...
        )
        sim_coef_mat <- sample_coefs(object, nsim)
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(
            get_sim_ci_cumu,
            object = object,
            nsim = nsim,
            interval_length = interval_length_name,
            sim_coef_mat = sim_coef_mat,
            ...
          )
      }
    }
  } else {
    newdata <-
      get_hazard(
        object,
        newdata,
        type = "response",
        ci = ci,
        ci_type = ci_type,
        time_var = time_var,
        se_mult = se_mult,
        ...
      )
  }
  newdata <- newdata %>%
    mutate(!!!mutate_args)

  vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
  if (length(vars_exclude) != 0) {
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
#' @details
#' When computing cumulative hazards or survival probabilities across groups,
#' the input data must be grouped via \code{group_by()} prior to calling
#' \code{add_cumu_hazard()} or \code{add_surv_prob()}. Omitting
#' \code{group_by()} will not produce an error or warning but will return
#' silently incorrect results, as the cumulative hazard will be accumulated
#' over the entire dataset rather than within each group.
#' See the \href{https://adibender.github.io/pammtools/articles/convenience.html#cumulative-hazard}{workflow vignette}
#' for a worked example.
#'
#' The returned data contains one boundary row per group at \code{time_var = 0}
#' for plotting cumulative quantities from the time origin. On this row,
#' \code{surv_prob = 1}; if confidence intervals are requested,
#' \code{surv_lower = surv_upper = 1}. If an interval-length column is present,
#' it is set to \code{0} on the boundary row.
#'
#' @inherit add_cumu_hazard
#' @examples
#' ped <- tumor[1:50,] %>% as_ped(Surv(days, status)~ age)
#' pam <- mgcv::gam(ped_status ~ s(tend)+age, data=ped, family=poisson(), offset=offset)
#' ped_info(ped) %>% add_surv_prob(pam, ci=TRUE)
#' @export
add_surv_prob <- function(newdata, object, ...) {
  UseMethod("add_surv_prob", object)
}

#' @rdname add_surv_prob
#' @export
add_surv_prob.default <- function(
  newdata,
  object,
  ci = TRUE,
  se_mult = 2,
  overwrite = FALSE,
  time_var = NULL,
  interval_length = "intlen",
  boundary = TRUE,
  ...
) {
  interval_length <- quo_name(enquo(interval_length))
  time_var <- resolve_time_var(time_var, object, newdata)
  # The boundary is a continuous-time row at time == 0 (survival 1); see
  # add_cumu_hazard(). Only added for models predicted on the continuous time
  # axis (gam/scam/pamm), not for interval-factor models (glm/PEM).
  boundary <- boundary &&
    (inherits(object, "gam") || inherits(object, "scam"))

  if (!overwrite) {
    if ("surv_prob" %in% names(newdata)) {
      stop(
        "Data set already contains 'surv_prob' column.
        Set `overwrite=TRUE` to overwrite"
      )
    }
  } else {
    rm.vars <- intersect(
      c("surv_prob", "surv_lower", "surv_upper"),
      names(newdata)
    )
    newdata <- newdata %>% select(-one_of(rm.vars))
  }

  if (boundary) {
    newdata <- drop_cumulative_boundary(newdata, time_var)
  }

  if (!interval_length %in% colnames(newdata)) {
    newdata <- reconstruct_intlen(
      newdata,
      time_var = time_var,
      interval_length = interval_length
    )
  }

  out <- get_surv_prob(
    newdata,
    object,
    ci = ci,
    se_mult = se_mult,
    time_var = time_var,
    interval_length = interval_length,
    ...
  )
  out <- restore_prediction_attrs(out, newdata)

  if (boundary) {
    out <- add_cumulative_boundary(
      out,
      time_var = time_var,
      values = c(surv_prob = 1, surv_lower = 1, surv_upper = 1),
      interval_length = interval_length
    )
  }

  out
}


#' Calculate survival probabilities
#'
#' @inheritParams add_surv_prob
#' @keywords internal
get_surv_prob <- function(
  newdata,
  object,
  ci = TRUE,
  ci_type = c("default", "delta", "sim"),
  se_mult = 2L,
  time_var = NULL,
  interval_length = "intlen",
  nsim = 100L,
  ...
) {
  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)
  assert_multi_class(object, classes = c("glm", "scam"))

  ci_type <- match.arg(ci_type)

  interval_length_name <- interval_length
  interval_length <- sym(interval_length)

  mutate_args <- list(
    surv_prob = quo(exp(
      -cumsum(
        .data[["hazard"]] *
          (!!interval_length)
      )
    ))
  )
  haz_vars_in_data <- map(
    c("hazard", "se", "ci_lower", "ci_upper"),
    ~ grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)
  ) %>%
    flatten_chr()
  vars_exclude <- c("hazard")

  if (ci) {
    if (ci_type == "default" | ci_type == "delta") {
      vars_exclude <- c(vars_exclude, "se", "ci_lower", "ci_upper")
      newdata <- get_hazard(
        object,
        newdata,
        type = "response",
        ci = ci,
        ci_type = ci_type,
        time_var = time_var,
        se_mult = se_mult,
        ...
      )
      if (ci_type == "default") {
        mutate_args <- mutate_args %>%
          append(list(
            surv_upper = quo(exp(
              -cumsum(.data[["ci_lower"]] * (!!interval_length))
            )),
            surv_lower = quo(exp(
              -cumsum(.data[["ci_upper"]] * (!!interval_length))
            ))
          ))
      } else {
        # ci delta rule
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(
            add_delta_ci_surv,
            object = object,
            se_mult = se_mult,
            interval_length = interval_length_name,
            ...
          )
      }
    } else {
      if (ci_type == "sim") {
        newdata <- get_hazard(
          object,
          newdata,
          type = "response",
          ci = FALSE,
          time_var = time_var,
          ...
        )
        sim_coef_mat <- sample_coefs(object, nsim)
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(
            get_sim_ci_surv,
            object = object,
            nsim = nsim,
            interval_length = interval_length_name,
            sim_coef_mat = sim_coef_mat,
            ...
          )
      }
    }
  } else {
    newdata <-
      get_hazard(
        object = object,
        newdata,
        type = "response",
        ci = FALSE,
        time_var = time_var,
        ...
      )
  }

  newdata <- newdata %>%
    mutate(!!!mutate_args)

  vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
  if (length(vars_exclude) != 0) {
    newdata <- newdata %>% select(-one_of(vars_exclude))
  }

  newdata
}

drop_cumulative_boundary <- function(newdata, time_var) {
  boundary_rows <- !is.na(newdata[[time_var]]) & newdata[[time_var]] == 0
  if (!any(boundary_rows)) {
    return(newdata)
  }

  out <- newdata %>%
    filter(is.na(.data[[time_var]]) | .data[[time_var]] != 0)

  restore_prediction_attrs(out, newdata)
}

add_cumulative_boundary <- function(
  newdata,
  time_var,
  values,
  interval_length = NULL
) {
  if (!nrow(newdata)) {
    return(newdata)
  }

  old_groups <- group_vars(newdata)
  out <- newdata %>%
    slice(1, .preserve = TRUE) %>%
    set_cumulative_boundary_values(
      time_var = time_var,
      values = values,
      interval_length = interval_length
    ) %>%
    bind_rows(newdata)

  if (length(old_groups)) {
    out <- out %>% group_by(across(all_of(old_groups)))
  }
  out <- out %>% arrange(.data[[time_var]], .by_group = TRUE)

  restore_prediction_attrs(out, newdata)
}

set_cumulative_boundary_values <- function(
  newdata,
  time_var,
  values,
  interval_length = NULL
) {
  newdata[[time_var]][] <- 0

  if (!is.null(interval_length) && interval_length %in% colnames(newdata)) {
    newdata[[interval_length]][] <- 0
  }
  if ("tstart" %in% colnames(newdata)) {
    newdata[["tstart"]] <- 0
  }
  if ("intmid" %in% colnames(newdata)) {
    newdata[["intmid"]] <- 0
  }
  if ("interval" %in% colnames(newdata)) {
    is.na(newdata[["interval"]]) <- TRUE
  }
  if ("offset" %in% colnames(newdata)) {
    newdata[["offset"]][] <- NA_real_
  }

  value_names <- intersect(names(values), colnames(newdata))
  for (value_name in value_names) {
    newdata[[value_name]] <- values[[value_name]]
  }

  newdata
}

restore_prediction_attrs <- function(newdata, template) {
  structural_attrs <- c("names", "row.names", "class", "groups")
  custom_attrs <- setdiff(names(attributes(template)), structural_attrs)
  for (attr_name in custom_attrs) {
    attr(newdata, attr_name) <- attr(template, attr_name)
  }
  if (inherits(template, "ped") && !inherits(newdata, "ped")) {
    class(newdata) <- c("ped", class(newdata))
  }

  newdata
}

add_ci <- function(
  newdata,
  object,
  X,
  type = c("response", "link"),
  se_mult = 2,
  ci_type = c("default", "delta", "sim"),
  nsim = 100,
  ...
) {
  ci_type <- match.arg(ci_type)

  V <- get_Vp(object)
  se <- unname(sqrt(rowSums((X %*% V) * X)))
  newdata$se <- se
  if (type == "link") {
    newdata <- newdata %>%
      mutate(
        ci_lower = .data[["hazard"]] - se_mult * .data[["se"]],
        ci_upper = .data[["hazard"]] + se_mult * .data[["se"]]
      )
  }

  if (type != "link") {
    if (ci_type == "default") {
      newdata <- newdata %>%
        mutate(
          ci_lower = exp(.data[["hazard"]] - se_mult * .data[["se"]]),
          ci_upper = exp(.data[["hazard"]] + se_mult * .data[["se"]])
        )
    } else {
      if (ci_type == "delta") {
        newdata <- split(newdata, group_indices(newdata)) %>%
          map_dfr(add_delta_ci, object = object, se_mult = se_mult, ...)
      } else {
        if (ci_type == "sim") {
          sim_coef_mat <- sample_coefs(object, nsim)
          newdata <- split(newdata, group_indices(newdata)) %>%
            map_dfr(
              get_sim_ci,
              object = object,
              nsim = nsim,
              sim_coef_mat = sim_coef_mat,
              ...
            )
        }
      }
    }
  }
  newdata
}

add_delta_ci <- function(newdata, object, se_mult = 2, ...) {
  X <- make_X(object, newdata, ...)
  V <- get_Vp(object)

  Jacobi <- diag(exp(newdata$hazard)) %*% X
  newdata %>%
    mutate(
      se = sqrt(rowSums((Jacobi %*% V) * Jacobi)),
      ci_lower = exp(.data[["hazard"]]) - .data[["se"]] * se_mult,
      ci_upper = exp(.data[["hazard"]]) + .data[["se"]] * se_mult
    )
}

add_delta_ci_cumu <- function(
  newdata,
  object,
  se_mult = 2,
  interval_length = "intlen",
  ...
) {
  X <- make_X(object, newdata, ...)
  V <- get_Vp(object)
  intlen <- newdata[[interval_length]]

  Delta <- lower.tri(diag(nrow(X)), diag = TRUE) %*% diag(intlen)
  Jacobi <- diag(newdata$hazard) %*% X
  LHS <- Delta %*% Jacobi
  newdata %>%
    mutate(
      se = sqrt(rowSums((LHS %*% V) * LHS)),
      cumu_lower = cumsum(intlen * .data[["hazard"]]) - .data[["se"]] * se_mult,
      cumu_upper = cumsum(intlen * .data[["hazard"]]) + .data[["se"]] * se_mult
    )
}

add_delta_ci_surv <- function(
  newdata,
  object,
  se_mult = 2,
  interval_length = "intlen",
  ...
) {
  X <- make_X(object, newdata, ...)
  V <- get_Vp(object)
  intlen <- newdata[[interval_length]]

  Delta <- lower.tri(diag(nrow(X)), diag = TRUE) %*% diag(intlen)
  Jacobi <- diag(newdata$hazard) %*% X
  LHS <- -diag(exp(-rowSums(Delta %*% diag(newdata$hazard)))) %*%
    (Delta %*% Jacobi)
  newdata %>%
    mutate(
      se = sqrt(rowSums((LHS %*% V) * LHS)),
      surv_lower = exp(-cumsum(.data[["hazard"]] * intlen)) -
        .data[["se"]] * se_mult,
      surv_upper = exp(-cumsum(.data[["hazard"]] * intlen)) +
        .data[["se"]] * se_mult
    )
}

#' Calculate simulation based confidence intervals
#'
#' @keywords internal
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef
get_sim_ci <- function(
  newdata,
  object,
  alpha = 0.05,
  nsim = 100L,
  sim_coef_mat = NULL,
  ...
) {
  X <- make_X(object, newdata, ...)

  if (is.null(sim_coef_mat)) {
    sim_coef_mat <- sample_coefs(object, nsim)
  }
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))

  newdata$ci_lower <- apply(
    sim_fit_mat,
    1,
    quantile,
    probs = alpha / 2,
    type = 6
  )
  newdata$ci_upper <- apply(
    sim_fit_mat,
    1,
    quantile,
    probs = 1 - alpha / 2,
    type = 6
  )

  newdata
}


get_sim_ci_cumu <- function(
  newdata,
  object,
  alpha = 0.05,
  nsim = 100L,
  interval_length = "intlen",
  sim_coef_mat = NULL,
  ...
) {
  X <- make_X(object, newdata, ...)
  intlen <- newdata[[interval_length]]

  if (is.null(sim_coef_mat)) {
    sim_coef_mat <- sample_coefs(object, nsim)
  }
  sim_fit_mat <- apply(
    sim_coef_mat,
    1,
    function(z) cumsum(intlen * exp(X %*% z))
  )

  newdata$cumu_lower <- apply(
    sim_fit_mat,
    1,
    quantile,
    probs = alpha / 2,
    type = 6
  )
  newdata$cumu_upper <- apply(
    sim_fit_mat,
    1,
    quantile,
    probs = 1 - alpha / 2,
    type = 6
  )

  newdata
}

get_sim_ci_surv <- function(
  newdata,
  object,
  alpha = 0.05,
  nsim = 100L,
  interval_length = "intlen",
  sim_coef_mat = NULL,
  ...
) {
  X <- make_X(object, newdata, ...)
  intlen <- newdata[[interval_length]]

  if (is.null(sim_coef_mat)) {
    sim_coef_mat <- sample_coefs(object, nsim)
  }
  sim_fit_mat <- apply(
    sim_coef_mat,
    1,
    function(z) exp(-cumsum(intlen * exp(X %*% z)))
  )

  newdata$surv_lower <- apply(
    sim_fit_mat,
    1,
    quantile,
    probs = alpha / 2,
    type = 6
  )
  newdata$surv_upper <- apply(
    sim_fit_mat,
    1,
    quantile,
    probs = 1 - alpha / 2,
    type = 6
  )

  newdata
}


## Cumulative Incidence Function (CIF) for competing risks data

#' Add cumulative incidence function to data
#'
#' @inheritParams add_hazard
#' @param alpha The alpha level for confidence/credible intervals.
#' @param nsim Number of simulations (draws from posterior of estimated coefficients)
#' on which estimation of CIFs and their confidence/credible intervals will be
#' based on. Interval bounds are empirical type-6 quantiles of the \code{nsim}
#' draws; larger values of \code{nsim} yield more stable interval bounds.
#' @param cause_var Character. Column name of the 'cause' variable.
#' @param interval_length \code{Character}, defaults to \code{"intlen"}.
#'   contains the interval length in `newdata`.
#'
#' @details
#' When computing cumulative incidence for multiple groups, the input data must
#' be grouped via \code{group_by()} before calling this function. Omitting
#' \code{group_by()} will not produce an error or warning but will return
#' silently incorrect results, as the cumulative incidence will be accumulated
#' over the entire dataset rather than within each group.
#'
#' The returned data contains one boundary row per group at \code{time_var = 0}
#' for plotting cumulative incidence from the time origin. On this row,
#' \code{cif = 0}; if confidence intervals are requested,
#' \code{cif_lower = cif_upper = 0}. If an interval-length column is present,
#' it is set to \code{0} on the boundary row. \code{add_cumu_hazard()} keeps the
#' original prediction grid and does not add this plotting boundary row.
#'
#' @examples
#' \donttest{
#' if (require("etm")) {
#'   data("fourD", package = "etm")
#'   ped_stacked <- fourD |>
#'     dplyr::select(-medication, -treated) |>
#'     as_ped(Surv(time, status) ~., id = "id") |>
#'     dplyr::mutate(cause = as.factor(cause))
#'   pam <- pamm(
#'     ped_status ~ s(tend, by = cause) + sex + sex:cause + age + age:cause,
#'     data = ped_stacked)
#'   ped_stacked |>
#'     make_newdata(tend = unique(tend), cause = unique(cause)) |>
#'     group_by(cause) |>
#'     add_cif(pam)
#' }
#' }
#'
#' @export
add_cif <- function(
  newdata,
  object,
  ...
) {
  UseMethod("add_cif", object)
}


#' @rdname add_cif
#' @export
add_cif.default <- function(
  newdata,
  object,
  ci = TRUE,
  overwrite = FALSE,
  alpha = 0.05,
  nsim = 500L,
  cause_var = "cause",
  time_var = NULL,
  interval_length = "intlen",
  ...
) {
  interval_length <- quo_name(enquo(interval_length))
  time_var <- resolve_time_var(time_var, object, newdata)

  if (!overwrite) {
    if ("cif" %in% names(newdata)) {
      stop(
        "Data set already contains 'cif' column.
        Set `overwrite=TRUE` to overwrite"
      )
    }
  } else {
    rm.vars <- intersect(
      c("cif", "cif_lower", "cif_upper"),
      names(newdata)
    )
    newdata <- newdata %>% select(-one_of(rm.vars))
  }

  newdata <- drop_cumulative_boundary(newdata, time_var)

  joindata <- reconstruct_cutpoints(newdata, object, time_var, interval_length)

  coefs <- get_coefs(object)
  V <- get_Vp(object)
  sim_coef_mat <- if (!ci) {
    matrix(coefs, nrow = 1)
  } else {
    sample_coefs(object, nsim)
  }

  joindata <- map_dfr(
    split(joindata, group_indices(joindata)),
    ~ get_cif(
      newdata = .x,
      object = object,
      ci = ci,
      alpha = alpha,
      nsim = nsim,
      cause_var = cause_var,
      coefs = coefs,
      V = V,
      sim_coef_mat = sim_coef_mat,
      time_var = time_var,
      interval_length = interval_length,
      ...
    )
  )

  out <- suppressMessages(
    newdata %>% left_join(joindata)
  )
  out <- restore_prediction_attrs(out, newdata)

  add_cumulative_boundary(
    out,
    time_var = time_var,
    values = c(cif = 0, cif_lower = 0, cif_upper = 0),
    interval_length = interval_length
  )
}

#' Calculate CIF for one cause
#'
#' Internal generic dispatching CIF calculation based on the model class.
#'
#' @param newdata A data frame of new observations, typically created via
#'   [make_newdata()].
#' @param object A fitted model object. The method is dispatched on this
#'   argument.
#' @param ... Additional arguments passed to the respective method.
#'
#' @return A data frame with CIF estimates appended.
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
  interval_length = "intlen",
  alpha,
  nsim,
  cause_var,
  coefs,
  V,
  sim_coef_mat,
  ...
) {
  time_var <- resolve_time_var(time_var, object, newdata)
  assert_string(interval_length)
  assert_choice(interval_length, colnames(newdata))

  # causes_model <- as.factor(object$attr_ped$risks)
  causes_model <- as.factor(levels(newdata[[cause_var]]))
  cause_data <- unique(newdata[[cause_var]])

  if (length(cause_data) > 1) {
    stop("Did you forget to group by cause?")
  }

  hazards <- map(
    causes_model,
    ~ {
      .df <- mutate(newdata, cause = .x) %>%
        arrange(.data[[time_var]], .by_group = TRUE)
      X <- make_X(object, .df)
      apply(sim_coef_mat, 1, function(z) exp(X %*% z))
    }
  )
  overall_survivals <- apply(
    Reduce("+", hazards),
    2,
    function(z) exp(-cumsum(z * newdata[[interval_length]]))
  )
  names(hazards) <- causes_model
  # calculate cif
  hazard <- hazards[[cause_data]]

  survival <- rbind(1, head(overall_survivals, -1))

  # total hazard h_j
  total_hazard <- Reduce("+", hazards)

  # interval lengths
  dt <- newdata[[interval_length]]

  # CIF increment using exact formula
  cif_increments <- (hazard / total_hazard) *
    survival *
    (1 - exp(-total_hazard * dt))

  # cumulative CIF
  cifs <- apply(cif_increments, 2, cumsum)

  newdata[["cif"]] <- pmin(pmax(rowMeans(cifs), 0), 1)
  if (ci) {
    newdata[["cif_lower"]] <- pmin(
      pmax(
        apply(cifs, 1, quantile, probs = alpha / 2, na.rm = TRUE, type = 6),
        0
      ),
      1
    )
    newdata[["cif_upper"]] <- pmin(
      pmax(
        apply(cifs, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 6),
        0
      ),
      1
    )
  }

  newdata
}

## Transition Probability Matrix for multi-state data
#' @keywords internal
get_trans_prob <- function(
  newdata,
  time_var = "tend",
  interval_length = "intlen",
  transition = "transition",
  tend = "tend",
  cumu_hazard = "cumu_hazard",
  sep = "->",
  ...
) {
  if (is.null(time_var)) {
    time_var <- tend
  }

  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_character(transition)
  assert_subset(transition, colnames(newdata))
  assert_character(time_var)
  assert_subset(time_var, colnames(newdata))
  assert_character(cumu_hazard)
  assert_subset(cumu_hazard, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)

  transition_var <- transition
  unique_transition <- transition_state_table(
    newdata[[transition_var]],
    sep = sep
  )
  names(unique_transition)[1] <- transition_var
  unique_tend <- sort(unique(newdata[[time_var]]))

  n_trans <- nrow(unique_transition)
  n_t <- length(unique_tend)

  m <- max(unique_transition$to_int)
  M <- array(0, dim = c(m, m, n_trans))

  idx <- cbind(
    unique_transition$from_int,
    unique_transition$to_int,
    seq_len(n_trans)
  )
  M[idx] <- 1
  M[cbind(
    unique_transition$from_int,
    unique_transition$from_int,
    seq_len(n_trans)
  )] <- -1

  newdata$delta_cumu_hazard <- 0
  split_idx <- split(
    seq_len(nrow(newdata)),
    as.character(newdata[[transition_var]])
  )
  for (idx_group in split_idx) {
    idx_group <- idx_group[order(newdata[[time_var]][idx_group])]
    ch <- newdata[[cumu_hazard]][idx_group]
    newdata$delta_cumu_hazard[idx_group] <- ch - c(0, ch[-length(ch)])
  }

  alpha <- matrix(0, nrow = n_t, ncol = n_trans)
  t_idx <- match(newdata[[time_var]], unique_tend)
  tr_idx <- match(
    as.character(newdata[[transition_var]]),
    as.character(unique_transition[[transition_var]])
  )
  alpha[cbind(t_idx, tr_idx)] <- newdata$delta_cumu_hazard

  M_mat <- matrix(M, nrow = m * m, ncol = n_trans)
  I <- diag(m)

  A_flat <- matrix(NA_real_, nrow = m * m, ncol = n_t)
  for (t in seq_len(n_t)) {
    A_flat[, t] <- as.vector(matrix(M_mat %*% alpha[t, ], m, m) + I)
  }

  mat_mult_flat <- function(x, y) {
    as.vector(matrix(x, nrow = m, ncol = m) %*% matrix(y, nrow = m, ncol = m))
  }
  cum_A_list <- Reduce(
    f = mat_mult_flat,
    x = split(A_flat, col(A_flat)),
    accumulate = TRUE
  )

  cum_A <- array(NA_real_, dim = c(m, m, n_t))
  for (t in seq_len(n_t)) {
    cum_A[,, t] <- matrix(cum_A_list[[t]], m, m)
  }

  prob_mat <- matrix(0, nrow = n_t, ncol = n_trans)
  for (k in seq_len(n_trans)) {
    prob_mat[, k] <- cum_A[
      unique_transition$from_int[k],
      unique_transition$to_int[k],
    ]
  }
  colnames(prob_mat) <- as.character(unique_transition[[transition_var]])

  t_idx <- match(newdata[[time_var]], unique_tend)
  tr_idx <- match(
    as.character(newdata[[transition_var]]),
    as.character(unique_transition[[transition_var]])
  )
  newdata$trans_prob <- pmin(pmax(prob_mat[cbind(t_idx, tr_idx)], 0), 1)

  newdata$delta_cumu_hazard <- NULL
  # `make_newdata()` averages numeric covariates, so any `from` / `to` columns
  # in `newdata` carry meaningless mean values. The integer codes used
  # internally come from `transition_state_table()`, so drop the noisy columns
  # from the user-facing output.
  newdata[["from"]] <- NULL
  newdata[["to"]] <- NULL

  list(
    data = newdata,
    matrix = cum_A
  )
}

## Build a per-transition state table from a vector of transition labels.
##
## Splits each unique transition label on `sep` into ("from", "to") halves and
## assigns integer codes for downstream matrix indexing. Two encoding modes:
##
##   - Numeric: if both halves of every label parse as integers, the integers
##     are used directly (with a +1 shift applied if the minimum is 0, so that
##     state codes start at 1 as required by the M / cum_A indexing).
##   - Categorical: otherwise, the union of unique "from"/"to" strings is
##     collected in sorted order and mapped to consecutive integers 1..K.
##
## Returns a data.frame with columns: transition (chr), from_int, to_int.
#' @keywords internal
transition_state_table <- function(transitions, sep = "->") {
  labs <- as.character(unique(transitions))
  if (any(is.na(labs)) || !length(labs)) {
    stop("`transition` contains no usable levels.")
  }
  pattern <- paste0("^([^", sep, "]+)", sep, "([^", sep, "]+)$")
  parts <- regmatches(labs, regexec(pattern, labs, perl = FALSE))
  bad <- vapply(parts, length, integer(1)) != 3L
  if (any(bad)) {
    stop(
      "Could not split transition label(s) ",
      paste(shQuote(labs[bad]), collapse = ", "),
      " on separator '",
      sep,
      "'. ",
      "Each transition must look like 'from",
      sep,
      "to'."
    )
  }
  from_str <- vapply(parts, `[[`, character(1), 2L)
  to_str <- vapply(parts, `[[`, character(1), 3L)

  from_int_try <- suppressWarnings(as.integer(from_str))
  to_int_try <- suppressWarnings(as.integer(to_str))
  numeric_mode <- !anyNA(from_int_try) &&
    !anyNA(to_int_try) &&
    all(from_str == as.character(from_int_try)) &&
    all(to_str == as.character(to_int_try))

  if (numeric_mode) {
    state_min <- min(c(from_int_try, to_int_try))
    shift <- if (state_min == 0L) 1L else 0L
    from_int <- from_int_try + shift
    to_int <- to_int_try + shift
    if (min(c(from_int, to_int)) < 1L) {
      stop(
        "Negative state indices are not supported in numeric transition labels."
      )
    }
  } else {
    states <- sort(unique(c(from_str, to_str)))
    from_int <- match(from_str, states)
    to_int <- match(to_str, states)
  }

  data.frame(
    transition = labs,
    from_int = from_int,
    to_int = to_int,
    stringsAsFactors = FALSE
  )
}

#' Add transition probabilities
#' @description
#' \code{add_trans_prob} adds transition probabilities on the provided data set and model.
#' Optionally, confidence intervals (CI) are added if \code{ci=TRUE}.
#' The function builds on cumulative hazards \code{cumu_hazard} and \code{mgcv::gam} models.
#'
#' @param newdata A data frame or list containing the values of the model
#' covariates at which predictions are required. If this is not provided then
#' predictions corresponding to the original data are returned. If newdata is
#' provided then it should contain all the variables needed for prediction:
#' a warning is generated if not. See details for use with
#' \link[mgcv]{linear.functional.terms}.
#' @param object A fitted \code{gam} object as produced by \code{mgcv::gam}
#' @param overwrite Should transition probability columns be overwritten if
#' already present in the data set? Defaults to \code{FALSE}.
#' If \code{TRUE}, columns with names \code{c("trans_prob", "trans_upper", "trans_lower")}
#' will be overwritten.
#' @param ci \code{Logical}, defaults to \code{TRUE}. Decides if confidence
#' intervals for transition probabilities are calculated.
#' @param alpha Sets the confidence intervals' \eqn{\alpha} level, Defaults to \code{0.05}
#' @param nsim Sets the number of iterations for simulated confidence intervals.
#' Defaults to \code{100L}. Interval bounds are empirical type-6 quantiles of
#' the \code{nsim} draws; larger values of \code{nsim} yield more stable
#' interval bounds.
#' @param time_var Name of the variable used for the baseline hazard. Defaults
#'   to \code{"tend"}.
#' @param interval_length \code{Character}, defaults to \code{"intlen"}.
#'   contains the interval length in `newdata`.
#' @param transition \code{Character}, defaults to \code{"transition"}.
#'   contains the transition labels in `newdata`.
#' @param ... Further arguments passed to underlying methods.
#'
#' @details
#' When computing transition probabilities for multiple groups, the input data must
#' be grouped via \code{group_by()} before calling this function. Omitting
#' \code{group_by()} will not produce an error or warning but will return
#' silently incorrect results, as the transition probability will be accumulated
#' over the entire dataset rather than within each group.
#'
#' The returned data contains one boundary row per group and transition at
#' \code{time_var = 0} for plotting transition probabilities from the time
#' origin. On this row, \code{trans_prob = 0}; if confidence intervals are
#' requested, \code{trans_lower = trans_upper = 0}. If an interval-length
#' column is present, it is set to \code{0} on the boundary row.
#'
#' @examplesIf require("mstate")
#'   data("prothr", package = "mstate")
#'   prothr <- prothr |>
#'     mutate(transition = as.factor(paste0(from, "->", to))
#'     , treat = as.factor(treat)) |>
#'     filter(Tstart != Tstop, id <= 100) |> select(-trans)
#'   ped <- as_ped(data= prothr, formula= Surv(Tstart, Tstop, status)~ .,
#'     transition = "transition", id= "id", timescale  = "calendar")
#'   pam <- mgcv::bam(ped_status ~ s(tend, by=transition) + transition * treat,
#'     data = ped, family = poisson(), offset = offset,
#'     method = "fREML", discrete = TRUE)
#'   ndf <- make_newdata(ped, tend  = unique(tend),
#'     treat  = unique(treat),
#'     transition = unique(transition)) |>
#'     group_by(treat, transition) |>  # important!
#'     add_trans_prob(pam)
#' @export
add_trans_prob <- function(
  newdata,
  object,
  overwrite = FALSE,
  ci = FALSE,
  alpha = 0.05,
  nsim = 100L,
  time_var = "tend",
  interval_length = "intlen",
  transition = "transition",
  ...
) {
  orig_names <- names(newdata)
  interval_length <- quo_name(enquo(interval_length))
  transition_var <- quo_name(enquo(transition))
  time_var <- resolve_time_var(time_var, object, newdata)
  assert_string(transition_var)
  assert_choice(transition_var, colnames(newdata))

  if (!interval_length %in% colnames(newdata)) {
    newdata <- reconstruct_intlen(
      newdata,
      time_var = time_var,
      interval_length = interval_length
    )
  }
  if (!overwrite) {
    if ("trans_prob" %in% names(newdata)) {
      stop(
        "Data set already contains 'trans_prob' column.
        Set `overwrite=TRUE` to overwrite"
      )
    }
  } else {
    rm.vars <- intersect(
      c("trans_prob", "trans_lower", "trans_upper"),
      names(newdata)
    )
    newdata <- newdata %>% select(-one_of(rm.vars))
  }
  newdata <- drop_cumulative_boundary(newdata, time_var)

  assert_subset(transition_var, names(newdata))
  assert_subset(time_var, names(newdata))
  assert_subset(interval_length, names(newdata))

  # add confidence intervals if wanted
  has_cumu <- "cumu_hazard" %in% colnames(newdata)
  old_groups <- group_vars(newdata)

  if (!(transition_var %in% old_groups)) {
    new_groups <- c(old_groups, transition_var)
    newdata <- newdata %>% group_by(across(all_of(new_groups)))
    old_groups <- group_vars(newdata)
  }

  split_groups <- setdiff(old_groups, transition_var)
  ordering_vars <- c(split_groups, transition_var, time_var)
  newdata <- newdata[
    do.call(order, newdata[, ordering_vars, drop = FALSE]),
    ,
    drop = FALSE
  ]

  if (!has_cumu) {
    newdata <- do.call(
      add_cumu_hazard,
      list(
        newdata = newdata,
        object = object,
        ci = FALSE,
        time_var = time_var,
        interval_length = interval_length,
        boundary = FALSE
      )
    )
  }
  if (ci) {
    newdata <- newdata |>
      add_trans_ci(
        object,
        nsim = nsim,
        alpha = alpha,
        time_var = time_var,
        interval_length = interval_length,
        transition = transition_var
      )
  }
  newdata <- add_cumulative_boundary(
    newdata,
    time_var = time_var,
    values = c(
      cumu_hazard = 0,
      trans_prob = 0,
      trans_lower = 0,
      trans_upper = 0
    ),
    interval_length = interval_length
  )

  if (length(split_groups) == 0L) {
    grp_index <- list(seq_len(nrow(newdata)))
  } else {
    grp <- interaction(
      newdata[, split_groups, drop = FALSE],
      drop = TRUE,
      lex.order = TRUE
    )
    grp_index <- split(seq_len(nrow(newdata)), grp)
  }

  n_grp <- length(grp_index)
  res_data <- vector("list", n_grp)
  res_matrix <- vector("list", n_grp)

  for (i in seq_len(n_grp)) {
    idx <- grp_index[[i]]
    tmp <- newdata[idx, , drop = FALSE] %>%
      group_by(across(all_of(transition_var)))

    res <- get_trans_prob(
      newdata = tmp,
      time_var = time_var,
      interval_length = interval_length,
      transition = transition_var
    )

    res_data[[i]] <- res$data
    res_matrix[[i]] <- res$matrix
  }

  out_df <- bind_rows(res_data)

  if (length(split_groups) == 0L) {
    group_keys <- data.frame(trans_prob_matrix = I(list(res_matrix[[1]])))
  } else {
    group_keys <- out_df %>%
      ungroup() %>%
      distinct(across(all_of(split_groups)))
    group_keys$trans_prob_matrix <- res_matrix
  }

  if (!has_cumu) {
    out_df$cumu_hazard <- NULL
  }

  if (length(old_groups) > 0L) {
    out_df <- out_df %>% group_by(across(all_of(old_groups)))
  }
  if (!"intlen" %in% orig_names) out_df[["intlen"]] <- NULL

  attr(out_df, "matrix") <- group_keys
  out_df
}

#' helper function for add_trans_ci
#' @keywords internal
get_sim_cumu <- function(newdata, interval_length = "intlen", ...) {
  newdata$cumu_hazard <- cumsum(newdata[[interval_length]] * newdata$hazard)

  newdata
}

#' Add transition probabilities confidence intervals
#' @keywords internal
add_trans_ci <- function(newdata, object, nsim = 100L, alpha = 0.05, ...) {
  dots <- list(...)
  time_var <- dots[["time_var"]]
  interval_length <- dots[["interval_length"]]
  transition <- dots[["transition"]]
  if (is.null(time_var)) time_var <- "tend"
  if (is.null(interval_length)) interval_length <- "intlen"
  if (is.null(transition)) transition <- "transition"

  assert_subset(time_var, names(newdata))
  assert_subset(interval_length, names(newdata))
  assert_subset(transition, names(newdata))

  group_vars_orig <- dplyr::group_vars(newdata)
  if (!(transition %in% group_vars_orig)) {
    new_groups <- c(group_vars_orig, transition)
    newdata <- newdata %>% group_by(across(all_of(new_groups)))
    group_vars_orig <- dplyr::group_vars(newdata)
  }
  group_vars_trans <- setdiff(group_vars_orig, transition)
  group_vars_ordered <- c(group_vars_trans, transition)

  df <- as.data.frame(newdata)
  df$.orig_row <- seq_len(nrow(df))
  ordering_vars <- c(group_vars_ordered, time_var)
  df <- df[do.call(order, df[, ordering_vars, drop = FALSE]), , drop = FALSE]

  X <- make_X(object, df)

  groups_array <- interaction(
    df[, group_vars_ordered, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  array_idx_list <- split(seq_len(nrow(df)), groups_array)

  if (length(group_vars_trans) > 0) {
    groups_trans <- interaction(
      df[, group_vars_trans, drop = FALSE],
      drop = TRUE,
      lex.order = TRUE
    )
    trans_idx_list <- split(seq_len(nrow(df)), groups_trans)
  } else {
    trans_idx_list <- list(seq_len(nrow(df)))
  }

  sim_coef_mat <- sample_coefs(object, nsim)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) exp(X %*% z))
  if (is.null(dim(sim_fit_mat))) {
    sim_fit_mat <- matrix(sim_fit_mat, ncol = 1L)
  }

  nrows <- nrow(df)
  sim_trans_probs <- matrix(NA_real_, nrow = nrows, ncol = nsim)
  time_vals <- df[[time_var]]
  intlen_vals <- df[[interval_length]]

  for (i in seq_len(nsim)) {
    haz_i <- sim_fit_mat[, i]
    cumu_hazard <- numeric(nrows)
    for (idx in array_idx_list) {
      idx_ord <- idx[order(time_vals[idx])]
      cumu_hazard[idx_ord] <- cumsum(haz_i[idx_ord] * intlen_vals[idx_ord])
    }

    trans_prob <- numeric(nrows)
    for (idx in trans_idx_list) {
      idx_ord <- idx[order(time_vals[idx], as.character(df[[transition]][idx]))]
      tmp <- df[idx_ord, , drop = FALSE]
      tmp$cumu_hazard <- cumu_hazard[idx_ord]
      res <- get_trans_prob(
        newdata = tmp,
        time_var = time_var,
        interval_length = interval_length,
        transition = transition
      )
      trans_prob[idx_ord] <- res$data$trans_prob
    }

    sim_trans_probs[, i] <- trans_prob
  }

  df$trans_lower <- apply(
    sim_trans_probs,
    1,
    quantile,
    probs = alpha / 2,
    na.rm = TRUE,
    type = 6
  )
  df$trans_upper <- apply(
    sim_trans_probs,
    1,
    quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE,
    type = 6
  )

  df <- df[order(df$.orig_row), , drop = FALSE]
  df$.orig_row <- NULL

  out <- as_tibble(df)
  if (length(group_vars_orig) > 0L) {
    out <- out %>% group_by(across(all_of(group_vars_orig)))
  }

  out
}
