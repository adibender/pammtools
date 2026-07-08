#' Extract cumulative coefficients (cumulative hazard differences)
#'
#' These functions are designed to extract (or mimic) the cumulative coefficients
#' usually used in additive hazards models (Aalen model) to depict (time-varying)
#' covariate effects. For PAMMs, these are the differences
#' between the cumulative hazard rates where all covariates except one have the
#' identical values. For a numeric covariate of interest, this calculates
#' \eqn{\Lambda(t|x+1) - \Lambda(t|x)}.  For non-numeric covariates
#' the cumulative hazard of the reference level is subtracted from
#' the cumulative hazards evaluated at all non reference levels. Standard
#' errors are calculated using the delta method.
#'
#' @rdname cumulative_coefficient
#' @param model Object from which to extract cumulative coefficients.
#' @param data Additional data if necessary.
#' @param terms A character vector of variables for which the cumulative
#' coefficient should be calculated.
#' @param ... Further arguments passed to methods.
#' @export
get_cumu_coef <- function(model, data = NULL, terms, ...) {
  UseMethod("get_cumu_coef", model)
}


#' @rdname cumulative_coefficient
#' @param time_var Name of the evaluation time variable in `data`. Defaults to
#'   \code{"tend"}.
#' @param interval_length Name of the interval-length variable in `data`.
#'   Defaults to \code{"intlen"}.
#' @export
get_cumu_coef.gam <- function(
  model,
  data,
  terms,
  time_var = "tend",
  interval_length = "intlen",
  ...
) {
  data <- ped_info(data)
  if (time_var != "tend" && "tend" %in% names(data)) {
    data[[time_var]] <- data[["tend"]]
  }
  if (interval_length != "intlen" && "intlen" %in% names(data)) {
    data[[interval_length]] <- data[["intlen"]]
  }
  if (!time_var %in% names(data)) {
    if ("tend" %in% names(data)) {
      data[[time_var]] <- data[["tend"]]
    } else {
      stop(
        "Column '",
        time_var,
        "' not found in `data`, and fallback column 'tend' is unavailable."
      )
    }
  }
  map(
    terms,
    \(term_i)
      cumu_coef(
        data,
        model,
        quo_name(sym(term_i)),
        time_var = time_var,
        interval_length = interval_length,
        ...
      )
  ) %>%
    bind_rows()
}

#' @rdname cumulative_coefficient
#' @export
get_cumu_coef.scam <- function(
  model,
  data,
  terms,
  time_var = "tend",
  interval_length = "intlen",
  ...
) {
  get_cumu_coef.gam(
    model = model,
    data = data,
    terms = terms,
    time_var = time_var,
    interval_length = interval_length,
    ...
  )
}

#' @rdname cumulative_coefficient
#' @param ci Logical. Indicates if confidence intervals should be returned as
#' well.
#' @export
get_cumu_coef.aalen <- function(model, data = NULL, terms, ci = TRUE, ...) {
  terms <- map(
    c("time", terms),
    ~ grep(.x, colnames(model$cum), value = TRUE)
  ) %>%
    reduce(union)
  cumu_coef <- model[["cum"]] %>%
    as_tibble() %>%
    select(one_of(terms)) %>%
    gather("variable", "cumu_hazard", -.data[["time"]])
  cumu_var <- model[["var.cum"]] %>%
    as_tibble() %>%
    select(terms) %>%
    gather("variable", "cumu_var", -.data[["time"]])

  suppressMessages(
    left_join(cumu_coef, cumu_var) %>%
      mutate(
        method = class(model)[1],
        cumu_lower = .data$cumu_hazard - 2 * .data$cumu_var**0.5,
        cumu_upper = .data$cumu_hazard + 2 * .data$cumu_var**0.5
      ) %>%
      select(
        one_of(c("method", "variable", "time")),
        everything(),
        -one_of("cumu_var")
      )
  )
}

#' @rdname cumulative_coefficient
#' @export
get_cumu_coef.cox.aalen <- function(model, data = NULL, terms, ci = TRUE, ...) {
  get_cumu_coef.aalen(model = model, data = data, terms = terms, ci = ci, ...)
}

get_cumu_diff <- function(
  d1,
  d2,
  model,
  nsim = 100L,
  alpha = 0.05,
  time_var = "tend",
  interval_length = "intlen"
) {
  lp <- compute_cumu_diff(
    d1,
    d2,
    model,
    alpha = alpha,
    nsim = nsim,
    time_var = time_var,
    interval_length = interval_length
  )
  d2 %>%
    mutate(
      cumu_hazard = lp[["cumu_diff"]],
      cumu_lower = lp[["cumu_lower"]],
      cumu_upper = lp[["cumu_upper"]]
    )
}

#' @import dplyr purrr
#' @importFrom rlang sym enquo quo_name
#' @keywords internal
cumu_coef <- function(
  data,
  model,
  term,
  nsim = 100L,
  alpha = 0.05,
  time_var = "tend",
  interval_length = "intlen",
  ...
) {
  if (quo_name(term) == "(Intercept)") {
    return(
      get_cumu_coef_baseline(
        data,
        model,
        time_var = time_var,
        interval_length = interval_length
      )
    )
  }

  if (is.character(term)) {
    term <- sym(term)
  } else {
    term <- enquo(term)
  }
  qname_term <- quo_name(term)

  if (!is.numeric(data[[qname_term]])) {
    x <- levels(as.factor(unique(data[[qname_term]])))
  } else {
    x <- mean(data[[qname_term]], na.rm = TRUE)
    x <- c(x, x + 1)
  }
  dat_list <- map(.x = x, function(z) {
    mutate_at(.tbl = data, .vars = qname_term, .funs = ~ identity(z)) %>%
      mutate(
        variable = paste0(
          qname_term,
          ifelse(is.numeric(z), "", paste0(" (", z, ")"))
        )
      )
  })

  map2(
    .x = dat_list[1],
    .y = dat_list[-1],
    .f = ~ get_cumu_diff(
      .x,
      .y,
      model,
      nsim = nsim,
      alpha = alpha,
      time_var = time_var,
      interval_length = interval_length
    )
  ) %>%
    map(
      ~ select(., one_of(c("variable", time_var)), contains("cumu")) %>%
        rename(time = all_of(time_var)) %>%
        mutate(method = class(model)[1])
    ) %>%
    bind_rows() %>%
    select(one_of(c("method", "variable", "time")), everything())
}

#' @keywords internal
get_cumu_coef_baseline <- function(
  data,
  model,
  time_var = "tend",
  interval_length = "intlen",
  ...
) {
  vars_modify <- colnames(data)[map_lgl(data, is.numeric)] %>%
    setdiff(c("tstart", interval_length, "intmid", time_var))

  data %>%
    mutate_at(
      .vars = vars(one_of(vars_modify)),
      .funs = ~ c(0)
    ) %>%
    add_cumu_hazard(
      model,
      time_var = time_var,
      interval_length = interval_length,
      boundary = FALSE,
      check_grouping = FALSE
    ) %>%
    mutate(
      method = class(model)[1],
      variable = "(Intercept)"
    ) %>%
    rename(time = all_of(time_var)) %>%
    select(one_of(c(
      "method",
      "variable",
      "time",
      "cumu_hazard",
      "cumu_lower",
      "cumu_upper"
    )))
}


#' Calculate difference in cumulative hazards and respective standard errors
#'
#' CIs are calculated by sampling coefficients from their posterior and
#' calculating the cumulative hazard difference \code{nsim} times. The CI
#' are obtained by the 2.5\% and 97.5\% empirical (type-6) quantiles.
#'
#' @param d1 A data set used as \code{newdata} in \code{\link{make_X}}
#' @param d2 See \code{d1}
#' @param model A model object for which a predict method is implemented which
#' returns the design matrix (e.g., \code{mgcv::gam}).
#' @importFrom stats coef
#' @importFrom mvtnorm rmvnorm
#' @keywords internal
compute_cumu_diff <- function(
  d1,
  d2,
  model,
  alpha = 0.05,
  nsim = 100L,
  time_var = "tend",
  interval_length = "intlen"
) {
  if (!interval_length %in% colnames(d1)) {
    d1 <- reconstruct_intlen(
      d1,
      time_var = time_var,
      interval_length = interval_length
    )
  }
  if (!interval_length %in% colnames(d2)) {
    d2 <- reconstruct_intlen(
      d2,
      time_var = time_var,
      interval_length = interval_length
    )
  }
  intlen1 <- d1[[interval_length]]
  intlen2 <- d2[[interval_length]]

  X1 <- make_X(model, d1)
  X2 <- make_X(model, d2)
  coefs <- get_coefs(model)
  sim_coef_mat <- sample_coefs(model, nsim)
  sim_fit_mat <- apply(sim_coef_mat, 1, function(z) {
    cumsum(intlen2 * exp(drop(X2 %*% z))) -
      cumsum(intlen1 * exp(drop(X1 %*% z)))
  })

  cumu_lower <- apply(sim_fit_mat, 1, quantile, probs = alpha / 2, type = 6)
  cumu_upper <- apply(sim_fit_mat, 1, quantile, probs = 1 - alpha / 2, type = 6)
  haz1 <- exp(drop(X1 %*% coefs))
  haz2 <- exp(drop(X2 %*% coefs))
  cumu_diff <- cumsum(haz2 * intlen2) - cumsum(haz1 * intlen1)

  list(cumu_diff = cumu_diff, cumu_lower = cumu_lower, cumu_upper = cumu_upper)
}
