#' Visualize effect estimates for specific covariate combinations
#'
#' Depending on the plot function and input, creates either a 1-dimensional slices,
#' bivariate surface or  (1D) cumulative effect.
#'
#' @import ggplot2
#' @importFrom rlang quos
#'
#' @inheritParams make_newdata
#' @param data Data used to fit the \code{model}.
#' @param model A suitable model object which will be used to estimate the
#' partial effect of \code{term}.
#' @param term A character string indicating the model term for which partial
#' effects should be plotted.
#' @param reference If specified, should be a list with covariate value pairs,
#' e.g. \code{list(x1 = 1, x2=50)}. The calculated partial effect will be relative
#' to an observation specified in \code{reference}.
#' @param ci Logical. Indicates if confidence intervals for the \code{term}
#' of interest should be calculated/plotted. Defaults to \code{TRUE}.
#' @export
gg_partial <- function(data, model, term, ..., reference = NULL, ci = TRUE) {

  expressions <- quos(...)
  vars        <- names(expressions)
  n_vars      <- length(expressions)
  ndf <- make_newdata(data, ...) %>%
    add_term(model, term, reference = reference, ci = ci)

  n_unique <- map_int(vars, ~length(unique(ndf[[.x]])))
  vars     <- vars[n_unique > 1]
  # vars     <- vars[n_unique[rev(order(n_unique))] > 1]
  # ndf      <- ndf %>% mutate_at(vars[-1], ~as.factor(.x))
  n_vars   <- length(vars)

  gg_base <- ggplot(ndf, aes(x = .data[[vars[1]]])) + xlab(vars[1])
  if (n_vars == 1) {
    gg_out <- gg_base +
      geom_ribbon(aes(ymin = .data[["ci_lower"]], ymax = .data[["ci_upper"]]),
        alpha = 0.3) +
      geom_line(aes(y = .data[["fit"]]))
  } else {
    # if (n_vars == 2) {
      gg_out <- gg_base + aes(y = .data[[vars[2]]], z = .data[["fit"]]) +
        geom_tile(aes(fill = .data[["fit"]])) +
        geom_contour(col = "grey30") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_fill_gradient2(high = "firebrick2", low = "steelblue")
    # }
  }

  gg_out + theme(legend.position = "bottom")

}

#' @rdname gg_partial
#' @inherit gg_partial
#' @importFrom tidyr complete
#' @param time_var The name of the variable that was used in \code{model} to
#' represent follow-up time.
#' @export
gg_partial_ll <- function(
  data,
  model,
  term,
  ...,
  reference = NULL,
  ci        = FALSE,
  time_var  = "tend") {

  ind_term <- which(map_lgl(attr(data, "func_mat_names"),
    ~any(grepl(term, .x))))
  tv_sym <- sym(time_var)
  tz_var <- attr(data, "tz_vars")[[ind_term]]
  tz_val <- attr(data, "tz")[[ind_term]]
  ll_var <- grep("LL", attr(data, "func_mat_names")[[ind_term]], value = TRUE)
  select_vars <- c(time_var, tz_var, ll_var, "fit")
  if (ci) {
    select_vars <- c(select_vars,  "ci_lower", "ci_upper")
  }
  ll_df <- get_partial_ll(data, model, term, ..., reference = reference,
      ci = ci, time_var = time_var) %>%
    select(one_of(select_vars)) %>%
    mutate(fit = ifelse(.data[[ll_var]] == 0, NA_real_, .data$fit)) %>%
    complete(!!tv_sym := unique(!!tv_sym), !!sym(tz_var) := tz_val) %>%
    left_join(int_info(attr(data, "breaks")) %>%
      rename(!!tv_sym := "tend"), by = time_var)
  if (ci) {
    ll_df <- ll_df %>%
      mutate(ci_lower = ifelse(is.na(.data$fit), NA, .data$ci_lower)) %>%
      mutate(ci_upper = ifelse(is.na(.data$fit), NA, .data$ci_upper)) %>%
      gather("type", "fit", one_of(c("fit", "ci_lower", "ci_upper"))) %>%
      mutate(type = factor(.data$type,
        levels = c("ci_lower", "fit", "ci_upper")))
  }

  gg_base <- ggplot(ll_df, aes(x = .data[["intmid"]], y =  tz_var)) +
    geom_tile(aes(fill = .data[["fit"]]), colour = "grey30") +
    scale_fill_gradient2(high = "firebrick2", low = "steelblue",
      na.value = "grey30") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("time") + ylab(expression(t[z]))
  if (ci) {
    gg_base + facet_wrap(~.data$type)
  } else {
    gg_base
  }

}

#' Plot 1D (smooth) effects
#'
#' Flexible, high-level plotting function for (non-linear) effects conditional
#' on further covariate specifications and potentially relative to
#' a comparison specification.
#'
#' @inheritParams gg_partial
#' @importFrom purrr map_int
#' @importFrom rlang quos
#' @examples
#' ped <- tumor[1:200, ] %>% as_ped(Surv(days, status) ~ . )
#' model <- mgcv::gam(ped_status~s(tend) + s(age, by = complications), data=ped,
#'   family = poisson(), offset=offset)
#' make_newdata(ped, age = seq_range(age, 20), complications = levels(complications))
#' gg_slice(ped, model, "age", age=seq_range(age, 20), complications=levels(complications))
#' gg_slice(ped, model, "age", age=seq_range(age, 20), complications=levels(complications),
#'  ci = FALSE)
#' gg_slice(ped, model, "age", age=seq_range(age, 20), complications=levels(complications),
#'   reference=list(age = 50))
#' @export
gg_slice <- function(data, model, term, ..., reference = NULL, ci = TRUE) {

  expressions <- quos(...)
  vars        <- names(expressions)
  ndf         <- make_newdata(data, ...) %>%
    add_term(model, term, reference = reference, ci = ci)

  n_unique <- map_int(vars, ~length(unique(ndf[[.x]])))
  vars     <- vars[rev(order(n_unique))]
  vars     <- vars[n_unique[rev(order(n_unique))] > 1]
  ndf      <- ndf %>% mutate_at(vars[-1], ~as.factor(.x))
  n_vars   <- length(vars)

  gg_out <- ggplot(ndf, aes(x = .data[[vars[1]]], y = .data[["fit"]]))
  if (ci) {
    gg_out <- gg_out +
      geom_ribbon(aes(ymin = .data[["ci_lower"]], ymax = .data[["ci_upper"]]), alpha = 0.3)
  }
  gg_out <- gg_out + geom_line()
  if (n_vars > 1) {
    if(ci) {
      gg_out <- gg_out + aes(group = .data[[vars[2]]], fill = .data[[vars[2]]]) +
      geom_line(aes(col = .data[[vars[2]]]))
    } else {
      gg_out <- gg_out + aes(group = .data[[vars[2]]]) +
        geom_line(aes(col = .data[[vars[2]]]))
    }
    if (n_vars > 2) {
      form   <- as.formula(paste0("~", vars[-1:-2], collapse = "+"))
      gg_out <- gg_out + facet_wrap(form, labeller = label_both)
    }
  }

  gg_out + theme(legend.position = "bottom")

}


#' @rdname get_cumu_eff
#' @inherit get_cumu_eff
#' @inheritParams get_cumu_eff
#' @inheritParams gg_partial
#' @export
gg_cumu_eff <- function(data, model, term, z1, z2=NULL, se_mult = 2, ci = TRUE) {

  cumu_eff_df <- get_cumu_eff(data, model, term, z1, z2, se_mult)

  gg_out <- ggplot(cumu_eff_df, aes(x = .data[["tend"]], y = .data[["cumu_eff"]]))
  if (ci) {
    gg_out <- gg_out +
      geom_ribbon(aes(ymin = .data[["cumu_eff_lower"]], ymax = .data[["cumu_eff_upper"]]),
        alpha = 0.3)
  }

  gg_out + geom_line() + xlab("time") + ylab("cumulative effect")

}




#' @inherit gg_partial_ll
#' @rdname gg_partial
#' @export
get_partial_ll <- function(
  data,
  model,
  term, ...,
  reference = NULL,
  ci        = FALSE,
  time_var  = "tend") {

  ind_term <- which(map_lgl(attr(data, "func_mat_names"), ~any(grepl(term, .x))))
  tz_var   <- attr(data, "tz_vars")[[ind_term]]
  tz_val   <- attr(data, "tz")[[ind_term]]

  ll_df <- get_ll(data, ind_term, ..., time_var = time_var)  %>%
    add_term(object = model, term = term, reference = reference, ci = ci)

}


#' @keywords internal
get_ll <- function(x, ind_term, ..., time_var = "tend") {

  n_func  <- length(attr(x, "ll_funs"))
  ll_fun  <- attr(x, "ll_funs")[[ind_term]]
  int_df  <- int_info(x)
  tz_vars <- attr(x, "tz_vars")
  func    <- attr(x, "func")[[ind_term]]
  tz_val  <- attr(x, "tz")[[ind_term]]
  tz_var  <- attr(x, "tz_vars")[[ind_term]]
  tz_var_mat <- make_mat_names(tz_var, func$latency_var, func$tz_var,
    func$suffix, n_func)
  ll_var_mat <- make_mat_names("LL", func$latency_var, func$tz_var, func$suffix,
    n_func)
  tz_mat_val <- x %>% pull(tz_var_mat) %>% as.numeric() %>% unique() %>% sort()
  nd <- make_newdata(x, ...)
  if (func$latency_var != "") {
    nd <- nd %>% mutate(!!sym(tz_var) := .data[[time_var]] - !!sym(tz_var_mat))
  }
  nd %>% filter(!!sym(tz_var) %in% tz_val) %>%
    mutate(
      !!sym(ll_var_mat) := ll_fun(.data[[time_var]], .data[[tz_var]]) * 1L) %>%
    arrange(.data[[time_var]], .data[[tz_var]]) %>%
    group_by(.data[[tz_var]]) %>%
    mutate(!!sym(ll_var_mat) := lag(!!sym(ll_var_mat), default = 0)) %>%
    ungroup()

}
