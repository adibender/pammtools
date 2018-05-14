#' Visualize effect estimates for specific covariate combinations
#'
#' Depending on the plot function and input, creates either a 1-dimensional slices,
#' bivariate surface or  (1D) cumulative effect.
#'
#' @import ggplot2
#' @importFrom rlang quos
#'
#' @inheritParams make_newdata
#' @param model A suitable model object which will be used to estimate the
#' partial effect of \code{term}.
#' @param reference If specified, should be a list with covariate value pairs,
#' e.g. \code{list(x1 = 1, x=50)}. The calculated partial effect will be relative
#' to an observation specified in \code{reference}.
#' @export
#' @examples
#' ped <- tumor[1:200, ] %>% as_ped(Surv(days, status)~.)
#' model <- mgcv::gam(ped_status~s(tend) + s(age, by = complications), data=ped,
#'   family = poisson(), offset=offset)
#' make_newdata(ped, age = seq_range(age, 20), complications = levels(complications))
#' gg_slice(ped, model, "age", age=seq_range(age, 20), complications=levels(complications))
#' gg_slice(ped, model, "age", age=seq_range(age, 20), complications=levels(complications),
#'   reference=list(age = 50))
#' @keywords internal
gg_partial <- function(data, model, term, ..., reference = NULL, ci = TRUE) {

  expressions <- quos(...)
  vars        <- names(expressions)
  n_vars      <- length(expressions)

  ndf <- make_newdata(data, ...) %>%
    add_term2(model, term, reference=reference, ci = ci)


  gg_base <- ggplot(ndf, aes_string(x = vars[1])) +
    xlab(vars[1])
  if(n_vars == 1) {
    gg_out <- gg_base +
      geom_ribbon(aes_string(ymin = "ci_lower", ymax = "ci_upper"), alpha = 0.3) +
      geom_line(aes_string(y = "fit"))
  } else {
    if(n_vars == 2) {
      gg_out <- gg_base + aes_string(y = vars[2], z = "fit") +
        geom_tile(aes_string(fill = "fit")) +
        geom_contour(col = "grey30") +
        scale_y_reverse(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        scale_fill_gradient2(high = "firebrick2", low = "steelblue")
    }
  }

  gg_out + theme(legend.position = "bottom")
}

#' @rdname gg_partial
#' @inherit gg_partial
#' @importFrom tidyr complete
#' @export
gg_partial_ll <- function(data, model, term, ..., reference=NULL, ci=FALSE) {

  ind_term <- which(map_lgl(attr(data, "func_mat_names"), ~any(grepl(term, .x))))
  tz_var   <- attr(data, "tz_vars")[[ind_term]]
  tz_val   <- attr(data, "tz")[[ind_term]]
  ll_fun <- attr(data, "ll_funs")[[ind_term]]
  ll_var <- grep("LL", attr(data, "func_mat_names")[[ind_term]], value=TRUE)
  select_vars <- c("tend", tz_var, ll_var, "fit")
  if (ci) {
    select_vars <- c(select_vars,  "ci_lower", "ci_upper")
  }
  ll_df <- get_partial_ll(data, model, term, ..., reference=reference, ci=ci) %>%
    select(one_of(select_vars)) %>%
    mutate(fit = ifelse(.data[[ll_var]] == 0, NA_real_, .data$fit)) %>%
    complete(tend=attr(data, "breaks")[-1], !!sym(tz_var):= tz_val) %>%
    left_join(int_info(attr(data, "breaks")), by = "tend")
  if(ci) {
    ll_df <- ll_df %>%
      mutate(ci_lower = ifelse(is.na(.data$fit), NA, .data$ci_lower)) %>%
      mutate(ci_upper = ifelse(is.na(.data$fit), NA, .data$ci_upper)) %>%
      gather(type, fit, one_of(c("fit", "ci_lower", "ci_upper"))) %>%
      mutate(type = factor(.data$type, levels = c("ci_lower", "fit", "ci_upper")))
  }

  gg_base <- ggplot(ll_df, aes_string(x = "intmid", y = "tz")) +
    geom_tile(aes_string(fill = "fit"), colour = "grey30") +
    scale_fill_gradient2(high = "firebrick2", low = "steelblue",
      na.value = "grey30") +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    xlab("time") + ylab(expression(t[z]))
  if(ci) {
    gg_base + facet_wrap(~type)
  } else {
    gg_base
  }

}

#' @rdname gg_partial
#' @importFrom purrr map_int
#' @importFrom rlang quos
#' @export
gg_slice <- function(data, model, term, ..., reference=NULL) {

  expressions <- quos(...)
  vars        <- names(expressions)
  ndf         <- make_newdata(data, ...) %>%
    add_term2(model, term, reference = reference)

  n_unique <- map_int(vars, ~length(unique(ndf[[.x]])))
  vars     <- vars[rev(order(n_unique))]
  vars     <- vars[n_unique[rev(order(n_unique))] > 1]
  ndf      <- ndf %>% mutate_at(vars[-1], ~as.factor(.x))
  n_vars   <- length(vars)

  gg_out <- ggplot(ndf, aes_string(x = vars[1], y = "fit")) +
    geom_ribbon(aes_string(ymin = "ci_lower", ymax = "ci_upper"), alpha = 0.3) +
    geom_line()
  if(n_vars > 1) {
    gg_out <- gg_out + aes_string(group=vars[2], fill=vars[2]) +
      geom_line(aes_string(col = vars[2]))
    if(n_vars > 2) {
      form   <- as.formula(paste0("~", vars[-1:-2], collapse="+"))
      gg_out <- gg_out + facet_wrap(form, labeller=label_both)
    }
  }

  gg_out + theme(legend.position = "bottom")

}


#' @rdname gg_partial
#' @inherit gg_partial
#'
#' @export
gg_cumu_eff <- function(data, model, term, z1, z2=NULL, se_mult = 2) {

  cumu_eff_df <- get_cumu_eff(data, model, term, z1, z2, se_mult)

  ggplot(cumu_eff_df, aes_string(x = "tend", y = "cumu_eff")) +
    geom_ribbon(aes_string(ymin = "cumu_eff_lower", ymax = "cumu_eff_upper"), alpha = 0.3) +
    geom_line() +
    xlab("time") + ylab("cumulative effect")

}




#' @inherit gg_partial
#' @export
#' @keywords internal
get_partial_ll <- function(x, model, term, ..., reference=NULL, ci = FALSE) {

  ind_term <- which(map_lgl(attr(x, "func_mat_names"), ~any(grepl(term, .x))))
  tz_var   <- attr(x, "tz_vars")[[ind_term]]
  tz_val   <- attr(x, "tz")[[ind_term]]

  ll_df <- get_ll(x, ind_term, ...)  %>%
    add_term2(object=model, term=term, reference=reference, se.fit = ci)

}


#' @keywords internal
get_ll <- function(x, ind_term, ...) {

  n_func  <- length(attr(x, "ll_funs"))
  ll_fun  <- attr(x, "ll_funs")[[ind_term]]
  int_df  <- int_info(x)
  tz_vars <- attr(x, "tz_vars")
  func    <- attr(x, "func")[[ind_term]]
  tz_val      <- attr(x, "tz")[[ind_term]]
  tz_var  <- attr(x, "tz_vars")[[ind_term]]
  tz_var_mat <- make_mat_names(tz_var, func$latency_var, func$tz_var, func$suffix, n_func)
  ll_var_mat <- make_mat_names("LL", func$latency_var, func$tz_var, func$suffix, n_func)
  tz_mat_val <- x %>% pull(tz_var_mat) %>% as.numeric() %>% unique() %>% sort()
  nd <- make_newdata(x, tend = int_df$tend, !!sym(tz_var_mat) := tz_mat_val,...)
  if(func$latency_var != "") {
    nd <- nd %>% mutate(!!sym(tz_var) := .data$tend - !!sym(tz_var_mat))
  }
  nd <- nd %>% filter(!!sym(tz_var) %in% tz_val) %>%
    mutate(!!sym(ll_var_mat) := ll_fun(.data$tend, .data[[tz_var]])*1L)

}
