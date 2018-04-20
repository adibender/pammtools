#' Extract cumulative coefficients (cumulative hazard differences)
#'
#' These functions are designed to extract (or mimick) the cumulative coefficients
#' usually used in additive hazards models (Aalen model) to depict (time-varying)
#' covariate effects. For PAMMs, these are the differences
#' between the cumulative hazard rates where all covariates except one have the
#' identical values. For a numeric covariate of interest, we calculate
#' \eqn{Lambda(t|x+1) - Lambda(t|x)}.  For non-numeric covariates
#' the cumulative hazard of the reference level is substracted from
#' the cumulative hazards eveluated at all non reference levels. Standard
#' errors are calculated using the delta method.
#'
#' @rdname cumulative_coefficient
#' @param model Object from which to extract cumulative coefficients.
#' @param data Additional data if necessary.
#' @param terms A character vector of variables for which the cumulative
#' coefficient should be calculated.
#' @param ... Further arguments passed to \code{\link[pammtools]{add_cumu_hazard}}.
#' @export
get_cumu_coef <- function(model, data=NULL, terms, ...) {
  UseMethod("get_cumu_coef", model)
}


#' @rdname cumulative_coefficient
#' @inherit get_cumu_coef
#' @export
get_cumu_coef.gam <- function(model, data, terms, ...) {

  data    <- ped_info(data)
  map(terms, ~cumu_coef(data, model, quo_name(sym(.)), ...)) %>%
    bind_rows()

}

#' @rdname cumulative_coefficient
#' @inherit get_cumu_coef
#' @param ci Logical. Indicates if confidence intervals should be returned as
#' well.
#' @export
get_cumu_coef.aalen <- function(model, data=NULL, terms, ci = TRUE, ...) {

    cumu_coef <- model$cum %>% as_tibble() %>%
      select(map_int(c("time", terms), ~which(grepl(., colnames(model$cum))))) %>%
      gather("variable", "cumu_hazard", -.data$time)
    cumu_var <- model$var %>% as_tibble() %>%
     select(map_int(c("time", terms), ~which(grepl(., colnames(model$cum))))) %>%
      gather("variable", "cumu_var", -.data$time)

    suppressMessages(
      left_join(cumu_coef, cumu_var) %>%
      mutate(
        method     = class(model)[1],
        cumu_lower = .data$cumu_hazard - 2 * .data$cumu_var^0.5,
        cumu_upper = .data$cumu_hazard + 2 * .data$cumu_var^0.5) %>%
      select(one_of(c("method", "variable", "time")), everything(), -one_of("cumu_var"))
      )

}

get_cumu_diff <- function(d1, d2, model) {
  lp <- compute_cum_diff(d1, d2, model)
  d2 %>%
    mutate(
      cumu_hazard = lp[["cumu_diff"]],
      cumu_lower = lp[["cumu_diff"]] - 2*lp[["se_cumu"]],
      cumu_upper = lp[["cumu_diff"]] + 2*lp[["se_cumu"]])
}

#' @inheritParams get_cumu_coef
#' @import dplyr purrr
#' @importFrom rlang sym
#' @keywords internal
cumu_coef <- function(data, model, term, ...) {

  if (quo_name(term) == "(Intercept)") {
    return(get_cumu_coef_baseline(data, model))
  }

  if(is.character(term)) {
    term <- sym(term)
  } else {
    term <- enquo(term)
  }
  qname_term   <- quo_name(term)
  if(!is.numeric(data[[qname_term]])) {
    x <- levels(as.factor(unique(data[[qname_term]])))
  } else {
    x <- mean(data[[qname_term]], na.rm =TRUE)
    x <- c(x, x + 1)
  }
  dat_list <- map(.x=x, function(z) {
    mutate_at(.tbl=data, .vars=qname_term, .funs = ~identity(z)) %>%
    mutate(variable = paste0(qname_term,
      ifelse(is.numeric(z), "", paste0(" (", z, ")"))))})

  map2(
    .x = dat_list[1],
    .y = dat_list[-1],
    .f = ~ get_cumu_diff(.x, .y, model)) %>%
    map(
    ~ select(., one_of(c("variable", "tend")), contains("cumu")) %>%
      rename("time"="tend") %>%
      mutate(method = class(model)[1]) ) %>%
  bind_rows() %>%
  select(one_of(c("method", "variable", "time")), everything())

}

#' @inherit get_cumu_coef
#' @keywords internal
get_cumu_coef_baseline <- function(data, model, ...) {
  ped_info(data) %>%
    mutate_at(
      .vars = vars(-one_of(c("tstart", "tend", "intlen", "intmid", "interval"))),
      .funs = ~0) %>% add_cumu_hazard(model) %>%
    mutate(
      method      = class(model)[1],
      variable    = "(Intercept)") %>%
    rename("time" = "tstart") %>%
    select(one_of(c("method", "variable", "time")), everything())
}
