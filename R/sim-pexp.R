#' Simulate survival times from the piece-wise exponential distribution
#'
#' @param formula An extended formula that specifies the linear predictor.
#' If you want to include a smooth baseline
#' or time-varying effects, use \code{t} within your formula as
#' if it was a covariate in the data, although it is not and should not
#' be included in the \code{data} provided to \code{sim_pexp}. See examples
#' below.
#'
#' @param data A data set with variables specified in \code{formula}.
#' @param cut A sequence of time-points starting with 0.
#' @import dplyr
#' @import Formula
#' @importFrom lazyeval f_eval
#' @importFrom tidyr replace_na
#' @examples
#' library(survival)
#' library(dplyr)
#' library(pammtools)
#'
#' # set number of observations/subjects
#' n <- 250
#' # create data set with variables which will affect the hazard rate.
#' df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
#'  as_tibble()
#' # the formula which specifies how covariates affet the hazard rate
#' f0 <- function(t) {
#'  dgamma(t, 8, 2) *6
#' }
#' form <- ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)
#' set.seed(24032018)
#' sim_df <- sim_pexp(form, df, 1:10)
#' head(sim_df)
#' plot(survfit(Surv(time, status)~1, data = sim_df ))
#'
#' # for control, estimate with Cox PH
#' mod <- coxph(Surv(time, status) ~ x1 + pspline(x2), data=sim_df)
#' coef(mod)[1]
#' layout(matrix(1:2, nrow=1))
#' termplot(mod, se = TRUE)
#'
#' # and using PAMs
#' layout(1)
#' ped <- sim_df %>% as_ped(Surv(time, status)~., max_time=10)
#' library(mgcv)
#' pam <- gam(ped_status ~ s(tend) + x1 + s(x2), data=ped, family=poisson, offset=offset)
#' coef(pam)[2]
#' plot(pam, page=1)
#'
#'\dontrun{
#' # Example 2: Functional covariates/cumulative coefficients
#' # function to generate one exposure profile, tz is a vector of time points
#' # at which TDC z was observed
#' rng_z = function(nz) {
#'   as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
#' }
#' # two different exposure times  for two different exposures
#' tz1 <- 1:10
#' tz2 <- -5:5
#' # generate exposures and add to data set
#' df <- df %>%
#'   add_tdc(tz1, rng_z) %>%
#'   add_tdc(tz2, rng_z)
#' df
#'
#' # define tri-variate function of time, exposure time and exposure z
#' ft <- function(t, tmax) {
#'   -1*cos(t/tmax*pi)
#' }
#' fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))
#' wpeak2 <- function(lag) 15*dnorm(lag,8,10)
#' wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))
#' f_xyz1 <- function(t, tz, z) {
#'   ft(t, tmax=10) * 0.8*fdnorm(z)* wpeak2(t - tz)
#' }
#' f_xyz2 <- function(t, tz, z) {
#'   wdnorm(t-tz) * z
#' }
#'
#' # define lag-lead window function
#' ll_fun <- function(t, tz) {t >= tz}
#' ll_fun2 <- function(t, tz) {t - 2 >= tz}
#' # simulate data with cumulative effect
#' sim_df <- sim_pexp(
#'   formula = ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)|
#'      fcumu(t, tz1, z.tz1, f_xyz=f_xyz1, ll_fun=ll_fun) +
#'      fcumu(t, tz2, z.tz2, f_xyz=f_xyz2, ll_fun=ll_fun2),
#'   data = df,
#'   cut = 0:10)
#'}
#' @export
sim_pexp <- function(formula, data, cut) {

  data <- data %>%
    mutate(
      id     = row_number(),
      time   = max(cut),
      status = 1)

  # extract formulas for different components
  Form <- Formula(formula)
  f1   <- formula(Form, rhs = 1)
  # later more sophisticated checks + could be used to map over all rhs
  # formulae, check what type of evaluation is needed and return ETAs for
  # each part of the formula separated by |, such that model estimation may
  # be checked for individuals terms/parts
  if (length(Form)[2] > 1) {
    f2  <- formula(Form, rhs = 2)
  } else {
    f2 <- NULL
  }

  # construct eta for time-constant part
  ped  <- split_data(
      formula = Surv(time, status)~.,
      data    = select_if(data, is_atomic),
      cut     = cut,
      id      = "id") %>%
    rename("t" = "tstart") %>%
    mutate(rate = exp(f_eval(f1, .)))

  # construct eta for time-dependent part
  if (!is.null(f2)) {
    terms_f2  <- terms(f2, specials = "fcumu")
    f2_ev     <- list()
    f2_tl <- attr(terms_f2, "term.labels")
    for (i in seq_along(f2_tl)) {
      f2_ev[[i]] <- eval(expr = parse(text = f2_tl[[i]]), envir = .GlobalEnv)
    }
    ll_funs   <- map(f2_ev, ~.x[["ll_fun"]])
    tz_vars   <- map_chr(f2_ev, ~.x[["vars"]][1])
    cumu_funs <- map(f2_ev, ~.x[["f_xyz"]])
    names(tz_vars) <- names(ll_funs) <- names(cumu_funs) <- tz_vars
    z_form <- list("eta_", map_chr(f2_ev, ~.x[["vars"]][2])) %>%
      reduce(paste0, collapse = "+") %>% paste0("~", .) %>% as.formula()

    df2 <- map(f2_ev, function(fc) eta_cumu(data = data, fc, cut = cut))
    suppressMessages(
      ped <- ped %>%
        left_join(reduce(df2, full_join))
    )
    ped <- ped %>%
      mutate_at(vars(contains("eta_")), replace_na, 0) %>%
      group_by(.data$id, .data$t) %>%
      mutate(eta_z = !!rlang::get_expr(z_form)) %>%
      mutate(rate = .data$rate * exp(.data$eta_z))
  } else {
    tz_vars <- NULL
  }

  sim_df <- ped %>%
    group_by(id) %>%
    summarize(time = rpexp(rate = .data$rate, t = .data$t)) %>%
    mutate(
      status = 1L * (.data$time <= max(cut)),
      time   = pmin(.data$time, max(cut)))

  suppressMessages(
    sim_df <- sim_df %>%
      left_join(select(data, -all_of(c("time", "status"))))
  )

  attr(sim_df, "id_var")     <- "id"
  attr(sim_df, "time_var")   <- "time"
  attr(sim_df, "status_var") <- "status"
  attr(sim_df, "tz_var")     <- tz_vars
  attr(sim_df, "cens_value") <- 0
  attr(sim_df, "breaks")     <- cut
  attr(sim_df, "tz")         <- imap(tz_vars, ~select(sim_df, all_of(.x)) %>%
    pull(.x) %>% unique()) %>% flatten()
  if (exists("ll_funs")) attr(sim_df, "ll_funs") <- ll_funs
  if (exists("cumu_funs")) attr(sim_df, "cumu_funs") <- cumu_funs
  attr(sim_df, "id_n") <- sim_df %>% pull("time") %>%
    pmin(max(cut)) %>%
    map_int(findInterval, vec = cut, left.open = TRUE, rightmost.closed = TRUE)
  attr(sim_df, "id_tseq") <- attr(sim_df, "id_n") %>%
    map(seq_len) %>% unlist()
  attr(sim_df, "id_tz_seq") <- rep(seq_along(pull(sim_df, id)),
    times = attr(sim_df, "id_n"))
  attr(sim_df, "sim_formula") <- formula

  class(sim_df) <- c("sim_df", class(unped(sim_df)))

  if (any(!map_lgl(sim_df, is_atomic))) {
    class(sim_df) <- c("nested_fdf", class(sim_df))
  }

  sim_df

}


#' Add time-dependent covariate to a data set
#'
#' Given a data set in standard format (with one row per subject/observation),
#' this function adds a column with the specified exposure time points
#' and a column with respective exposures, created from \code{rng_fun}.
#' This function should usually only be used to create data sets passed
#' to \code{\link[pammtools]{sim_pexp}}.
#'
#' @inheritParams sim_pexp
#' @param tz A numeric vector of exposure times (relative to the
#' beginning of the follow-up time \code{t})
#' @param rng_fun A random number generating function that creates
#' the time-dependent covariates at time points \code{tz}.
#' First argument of the function should be \code{n}, the number of
#' random numbers to generate. Within \code{add_tdc}, \code{n} will be set
#' to \code{length(tz)}.
#' @param ... Currently not used.
#' @import dplyr
#' @importFrom rlang eval_tidy :=
#' @importFrom purrr map
#' @export
add_tdc <- function(data, tz, rng_fun, ...) {

  tz      <- enquo(tz)
  nz      <- length(eval_tidy(tz))
  name_tz <- quo_name(tz)
  z_var   <- paste0("z.", name_tz)

  data %>%
    mutate(
      !!name_tz := map(seq_len(n()), ~ !!tz),
      !!z_var   := map(seq_len(n()), ~ rng_fun(nz = nz))) %>%
    as_tibble()

}



#' A formula special used to handle cumulative effect specifications
#'
#' Can be used in the second part of the formula specification provided
#' to \code{\link[pammtools]{sim_pexp}} and should only be used in this
#' context.
#'
#' @importFrom purrr map
#' @export
#' @keywords internal
fcumu <- function(..., by = NULL, f_xyz, ll_fun) {

  vars   <- as.list(substitute(list(...)))[-1] %>%
    map(~as.character(.x)) %>%
    unlist()
  vars <- vars[vars != "t"]

  list(
    vars   = vars,
    f_xyz  = f_xyz,
    ll_fun = ll_fun)

}

#' @import dplyr
#' @importFrom tidyr unnest
#' @importFrom rlang sym :=
#' @keywords internal
eta_cumu <- function(data, fcumu, cut, ...) {

  vars   <- fcumu$vars
  f_xyz  <- fcumu$f_xyz
  ll_fun <- fcumu$ll_fun
  eta_name <- paste0("eta_", vars[2])
  comb_df <- combine_df(
    data.frame(t = cut),
    select(data, one_of("id", vars)))
  comb_df <- comb_df %>% unnest(cols = -one_of("id"))
  comb_df %>%
    group_by(.data$id, .data$t) %>%
    mutate(
      LL = ll_fun(t, !!sym(vars[1])) * 1,
      delta = c(mean(abs(diff(!!sym(vars[1])))), abs(diff(!!sym(vars[1]))))) %>%
    ungroup() %>%
    filter(.data$LL != 0) %>%
    group_by(.data$id, .data$t) %>%
    summarize(!!eta_name :=
      sum(.data$delta * f_xyz(.data$t, .data[[vars[1]]], .data[[vars[2]]])))

}

#' Simulate data for competing risks scenario
#'
#'
#' @keywords internal
sim_pexp_cr <- function(formula, data, cut) {

  # Formula extends the base class formula by allowing for multiple responses and multiple parts of regressors
  Form    <- Formula(formula)
  # Extract the right handside of the Formula
  F_rhs   <- attr(Form, "rhs")
  l_rhs   <- length(F_rhs)
  seq_rhs <- seq_len(l_rhs)

  if (!("id" %in% names(data))) {
    data$id <- 1:(nrow(data))
  }

  if (!("t" %in% names(data))) {
    data$t <- 0
  }

  data <- data %>%
    mutate(
      time   = max(cut),
      status = 1
    )

  # construct eta for time-constant part
  # offset (the log of the duration during which the subject was under risk in that interval)

  ped  <- split_data(
    formula = Surv(time, status)~.,
    data    = select_if(data, is_atomic),
    cut     = cut,
    id      = "id") %>%
    mutate(
      t = t + .data$tstart
    )

  # calculate cause specific hazards

  for (i in seq_rhs) {
    ped[[paste0("hazard", i)]] <-  exp(eval(F_rhs[[i]], ped))
  }
  ped[["rate"]] <- reduce(ped[paste0("hazard", seq_rhs)], `+`)

  # simulate survival times

  sim_df <- ped %>%
    group_by(id) %>%
    mutate(
      time   = rpexp(rate = .data$rate, t = .data$tstart),
      status = 1L * (.data$time <= max(cut)),
      time   = pmin(.data$time, max(cut)),
      # t wieder ins "Original" zurückrechnen, muss später auf die Waitingtime drauf gerechnet werden
      t = .data$t - .data$tstart
    ) %>%
    filter(.data$tstart < .data$time & .data$time <= .data$tend)



  # Ziehe aus den möglichen hazards eins mit den entsprechenden Wahrscheinlichkeiten
  sim_df$type <- apply(sim_df[paste0("hazard", seq_rhs)], 1,
    function(probs)
      sample(seq_rhs, 1, prob = probs))

  sim_df %>%
    select(-one_of(c("tstart", "tend", "interval", "offset", "ped_status", "rate")))

}
