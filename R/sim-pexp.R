#' Simulate survival times from the piece-wise exponential distribution
#'
#' @param formula An extended formula that denotes how the linear predictor
#' should be comprised. If you want to include a smooth baseline
#' or time-varying effects, use \code{t} within your formula as
#' if it was a covariate in the data, although it is not and should not
#' be included in the \code{data} provided to \code{sim_pexp}. See examples
#' below.
#'
#' @param data Data in which to look for variables in formula.
#' @param cut A sequence of time-points starting with 0
#' @import dplyr
#' @import Formula
#' @importFrom msm rpexp
#' @importFrom lazyeval f_eval
#' @examples
#' library(survival)
#' library(dplyr)
#' library(pammtools)
#'
#' # set number of observations/subjects
#' n <- 250
#' # create data set with variables which will affect the hazard rate.
#' df <- cbind.data.frame(x1 = runif(n, -3, 3), x2 = runif(n, 0, 6)) %>%
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
#' ped <- split_data(Surv(time, status)~., data=sim_df, max_time=10)
#' library(mgcv)
#' pam <- gam(ped_status ~ s(tend) + x1 + s(x2), data=ped, family=poisson, offset=offset)
#' coef(pam)[2]
#' plot(pam, page=1)
#'
#' ## Example 2: Functional covariates/cumulative coefficients
#' # function to generate one exposure profile, te is a vector of exposure time points
#' rng_z = function(nz) {
#'   as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
#' }
#' # two different exposure times  for two different exposures
#' te1 <- 1:10
#' te2 <- -5:5
#' # generate exposures and add to data set
#' df <- df %>%
#'   add_tdc(te1, rng_z) %>%
#'   add_tdc(te2, rng_z)
#' df
#'
#' # define tri-variate function of time, exposure time and exposure z
#' ft <- function(t, tmax) {
#'   -1*cos(t/tmax*pi)
#' }
#' fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))
#' wpeak2 <- function(lag) 15*dnorm(lag,8,10)
#' f_ttez <- function(t, te, z) {
#'   ft(t, tmax=10) * 0.8*fdnorm(z)* wpeak2(t-te)
#' }
#' # define lag-lead window function
#' ll_fun <- function(t, te) {t <= te}
#' # simulate data with cumulative effect
#' sim_df <- sim_pexp(
#'   formula = ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)|
#'      fcumu(t, te1, z.te1, f_xyz=f_ttez, ll_fun=ll_fun),
#'   data = df,
#'   cut = 0:10)
#' plot(survfit(Surv(time, status)~1, data = sim_df ))
#' @export
sim_pexp <- function(formula, data, cut) {

  n <- nrow(data)
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
  if(length(Form)[2] > 1) {
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
  if(!is.null(f2)) {
    terms_f2 <- terms(f2, specials="fcumu")
    f2_ev    <- map(attr(terms_f2, "term.labels"), ~eval(expr=parse(text=.x)))
    df2      <- map_dfc(f2_ev, ~eta_cumu(data, ., cut))
    ped <- ped %>%
      left_join(df2) %>%
      mutate(rate = rate*exp(eta_z))
  }

  sim_df <- ped %>%
    group_by(id) %>%
    summarize(time = rpexp(rate = .data$rate, t = .data$t)) %>%
    mutate(
      status = 1L*(.data$time <= max(cut)),
      time = pmin(.data$time, max(cut))) %>%
  left_join(select(data, -.data$time, -.data$status))

  attr(sim_df, "formula") <- formula

  sim_df

}


#' @import dplyr
#' @importFrom rlang eval_tidy
#' @importFrom purrr map
#' @export
add_tdc <- function(data, te, rng_fun, ...) {
  te      <- enquo(te)
  nz      <- length(eval_tidy(te))
  name_te <- quo_name(te)
  z_var   <- paste0("z.", name_te)

  data %>% mutate(
    !!name_te := map(seq_len(n()), ~ !!te),
    !!z_var   := map(seq_len(n()), ~ rng_fun(nz)))

}



#' @importFrom purrr map
#' @export
#' @keywords internal
fcumu <- function(..., by = NULL, f_xyz, ll_fun) {

  vars   <- as.list(substitute(list(...)))[-1] %>%
    map(~as.character(.)) %>%
    unlist()
  vars <- vars[vars != "t"]

  list(vars = vars, f_xyz = f_xyz, ll_fun = ll_fun)

}


#' @import dplyr
#' @importFrom tidyr unnest
#' @importFrom rlang sym
#' @keywords internal
eta_cumu <- function(data, fcumu, cut, ...) {

  vars <- fcumu$vars
  # func_var <- grep(".", vars, fixed=T, value=T)
  f_xyz <- fcumu$f_xyz
  ll_fun <- fcumu$ll_fun
  combine_df(
      data.frame(t=cut),
      select(data, one_of("id", vars))) %>%
    unnest() %>%
    filter(ll_fun(t, !!sym(vars[1]))) %>%
    group_by(.data$id, .data$t) %>%
    summarize(eta_z = sum(f_xyz(.data$t, .data[[vars[1]]], .data[[vars[2]]])))
}
