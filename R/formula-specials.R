#' Formula specials for defining  time-dependent covariates
#'
#' So far, two specials are implemented. \code{concurrent} is used when
#' the goal is to estimate a concurrent effect of the TDC. \code{func}
#' is used when the goal is to estimate a cumulative effect of the TDC.
#'
#' @rdname specials
#' @importFrom purrr map
#' @export
#' @keywords internal
func <- function(...,
  te_var,
  ll_fun = function(t, te) {t >= te},
  suffix = NULL) {

  vars        <- as.list(substitute(list(...)))[-1]
  vars_chr    <- vars %>% map(~as.character(.))
  lgl_latency <- map_lgl(vars_chr, ~any(. %in% "latency"))

  if (any(lgl_latency)) {
    latency_var <- unlist(vars_chr)[unlist(vars_chr) != "latency"][lgl_latency]
    col_vars    <- unlist(vars_chr)[unlist(vars_chr) != "latency"]
  } else {
    latency_var <- ""
    col_vars    <- unlist(vars_chr)
  }

  list(
    col_vars    = col_vars,
    latency_var = latency_var,
    te_var      = te_var,
    suffix      = suffix,
    ll_fun      = ll_fun)

}


#' @rdname specials
#' @inherit func
#' @keywords internal
concurrent <- function(...,
  te_var,
  ll_fun = function(t) {t == t},
  suffix = NULL) {

  vars        <- as.list(substitute(list(...)))[-1]
  vars_chr    <- vars %>% map(~as.character(.)) %>% unlist()


  list(
    col_vars    = vars_chr,
    te_var      = te_var,
    suffix      = suffix,
    ll_fun      = ll_fun)

}



#' Expand time-dependent covariates to functionals
#'
#' Given formula specification on how time-dependent covariates affect the
#' outcome, creates respective functional covariate as well as auxiliary
#' matrices for time/latency etc.
#'
#' @param data Data frame (or similar) in which variables specified in ...
#' will be looked for
#' @param ... One sided formula of the form \code{~func(t,te,x) + func(t-te, z)}
#' that specifies the type of cumulative effect desired (see examples).
#' Currently it is assumed that observations of covariates \code{x}, \code{z}
#' happened on the same time scale \code{te}.
#' Possible specifications are
#' \itemize{
#' \item \code{func(t,te,x)}: The most general specification.
#' \item \code{func(t, t-te,x)}: A time-varying DLNM.
#' \item \code{func(t-te,x}: A DLNM.
#' \item \code{func(t-te,by=x) A WCE}}
#' @importFrom purrr flatten map
#' @importFrom stats terms
#' @export
#' @keywords internal
get_func <- function(data, formula) {

  func_list <- eval_func(formula)

  n_func <- length(func_list)
  ll_funs <- map(func_list, ~.x[["ll_fun"]])
  te_vars <- map(func_list, ~.x[["te_var"]])
  te <- map(te_vars, ~pull(data, .x) %>% unlist() %>% unique() %>% sort())

  names(te) <- names(te_vars) <- names(ll_funs) <- te_vars

  ## create matrices
  func_mats <- map(func_list, ~expand_func(data=data, ., n_func=n_func)) %>% flatten()

  list(
    func_mats = func_mats,
    ll_funs   = ll_funs,
    te_vars   = te_vars,
    te        = te)

}

#' @keywords internal
eval_func <- function(formula) {

  tf  <- terms(get_tdc_form(formula), specials="func")
  # extract components
  terms_vec <- attr(tf, "term.labels")

  map(terms_vec, ~eval(expr=parse(text=.)))

}


#' @inherit get_func
has_special <- function(formula, special = "func") {
  if(!has_tdc_form) {
    return(FALSE)
  } else {
    formula <- formula(Formula(formula), lhs=FALSE, rhs = 2)
    terms <- terms(formula, specials = specials)
    if(is.null(attr(terms, "specials")[[special]])) {
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}


get_te_from_concurrent <- function(concurrent) {

  concurrent$te %>% unlist() %>% unique() %>% sort()

}

#' @rdname get_func
#' @inheritParams get_func
#' @param func Single evaluated \code{\link{func}} term.
#' @importFrom purrr map invoke_map
#' @keywords internal
expand_func <- function(data, func, n_func) {

  col_vars <- func$col_vars
  te_var <- func$te_var
  te <- pull(data, te_var) %>% unlist() %>% unique() %>% sort()
  time_var <- attr(data, "time_var")
  id_var <- attr(data, "id_var")
  lgl_var_in_data <- map_lgl(col_vars, ~ . %in% colnames(data))
  if (!all(lgl_var_in_data)) {
    stop(paste0("The following variables provided to 'formula' are not contained
      in 'data': ", col_vars[!lgl_var_in_data]))
  }
  ncols_vars <- get_ncols(data, col_vars[!(col_vars == time_var)])
  if (!all(diff(ncols_vars) == 0)) {
    stop(paste0("The following variables have unequal maximum number of elements per ",
      id_var, ": ", paste0(col_vars[!(col_vars == time_var)], sep="; ")))
  } else {
    nz <- ncols_vars[1]
  }

  # create list of matrices for covariates/time matrices provided in func
  hist_mats <- list()
  for(i in seq_along(col_vars)) {
    hist_mats[[i]] <- if(col_vars[i] == attr(data, "time_var")) {
      make_time_mat(data, nz)
    } else if (col_vars[i] == func$latency_var) {
      make_latency_mat(data, te)
    } else {
      make_z_mat(data, col_vars[i], nz)
    }
  }
  if (is.null(func$suffix) & n_func == 1) {
    suffix <- ""
  } else {
    suffix <- func$suffix
  }
  if(any(c(time_var, te_var) %in% col_vars)) {
      hist_mats <- c(hist_mats, list(make_lag_lead_mat(data, te, func$ll_fun)))
      names(hist_mats) <- make_mat_names(c(col_vars, "LL"), func$latency_var, te_var, suffix)
  } else {
       names(hist_mats) <- make_mat_names(col_vars, func$latency_var, te_var, suffix)
  }

  hist_mats

}
