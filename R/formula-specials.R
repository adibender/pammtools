#' Formula specials for defining  time-dependent covariates
#'
#' So far, two specials are implemented. \code{concurrent} is used when
#' the goal is to estimate a concurrent effect of the TDC. \code{cumulative}
#' is used when the goal is to estimate a cumulative effect of the TDC.
#'
#' @rdname specials
#' @importFrom purrr map
#' @export
#' @keywords internal
cumulative <- function(...,
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
#' @inherit cumulative
#' @keywords internal
concurrent <- function(...,
  te_var,
  ll_fun = function(t) {t == t},
  suffix = NULL) {

  vars     <- as.list(substitute(list(...)))[-1]
  vars_chr <- vars %>% map(~as.character(.)) %>% unlist()


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
#' @param formula  A formula containing \code{cumulative} specials,
#' that specify the type of cumulative effect one wants to estimate. For details
#' see the vignettes on data transformation and time-dependent covariates.
#' @importFrom purrr flatten map
#' @importFrom stats terms
#' @keywords internal
get_cumulative <- function(data, formula) {

  stopifnot(has_tdc_form(formula))

  func_list <- eval_special(formula)

  n_func <- length(func_list)
  ll_funs <- map(func_list, ~.x[["ll_fun"]])
  te_vars <- map(func_list, ~.x[["te_var"]])
  te <- map(te_vars, ~pull(data, .x) %>% unlist() %>% unique() %>% sort())

  names(te) <- names(te_vars) <- names(ll_funs) <- te_vars

  ## create matrices
  func_mats <- map(func_list, ~expand_cumulative(data=data, ., n_func=n_func)) %>%
    flatten()

  list(
    func_list = func_list,
    func_mats = func_mats,
    ll_funs   = ll_funs,
    te_vars   = te_vars,
    te        = te)

}

#' @keywords internal
eval_special <- function(formula, special="cumulative") {

  tf  <- terms(get_tdc_form(formula), specials = special)
  ind_special <- attr(tf, "specials")[[special]]
  # extract components
  if (!is.null(ind_special)) {
    terms_vec <- attr(tf, "term.labels")[ind_special]
    map(terms_vec, ~eval(expr=parse(text=.)))
  } else {
    NULL
  }

}


#' Querries on formula objects
#'
#' @rdname specials
#' @inheritParams as_ped
#' @param special The name of the special whose existence in the
#' \code{formula} should be checked
#' @keywords internal
has_special <- function(formula, special = "cumulative") {
  if(!has_tdc_form(formula)) {
    return(FALSE)
  } else {
    formula <- formula(Formula(formula), lhs=FALSE, rhs = 2)
    terms <- terms(formula, specials = special)
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

#' @rdname get_cumulative
#' @inheritParams get_cumulative
#' @param func Single evaluated \code{\link{cumulative}} term.
#' @importFrom purrr map invoke_map
#' @keywords internal
expand_cumulative <- function(data, func, n_func) {

  col_vars <- func$col_vars
  te_var   <- func$te_var
  te       <- pull(data, te_var) %>% unlist() %>% unique() %>% sort()
  time_var <- attr(data, "time_var")
  id_var   <- attr(data, "id_var")
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

  if(any(c(time_var, te_var) %in% col_vars)) {
    hist_mats <- c(hist_mats, list(make_lag_lead_mat(data, te, func$ll_fun)))
    names(hist_mats) <- make_mat_names(c(col_vars, "LL"), func$latency_var,
      te_var, func$suffix, n_func)
  } else {
    names(hist_mats) <- make_mat_names(col_vars, func$latency_var, te_var,
      func$suffix, n_func)
  }

  hist_mats

}

#' Extract information on concurrent effects
#'
#' @keywords internal
#' @param x A suitable object from which variables contained in
#' \code{formula} can be extracted.
#' @param ... Further arguments passed to methods.
prep_concurrent <- function(x, formula, ...) {
  UseMethod("prep_concurrent", x)
}

#' @rdname prep_concurrent
#' @inherit prep_concurrent
#' @keywords internal
prep_concurrent.list <- function(x, formula, ...) {

  lgl_concurrent <- has_special(formula, "concurrent")

  if(lgl_concurrent) {
    ccr_list    <- eval_special(formula, special="concurrent")
    ccr_te_vars <- map_chr(ccr_list, ~.x[["te_var"]]) %>% unique()
    ccr_time    <- map2(ccr_te_vars, x, ~get_te(.y, .x)) %>%
      reduce(union) %>% sort()
  }

  list(
    ccr_list = ccr_list,
    ccr_time = ccr_time)

}


#' @keywords internal
get_te <- function(data, te_var) {
  if (te_var %in% colnames(data)) {
    te <- pull(data, te_var) %>% unique()
  } else {
    te <- NULL
  }
  te
}

#' @keywords internal
add_concurrent <- function(ped, data, id_var) {

  ccr <- attr(data, "ccr")

  for(ccr_i in ccr[["ccr_list"]]) {
    tdc_vars_i <- ccr_i[["col_vars"]]
    te_var_i   <- ccr_i[["te_var"]]
    ccr_vars_i <- c(te_var_i, tdc_vars_i)
    ccr_i_df   <- data %>% select(one_of(c(id_var, ccr_vars_i))) %>%
      unnest()
    ped <- ped %>%
      left_join(ccr_i_df, by = c(id_var, "tstart"=te_var_i)) %>%
      group_by(!!sym(id_var)) %>%
      fill(tdc_vars_i)

    attr(ped, "ccr") <- ccr

  }

  ped


}

#' @keywords internal
add_cumulative <- function(ped, data, formula) {

  func_components <- get_cumulative(data, formula)
  func_matrices <- func_components$func_mats
  for(i in seq_along(func_matrices)) {
    ped[[names(func_matrices)[i]]] <- func_matrices[[i]]
  }
  attr(ped, "func")           <- func_components$func_list
  attr(ped, "func_mat_names") <- names(func_matrices)
  attr(ped, "ll_funs")        <- func_components$ll_funs
  attr(ped, "te")             <- func_components$te
  attr(ped, "te_vars")        <- func_components$te_vars

  ped

}

#' @keywords internal
make_mat_names <- function(col_vars, latency_var=NULL, te_var=NULL, suffix=NULL,
  nfunc = 1) {

  if (!is.null(suffix)) {
    return(paste(col_vars, suffix, sep="_"))
  } else {
    if (!is.null(te_var) & nfunc > 1)  {
      te_ind <- col_vars == te_var
      col_vars[!te_ind] <- paste(col_vars[!te_ind], te_var,  sep="_")
    }
    if (!is.null(latency_var)) {
      latency_ind <- col_vars == latency_var
      col_vars[latency_ind] <- paste(col_vars[latency_ind], "latency", sep="_")
    }
  }

  return(col_vars)

}

#' Create matrix components for cumulative effects
#'
#' These functions are called internally by \code{\link{get_cumulative}} and
#' should usually not be called directly.
#' @rdname elra_matrix
#' @param data A data set (or similar) from which meta information on cut-points,
#' interval-specific time, covariates etc. can be obtained.
#'
#' @keywords internal
make_time_mat <- function(data, nz) {

  brks    <- attr(data, "breaks")
  id_tseq <- attr(data, "id_tseq")
  Tmat    <- matrix(brks[id_tseq], nrow = length(id_tseq), ncol=nz)
  Tmat

}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @keywords internal
make_latency_mat <- function(data, te) {

  time        <- attr(data, "breaks")
  id_tseq     <- attr(data, "id_tseq")
  Latency_mat <- outer(time, te, FUN = "-")
  Latency_mat[Latency_mat < 0] <- 0
  Latency_mat[id_tseq, , drop=FALSE]

}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @keywords internal
make_lag_lead_mat <- function(
  data,
  te,
  ll_fun = function(t, te) {(t >= te)}) {

  LL <- outer(attr(data, "breaks"), te, FUN=ll_fun)*1L
  LL[attr(data, "id_tseq"), , drop=FALSE]

}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @param z_var Which should be transformed into functional covariate format
#' suitable to fit cumulative effects in \code{mgcv::gam}.
#' @importFrom purrr map
#' @importFrom dplyr pull
#' @keywords internal
make_z_mat <- function(data, z_var, nz, ...) {

  te_ind <- seq_len(nz)
  Z <- map(data[[z_var]], .f=~unlist(.x)[te_ind])
  Z <- do.call(rbind, Z)
  Z[attr(data, "id_teseq"), , drop=FALSE]

 }

get_ncols <- function(data, col_vars) {
  map(col_vars, function(var) pull(data, var)) %>%
    map_int(function(z) max(map_int(z, ~ifelse(is_atomic(.), length(.), nrow(.)))))

}
