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

  tf  <- terms(formula, specials="func")
  # extract components
  terms_vec <- attr(tf, "term.labels")
  func_list <- map(terms_vec, ~eval(expr=parse(text=.)))

  ## create matrices
  map(func_list, ~expand_func(data=data, .)) %>% flatten()

}

#' @rdname get_func
#' @inheritParams get_func
#' @param func Single evaluated \code{\link{func}} term.
#' @importFrom purrr map invoke_map
#' @keywords internal
expand_func <- function(data, func) {

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
  if(any(c(time_var, te_var) %in% col_vars)) {
      hist_mats <- c(hist_mats, list(make_lag_lead_mat(data, te, func$ll_fun)))
      names(hist_mats) <- make_mat_names(c(col_vars, "LL"), func$latency_var, te_var, func$suffix)
  } else {
       names(hist_mats) <- make_mat_names(col_vars, func$latency_var, te_var, func$suffix)
  }

  hist_mats

}

make_mat_names <- function(col_vars, latency_var=NULL, te_var=NULL, suffix="") {

  if (suffix != "") {
    return(paste(col_vars, suffix, sep="_"))
  } else {
    if (!is.null(te_var))  {
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
#' These functions are called internally by \code{\link{get_func}} and
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

  time        <- attr(data, "breaks")[-1]
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

  LL <- outer(attr(data, "breaks"), te, ll_fun)*1L
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
