# #' @keywords internal
# as_historical <- function(ntdc_df, tdc_var, idx_rep) {

#   tdc_col <- attr(ntdc_df, "tdc_col")
#   te_var  <- attr(ntdc_df, "te_var")
#   tdc_df  <- ntdc_df %>%
#     unnest() %>%
#     select(!!id_var, !!te_var, !!tdc_var)

# }


#' Transform nested to functional piece-wise exponential data format
#'
#' @param ped Piece-wise Exponential Data (PED) with time-dependet covariates.
#' @inheritParams split_tdc
#' @keywords internal
as_fped <- function(ped, formula) {

  # func_mat_list <- as_func(ped, formula)
  # for(i in seq_along(func_mat_list)) {
  #   ped[[names(func_mat_list)[i]]] <- I(func_mat_list[[i]])
  # }

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

  tf        <- terms(formula, specials="func")
  # extract components
  terms_vec <- attr(tf, "term.labels")
  func_list <- map(terms_vec, ~eval(expr=parse(text=.)))
  ll_funs <- map(func_list, ~.[["ll_fun"]]) # hier mit unique und duplicated arbeiten.
  time_components <- map(func_list, get_components) %>% transpose() %>%
    map_lgl(~any(unlist(.)))

  ## create matrices
  # time components
  time_mats <- invoke_map(
    .f   = list(make_time_mat, make_te_mat, make_latency_mat)[time_components],
    data = data)
  names(time_mats) <- c("Tmat", "TEmat", "Latency")[time_components]
  # lag-lead
  LL_mats <- map(unique(ll_funs), make_lag_lead_mat, data=data)
  names(LL_mats) <- get_ll_names(ll_funs)
  # covariates
  func_covars <- map(func_list, expand_func, data=data) %>% flatten()

  list(time_mats, LL_mats, func_covars) %>% flatten()

}

#' @rdname get_func
#' @inheritParams get_func
#' @param func Single \code{\link{func}} term \code{\link{get_func}} \code{formula}
#' evaluated.
#' @importFrom purrr map invoke_map
#' @keywords internal
expand_func <- function(data, func) {

  time_components  <- get_components(func)
  fun_covars       <- func %>% get_fun_covars()
  func_mats        <- map(fun_covars, ~make_z_mat(data=data, z_var = .))
  names(func_mats) <- fun_covars

  return(func_mats)

}

#' Extract time-scales information from \code{func} object
#'
#' Given \code{func} object, looks up which of the three time-components
#' (follow-up time, time of exposure, latency) are needed
#'
#' @inheritParams expand_func
#' @importFrom purrr map_lgl
#' @keywords internal
get_components <- function(func) {
  map_lgl(func[c("lgl_time", "lgl_te", "lgl_latency")], any)
}

#' @inherit get_components
#' @importFrom purrr flatten_chr
#' @keywords internal
get_fun_covars <- function(func) {
   func[c("fun_covar", "by_var")] %>% flatten_chr()
}

#' @inherit get_components
#' @keywords internal
get_vname <- function(func) {
  get_fun_covars(func) %>% paste0(collapse=".")
}

#' @keywords internal
get_ll_names <- function(ll_funs) {

  if (sum(!duplicated(ll_funs)) > 1) {
    suffix <- seq_len(sum(!duplicated(ll_funs)))
  } else {
    suffix <- ""
  }
  paste0("LL", suffix)

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
make_time_mat <- function(data) {
  brks    <- attr(data, "breaks")
  id_tseq <- attr(data, "id_tseq")
  te      <- attr(data, "te")
  Tmat    <- matrix(brks[id_tseq], nrow = length(id_tseq), ncol=length(te))
  Tmat
}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @keywords internal
make_te_mat <- function(data) {

    te      <- attr(data, "te")
    id_tseq <- attr(data, "id_tseq")
    TEmat   <- matrix(te, nrow=length(id_tseq), ncol=length(te), byrow=TRUE)
    TEmat
}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @keywords internal
make_latency_mat <- function(data) {

  time        <- attr(data, "breaks")[-1]
  id_tseq     <- attr(data, "id_tseq")
  te          <- attr(data, "te")
  Latency_mat <- outer(time, te, "-")
  Latency_mat[Latency_mat < 0] <- 0
  Latency_mat[id_tseq, ]
}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @inheritParams as_fped
#' @keywords internal
make_lag_lead_mat <- function(
  data,
  ll_fun = function(t, te) {(t >= te)}) {

  LL <- outer(attr(data, "breaks"), attr(data, "te"), ll_fun)*1L
  LL[attr(data, "id_tseq"), ]

}

#' @rdname elra_matrix
#' @inherit make_time_mat
#' @param z_var Which should be transformed into functional covariate format
#' suitable to fit cumulative effects in \code{mgcv::gam}.
#' @importFrom purrr map
#' @importFrom dplyr pull
#' @keywords internal
make_z_mat <- function(data, z_var) {

  te_ind <- seq_along(attr(data, "te"))
  Z <- map(data[[attr(data, "tdc_col")]], ~pull(.x, z_var)) %>%
    map(.f=~.x[te_ind])
  Z <- do.call(rbind, Z)
  Z[attr(data, "id_teseq"), ]

 }
