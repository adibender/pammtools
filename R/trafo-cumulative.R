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
