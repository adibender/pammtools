#' Calculate (or plot) cumulative effect for all time-points of the follow-up
#'
#' @inheritParams gg_partial
#' @param z1 The exposure profile for which to calculate the cumulative effect.
#' Can be either a single number or a vector of same length as unique observation
#' time points.
#' @param z2 If provided, calculated cumulative effect is for the difference
#' between the two exposure profiles (g(z1,t)-g(z2,t)).
#' @param se_mult Multiplicative factor used to calculate confidence intervals
#' (e.g., lower = fit - 2*se).
#' @export
get_cumu_eff <- function(data, model, term, z1, z2 = NULL, se_mult = 2) {

  assert_class(data, "fped")
  ped     <- make_ped_dat(data, term, z1)
  coefs   <- coef(model)
  col_ind <- grep(term, names(coefs))
  coefs   <- coefs[col_ind]
  Vp      <- model$Vp[col_ind, col_ind]
  X       <- predict(model, ped, type = "lpmatrix")[, col_ind]
  if (!is.null(z2)) {
    X2 <- predict(model, make_ped_dat(data, term, z2),
      type = "lpmatrix")[, col_ind]
    X <- X - X2
  }
  ped$cumu_eff <- drop(X %*% coefs)
  ped$se_cumu_eff <- drop(sqrt(rowSums( (X %*% Vp) * X) ))
  ped$cumu_eff_lower <- ped$cumu_eff - se_mult * ped$se_cumu_eff
  ped$cumu_eff_upper <- ped$cumu_eff + se_mult * ped$se_cumu_eff

  ped

}




#' @keywords internal
make_ped_dat <- function(x, term, z_vec) {

  nfunc <- length(attr(x, "ll_funs"))
  ind_term <- get_term_ind(x, term)
  nz       <- length(attr(x, "tz")[[ind_term]])
  tz_var   <- attr(x, "tz_vars")[[ind_term]]
  tz       <- attr(x, "tz")[[ind_term]]
  func <- attr(x, "func")[[ind_term]]
  ll_fun <- attr(x, "ll_funs")[[ind_term]]
  func_mat_names <- attr(x, "func_mat_names")[[ind_term]]
  LL_name <- grep("LL", func_mat_names, value = TRUE)
  tz_var_mat <- make_mat_names(tz_var, func$latency_var, func$tz_var,
    func$suffix, nfunc)
  q_weights <- attr(x, "ll_weights")[[ind_term]]
  stopifnot(length(z_vec) == nz | length(z_vec) == 1)

  z_vec <- if (length(z_vec) == 1)  {
    rep(z_vec, nz)
  } else {
    z_vec
  }

  ped_df <- make_newdata(x, tend = unique(.data$tend))
  ped_df[[LL_name]] <- outer(ped_df$tend, tz, FUN = ll_fun) * 1L *
    matrix(q_weights$ll_weight, nrow = nrow(ped_df), ncol = nz, byrow = TRUE)
  if (func$latency_var != "") {
    ped_df[[tz_var_mat]] <- outer(ped_df$tend, tz, FUN = "-")
    ped_df[[tz_var_mat]] * (ped_df[[LL_name]] != 0)
  } else {
    ped_df[[tz_var]] <- matrix(tz, nrow = nrow(ped_df), ncol = nz, byrow = TRUE)
    ped_df[[tz_var]] <- ped_df[[tz_var]] * (ped_df[[LL_name]] != 0)
  }
  ped_df[[term]] <- matrix(z_vec, nrow = nrow(ped_df), ncol = nz, byrow = TRUE)
  t_mat_var <- grep(attr(x, "time_var"), func_mat_names, value = TRUE)
  if (length(t_mat_var) != 0) {
    ped_df[[t_mat_var]] <- matrix(unique(x[[t_mat_var]][, 1]),
      nrow = nrow(ped_df), ncol = nz)
  }

  ped_df

}


get_term_ind <- function(x, term) {
  which(map_lgl(attr(x, "func_mat_names"), ~any(grepl(term, .x))))
}
