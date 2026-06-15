#' Pooling of multiple-imputation PAMM fits
#'
#' Inference for interval-censored PAMMs (\code{\link{pamm_ic}}) pools the
#' \code{m} re-fits by drawing from each fit's empirical-Bayes posterior
#' \eqn{N(\hat\beta^{(m)}, V_\beta^{(m)})} and propagating every draw through the
#' quantity of interest using that fit's \emph{own} design matrix, then taking
#' empirical quantiles of the combined draws. Because mgcv's identifiability
#' constraints make the (centered) spline basis depend on each imputed data set,
#' the design matrix is \emph{not} shared across fits, so each fit must be
#' evaluated with its own \code{lpmatrix}. Before empirical quantiles are taken,
#' the per-fit prediction draws are shifted on the quantity-of-interest scale so
#' their between-imputation component has Rubin's finite-\code{m} variance
#' \eqn{(1 + 1/M)B} rather than the raw mixture variance
#' \eqn{(M - 1)B/M}. Point estimates are the average of the per-fit point
#' estimates (the MI estimate).
#'
#' These methods are dispatched automatically by \code{\link{add_hazard}},
#' \code{\link{add_cumu_hazard}}, \code{\link{add_surv_prob}} and
#' \code{\link{add_cif}} when given a \code{pamm_ic} object.
#'
#' @name pamm_ic_pooling
#' @keywords internal
NULL

# Average a per-fit point estimate (the value column produced by an `add_*`
# default with ci = FALSE) across the imputation fits -> MI point estimate.
pooled_point <- function(object, newdata, adder, value_col, ...) {
  preds <- lapply(object[["fits"]], function(f) {
    adder(newdata, f, ci = FALSE, ...)[[value_col]]
  })
  rowMeans(do.call(cbind, preds))
}

ic_prediction_grid <- function(object, newdata, time_var, interval_length) {
  fit1 <- object[["fits"]][[1]]
  tv <- resolve_time_var(time_var, fit1, newdata)
  list(
    data = reconstruct_cutpoints(newdata, fit1, tv, interval_length),
    time_var = tv
  )
}

# Per-group cumulative sum of intlen * hazard, applied column-wise to a draw
# matrix (one column per posterior draw). Rows are assumed time-ordered within
# each group, as produced by make_newdata().
ic_group_cumsum <- function(mat, intlen, grp) {
  out <- mat
  for (g in unique(grp)) {
    ix <- which(grp == g)
    out[ix, ] <- apply(mat[ix, , drop = FALSE] * intlen[ix], 2, cumsum)
  }
  out
}

# Inflate QOI-level draws so their between-imputation variance matches Rubin's
# finite-m total variance. The raw equal-weight mixture has (m - 1) / m * B; the
# scale below turns that into (1 + 1 / m) * B without assuming a common GAM basis.
rubin_inflate_qoi_draws <- function(draws, estimates) {
  m <- length(draws)
  if (m <= 1L) {
    return(draws)
  }
  qhat <- do.call(cbind, estimates)
  qbar <- rowMeans(qhat)
  between_scale <- sqrt((m + 1) / (m - 1))

  lapply(seq_len(m), function(i) {
    shift <- (between_scale - 1) * (qhat[, i] - qbar)
    sweep(draws[[i]], 1L, shift, `+`)
  })
}

# Pooled posterior draws of hazard / cumulative hazard / survival, evaluated with
# each fit's own design matrix, returning lower/upper quantiles.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef quantile
ic_ci_draws <- function(
  object,
  newdata,
  nsim,
  kind,
  alpha,
  time_var,
  interval_length
) {
  fits <- object[["fits"]]
  m <- length(fits)
  per <- ceiling(nsim / m)
  tv <- resolve_time_var(time_var, fits[[1]], newdata)

  nd <- newdata
  intlen <- NULL
  grp <- NULL
  if (kind != "hazard") {
    nd <- reconstruct_cutpoints(nd, fits[[1]], tv, interval_length)
    intlen <- nd[[interval_length]]
    grp <- group_indices(nd)
  }

  pieces <- lapply(fits, function(f) {
    X <- make_X(f, newdata = nd)
    B <- sample_coefs(f, per)
    H <- exp(X %*% t(B)) # nrow x per hazard draws
    coefs <- get_coefs(f)
    h0 <- as.numeric(exp(X %*% coefs))
      return(list(draws = H, estimate = h0))
    }
    C <- ic_group_cumsum(H, intlen, grp)
    c0 <- ic_group_cumsum(matrix(h0, ncol = 1L), intlen, grp)[, 1L]
    if (kind == "cumu") {
      list(draws = C, estimate = c0)
    } else {
      list(draws = exp(-C), estimate = exp(-c0))
    }
  })
  M <- do.call(
    cbind,
    rubin_inflate_qoi_draws(
      lapply(pieces, `[[`, "draws"),
      lapply(pieces, `[[`, "estimate")
    )
  )

  lower <- apply(M, 1, quantile, probs = alpha / 2, na.rm = TRUE, type = 6)
  upper <- apply(M, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 6)
  if (kind %in% c("hazard", "cumu")) {
    lower <- pmax(lower, 0)
    upper <- pmax(upper, 0)
  } else if (kind == "surv") {
    lower <- pmin(pmax(lower, 0), 1)
    upper <- pmin(pmax(upper, 0), 1)
  }

  list(
    lower = lower,
    upper = upper,
    data = nd
  )
}

#' @rdname add_hazard
#' @param alpha Significance level for pooled confidence intervals (a
#'   \eqn{(1-\alpha)} interval).
#' @param nsim Total number of pooled posterior draws used for the interval.
#' @export
add_hazard.pamm_ic <- function(
  newdata,
  object,
  ci = TRUE,
  alpha = 0.05,
  nsim = 500L,
  time_var = NULL,
  ...
) {
  newdata[["hazard"]] <- pooled_point(
    object,
    newdata,
    add_hazard,
    "hazard",
    time_var = time_var,
    ...
  )
  if (ci) {
    d <- ic_ci_draws(object, newdata, nsim, "hazard", alpha, time_var, "intlen")
    newdata[["ci_lower"]] <- d[["lower"]]
    newdata[["ci_upper"]] <- d[["upper"]]
  }
  newdata
}

#' @rdname add_hazard
#' @export
add_cumu_hazard.pamm_ic <- function(
  newdata,
  object,
  ci = TRUE,
  alpha = 0.05,
  nsim = 500L,
  time_var = NULL,
  interval_length = "intlen",
  ...
) {
  grid <- ic_prediction_grid(object, newdata, time_var, interval_length)
  joindata <- grid[["data"]]
  time_var <- grid[["time_var"]]
  joindata[["cumu_hazard"]] <- pooled_point(
    object,
    joindata,
    add_cumu_hazard,
    "cumu_hazard",
    time_var = time_var,
    interval_length = interval_length,
    boundary = FALSE,
    ...
  )
  if (ci) {
    d <- ic_ci_draws(
      object,
      joindata,
      nsim,
      "cumu",
      alpha,
      time_var,
      interval_length
    )
    joindata <- d[["data"]]
    joindata[["cumu_lower"]] <- d[["lower"]]
    joindata[["cumu_upper"]] <- d[["upper"]]
  }
  suppressMessages(newdata %>% left_join(joindata))
}

#' @rdname add_surv_prob
#' @param alpha Significance level for pooled confidence intervals.
#' @param nsim Total number of pooled posterior draws used for the interval.
#' @export
add_surv_prob.pamm_ic <- function(
  newdata,
  object,
  ci = TRUE,
  alpha = 0.05,
  nsim = 500L,
  time_var = NULL,
  interval_length = "intlen",
  ...
) {
  grid <- ic_prediction_grid(object, newdata, time_var, interval_length)
  joindata <- grid[["data"]]
  time_var <- grid[["time_var"]]
  joindata[["surv_prob"]] <- pooled_point(
    object,
    joindata,
    add_surv_prob,
    "surv_prob",
    time_var = time_var,
    interval_length = interval_length,
    boundary = FALSE,
    ...
  )
  if (ci) {
    d <- ic_ci_draws(
      object,
      joindata,
      nsim,
      "surv",
      alpha,
      time_var,
      interval_length
    )
    joindata <- d[["data"]]
    joindata[["surv_lower"]] <- d[["lower"]]
    joindata[["surv_upper"]] <- d[["upper"]]
  }
  suppressMessages(newdata %>% left_join(joindata))
}

# CIF values for one cause-by-covariate group and one fit. `coef_mat` has one
# row per coefficient vector, so posterior draws and plug-in estimates share the
# same computation.
#' @importFrom stats predict
ic_cif_fit_group <- function(
  group_df,
  fit,
  coef_mat,
  cause_levels,
  cause_var,
  interval_length
) {
  cause_data <- unique(as.character(group_df[[cause_var]]))
  if (length(cause_data) > 1) {
    stop("Did you forget to group by cause?", call. = FALSE)
  }
  dt <- group_df[[interval_length]]

  hazards <- lapply(cause_levels, function(cl) {
    dfc <- group_df
    dfc[[cause_var]] <- factor(cl, levels = cause_levels)
    X <- predict(fit, dfc, type = "lpmatrix")
    exp(X %*% t(coef_mat))
  })
  names(hazards) <- as.character(cause_levels)

  total_hazard <- Reduce(`+`, hazards)
  overall_surv <- apply(total_hazard, 2, function(z) exp(-cumsum(z * dt)))
  if (!is.matrix(overall_surv)) {
    overall_surv <- matrix(overall_surv, nrow = nrow(total_hazard))
  }
  survival <- rbind(1, overall_surv[-nrow(overall_surv), , drop = FALSE])
  hazard <- hazards[[cause_data]]
  cif_inc <- (hazard / total_hazard) *
    survival *
    (1 - exp(-total_hazard * dt))
  cif <- apply(cif_inc, 2, cumsum)
  if (is.matrix(cif)) {
    cif
  } else {
    matrix(cif, ncol = nrow(coef_mat))
  }
}

# Deterministic MI point estimate: average per-imputation plug-in CIFs.
#' @importFrom stats coef
ic_cif_point_group <- function(group_df, object, cause_var, interval_length) {
  cols <- lapply(object[["fits"]], function(f) {
    ic_cif_fit_group(
      group_df,
      f,
      matrix(get_coefs(f), nrow = 1),
      object[["cause_levels"]],
      cause_var,
      interval_length
    )
  })
  rowMeans(do.call(cbind, cols))
}

# Pooled CIF draws for a single cause-by-covariate group, evaluated with each
# fit's own design matrices. Returns the [nrow x nsim] matrix of CIF draws.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef
ic_cif_draws_group <- function(
  group_df,
  object,
  per,
  cause_var,
  interval_length
) {
  pieces <- lapply(object[["fits"]], function(f) {
    B <- rmvnorm(per, mean = coef(f), sigma = f[["Vp"]])
    list(
      draws = ic_cif_fit_group(
        group_df,
        f,
        B,
        object[["cause_levels"]],
        cause_var,
        interval_length
      ),
      estimate = as.numeric(ic_cif_fit_group(
        group_df,
        f,
        matrix(coef(f), nrow = 1L),
        object[["cause_levels"]],
        cause_var,
        interval_length
      ))
    )
  })
  do.call(
    cbind,
    rubin_inflate_qoi_draws(
      lapply(pieces, `[[`, "draws"),
      lapply(pieces, `[[`, "estimate")
    )
  )
}

#' @rdname add_cif
#' @param alpha Significance level for pooled confidence intervals.
#' @param nsim Total number of pooled posterior draws used for the interval.
#' @export
add_cif.pamm_ic <- function(
  newdata,
  object,
  ci = TRUE,
  alpha = 0.05,
  nsim = 500L,
  cause_var = "cause",
  time_var = NULL,
  interval_length = "intlen",
  ...
) {
  if (!identical(object[["type"]], "cr")) {
    stop(
      "add_cif() requires a competing-risks `pamm_ic` (see pamm_ic_cr()).",
      call. = FALSE
    )
  }
  fit1 <- object[["fits"]][[1]]
  m <- length(object[["fits"]])
  per <- ceiling(nsim / m)
  time_var <- resolve_time_var(time_var, fit1, newdata)
  joindata <- reconstruct_cutpoints(newdata, fit1, time_var, interval_length)

  joindata <- map_dfr(
    split(joindata, group_indices(joindata)),
    function(.x) {
      .x <- arrange(.x, .data[[time_var]])
      .x[["cif"]] <- pmin(
        pmax(
          ic_cif_point_group(.x, object, cause_var, interval_length),
          0
        ),
        1
      )
      if (ci) {
        cif <- ic_cif_draws_group(.x, object, per, cause_var, interval_length)
        .x[["cif_lower"]] <- pmin(
          pmax(
            apply(cif, 1, quantile, alpha / 2, na.rm = TRUE, type = 6),
            0
          ),
          1
        )
        .x[["cif_upper"]] <- pmin(
          pmax(
            apply(cif, 1, quantile, 1 - alpha / 2, na.rm = TRUE, type = 6),
            0
          ),
          1
        )
      }
      .x
    }
  )

  suppressMessages(newdata %>% left_join(joindata))
}
