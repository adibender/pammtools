#' Build the subject-by-interval prediction grid used for IC imputation
#'
#' Constructs the (subjects \eqn{\times} intervals) grid on the \emph{fixed}
#' cut-points and evaluates the \code{lpmatrix} of the fitted PAMM once, so that
#' across imputations only the linear predictor (and hence the hazard) needs to
#' be recomputed for a new coefficient draw. Rows are subject-major: the first
#' \code{n_int} rows belong to subject 1, the next \code{n_int} to subject 2,
#' etc., so that \code{matrix(h, nrow = n_int)} has one column per subject.
#'
#' @param object A fitted \code{\link{pamm}} model used as imputation model.
#' @param ic A data frame as returned by \code{\link{parse_ic_surv}} (subject-
#'   level, with covariates and \code{ic_L}/\code{ic_R}/\code{ic_kind}).
#' @param cut The fixed vector of interval cut-points (shared across imputations).
#' @param cause_levels Optional character vector of competing-risk cause levels.
#'   When supplied, one \code{lpmatrix} per cause is built and returned in
#'   \code{X_list}.
#' @param cause_var Name of the cause column expected by the model.
#' @return A list with the interval information \code{ii}, \code{n_int},
#'   \code{n_sub}, and either a single design matrix \code{X} or a list
#'   \code{X_list} (competing risks).
#' @importFrom mgcv predict.gam
#' @keywords internal
ic_pred_cache <- function(
  object,
  ic,
  cut,
  cause_levels = NULL,
  cause_var = "cause"
) {
  ii <- int_info(cut)
  n_int <- nrow(ii)
  n_sub <- nrow(ic)

  grid <- ic[rep(seq_len(n_sub), each = n_int), , drop = FALSE]
  grid[["tstart"]] <- rep(ii[["tstart"]], times = n_sub)
  grid[["tend"]] <- rep(ii[["tend"]], times = n_sub)
  grid[["intlen"]] <- rep(ii[["intlen"]], times = n_sub)
  grid[["intmid"]] <- rep(ii[["intmid"]], times = n_sub)
  grid[["interval"]] <- rep(ii[["interval"]], times = n_sub)

  if (is.null(cause_levels)) {
    X <- make_X(object, newdata = grid)
    return(list(ii = ii, n_int = n_int, n_sub = n_sub, X = X, grid = grid))
  }

  X_list <- lapply(cause_levels, function(cl) {
    g <- grid
    g[[cause_var]] <- factor(cl, levels = cause_levels)
    make_X(object, newdata = g)
  })
  names(X_list) <- cause_levels
  list(
    ii = ii,
    n_int = n_int,
    n_sub = n_sub,
    X_list = X_list,
    grid = grid,
    cause_levels = cause_levels,
    cause_var = cause_var
  )
}

# Invert a (per-subject) piecewise-linear cumulative hazard at a target value.
# Hcut: (n_int + 1) x n_sub matrix of H at cut-points c_0, ..., c_J.
# cutv: c_0, ..., c_J. hm: n_int x n_sub interval hazards. Returns event times.
.invert_cumu_hazard <- function(target, s, Hcut, cutv, hm, n_int) {
  out <- numeric(length(target))
  for (k in seq_along(target)) {
    si <- s[k]
    jc <- findInterval(
      target[k],
      Hcut[, si],
      left.open = TRUE,
      rightmost.closed = TRUE
    )
    jc <- min(max(jc, 1L), n_int)
    out[k] <- cutv[jc] + (target[k] - Hcut[jc, si]) / hm[jc, si]
  }
  out
}

#' Draw event times for interval-censored subjects from the conditional hazard
#'
#' For a fitted PAMM with piecewise-constant hazard, draws
#' \eqn{T_i \sim p(T \mid L_i < T \le R_i, x_i, \theta)} by inverting the
#' cumulative-hazard increment between \eqn{L_i} and \eqn{R_i}. Exact and
#' right-censored observations are returned unchanged (right-censored subjects
#' are \emph{not} imputed: they contribute correctly as censored at \code{ic_L}).
#'
#' @inheritParams ic_pred_cache
#' @param beta Coefficient vector to evaluate the hazard at. Defaults to
#'   \code{coef(object)}; pass a posterior draw for proper multiple imputation.
#' @param cache Optional pre-built cache from \code{\link{ic_pred_cache}}; avoids
#'   recomputing the (expensive) design matrix across imputations.
#' @return Numeric vector of (possibly imputed) event times, length
#'   \code{nrow(ic)}.
#' @importFrom stats coef runif
#' @keywords internal
impute_ic_times <- function(object, ic, cut, beta = NULL, cache = NULL) {
  if (is.null(cache)) {
    cache <- ic_pred_cache(object, ic, cut)
  }
  ii <- cache[["ii"]]
  X <- cache[["X"]]
  n_int <- cache[["n_int"]]
  n_sub <- nrow(ic)
  if (is.null(beta)) {
    beta <- get_coefs(object)
  }

  h <- pmax(as.numeric(exp(X %*% beta)), .Machine$double.xmin)
  hm <- matrix(h, nrow = n_int, ncol = n_sub)
  Hend <- matrix(apply(hm * ii[["intlen"]], 2, cumsum), nrow = n_int)
  Hcut <- rbind(0, Hend) # H at c_0, ..., c_J
  cutv <- c(ii[["tstart"]][1], ii[["tend"]]) # c_0, ..., c_J

  kind <- as.character(ic[["ic_kind"]])
  t_imp <- ifelse(kind %in% c("exact", "right"), ic[["ic_L"]], NA_real_)
  idx <- which(kind %in% c("interval", "left"))

  if (length(idx)) {
    L <- ic[["ic_L"]][idx]
    R <- pmin(ic[["ic_R"]][idx], max(cut))

    eval_H <- function(t, s) {
      j <- findInterval(t, cutv, left.open = TRUE, rightmost.closed = TRUE)
      j <- pmin(pmax(j, 1L), n_int)
      Hcut[cbind(j, s)] + hm[cbind(j, s)] * (t - cutv[j])
    }

    HL <- eval_H(L, idx)
    HR <- eval_H(R, idx)
    U <- runif(length(idx))
    delta <- HR - HL
    small <- delta < 1e-8

    ti <- numeric(length(idx))
    # (near-)flat cumulative hazard: conditional distribution ~ Uniform(L, R]
    ti[small] <- L[small] + U[small] * (R[small] - L[small])

    big <- which(!small)
    if (length(big)) {
      target <- HL[big] - log1p(-U[big] * (1 - exp(-delta[big])))
      ti[big] <- .invert_cumu_hazard(target, idx[big], Hcut, cutv, hm, n_int)
    }

    # numerical safety: keep strictly within (L, R]
    ti <- pmin(pmax(ti, L + .Machine$double.eps), R)
    t_imp[idx] <- ti
  }

  t_imp
}

#' Draw event times and causes for interval-censored competing-risks subjects
#'
#' Draws the event time from the all-cause conditional hazard within
#' \eqn{(L, R]} (as in \code{\link{impute_ic_times}}) and assigns a cause. If the
#' cause is observed it is retained (the time is then drawn by a rejection step
#' so that it follows the cause-specific conditional density); if the cause is
#' unknown it is sampled with probability \eqn{h_k(T)/h_\bullet(T)} at the imputed
#' time, mirroring the cause-assignment in \code{\link{sim_pexp}} and the CIF
#' increment in \code{get_cif}.
#'
#' @inheritParams impute_ic_times
#' @param cause_known Optional vector (length \code{nrow(ic)}) of observed causes
#'   (as levels of \code{cache$cause_levels}); \code{NA} marks unknown cause.
#'   Censored and exact rows are ignored.
#' @return A list with numeric \code{time} and character \code{cause} (both
#'   length \code{nrow(ic)}; \code{cause} is \code{NA} for censored rows).
#' @importFrom stats coef runif
#' @keywords internal
impute_ic_cr <- function(
  object,
  ic,
  cut,
  beta = NULL,
  cache = NULL,
  cause_known = NULL
) {
  if (is.null(cache) || is.null(cache[["X_list"]])) {
    stop(
      "`cache` with per-cause design matrices (`X_list`) is required for ",
      "competing-risks imputation.",
      call. = FALSE
    )
  }
  ii <- cache[["ii"]]
  n_int <- cache[["n_int"]]
  n_sub <- nrow(ic)
  cause_levels <- cache[["cause_levels"]]
  if (is.null(beta)) {
    beta <- coef(object)
  }

  # cause-specific hazard matrices [interval x subject]
  hk <- lapply(cache[["X_list"]], function(Xc) {
    h <- pmax(as.numeric(exp(Xc %*% beta)), .Machine$double.xmin)
    matrix(h, nrow = n_int, ncol = n_sub)
  })
  htot <- Reduce(`+`, hk)
  Hend <- matrix(apply(htot * ii[["intlen"]], 2, cumsum), nrow = n_int)
  Hcut <- rbind(0, Hend)
  cutv <- c(ii[["tstart"]][1], ii[["tend"]])

  eval_H <- function(t, s) {
    j <- findInterval(t, cutv, left.open = TRUE, rightmost.closed = TRUE)
    j <- pmin(pmax(j, 1L), n_int)
    Hcut[cbind(j, s)] + htot[cbind(j, s)] * (t - cutv[j])
  }
  draw_time <- function(s) {
    L <- ic[["ic_L"]][s]
    R <- min(ic[["ic_R"]][s], max(cut))
    HL <- eval_H(L, s)
    HR <- eval_H(R, s)
    U <- runif(1)
    delta <- HR - HL
    if (delta < 1e-8) {
      return(L + U * (R - L))
    }
    target <- HL - log1p(-U * (1 - exp(-delta)))
    t <- .invert_cumu_hazard(target, s, Hcut, cutv, htot, n_int)
    min(max(t, L + .Machine$double.eps), R)
  }
  interval_of <- function(t) {
    j <- findInterval(t, cutv, left.open = TRUE, rightmost.closed = TRUE)
    min(max(j, 1L), n_int)
  }
  cause_probs <- function(t, s) {
    j <- interval_of(t)
    p <- vapply(hk, function(m) m[j, s], numeric(1))
    p / sum(p)
  }
  max_rejection_tries <- 1000L

  time <- rep(NA_real_, n_sub)
  cause <- rep(NA_character_, n_sub)
  to_imp <- which(as.character(ic[["ic_kind"]]) %in% c("interval", "left"))
  exact <- which(as.character(ic[["ic_kind"]]) == "exact")

  for (s in exact) {
    time[s] <- ic[["ic_L"]][s]
    known <- if (!is.null(cause_known)) cause_known[s] else NA
    cause[s] <- if (!is.na(known)) {
      as.character(known)
    } else {
      sample(cause_levels, 1L, prob = cause_probs(time[s], s))
    }
  }
  for (s in to_imp) {
    known <- if (!is.null(cause_known)) cause_known[s] else NA
    if (!is.na(known)) {
      known <- as.character(known)
      # rejection: propose from all-cause conditional, accept w.p. h_k/h_tot
      accepted <- FALSE
      for (try in seq_len(max_rejection_tries)) {
        t <- draw_time(s)
        p <- cause_probs(t, s)[known]
        if (!is.finite(p)) {
          stop(
            "Could not evaluate the cause-specific acceptance probability ",
            "for known cause `",
            known,
            "`.",
            call. = FALSE
          )
        }
        if (runif(1) <= p) {
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        stop(
          "Failed to draw an interval-censored event time for known cause `",
          known,
          "` after ",
          max_rejection_tries,
          " rejection-sampling proposals. The fitted cause-specific hazard ",
          "is too small relative to the all-cause hazard in the censoring ",
          "interval.",
          call. = FALSE
        )
      }
      time[s] <- t
      cause[s] <- known
    } else {
      t <- draw_time(s)
      time[s] <- t
      cause[s] <- sample(cause_levels, 1L, prob = cause_probs(t, s))
    }
  }

  list(time = time, cause = cause)
}
