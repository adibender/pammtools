#' Fit a PAMM to interval-censored data via multiple imputation
#'
#' Fits a piecewise exponential additive (mixed) model to interval-censored
#' time-to-event data using a multiple-imputation (MI) and re-fit strategy: exact
#' event times are repeatedly drawn from the model-based conditional distribution
#' \eqn{p(T \mid L < T \le R, x, \theta)} (see \code{\link{impute_ic_times}}),
#' with \eqn{\theta} drawn from the imputation model's asymptotic posterior
#' before each imputation ("proper" MI -- this is what makes the pooled
#' intervals calibrated),
#' each completed data set is transformed to PED format with the standard
#' (right-censored) pipeline and re-fit, and the resulting fits are pooled for
#' inference with the existing \code{add_*} family (see \code{\link{add_surv_prob}}
#' and the \code{pamm_ic} methods).
#'
#' An imputed event time is an exact event time, so once imputation has produced
#' it, the entire downstream pipeline (\code{\link{split_data}} -> \code{\link{pamm}}
#' -> \code{add_*}) is reused unchanged. The interval cut-points are resolved once
#' and shared across all imputations, but \code{mgcv}'s smooth bases and
#' centering constraints can still differ across completed data sets. Pooled
#' predictions therefore evaluate each fitted imputation model with its own
#' design matrix; \code{object$pooled} is a summary container, not a
#' \code{gam}-like model for direct \code{predict()} or \code{plot()} calls.
#'
#' @param formula A two-sided formula whose left-hand side is an interval-censored
#'   response \code{Surv(L, R, type = "interval2")} and whose right-hand side lists
#'   the covariates to retain (as in \code{\link{as_ped}}).
#' @param data A data frame in standard (one row per subject) format.
#' @param model_formula Optional model formula passed to \code{\link{pamm}}
#'   (e.g.\ \code{ped_status ~ s(tend) + x}). If \code{NULL}, a default
#'   \code{ped_status ~ s(tend) + <covariates>} formula is constructed.
#' @param cut Optional fixed vector of interval cut-points shared across all
#'   imputations. If \code{NULL}, the finite interval endpoints are used.
#' @param max_time Optional cap on the cut-points.
#' @param m Number of imputations (default 10).
#' @param iter Number of impute-refit iterations per imputation chain (default
#'   \code{1} = classic one-step MI: all \code{m} imputations are drawn from the
#'   single initialiser fit). For \code{iter = k > 1}, each chain alternates
#'   imputation and re-fitting on its own completed data set \code{k} times, so
#'   later imputations are drawn from fits whose dependence on the midpoint
#'   initialiser is progressively attenuated -- a sequential ("chained") MI
#'   scheme that
#'   progressively removes initialiser bias under sparse inspection, at roughly
#'   \code{iter}-fold fitting cost. Simulation evidence (see the package's
#'   interval-censoring benchmark): with inspection gaps that are small
#'   relative to the time scale, \code{iter = 1} is unbiased; with wide gaps
#'   (mean gap of order 1/3 of the follow-up), early-time survival estimates
#'   from \code{iter = 1} are biased upward and \code{iter = 3} removes most
#'   of that bias (\code{iter = 5} essentially all of it), with bias shrinking
#'   roughly geometrically in \code{iter}. Caveat: with flexible time-varying
#'   effect terms and small samples, iterating can occasionally amplify a
#'   weakly identified imputation chain into divergent estimates with very
#'   wide intervals (without \code{mgcv} warnings) -- inspect pooled smooth
#'   effects for plausibility when iterating such models.
#' @param init Initialiser for the first fit: \code{"midpoint"} (default) or
#'   \code{"uniform"} imputation within each interval.
#' @param id Name of the subject identifier column.
#' @param engine Estimation engine passed to \code{\link{pamm}} (\code{"gam"} or
#'   \code{"bam"}).
#' @param ... Further arguments passed to \code{\link{pamm}} / \code{mgcv}.
#' @return An object of class \code{pamm_ic}: a list with
#'   \describe{
#'     \item{\code{fits}}{the \code{m} imputation fits, each \emph{slimmed} (via
#'       \code{\link{strip_pamm_fit}}) to drop per-observation slots so memory
#'       does not scale with the number of imputations; they still support
#'       \code{coef}, \code{vcov} and \code{predict(type = "lpmatrix")}, which is
#'       all the pooled \code{add_*} methods need.}
#'     \item{\code{pooled}}{a pooled summary container with Rubin-pooled
#'       parametric coefficients and covariance, pooled parametric/smooth tables
#'       with median-p values (\code{$p.table}, \code{$s.table}), parametric
#'       coefficient FMI diagnostics (\code{$fmi.table}) and smooth-term FMI
#'       five-number summaries over the training grid (\code{$smooth.fmi}).}
#'     \item{\code{init_fit}}{the (slimmed) initialiser/imputation model.}
#'     \item{\code{unstable_chains}}{indices of imputation chains flagged as
#'       numerically unstable (extreme coefficients or coefficient SEs on the
#'       log-hazard scale; also raised as a \code{warning}). Degenerate chains
#'       can arise -- silently, without \code{mgcv} warnings -- when iterating
#'       flexible time-varying models on small samples.}
#'     \item{others}{the parsed bounds \code{ic}, the shared \code{cut}, and
#'       metadata.}
#'   }
#'   \code{print}/\code{summary} report the pooled summary; \code{add_*} compute
#'   pooled quantities of interest from \code{fits}.
#' @seealso \code{\link{impute_ic_times}}, \code{\link{add_surv_prob}},
#'   \code{\link{strip_pamm_fit}}
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef as.formula
#' @export
pamm_ic <- function(
  formula,
  data,
  model_formula = NULL,
  cut = NULL,
  max_time = NULL,
  m = 10L,
  iter = 1L,
  init = c("midpoint", "uniform"),
  id = "id",
  engine = "gam",
  ...
) {
  init <- match.arg(init)
  assert_count(m, positive = TRUE)
  assert_count(iter, positive = TRUE)

  ped0 <- as_ped_ic(data, formula, cut = cut, max_time = max_time, id = id)
  ic <- attr(ped0, "ic")
  cut <- attr(ped0, "breaks")

  if (init == "uniform") {
    # replace midpoint initialiser by a single uniform draw
    t_unif <- draw_uniform_ic(ic, cut)
    ped0 <- build_ic_ped(ic, t_unif, cut, formula, id)
  }

  if (is.null(model_formula)) {
    model_formula <- default_pamm_formula(formula, data = data, id = id)
  }
  fit0 <- pamm(model_formula, data = ped0, engine = engine, ...)
  cache <- ic_pred_cache(fit0, ic, cut)

  fits <- vector("list", m)
  smry <- vector("list", m)
  skeleton <- NULL
  n_obs <- NA_integer_
  for (mm in seq_len(m)) {
    # chained MI: iteration k > 1 re-imputes from this chain's own re-fit,
    # progressively attenuating the midpoint initialiser's influence.
    # iter = 1 reproduces classic one-step MI draw-for-draw.
    fit_mm <- fit0
    cache_mm <- cache
    ped_m <- NULL
    for (k in seq_len(iter)) {
      # proper MI: draw the imputation-model coefficients from their posterior
      beta_mm <- as.numeric(
        rmvnorm(1, mean = coef(fit_mm), sigma = fit_mm[["Vp"]])
      )
      t_imp <- impute_ic_times(
        fit_mm,
        ic,
        cut,
        beta = beta_mm,
        cache = cache_mm
      )
      ped_m <- build_ic_ped(ic, t_imp, cut, formula, id)
      if (k < iter) {
        fit_mm <- pamm(model_formula, data = ped_m, engine = engine, ...)
        cache_mm <- ic_pred_cache(fit_mm, ic, cut)
      }
    }
    fs <- fit_strip_summarise(model_formula, ped_m, engine, ...)
    fits[[mm]] <- fs[["fit"]]
    smry[[mm]] <- fs[["summary"]]
    if (mm == 1L) skeleton <- fs[["full"]]
    n_obs <- fs[["n"]]
  }

  structure(
    list(
      fits = fits,
      pooled = pool_pamm_fits(fits, smry, skeleton = skeleton),
      init_fit = strip_pamm_fit(fit0),
      ic = ic,
      cut = cut,
      formula = formula,
      model_formula = model_formula,
      m = m,
      iter = iter,
      id_var = id,
      n_obs = n_obs,
      unstable_chains = warn_unstable_chains(fits),
      type = "single"
    ),
    class = c("pamm_ic", "list")
  )
}

#' Fit a competing-risks PAMM to interval-censored data via multiple imputation
#'
#' Competing-risks extension of \code{\link{pamm_ic}}. The event time is drawn
#' from the all-cause conditional hazard within \eqn{(L, R]} and a cause is
#' assigned: observed causes are retained (with the time drawn so that it follows
#' the cause-specific conditional density, via rejection), unknown causes are
#' sampled with probability proportional to the cause-specific hazards at the
#' imputed time (see \code{\link{impute_ic_cr}}). Each completed data set is
#' transformed with \code{\link{as_ped_cr}} (cause-specific hazards) and re-fit.
#' Cf. Delord & Genin (2016) for MI of interval-censored competing-risks data.
#'
#' @inheritParams pamm_ic
#' @param cause Name of the column in \code{data} giving the observed cause for
#'   events (any factor/character coding). Rows with the censoring code are
#'   treated as right-censored; \code{NA} marks an event with unknown cause.
#' @param censor_code Value of \code{cause} that encodes censoring (default 0).
#' @return An object of class \code{pamm_ic} with \code{type = "cr"}; \code{fits}
#'   are cause-specific (stacked \code{ped_cr}) \code{pamm} objects and
#'   \code{cause_levels} records the competing causes.
#' @seealso \code{\link{pamm_ic}}, \code{\link{add_cif}}
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef as.formula
#' @export
pamm_ic_cr <- function(
  formula,
  data,
  cause,
  model_formula = NULL,
  cut = NULL,
  max_time = NULL,
  m = 10L,
  iter = 1L,
  censor_code = 0L,
  id = "id",
  engine = "gam",
  ...
) {
  assert_count(m, positive = TRUE)
  assert_count(iter, positive = TRUE)
  assert_string(cause)
  assert_subset(cause, names(data))

  ic <- parse_ic_surv(formula, data, id = id)
  cut <- resolve_ic_cut(ic, cut = cut, max_time = max_time)

  cause_raw <- data[[cause]]
  is_cens <- !is.na(cause_raw) &
    as.character(cause_raw) == as.character(censor_code)
  # for right-censored survival rows, force censoring regardless of cause column
  is_cens <- is_cens | as.character(ic[["ic_kind"]]) == "right"
  cause_known <- ifelse(is_cens, NA, as.character(cause_raw))
  cause_levels <- sort(unique(stats::na.omit(cause_known)))
  if (length(cause_levels) < 2) {
    stop(
      "Fewer than two competing causes found in `",
      cause,
      "`.",
      call. = FALSE
    )
  }

  if (is.null(model_formula)) {
    model_formula <- default_pamm_formula(
      formula,
      data = data,
      id = id,
      by_cause = TRUE,
      exclude = cause
    )
  }

  # initialiser: midpoint times + observed causes (unknown causes drawn from the
  # marginal cause distribution) -> stacked cause-specific PED
  L <- ic[["ic_L"]]
  R <- ic[["ic_R"]]
  t_mid <- pmin(ifelse(ic[["ic_kind"]] == "left", R / 2, (L + R) / 2), max(cut))
  cause0 <- cause_known
  unknown <- which(!is_cens & is.na(cause0))
  if (length(unknown)) {
    cause0[unknown] <- sample(cause_levels, length(unknown), replace = TRUE)
  }
  ped0 <- build_ic_ped_cr(
    ic,
    t_mid,
    cause0,
    is_cens,
    cut,
    formula,
    id,
    cause_levels,
    censor_code,
    cause_var = cause
  )
  fit0 <- pamm(model_formula, data = ped0, engine = engine, ...)
  cache <- ic_pred_cache(fit0, ic, cut, cause_levels = cause_levels)

  fits <- vector("list", m)
  smry <- vector("list", m)
  skeleton <- NULL
  n_obs <- NA_integer_
  for (mm in seq_len(m)) {
    # chained MI (see pamm_ic): iter = 1 reproduces one-step MI draw-for-draw
    fit_mm <- fit0
    cache_mm <- cache
    ped_m <- NULL
    for (k in seq_len(iter)) {
      # proper MI: draw the imputation-model coefficients from their posterior
      beta_mm <- as.numeric(
        rmvnorm(1, mean = coef(fit_mm), sigma = fit_mm[["Vp"]])
      )
      imp <- impute_ic_cr(
        fit_mm,
        ic,
        cut,
        beta = beta_mm,
        cache = cache_mm,
        cause_known = cause_known
      )
      ped_m <- build_ic_ped_cr(
        ic,
        imp[["time"]],
        imp[["cause"]],
        is_cens,
        cut,
        formula,
        id,
        cause_levels,
        censor_code,
        cause_var = cause
      )
      if (k < iter) {
        fit_mm <- pamm(model_formula, data = ped_m, engine = engine, ...)
        cache_mm <- ic_pred_cache(fit_mm, ic, cut, cause_levels = cause_levels)
      }
    }
    fs <- fit_strip_summarise(model_formula, ped_m, engine, ...)
    fits[[mm]] <- fs[["fit"]]
    smry[[mm]] <- fs[["summary"]]
    if (mm == 1L) skeleton <- fs[["full"]]
    n_obs <- fs[["n"]]
  }

  structure(
    list(
      fits = fits,
      pooled = pool_pamm_fits(fits, smry, skeleton = skeleton),
      init_fit = strip_pamm_fit(fit0),
      ic = ic,
      cut = cut,
      formula = formula,
      model_formula = model_formula,
      m = m,
      iter = iter,
      id_var = id,
      cause_levels = cause_levels,
      n_obs = n_obs,
      unstable_chains = warn_unstable_chains(fits),
      type = "cr"
    ),
    class = c("pamm_ic", "list")
  )
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# Stability guard for imputation chains. Iterating the impute-refit cycle can
# occasionally amplify a weakly identified chain (flexible time-varying terms,
# small n) into a degenerate fit with extreme coefficients and vacuous CIs --
# WITHOUT any mgcv warning (found in the package's IC benchmark, Gate-R3-ext).
# Flag chains whose log-hazard-scale coefficients or coefficient SEs are
# beyond any plausible magnitude. Returns the offending chain indices
# (integer(0) if none) and warns once.
warn_unstable_chains <- function(fits, coef_limit = 20, se_limit = 10) {
  unstable <- which(vapply(
    fits,
    function(f) {
      cf <- stats::coef(f)
      vp <- f[["Vp"]]
      max(abs(cf), na.rm = TRUE) > coef_limit ||
        (!is.null(vp) && sqrt(max(diag(vp), na.rm = TRUE)) > se_limit)
    },
    logical(1)
  ))
  if (length(unstable)) {
    warning(
      "Imputation chain(s) ",
      paste(unstable, collapse = ", "),
      " look numerically unstable (|coefficient| > ",
      coef_limit,
      " or coefficient SE > ",
      se_limit,
      " on the log-hazard scale). Pooled estimates may be dominated by ",
      "degenerate chains with vacuous intervals; inspect the affected fits ",
      "and consider a less flexible model_formula or fewer `iter`ations.",
      call. = FALSE
    )
  }
  unstable
}

# Construct a default model formula `ped_status ~ s(tend) [+ by-cause] + covars`.
default_pamm_formula <- function(
  formula,
  data = NULL,
  id = "id",
  by_cause = FALSE,
  exclude = character()
) {
  rhs <- resolve_rhs_vars(formula, data = data, exclude = c(id, exclude))
  base <- if (by_cause) "s(tend, by = cause) + cause" else "s(tend)"
  terms <- c(base, if (length(rhs)) paste(rhs, collapse = " + "))
  stats::as.formula(paste("ped_status ~", paste(terms, collapse = " + ")))
}

# Uniform initial draw within (L, R] (R/2 for left-censored).
draw_uniform_ic <- function(ic, cut) {
  L <- ic[["ic_L"]]
  R <- pmin(ic[["ic_R"]], max(cut))
  u <- stats::runif(nrow(ic))
  ifelse(
    as.character(ic[["ic_kind"]]) %in% c("exact", "right"),
    L,
    L + u * (R - L)
  )
}

# Build a single-event PED from imputed exact times via the standard pipeline.
build_ic_ped <- function(ic, t_imp, cut, formula, id) {
  rhs_vars <- resolve_rhs_vars(formula, data = ic, exclude = id)
  evd <- drop_zero_followup(ic_event_data(ic, t_imp), warn = FALSE)
  ped_form <- stats::as.formula(
    paste0(
      "Surv(.ped_time, .ped_status) ~ ",
      paste0(unique(c(rhs_vars, id)), collapse = " + ")
    )
  )
  split_data(ped_form, data = evd, cut = cut, id = id)
}

# Build a competing-risks PED (stacked ped_cr) from imputed times and causes.
build_ic_ped_cr <- function(
  ic,
  time,
  cause,
  is_cens,
  cut,
  formula,
  id,
  cause_levels,
  censor_code,
  cause_var = NULL
) {
  rhs_vars <- resolve_rhs_vars(formula, data = ic, exclude = c(id, cause_var))
  dat <- ic
  dat[["ic_L"]] <- dat[["ic_R"]] <- dat[["ic_kind"]] <- NULL
  dat[[".ped_time"]] <- ifelse(
    as.character(ic[["ic_kind"]]) %in%
      c("exact", "right"),
    ic[["ic_L"]],
    time
  )
  status_cr <- ifelse(is_cens, censor_code, cause)
  dat[[".status_cr"]] <- factor(
    status_cr,
    levels = c(censor_code, cause_levels)
  )

  keep <- dat[[".ped_time"]] > 0
  dat <- dat[keep, , drop = FALSE]

  cr_form <- stats::as.formula(
    paste0(
      "Surv(.ped_time, .status_cr) ~ ",
      paste0(unique(c(rhs_vars, id)), collapse = " + ")
    )
  )
  as_ped(dat, formula = cr_form, cut = cut, censor_code = censor_code, id = id)
}
