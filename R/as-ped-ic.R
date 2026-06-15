#' Detect, parse and transform interval-censored survival data
#'
#' Interval-censored (IC) data record the event time of subject \eqn{i} only up
#' to an interval \eqn{(L_i, R_i]}. \pkg{pammtools} handles such data via
#' multiple imputation (MI): exact event times are repeatedly drawn from the
#' model-based conditional distribution and the resulting (exact) data sets are
#' transformed and re-fit using the standard right-censored PAMM pipeline (see
#' \code{\link{pamm_ic}}). The functions documented here implement the
#' preprocessing building blocks of that workflow.
#'
#' IC data are specified through the standard \pkg{survival} interface, i.e. a
#' three-argument response of the form \code{Surv(L, R, type = "interval2")}.
#' The four observation types are encoded as usual:
#' \describe{
#'   \item{exact}{\eqn{L = R} (known event time).}
#'   \item{right-censored}{\eqn{R = \infty} (event after \eqn{L}).}
#'   \item{left-censored}{\eqn{L = 0} (event in \eqn{(0, R]}).}
#'   \item{interval-censored}{\eqn{0 < L < R < \infty} (event in \eqn{(L, R]}).}
#' }
#'
#' @param formula A two-sided formula whose left-hand side is an interval-
#'   censored \code{\link[survival]{Surv}} object.
#' @param data A data frame containing the variables referenced in
#'   \code{formula}.
#' @name interval_censored
#' @keywords internal
NULL

# Evaluation environment for the Surv() LHS that guarantees `Surv` resolves even
# when survival is imported but not attached (the usual case for a package user).
#' @importFrom survival Surv
ic_eval_env <- function(formula) {
  parent <- environment(formula)
  if (is.null(parent)) {
    parent <- parent.frame(2L)
  }
  e <- new.env(parent = parent)
  e[["Surv"]] <- survival::Surv
  e
}

#' @describeIn interval_censored Detect whether \code{formula} specifies
#'   interval-censored data. Returns \code{"interval2"} for interval-censored
#'   responses and \code{"none"} otherwise (right-censored and left-truncated
#'   counting-process responses both return \code{"none"} and are handled by the
#'   standard pipeline).
#' @keywords internal
detect_ic <- function(formula, data) {
  lhs <- tryCatch(rlang::f_lhs(formula), error = function(e) NULL)
  if (is.null(lhs)) {
    return("none")
  }
  so <- tryCatch(
    suppressWarnings(eval(lhs, envir = data, enclos = ic_eval_env(formula))),
    error = function(e) NULL
  )
  if (is.null(so) || !inherits(so, "Surv")) {
    return("none")
  }
  if (isTRUE(attr(so, "type") == "interval")) {
    return("interval2")
  }

  "none"
}

#' @describeIn interval_censored Parse the interval-censored response into a
#'   tibble of lower/upper bounds and observation type, augmenting \code{data}
#'   with the columns \code{ic_L}, \code{ic_R} and \code{ic_kind} (a factor with
#'   levels \code{exact}, \code{right}, \code{left}, \code{interval}) and, if
#'   absent, an \code{id} column.
#' @param id Name of the subject identifier column. If it does not exist in
#'   \code{data} it is created as a row index.
#' @import checkmate
#' @keywords internal
parse_ic_surv <- function(formula, data, id = "id") {
  lhs <- rlang::f_lhs(formula)
  so <- eval(lhs, envir = data, enclos = ic_eval_env(formula))
  if (!inherits(so, "Surv") || !isTRUE(attr(so, "type") == "interval")) {
    stop(
      "`formula` does not specify an interval-censored `Surv` response. ",
      "Use `Surv(L, R, type = \"interval2\")`.",
      call. = FALSE
    )
  }

  m <- as.matrix(so)
  status <- as.integer(m[, "status"])
  t1 <- as.numeric(m[, "time1"])
  t2 <- as.numeric(m[, "time2"])

  # survival status codes: 0 right, 1 exact, 2 left, 3 interval
  ic_L <- ifelse(status == 2L, 0, t1)
  ic_R <- ifelse(
    status == 0L,
    Inf,
    ifelse(status == 1L, t1, ifelse(status == 2L, t1, t2))
  )
  ic_kind <- factor(
    c("right", "exact", "left", "interval")[status + 1L],
    levels = c("exact", "right", "left", "interval")
  )

  if (any(is.na(ic_L)) || any(is.na(ic_R[status != 0L]))) {
    stop(
      "Interval bounds contain missing values after parsing the response.",
      call. = FALSE
    )
  }
  if (any(ic_R <= ic_L & ic_kind == "interval")) {
    stop("Found interval-censored observations with R <= L.", call. = FALSE)
  }

  out <- data
  if (!id %in% names(out)) {
    out[[id]] <- seq_len(nrow(out))
  }
  out[["ic_L"]] <- ic_L
  out[["ic_R"]] <- ic_R
  out[["ic_kind"]] <- ic_kind

  out
}

#' @describeIn interval_censored Resolve a fixed vector of interval cut-points
#'   for the IC transformation. When \code{cut} is supplied it is sorted and
#'   de-duplicated; otherwise the unique finite interval endpoints (the
#'   inspection times) are used, capped at \code{max_time}. The resolved
#'   \code{cut} \emph{must} be shared across all imputations so that the PED
#'   interval structure is consistent across refits. Note that \code{mgcv}'s
#'   centering constraints can still make the \code{lpmatrix} differ by fit.
#' @param cut Optional numeric vector of interval cut-points. If \code{NULL}
#'   (default) the finite interval endpoints are used.
#' @param max_time Optional numeric scalar; cut-points are capped at this value.
#' @keywords internal
resolve_ic_cut <- function(ic, cut = NULL, max_time = NULL) {
  if (!is.null(cut)) {
    cut <- sort(unique(cut))
    validate_ic_cut_covers_events(ic, cut)
    return(cut)
  }
  ep <- c(ic[["ic_L"]], ic[["ic_R"]])
  ep <- ep[is.finite(ep) & ep > 0]
  cut <- sort(unique(c(0, ep)))
  if (!is.null(max_time)) {
    cut <- cut[cut <= max_time]
  }
  if (length(cut) < 2) {
    stop(
      "Could not derive interval cut-points from the data; please supply ",
      "`cut` explicitly.",
      call. = FALSE
    )
  }
  validate_ic_cut_covers_events(ic, cut)

  cut
}

validate_ic_cut_covers_events <- function(ic, cut) {
  event <- as.character(ic[["ic_kind"]]) != "right"
  upper <- ic[["ic_R"]][event & is.finite(ic[["ic_R"]])]
  if (length(upper) && any(upper > max(cut))) {
    stop(
      "The IC cut-points must cover all finite event upper bounds. ",
      "Increase `cut`/`max_time` or recode events beyond the analysis horizon ",
      "as right-censored.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @describeIn interval_censored Build the subject-level data frame of (exact)
#'   event times implied by an imputation. Exact observations keep their event
#'   time; right-censored observations are censored at \code{ic_L}; left- and
#'   interval-censored observations take the imputed time \code{t_imp}. Returns a
#'   data frame with the response columns \code{.ped_time} and \code{.ped_status}
#'   ready for \code{\link{split_data}} via a two-argument \code{Surv}.
#' @param ic A data frame as returned by \code{parse_ic_surv}.
#' @param t_imp Numeric vector of imputed event times (length \code{nrow(ic)}).
#'   Ignored for exact and right-censored rows.
#' @keywords internal
ic_event_data <- function(ic, t_imp) {
  kind <- ic[["ic_kind"]]
  ped_time <- ifelse(kind %in% c("exact", "right"), ic[["ic_L"]], t_imp)
  ped_status <- ifelse(kind == "right", 0L, 1L)

  out <- ic
  out[["ic_L"]] <- NULL
  out[["ic_R"]] <- NULL
  out[["ic_kind"]] <- NULL
  out[[".ped_time"]] <- ped_time
  out[[".ped_status"]] <- as.integer(ped_status)

  out
}

#' @describeIn interval_censored Drop subjects with non-positive follow-up time
#'   (e.g.\ right-censored at time 0 with no observed inspection), which carry no
#'   information and would break the interval split. Returns the filtered data.
#' @param evd A data frame with a \code{.ped_time} column (output of
#'   \code{ic_event_data}).
#' @param warn Logical; emit a one-time warning when rows are dropped.
#' @keywords internal
drop_zero_followup <- function(evd, warn = TRUE) {
  keep <- evd[[".ped_time"]] > 0
  if (any(!keep)) {
    if (warn) {
      warning(
        sum(!keep),
        " subject(s) with non-positive follow-up time ",
        "(no observed inspection) were dropped.",
        call. = FALSE
      )
    }
    evd <- evd[keep, , drop = FALSE]
  }
  evd
}

#' @describeIn interval_censored Transform interval-censored data into an
#'   initial (midpoint-imputed) PED object. Left- and interval-censored event
#'   times are initialised at the interval midpoint (\eqn{(L+R)/2}, and
#'   \eqn{R/2} for left-censored observations); this object is only an
#'   initialiser for \code{\link{pamm_ic}} and should not be used for inference
#'   on its own. The parsed interval bounds and the resolved cut-points are
#'   attached as the \code{"ic"} and \code{"breaks"} attributes.
#' @inheritParams as_ped
#' @keywords internal
as_ped_ic <- function(
  data,
  formula,
  cut = NULL,
  max_time = NULL,
  id = "id",
  ...
) {
  rhs_vars <- resolve_rhs_vars(formula, data = data, exclude = id)
  ic <- parse_ic_surv(formula, data, id = id)
  cut <- resolve_ic_cut(ic, cut = cut, max_time = max_time)

  # midpoint initialiser: left-censored -> R/2, interval -> (L+R)/2
  L <- ic[["ic_L"]]
  R <- ic[["ic_R"]]
  t_mid <- ifelse(ic[["ic_kind"]] == "left", R / 2, (L + R) / 2)
  t_mid <- pmin(t_mid, max(cut))

  ev_data <- drop_zero_followup(ic_event_data(ic, t_mid))
  ped_form <- stats::as.formula(
    paste0(
      "Surv(.ped_time, .ped_status) ~ ",
      paste0(unique(c(rhs_vars, id)), collapse = " + ")
    )
  )

  ped <- split_data(ped_form, data = ev_data, cut = cut, id = id, ...)

  attr(ped, "ic") <- ic
  attr(ped, "ic_formula") <- formula
  attr(ped, "breaks") <- cut
  class(ped) <- unique(c("ped_ic_init", class(ped)))

  ped
}
