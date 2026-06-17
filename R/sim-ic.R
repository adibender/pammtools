#' Turn exact event times into interval-censored observations
#'
#' Convenience helper to manufacture interval-censored (panel) data from exact
#' simulated survival times (e.g.\ the output of \code{\link{sim_pexp}}), for
#' coverage studies and examples. Each subject is "inspected" at a sequence of
#' times; the true event time is then only known to lie between the last clean
#' and the first positive inspection. The exact time is retained (by default in
#' column \code{true_time}) so that coverage can be scored against the truth.
#'
#' @param data A data frame with one row per subject containing an exact event
#'   time and a status indicator (as produced by \code{\link{sim_pexp}}).
#' @param time_var,status_var Names of the (exact) event-time and status columns.
#'   \code{status_var} may be missing, in which case all rows are treated as
#'   events.
#' @param mechanism Inspection mechanism: \code{"random"} (default) draws
#'   inter-inspection gaps from an \code{Exp(rate)} distribution; \code{"fixed"}
#'   uses the common grid given in \code{schedule}; \code{"mixed"} jitters the
#'   fixed grid by a random offset per subject.
#' @param rate Inspection rate for \code{mechanism = "random"} / \code{"mixed"}
#'   (expected gap \eqn{1/\mathrm{rate}}).
#' @param schedule Numeric vector of inspection times for
#'   \code{mechanism = "fixed"}/\code{"mixed"}.
#' @param max_time Inspection horizon. Defaults to \code{max(data[[time_var]])}.
#' @param terminal_exam Logical; if \code{TRUE} (default), every subject is
#'   additionally examined at \code{max_time} (an end-of-study examination), so
#'   events before \code{max_time} always have a finite upper bound and only
#'   subjects event-free at \code{max_time} are right-censored. If \code{FALSE},
#'   there is no closing examination: events after a subject's last inspection
#'   are right-censored at that inspection, and \emph{subjects that exit
#'   event-free} (\code{status == 0}) are likewise right-censored at their last
#'   inspection before exit (not at their exact exit time). Both conventions
#'   yield coarsening-at-random data; mixing them (exact exit times for
#'   survivors but open intervals for undetected events) would make the
#'   right-censoring informative and bias every interval-censoring likelihood.
#' @param keep_truth Logical; keep the exact event time in \code{true_time}.
#' @param L,R Names of the created lower/upper bound columns.
#' @return \code{data} augmented with interval bounds in columns \code{L} and
#'   \code{R} (and \code{true_time}). Use
#'   \code{Surv(L, R, type = "interval2")} on the result.
#' @seealso \code{\link{pamm_ic}}, \code{\link{sim_pexp}}
#' @importFrom stats rexp runif
#' @examples
#' \donttest{
#' set.seed(1)
#' df <- data.frame(x = runif(100, -1, 1))
#' sdf <- sim_pexp(~ -2 + 0.4 * x, df, cut = seq(0, 10, by = 0.5))
#' icd <- add_inspections(sdf, rate = 1)
#' fit <- pamm_ic(Surv(L, R, type = "interval2") ~ x, icd, m = 5)
#' }
#' @export
add_inspections <- function(
  data,
  time_var = "time",
  status_var = "status",
  mechanism = c("random", "fixed", "mixed"),
  rate = 1,
  schedule = NULL,
  max_time = NULL,
  terminal_exam = TRUE,
  keep_truth = TRUE,
  L = "L",
  R = "R"
) {
  mechanism <- match.arg(mechanism)
  assert_data_frame(data, min.rows = 1)
  assert_choice(time_var, names(data))
  assert_flag(terminal_exam)

  tt <- data[[time_var]]
  st <- if (status_var %in% names(data)) data[[status_var]] else
    rep(1L, nrow(data))
  # NA status would silently be treated as an event below
  assert_integerish(st, any.missing = FALSE)
  hor <- if (is.null(max_time)) max(tt) else max_time

  if (mechanism != "random") {
    assert_numeric(schedule, min.len = 1, any.missing = FALSE)
  }

  gen_grid <- function() {
    g <- if (mechanism == "random") {
      g <- numeric(0)
      cur <- 0
      repeat {
        cur <- cur + rexp(1, rate)
        if (cur > hor) break
        g <- c(g, cur)
      }
      g
    } else if (mechanism == "fixed") {
      schedule[schedule <= hor]
    } else {
      off <- runif(1, 0, 1 / rate)
      s <- (schedule + off)
      s[s <= hor]
    }
    # the closing examination is deterministic, so adding it never changes
    # the RNG stream relative to terminal_exam = FALSE
    if (terminal_exam) g <- c(g, hor)
    sort(unique(g))
  }

  Lv <- numeric(nrow(data))
  Rv <- numeric(nrow(data))
  for (i in seq_len(nrow(data))) {
    if (isTRUE(st[i] == 0)) {
      if (terminal_exam) {
        # subject exits event-free at tt[i]: event known only to be > tt[i]
        Lv[i] <- tt[i]
        Rv[i] <- Inf
      } else {
        # no closing exam: last seen event-free at the final inspection
        # before exit (using the exact exit time here while undetected
        # events keep open intervals would be informative censoring)
        grid <- gen_grid()
        below <- grid[grid <= tt[i]]
        Lv[i] <- if (length(below)) max(below) else 0
        Rv[i] <- Inf
      }
      next
    }
    grid <- gen_grid()
    below <- grid[grid < tt[i]]
    above <- grid[grid >= tt[i]]
    Lv[i] <- if (length(below)) max(below) else 0
    Rv[i] <- if (length(above)) min(above) else Inf
  }

  out <- data
  out[[L]] <- Lv
  out[[R]] <- Rv
  if (keep_truth) {
    out[["true_time"]] <- tt
  }

  out
}
