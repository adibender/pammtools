# Data-generating process and scenario definitions for the CI coverage study
# (pammtools issue #285). Sourced by run-sim.R / sanity-checks.R.

# ------------------------------------------------------------------ constants
TMAX        <- 10
DGP_CUT     <- seq(0, TMAX, length.out = 101) # DGP grid: 100 intervals
BETA_X1     <- 0.7
PROFILES    <- c(0, 1.5)                      # covariate profiles for S(t|x)
S_QUANTILES <- c(S90 = 0.90, S75 = 0.75, S50 = 0.50, S25 = 0.25, S10 = 0.10)
S_TMAX_REF  <- 0.07                           # target true S(TMAX | x1 = 0)

# baseline log-hazard shapes f0(t)
baseline_funs <- list(
  constant    = function(t) 0 * t,
  increasing  = function(t) 0.25 * t,            # Gompertz-type, HR exp(2.5) over follow-up
  nonmonotone = function(t) dgamma(t, 8, 2) * 6  # as in vignettes/simulations.Rmd
)

# ------------------------------------------------------------------ truth
# sim_pexp() evaluates the log-hazard formula at interval *starts* of `cut`
# (R/sim-pexp.R: rename("t" = "tstart")), so the DGP hazard is the step
# function with rate exp(b0 + f0(tstart_l) + beta * x1) on (tstart_l, tend_l].
# All "truth" below replicates exactly this convention.
dgp_rates <- function(baseline, b0, x1) {
  f0 <- baseline_funs[[baseline]]
  exp(b0 + f0(DGP_CUT[-length(DGP_CUT)]) + BETA_X1 * x1)
}

# exact true S(t | x1) at arbitrary times under the step-hazard DGP
true_surv <- function(times, baseline, b0, x1) {
  rates  <- dgp_rates(baseline, b0, x1)
  starts <- DGP_CUT[-length(DGP_CUT)]
  intlen <- diff(DGP_CUT)
  cumhaz_ends <- cumsum(rates * intlen)
  vapply(times, function(t) {
    if (t <= 0) return(1)
    j  <- findInterval(t, DGP_CUT, left.open = TRUE, rightmost.closed = TRUE)
    ch <- if (j > 1) cumhaz_ends[j - 1] else 0
    exp(-(ch + rates[j] * (t - starts[j])))
  }, numeric(1))
}

# intercept b0 such that true S(TMAX | x1 = 0) equals S_TMAX_REF
solve_b0 <- function(baseline) {
  f0     <- baseline_funs[[baseline]]
  starts <- DGP_CUT[-length(DGP_CUT)]
  log(-log(S_TMAX_REF) / sum(exp(f0(starts)) * diff(DGP_CUT)))
}

# ------------------------------------------------------------------ censoring
# Censoring: C ~ Exp(c_rate), independent of T | x1, plus administrative
# censoring at TMAX. P(event observed | x1) has a closed form under the
# step-hazard DGP; the marginal over x1 ~ N(0,1) is computed on a quantile grid.
prob_event_observed <- function(c_rate, baseline, b0,
                                x_nodes = qnorm(ppoints(401))) {
  starts <- DGP_CUT[-length(DGP_CUT)]
  intlen <- diff(DGP_CUT)
  mean(vapply(x_nodes, function(x) {
    lam          <- dgp_rates(baseline, b0, x)
    cumhaz_start <- c(0, cumsum(lam * intlen))[seq_along(lam)]
    sum(lam / (lam + c_rate) * exp(-cumhaz_start - c_rate * starts) *
          (1 - exp(-(lam + c_rate) * intlen)))
  }, numeric(1)))
}

# Exp-rate of the censoring distribution s.t. overall P(censored) = target
# (administrative censoring at TMAX included in the target).
calibrate_censoring <- function(baseline, b0, target_cens) {
  f <- function(log_c) {
    (1 - prob_event_observed(exp(log_c), baseline, b0)) - target_cens
  }
  if (f(log(1e-8)) >= 0) return(0) # admin censoring alone reaches target
  exp(uniroot(f, lower = log(1e-8), upper = log(50), tol = 1e-10)$root)
}

# ------------------------------------------------------------------ scenarios
scenario_grid <- function() {
  grid <- expand.grid(
    baseline  = names(baseline_funs),
    n         = c(100L, 250L, 1000L),
    cens      = c(0.2, 0.6),
    grid_type = c("equidistant", "quantile"),
    stringsAsFactors = FALSE
  )
  grid$J <- ifelse(grid$grid_type == "equidistant", 10L, 40L)
  # k must stay safely below the number of *populated* intervals: under heavy
  # censoring the trailing equidistant intervals are often empty, and mgcv
  # errors when s(tend, k) has fewer unique tend values than k
  grid$k <- ifelse(grid$grid_type == "equidistant", 7L, 10L)
  # sensitivity arm: one mid scenario refit with larger basis dimension
  sens <- grid[grid$baseline == "nonmonotone" & grid$n == 250 &
                 grid$cens == 0.2 & grid$grid_type == "quantile", ]
  sens$k <- 20L
  grid <- rbind(grid, sens)
  grid$scenario_id <- seq_len(nrow(grid))
  b0s <- vapply(names(baseline_funs), solve_b0, numeric(1))
  grid$b0 <- b0s[grid$baseline]
  # censoring rates only depend on (baseline, cens)
  cc <- unique(grid[, c("baseline", "cens", "b0")])
  cc$cens_rate <- mapply(calibrate_censoring, cc$baseline, cc$b0, cc$cens)
  grid <- merge(grid, cc, by = c("baseline", "cens", "b0"), sort = FALSE)
  grid[order(grid$scenario_id), ]
}

# ------------------------------------------------------------------ fit grid
# Fitting cut points: either equidistant or event-time-quantile-based, the
# latter computed once per scenario from a fixed, seeded reference simulation
# (n = 20000) so that the grid -- and hence the evaluation times -- are
# constant across replications within a scenario.
make_fit_cuts <- function(scen, prep_seed = 20260610) {
  if (scen$grid_type == "equidistant") {
    return(seq(0, TMAX, length.out = scen$J + 1))
  }
  rng_old <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  on.exit(if (!is.null(rng_old)) {
    assign(".Random.seed", rng_old, envir = globalenv())
  })
  set.seed(prep_seed + scen$scenario_id)
  dat <- simulate_one_dataset(scen, n_override = 20000L)
  ev  <- dat$time[dat$status == 1L]
  cuts <- c(0, quantile(ev, probs = seq_len(scen$J - 1) / scen$J), TMAX)
  unname(unique(cuts))
}

# Evaluation times: per profile, the fitting-grid endpoints closest to the
# times where true S(t | x1) hits S_QUANTILES; NA when not reached by TMAX.
make_eval_times <- function(scen, cuts) {
  endpoints <- cuts[-1]
  do.call(rbind, lapply(PROFILES, function(x) {
    St <- true_surv(endpoints, scen$baseline, scen$b0, x)
    do.call(rbind, lapply(names(S_QUANTILES), function(lab) {
      q <- S_QUANTILES[[lab]]
      if (min(St) > q) {
        return(data.frame(x1 = x, t_label = lab, tend = NA_real_,
                          true_surv = NA_real_))
      }
      i <- which.min(abs(St - q))
      data.frame(x1 = x, t_label = lab, tend = endpoints[i], true_surv = St[i])
    }))
  }))
}

# all per-scenario preparation that is fixed across replications
prep_scenario <- function(scen) {
  cuts <- make_fit_cuts(scen)
  list(cuts = cuts, eval_tbl = make_eval_times(scen, cuts))
}
