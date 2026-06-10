# One replication: simulate -> fit -> CIs for all methods.
# `rng_stream` is a pre-generated L'Ecuyer-CMRG .Random.seed vector, making
# every (scenario, rep) individually reproducible regardless of scheduling.
# Fit/CI failures are recorded, not dropped.

run_one_rep <- function(scen, prep, rep_id, rng_stream) {
  assign(".Random.seed", rng_stream, envir = globalenv())
  t0 <- proc.time()[3]

  dat      <- simulate_one_dataset(scen)
  n_events <- sum(dat$status)

  failed_row <- function(stage) {
    data.frame(
      scenario_id = scen$scenario_id, rep_id = rep_id,
      method = NA_character_, x1 = NA_real_, t_label = NA_character_,
      tend = NA_real_, true_surv = NA_real_, surv_prob = NA_real_,
      surv_lower = NA_real_, surv_upper = NA_real_,
      fit_failed = TRUE, fail_stage = stage, converged = NA,
      n_events = n_events, time_sec = proc.time()[3] - t0
    )
  }

  fit <- tryCatch(fit_pamm(dat, prep$cuts, scen$k), error = identity)
  if (inherits(fit, "condition")) return(failed_row("fit"))

  cis <- tryCatch(compute_cis(fit, prep$eval_tbl), error = identity)
  if (inherits(cis, "condition")) return(failed_row("ci"))

  cbind(
    data.frame(scenario_id = scen$scenario_id, rep_id = rep_id),
    cis,
    data.frame(
      fit_failed = FALSE, fail_stage = NA_character_,
      converged = isTRUE(fit$mod$converged),
      n_events = n_events, time_sec = proc.time()[3] - t0
    )
  )
}
