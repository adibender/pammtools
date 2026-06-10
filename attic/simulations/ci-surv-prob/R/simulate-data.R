# Simulate one data set for a given scenario: piecewise-exponential event
# times via pammtools::sim_pexp() plus independent Exp censoring and
# administrative censoring at TMAX (the latter applied by sim_pexp itself).

simulate_one_dataset <- function(scen, n_override = NULL) {
  n  <- if (is.null(n_override)) scen$n else n_override
  df <- data.frame(x1 = rnorm(n))

  f0   <- baseline_funs[[scen$baseline]]
  b0   <- scen$b0
  beta <- BETA_X1
  form <- ~ b0 + f0(t) + beta * x1
  environment(form) <- environment() # sim_pexp() finds b0, f0, beta here

  sim_df <- pammtools::sim_pexp(form, df, cut = DGP_CUT)
  sim_df <- sim_df[order(sim_df$id), ]

  cens_time <- if (scen$cens_rate > 0) rexp(n, scen$cens_rate) else rep(Inf, n)
  data.frame(
    time   = pmin(sim_df$time, cens_time),
    status = as.integer(sim_df$status == 1L & sim_df$time <= cens_time),
    x1     = sim_df$x1
  )
}
