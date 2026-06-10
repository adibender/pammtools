# Fit a PAMM and compute S(t|x) confidence intervals with all methods under
# comparison: "default", "delta", and "sim" (with two values of nsim).

fit_pamm <- function(dat, cuts, k) {
  ped <- pammtools::as_ped(
    survival::Surv(time, status) ~ x1, data = dat, cut = cuts
  )
  mod <- pammtools::pamm(
    ped_status ~ s(tend, k = k) + x1, data = ped, method = "REML"
  )
  list(ped = ped, mod = mod)
}

# Returns one row per (method, profile, evaluation time).
# Note: newdata passed to add_surv_prob() must contain *all* interval
# endpoints (add_surv_prob only reconstructs intlen from the rows it is
# given), and must be grouped by the covariate profile; results are filtered
# to the evaluation times afterwards.
compute_cis <- function(fit, eval_tbl, se_mult = qnorm(0.975),
                        nsims = c(100L, 500L)) {
  nd <- pammtools::make_newdata(fit$ped, tend = unique(tend), x1 = PROFILES)
  nd <- dplyr::arrange(dplyr::group_by(nd, x1), tend, .by_group = TRUE)

  res <- list(
    default = pammtools::add_surv_prob(
      nd, fit$mod, ci = TRUE, ci_type = "default", se_mult = se_mult),
    delta = pammtools::add_surv_prob(
      nd, fit$mod, ci = TRUE, ci_type = "delta", se_mult = se_mult)
  )
  for (ns in nsims) {
    res[[paste0("sim", ns)]] <- pammtools::add_surv_prob(
      nd, fit$mod, ci = TRUE, ci_type = "sim", nsim = as.integer(ns))
  }

  res <- lapply(res, function(d) {
    d <- dplyr::ungroup(d)
    d[, c("x1", "tend", "surv_prob", "surv_lower", "surv_upper")]
  })
  res <- dplyr::bind_rows(res, .id = "method")
  stopifnot(all(res$surv_lower <= res$surv_upper + 1e-12, na.rm = TRUE))

  # evaluation times are fitting-grid endpoints, i.e. exactly the tend values
  # in nd; rounding guards against floating-point noise in the join
  res$.tkey <- round(res$tend, 8)
  et <- eval_tbl[!is.na(eval_tbl$tend), ]
  et$.tkey <- round(et$tend, 8)
  out <- merge(
    et[, c("x1", "t_label", ".tkey", "true_surv")],
    res[, setdiff(names(res), "tend")],
    by = c("x1", ".tkey")
  )
  out$tend <- out$.tkey
  out$.tkey <- NULL
  out
}
