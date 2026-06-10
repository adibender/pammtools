#!/usr/bin/env Rscript

# Oracle/sanity checks for the CI coverage study. Run before any pilot:
#   Rscript sanity-checks.R          # all checks (incl. micro coverage, ~mins)
#   Rscript sanity-checks.R --quick  # skip the micro coverage check

suppressPackageStartupMessages({
  library(pammtools)
  library(dplyr)
})

script_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
for (f in list.files(file.path(script_dir, "R"), full.names = TRUE)) source(f)

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
ok <- function(label) message("PASS: ", label)

grid <- scenario_grid()

## 1. truth: constant baseline has closed form S(t|x) = exp(-exp(b0+beta*x)*t)
b0c <- solve_b0("constant")
tt  <- c(0.3, 1, 2.5, 7.7, 10)
for (x in PROFILES) {
  stopifnot(isTRUE(all.equal(
    true_surv(tt, "constant", b0c, x),
    exp(-exp(b0c + BETA_X1 * x) * tt),
    tolerance = 1e-12
  )))
}
stopifnot(isTRUE(all.equal(true_surv(TMAX, "constant", b0c, 0), S_TMAX_REF)))
ok("true_surv matches closed form (constant baseline), S(TMAX|0) calibrated")

## 1b. b0 calibration holds for all baselines
for (b in names(baseline_funs)) {
  stopifnot(isTRUE(all.equal(true_surv(TMAX, b, solve_b0(b), 0), S_TMAX_REF)))
}
ok("b0 calibration for all baselines")

## 1c. censoring calibration: simulate and compare empirical censoring rate
set.seed(1)
for (sid in c(1, 36)) {
  scen <- grid[grid$scenario_id == sid, ]
  dat  <- simulate_one_dataset(scen, n_override = 50000L)
  emp  <- mean(dat$status == 0)
  stopifnot(abs(emp - scen$cens) < 0.02)
}
ok("calibrated censoring rates match empirical rates (within 2pp)")

## 2. estimand alignment: huge n, no censoring => surv_prob close to truth
scen_big <- grid[grid$baseline == "nonmonotone" & grid$n == 1000 &
                   grid$cens == 0.2 & grid$grid_type == "quantile", ][1, ]
scen_big$cens_rate <- 0
set.seed(2)
dat_big  <- simulate_one_dataset(scen_big, n_override = 20000L)
prep_big <- prep_scenario(scen_big)
fit_big  <- fit_pamm(dat_big, prep_big$cuts, k = 10)
cis_big  <- compute_cis(fit_big, prep_big$eval_tbl, nsims = 5000L)
dev <- with(cis_big, abs(surv_prob - true_surv))
stopifnot(max(dev) < 0.01)
ok(sprintf("large-n estimand alignment (max |S_hat - S_true| = %.4f)", max(dev)))

## 3a. large n: delta and sim(5000) intervals agree closely
cmp <- merge(
  subset(cis_big, method == "delta",
         c("x1", "t_label", "surv_lower", "surv_upper")),
  subset(cis_big, method == "sim5000",
         c("x1", "t_label", "surv_lower", "surv_upper")),
  by = c("x1", "t_label"), suffixes = c("_delta", "_sim")
)
dl <- with(cmp, max(abs(surv_lower_delta - surv_lower_sim),
                    abs(surv_upper_delta - surv_upper_sim)))
stopifnot(dl < 0.005)
ok(sprintf("delta ~ sim(5000) at large n (max bound diff = %.4f)", dl))

## 3b. default ~ delta at the *first* interval endpoint (single cumsum term)
nd1 <- make_newdata(fit_big$ped, tend = unique(tend), x1 = c(0)) %>%
  group_by(x1)
d_def <- add_surv_prob(nd1, fit_big$mod, ci_type = "default",
                       se_mult = qnorm(0.975)) %>% ungroup() %>% slice(1)
d_del <- add_surv_prob(nd1, fit_big$mod, ci_type = "delta",
                       se_mult = qnorm(0.975)) %>% ungroup() %>% slice(1)
stopifnot(abs(d_def$surv_lower - d_del$surv_lower) < 5e-4,
          abs(d_def$surv_upper - d_del$surv_upper) < 5e-4)
ok("default ~ delta at first interval (first-order agreement)")

## 4. mechanics on a small fit
scen_sm <- grid[grid$scenario_id == 1, ]
prep_sm <- prep_scenario(scen_sm)
set.seed(3)
dat_sm <- simulate_one_dataset(scen_sm, n_override = 300L)
fit_sm <- fit_pamm(dat_sm, prep_sm$cuts, k = 10)
nd <- make_newdata(fit_sm$ped, tend = unique(tend), x1 = PROFILES) %>%
  group_by(x1) %>% arrange(tend, .by_group = TRUE)
for (ct in c("default", "delta", "sim")) {
  r <- add_surv_prob(nd, fit_sm$mod, ci_type = ct, se_mult = qnorm(0.975)) %>%
    ungroup()
  stopifnot(all(r$surv_lower <= r$surv_prob + 1e-12),
            all(r$surv_prob <= r$surv_upper + 1e-12))
  if (ct != "delta") {
    stopifnot(all(r$surv_lower >= -1e-12), all(r$surv_upper <= 1 + 1e-12))
  }
  mono <- r %>% group_by(x1) %>%
    summarize(mono = all(diff(surv_prob) <= 1e-12), .groups = "drop")
  stopifnot(all(mono$mono))
}
ok("bounds ordered; default/sim within [0,1]; S(t) non-increasing")

## 4b. grouped two-profile call equals two single-profile calls (delta)
r2 <- add_surv_prob(nd, fit_sm$mod, ci_type = "delta",
                    se_mult = qnorm(0.975)) %>% ungroup()
for (x in PROFILES) {
  nd1 <- make_newdata(fit_sm$ped, tend = unique(tend), x1 = c(x)) %>%
    group_by(x1)
  r1 <- add_surv_prob(nd1, fit_sm$mod, ci_type = "delta",
                      se_mult = qnorm(0.975)) %>% ungroup()
  rx <- r2 %>% filter(x1 == x)
  stopifnot(isTRUE(all.equal(r1$surv_lower, rx$surv_lower, tolerance = 1e-10)),
            isTRUE(all.equal(r1$surv_upper, rx$surv_upper, tolerance = 1e-10)))
}
ok("grouped two-profile call equals single-profile calls")

## 5. reproducibility of run_one_rep given identical RNG stream
RNGkind("L'Ecuyer-CMRG")
set.seed(285)
stream1 <- parallel::nextRNGStream(.Random.seed)
stream2 <- parallel::nextRNGStream(stream1)
ra <- run_one_rep(scen_sm, prep_sm, 1L, stream1)
rb <- run_one_rep(scen_sm, prep_sm, 1L, stream1)
rc <- run_one_rep(scen_sm, prep_sm, 2L, stream2)
stopifnot(identical(ra[setdiff(names(ra), "time_sec")],
                    rb[setdiff(names(rb), "time_sec")]))
stopifnot(!isTRUE(all.equal(ra$surv_prob, rc$surv_prob)))
ok("run_one_rep reproducible per stream; different streams differ")

## 6. micro coverage check: easiest scenario, delta at S50, 200 reps
if (!quick) {
  scen_easy <- grid[grid$baseline == "constant" & grid$n == 1000 &
                      grid$cens == 0.2 & grid$grid_type == "equidistant", ]
  prep_easy <- prep_scenario(scen_easy)
  set.seed(285)
  base_stream <- .Random.seed
  covered <- vapply(seq_len(200), function(r) {
    base_stream <<- parallel::nextRNGStream(base_stream)
    res <- run_one_rep(scen_easy, prep_easy, r, base_stream)
    row <- subset(res, method == "delta" & t_label == "S50" & x1 == 0)
    row$surv_lower <= row$true_surv && row$true_surv <= row$surv_upper
  }, logical(1))
  cov <- mean(covered)
  stopifnot(cov >= 0.90, cov <= 0.99)
  ok(sprintf("micro coverage check: delta @ S50 = %.3f (200 reps)", cov))
} else {
  message("SKIP: micro coverage check (--quick)")
}

message("All sanity checks passed.")
