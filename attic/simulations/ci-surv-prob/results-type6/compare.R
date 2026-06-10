# Paired comparison of sim-CI coverage: type-7 (stored full run, reps 1-100)
# vs type-6 (rerun with patched pammtools, issue #288) on identical RNG
# streams. delta/default are unaffected by the patch and serve as a pairing
# check (they must agree if data + fits reproduced exactly).

library(dplyr)

scens <- c(1L, 30L, 33L)
reps_used <- 100L

read_run <- function(dir, sid, max_rep) {
  res <- readRDS(file.path(dir, sprintf("scen-%02d.rds", sid)))
  bind_rows(res[seq_len(min(max_rep, length(res)))])
}

old <- bind_rows(lapply(scens, read_run, dir = "results/raw", max_rep = reps_used))
new <- bind_rows(lapply(scens, read_run, dir = "results-type6/raw", max_rep = reps_used))

key <- c("scenario_id", "rep_id", "x1", "t_label", "method")
both <- inner_join(
  old, new, by = key, suffix = c("_t7", "_t6"), relationship = "one-to-one"
)

# ---- pairing checks -------------------------------------------------------
cat("== Pairing checks (reps 1-", reps_used, ", scenarios ",
    paste(scens, collapse = ", "), ") ==\n", sep = "")
cat("rows matched:", nrow(both), "of", nrow(old), "(old) /", nrow(new), "(new)\n")
pair <- both %>%
  filter(!fit_failed_t7, !fit_failed_t6, !is.na(surv_prob_t7)) %>%
  group_by(method) %>%
  summarize(
    max_dprob  = max(abs(surv_prob_t7 - surv_prob_t6)),
    max_dlower = max(abs(surv_lower_t7 - surv_lower_t6)),
    max_dupper = max(abs(surv_upper_t7 - surv_upper_t6)),
    .groups = "drop"
  )
print(as.data.frame(pair), digits = 3)
cat("(delta/default rows must be ~0 throughout; sim rows must differ only",
    "in the bounds)\n\n")

# ---- coverage / width, sim methods ----------------------------------------
# deduplicate collapsed labels (several S-labels can snap to the same
# endpoint on coarse grids), as in eval-sim.R: keep one row per realized
# endpoint
sim <- both %>%
  filter(method %in% c("sim100", "sim500"), !fit_failed_t7,
         !is.na(surv_prob_t7), !is.na(surv_prob_t6)) %>%
  distinct(scenario_id, rep_id, x1, method, tend_t7, .keep_all = TRUE) %>%
  mutate(
    cov_t7   = surv_lower_t7 <= true_surv_t7 & true_surv_t7 <= surv_upper_t7,
    cov_t6   = surv_lower_t6 <= true_surv_t6 & true_surv_t6 <= surv_upper_t6,
    width_t7 = surv_upper_t7 - surv_lower_t7,
    width_t6 = surv_upper_t6 - surv_lower_t6
  )

cat("== Coverage (nominal 0.95) and width, pooled over",
    "scenario x profile x time cells ==\n")
summ <- sim %>%
  group_by(method) %>%
  summarize(
    n_pairs     = n(),
    coverage_t7 = mean(cov_t7),
    coverage_t6 = mean(cov_t6),
    gain_pp     = 100 * (mean(cov_t6) - mean(cov_t7)),
    win_t6      = sum(cov_t6 & !cov_t7),   # type-6 covers, type-7 misses
    loss_t6     = sum(!cov_t6 & cov_t7),   # reverse (0 if draws identical)
    width_ratio = median(width_t6 / width_t7),
    .groups = "drop"
  )
print(as.data.frame(summ), digits = 4)

cat("\n== Per-scenario mean cell coverage ==\n")
per_scen <- sim %>%
  group_by(scenario_id, method, x1, t_label) %>%
  summarize(c7 = mean(cov_t7), c6 = mean(cov_t6), .groups = "drop") %>%
  group_by(scenario_id, method) %>%
  summarize(
    mean_cell_cov_t7 = mean(c7),
    mean_cell_cov_t6 = mean(c6),
    cells_below_93_t7 = mean(c7 < 0.93),
    cells_below_93_t6 = mean(c6 < 0.93),
    .groups = "drop"
  )
print(as.data.frame(per_scen), digits = 3)
