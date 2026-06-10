#!/usr/bin/env Rscript

# Aggregate simulation results: coverage, widths, width ratios, failure rates.
#   Rscript eval-sim.R [--raw results/raw] [--out output] [--tag pilot]
# Writes CSV summaries and figures to --out.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

script_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
for (f in list.files(file.path(script_dir, "R"), full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(key, default) {
  i <- which(args == paste0("--", key))
  if (length(i)) args[i + 1] else default
}
raw_dir <- get_opt("raw", file.path(script_dir, "results", "raw"))
out_dir <- get_opt("out", file.path(script_dir, "output"))
tag     <- get_opt("tag", "")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(raw_dir, pattern = "^scen-\\d+\\.rds$", full.names = TRUE)
if (!length(files)) stop("No result files in ", raw_dir)
res <- bind_rows(lapply(files, function(f) bind_rows(readRDS(f))))

grid <- scenario_grid() %>%
  mutate(
    baseline = factor(baseline,
                      levels = c("constant", "increasing", "nonmonotone")),
    cens_lab = sprintf("%.0f%% cens", 100 * cens),
    grid_lab = sprintf("%s J=%d", grid_type, J),
    arm      = ifelse(k == 20L, "sensitivity (k=20)", "main")
  )
res <- left_join(res, grid, by = "scenario_id")

# ----------------------------------------------------------------- failures
fail_tbl <- res %>%
  group_by(scenario_id, baseline, n, cens, grid_type, J, k) %>%
  summarize(
    n_reps      = n_distinct(rep_id),
    fail_rate   = mean(tapply(fit_failed, rep_id, any)),
    nonconv     = mean(tapply(!converged & !fit_failed, rep_id, any),
                       na.rm = TRUE),
    mean_events = mean(tapply(n_events, rep_id, max)),
    .groups = "drop"
  )
write.csv(fail_tbl, file.path(out_dir, paste0("failures", tag, ".csv")),
          row.names = FALSE)

ok <- res %>% filter(!fit_failed, !is.na(true_surv))
ok$method <- factor(ok$method,
                    levels = c("default", "delta", "sim100", "sim500"))
ok$t_label <- factor(ok$t_label, levels = names(S_QUANTILES))

wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n
  c2 <- z^2 / n
  lo <- (p + c2 / 2 - z * sqrt(p * (1 - p) / n + c2 / (4 * n))) / (1 + c2)
  hi <- (p + c2 / 2 + z * sqrt(p * (1 - p) / n + c2 / (4 * n))) / (1 + c2)
  cbind(lo, hi)
}

# ----------------------------------------------------------------- coverage
cov_tbl <- ok %>%
  mutate(
    covered    = surv_lower <= true_surv & true_surv <= surv_upper,
    miss_below = true_surv < surv_lower,
    miss_above = true_surv > surv_upper,
    width      = surv_upper - surv_lower,
    violates   = surv_lower < 0 | surv_upper > 1,
    cov_clip   = pmax(surv_lower, 0) <= true_surv &
                   true_surv <= pmin(surv_upper, 1)
  ) %>%
  group_by(scenario_id, baseline, n, cens, cens_lab, grid_type, grid_lab,
           J, k, arm, method, x1, t_label) %>%
  summarize(
    n_reps     = n(),
    coverage   = mean(covered),
    miss_below = mean(miss_below),
    miss_above = mean(miss_above),
    mean_width = mean(width),
    med_width  = median(width),
    bias       = mean(surv_prob - true_surv),
    rmse       = sqrt(mean((surv_prob - true_surv)^2)),
    viol_01    = mean(violates),
    cov_clip   = mean(cov_clip),
    true_surv  = first(true_surv),
    .groups = "drop"
  )
ci <- wilson(cov_tbl$coverage * cov_tbl$n_reps, cov_tbl$n_reps)
cov_tbl$cov_lo <- ci[, 1]
cov_tbl$cov_hi <- ci[, 2]
write.csv(cov_tbl, file.path(out_dir, paste0("coverage", tag, ".csv")),
          row.names = FALSE)

# --------------------------------------------------- paired width ratios
ratio_tbl <- ok %>%
  mutate(width = surv_upper - surv_lower) %>%
  select(scenario_id, rep_id, method, x1, t_label, width) %>%
  pivot_wider(names_from = method, values_from = width) %>%
  pivot_longer(cols = c("default", "sim100", "sim500"),
               names_to = "method", values_to = "width") %>%
  mutate(ratio = width / delta) %>%
  group_by(scenario_id, method, x1, t_label) %>%
  summarize(
    med_ratio = median(ratio, na.rm = TRUE),
    q25_ratio = quantile(ratio, 0.25, na.rm = TRUE),
    q75_ratio = quantile(ratio, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(distinct(select(grid, scenario_id, baseline, n, cens, cens_lab,
                            grid_type, grid_lab, J, k, arm)),
            by = "scenario_id")
write.csv(ratio_tbl, file.path(out_dir, paste0("width-ratio", tag, ".csv")),
          row.names = FALSE)

# ----------------------------------------------------------------- figures
main_cov <- cov_tbl %>% filter(arm == "main")
for (x_prof in unique(main_cov$x1)) {
  p <- main_cov %>%
    filter(x1 == x_prof) %>%
    ggplot(aes(x = t_label, y = coverage, color = method, group = method)) +
    geom_hline(yintercept = 0.95, linetype = 2, color = "grey40") +
    geom_line() +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = cov_lo, ymax = cov_hi), width = 0.15,
                  alpha = 0.5) +
    facet_grid(baseline + grid_lab ~ n + cens_lab,
               labeller = labeller(n = function(x) paste0("n=", x))) +
    labs(
      title = sprintf("Empirical coverage of 95%% CIs for S(t | x1 = %g)",
                      x_prof),
      subtitle = paste0("Error bars: Wilson 95% CI for the coverage ",
                        "estimate; x-axis: time points where true S(t) ",
                        "crosses 0.9, ..., 0.1"),
      x = NULL, y = "coverage"
    ) +
    theme_bw()
  ggsave(file.path(out_dir, sprintf("fig-coverage-x%s%s.png",
                                    gsub("[.]", "", x_prof), tag)),
         p, width = 14, height = 10, dpi = 150)
}

p_ratio <- ratio_tbl %>%
  filter(arm == "main") %>%
  ggplot(aes(x = t_label, y = med_ratio, color = method, group = method)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(baseline + grid_lab ~ n + cens_lab + x1,
             labeller = labeller(n = function(x) paste0("n=", x),
                                 x1 = function(x) paste0("x1=", x))) +
  labs(
    title = "Median CI width relative to the delta method (paired by rep)",
    x = NULL, y = "width ratio vs delta"
  ) +
  theme_bw()
ggsave(file.path(out_dir, paste0("fig-width-ratio", tag, ".png")),
       p_ratio, width = 16, height = 10, dpi = 150)

# ------------------------------------------------------------- text summary
overall <- cov_tbl %>%
  filter(arm == "main") %>%
  group_by(method) %>%
  summarize(
    mean_coverage = mean(coverage),
    min_coverage  = min(coverage),
    max_coverage  = max(coverage),
    mean_width    = mean(mean_width),
    viol_01       = mean(viol_01),
    .groups = "drop"
  )
by_factor <- lapply(
  c("baseline", "n", "cens", "grid_type", "t_label"),
  function(v) {
    cov_tbl %>%
      filter(arm == "main") %>%
      group_by(across(all_of(c(v, "method")))) %>%
      summarize(coverage = mean(coverage), mean_width = mean(mean_width),
                .groups = "drop")
  }
)
sink(file.path(out_dir, paste0("summary", tag, ".txt")))
cat("== Overall (averaged over scenarios, profiles, time points) ==\n")
print(as.data.frame(overall), digits = 3)
cat("\n== Marginal coverage by design factor ==\n")
for (tb in by_factor) {
  print(as.data.frame(tb), digits = 3)
  cat("\n")
}
sink()

message("Summaries and figures written to ", out_dir)
