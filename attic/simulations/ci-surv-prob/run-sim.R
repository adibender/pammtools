#!/usr/bin/env Rscript

# Run the CI coverage simulation study (pammtools issue #285).
#
# Usage:
#   Rscript run-sim.R [--scenarios 1-37] [--reps 25] [--cores 4]
#                     [--seed 285] [--out results/raw]
#
# Results are checkpointed per scenario (results/raw/scen-XX.rds); re-running
# with a higher --reps computes only the missing replications. Reproducible:
# every (scenario, rep) gets its own pre-generated L'Ecuyer-CMRG stream
# derived from --seed, independent of --cores and scheduling.

suppressPackageStartupMessages({
  library(pammtools)
  library(dplyr)
})

# ------------------------------------------------------------------ CLI args
parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) stop("Arguments must come in --key value pairs.")
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(defaults)) stop("Unknown argument: --", key)
    defaults[[key]] <- args[i + 1]
  }
  defaults
}

parse_scenarios <- function(x, n_scen) {
  out <- unlist(lapply(strsplit(x, ",")[[1]], function(part) {
    rng <- as.integer(strsplit(part, "-")[[1]])
    if (length(rng) == 2) seq(rng[1], rng[2]) else rng
  }))
  out <- sort(unique(out))
  stopifnot(all(out >= 1 & out <= n_scen))
  out
}

script_dir <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})

opts <- parse_args(list(
  scenarios = "all", reps = "25", cores = "4", seed = "285",
  out = file.path(script_dir, "results", "raw")
))
n_reps  <- as.integer(opts$reps)
n_cores <- as.integer(opts$cores)
seed    <- as.integer(opts$seed)
out_dir <- opts$out

for (f in list.files(file.path(script_dir, "R"), full.names = TRUE)) source(f)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------- scenario preparation
# deterministic (censoring calibration is analytic, quantile cuts use a fixed
# internal seed), cached because calibration takes a few seconds per scenario
grid <- scenario_grid()
scen_ids <- if (opts$scenarios == "all") {
  grid$scenario_id
} else {
  parse_scenarios(opts$scenarios, nrow(grid))
}

prep_file <- file.path(dirname(out_dir), "prep.rds")
prep <- if (file.exists(prep_file)) readRDS(prep_file) else list()
for (s in scen_ids) {
  key <- as.character(s)
  if (is.null(prep[[key]])) {
    prep[[key]] <- prep_scenario(grid[grid$scenario_id == s, ])
  }
}
dir.create(dirname(prep_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(prep, prep_file)

# ------------------------------------------------------------- RNG streams
# stream index of (scenario s, rep r) is (s-1) * RESERVED_REPS + r, so results
# are extendable up to RESERVED_REPS replications without changing streams
RESERVED_REPS <- 1000L
stopifnot(n_reps <= RESERVED_REPS)
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
n_streams <- max(scen_ids) * RESERVED_REPS
streams <- vector("list", n_streams)
s <- .Random.seed
for (m in seq_len(n_streams)) {
  s <- parallel::nextRNGStream(s)
  streams[[m]] <- s
}

# ------------------------------------------------------------------- run
for (sid in scen_ids) {
  scen      <- grid[grid$scenario_id == sid, ]
  scen_prep <- prep[[as.character(sid)]]
  res_file  <- file.path(out_dir, sprintf("scen-%02d.rds", sid))

  existing  <- if (file.exists(res_file)) readRDS(res_file) else list()
  done_reps <- as.integer(names(existing))
  todo      <- setdiff(seq_len(n_reps), done_reps)
  if (!length(todo)) {
    message(sprintf("scenario %02d: %d reps already done, skipping.",
                    sid, length(done_reps)))
    next
  }

  t0 <- proc.time()[3]
  new_res <- parallel::mclapply(todo, function(r) {
    tryCatch(
      run_one_rep(scen, scen_prep, r,
                  streams[[(sid - 1L) * RESERVED_REPS + r]]),
      error = function(e) e
    )
  }, mc.cores = n_cores, mc.preschedule = TRUE)
  names(new_res) <- as.character(todo)

  errs <- vapply(new_res, inherits, logical(1), what = "condition")
  if (any(errs)) {
    warning(sprintf("scenario %02d: %d replication(s) errored: %s",
                    sid, sum(errs),
                    conditionMessage(new_res[[which(errs)[1]]])))
    new_res <- new_res[!errs]
  }

  existing <- c(existing, new_res)
  existing <- existing[order(as.integer(names(existing)))]
  saveRDS(existing, res_file)
  message(sprintf(
    "scenario %02d (%s, n=%d, cens=%.0f%%, %s J=%d, k=%d): +%d reps in %.1fs",
    sid, scen$baseline, scen$n, 100 * scen$cens, scen$grid_type, scen$J,
    scen$k, sum(!errs), proc.time()[3] - t0
  ))
}

message("Done. Results in ", out_dir)
