# Precompute the Bayesian baseline-hazard comparison for the `bayesian.Rmd`
# article.
#
# The JAGS and Stan fits run on an ~11.6k-row PED and take ~15-20 minutes
# combined, which is far too slow to evaluate on every website/vignette build.
# We therefore fit the models and run the pammtools prediction pipeline ONCE
# here, and save only the small resulting `pred_df` (a few kB). The article
# shows the full fitting + prediction code (with `eval = FALSE`) and renders the
# comparison plot and tables from this saved data frame.
#
# This file and `bayesian_results.rds` live in `vignettes/`, which is listed in
# `.Rbuildignore`, so neither enters the package tarball: they only support the
# pkgdown article and are kept under version control for reproducible rebuilds.
# (Saving the fitted `brmsfit` itself was rejected: it serialises to ~28 MB,
# unsuitable for git; the derived `pred_df` is a few kB.)
#
# Re-run after changing the data or model specifications:
#   Rscript vignettes/precompute_bayesian.R

suppressMessages({
  library(dplyr)
  library(mgcv)
  library(pammtools)
  library(survival)
  library(rjags)
  library(brms)
})

data("tumor")
tumor <- tumor %>%
  slice(1:200) %>%
  mutate(id = row_number())

ped <- as_ped(tumor, Surv(days, status) ~ 1, id = "id")
message("PED rows: ", nrow(ped))

## 1) frequentist PAMM (fast; included so the saved bundle is self-contained) ---
pam_mgcv <- mgcv::gam(
  formula = ped_status ~ s(tend, k = 20),
  data = ped,
  family = poisson(),
  offset = offset,
  method = "REML"
)

## 2) Bayesian PAMM via jagam + rjags ------------------------------------------
set.seed(101)
jagam_file <- tempfile(fileext = ".jags")
pam_jagam_prep <- mgcv::jagam(
  formula = ped_status ~ s(tend, k = 20),
  data = ped,
  family = poisson(),
  offset = offset,
  file = jagam_file,
  sp.prior = "gamma",
  diagonalize = TRUE
)
rjags::load.module("glm")
jagam_sampler <- rjags::jags.model(
  file = jagam_file,
  data = pam_jagam_prep$jags.data,
  inits = pam_jagam_prep$jags.ini,
  n.chains = 2,
  quiet = TRUE
)
jagam_post <- rjags::jags.samples(
  model = jagam_sampler,
  variable.names = c("b", "lambda"),
  n.iter = 4000,
  thin = 10
)
pam_jagam <- mgcv::sim2jam(jagam_post, pam_jagam_prep$pregam)

## 3) Bayesian PAMM via brms ---------------------------------------------------
set.seed(101)
pam_brms <- brms::brm(
  formula = ped_status ~ s(tend, k = 20) + offset(offset),
  data = ped,
  family = poisson(),
  chains = 2,
  cores = 2,
  iter = 3000,
  warmup = 750,
  init = 0,
  seed = 101,
  refresh = 0,
  silent = 2
)

## prediction via the unified pammtools inference interface -------------------
## (identical to the code shown in the article).
conf_level <- 0.95
newdata_base <- make_newdata(ped, tend = sort(unique(ped$tend))) %>%
  arrange(tend)
method_labels <- c(
  mgcv = "PAMM (mgcv)",
  jagam = "Bayesian PAMM (jagam)",
  brms = "Bayesian PAMM (brms)"
)

# brms as a pammtools backend: the two primitives of the unified interface.
sim_hazard.brmsfit <- function(object, newdata, nsim = NULL, ...) {
  newdata$offset <- 0
  ep <- brms::posterior_epred(
    object,
    newdata = newdata,
    summary = FALSE,
    re_formula = NA
  )
  t(as.matrix(ep))
}
get_hazard.brmsfit <- function(object, newdata, ...) {
  rowMeans(sim_hazard.brmsfit(object, newdata, ...))
}

get_cumu <- function(model, method_id, newdata, conf_level = 0.95) {
  if (inherits(model, "jam")) {
    class(model) <- unique(c(class(model), "gam", "glm", "lm"))
  }
  ci_type <- if (inherits(model, "brmsfit")) "sim" else "default"
  pred <- add_cumu_hazard(
    newdata,
    model,
    ci = TRUE,
    ci_type = ci_type,
    se_mult = stats::qnorm((1 + conf_level) / 2),
    alpha = 1 - conf_level,
    boundary = FALSE
  )
  tibble::tibble(
    tend = pred$tend,
    method_id = method_id,
    method = unname(method_labels[method_id]),
    cumu_hazard = pred$cumu_hazard,
    cumu_lower = pred$cumu_lower,
    cumu_upper = pred$cumu_upper
  )
}

pred_df <- bind_rows(
  get_cumu(pam_mgcv, "mgcv", newdata_base, conf_level = conf_level),
  get_cumu(pam_jagam, "jagam", newdata_base, conf_level = conf_level),
  get_cumu(pam_brms, "brms", newdata_base, conf_level = conf_level)
)

saveRDS(pred_df, "vignettes/bayesian_results.rds")
message(sprintf(
  "saved vignettes/bayesian_results.rds (%.1f kB, %d rows)",
  file.size("vignettes/bayesian_results.rds") / 1e3,
  nrow(pred_df)
))
