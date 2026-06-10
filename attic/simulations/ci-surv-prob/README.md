# Simulation study: CIs for S(t|x) from PAMMs — coverage and width (issue #285)

[Issue #285](https://github.com/adibender/pammtools/issues/285) claims that
`add_surv_prob(..., ci_type = "default")` — which transforms pointwise hazard
confidence limits via `exp(-cumsum(ci_bound_hazard * intlen))` — is overly
conservative, because cumulating pointwise limits implicitly assumes perfect
correlation of the hazard estimates across intervals. That claim was based on a
single (interval-censored) data example and has been disputed. This study
estimates **pointwise coverage and width** of nominal-95% CIs for S(t|x) for
the methods implemented in `add_surv_prob()`:

| method | construction | a priori weakness |
|---|---|---|
| `default` | transform hazard CI bounds, cumulate | ignores between-interval covariance |
| `delta` | delta method on the probability scale (`Vp`) | symmetric, linearization error, can leave [0,1] |
| `sim100` | posterior draws of coefficients, empirical quantiles, `nsim = 100` (package default) | Monte Carlo noise in 2.5% quantiles |
| `sim500` | same with `nsim = 500` | — |

Nominal levels are aligned across methods: `se_mult = qnorm(0.975)` for
default/delta (the package default `se_mult = 2` corresponds to ~95.4% and
would inflate apparent coverage), equal-tail 2.5%/97.5% quantiles for sim.

## Design

Data are generated with `pammtools::sim_pexp()` on a fine grid
(`seq(0, 10, length.out = 101)`), log-hazard
`eta = b0 + f0(t) + 0.7 * x1`, `x1 ~ N(0,1)`; `b0` is calibrated per baseline
so that true S(10 | x1 = 0) = 0.07. Because `sim_pexp()` evaluates the formula
at interval starts, the DGP is *exactly* piecewise-constant on the DGP grid and
the true S(t|x) is computed exactly (`R/dgp.R::true_surv()`); there is no
numerical-integration error in the truth.

Crossed factors (36 scenarios + 1 sensitivity arm), chosen because they should
affect the methods *differently*:

| factor | levels | rationale |
|---|---|---|
| baseline log-hazard f0(t) | constant; increasing (`0.25 t`, Gompertz-type); non-monotone (`dgamma(t, 8, 2) * 6`) | wiggliness drives the between-interval correlation of fitted hazards — exactly what `default` ignores |
| sample size n | 100, 250, 1000 | SE magnitude → delta's linearization error and [0,1] violations |
| censoring | ~20% vs ~60% (`Exp(c)` + administrative at t = 10; `c` calibrated analytically) | information loss at late times |
| fitting cut grid | 10 equidistant vs 40 event-time-quantile-based intervals | number of cumulated terms governs `default`'s behavior; quantile cuts mimic `as_ped()` defaults |

Sensitivity arm (scenario 37): non-monotone baseline, n = 250, 20% censoring,
quantile grid, refit with `s(tend, k = 20)` instead of k = 10.

Fitted model (correctly specified): `pamm(ped_status ~ s(tend, k = k) + x1,
method = "REML")`, with k = 7 for the J = 10 equidistant grids and k = 10 for
the J = 40 quantile grids. k must stay safely below the number of *populated*
intervals: under heavy censoring the trailing equidistant intervals are often
empty and mgcv errors when `s(tend, k)` has fewer unique `tend` values than k
(with k = J = 10, up to 84% of fits failed in the constant-baseline/60%-
censoring/n=100 cell). Residual fit failures (~1–2% worst case) are recorded
and reported; note they coincide with datasets lacking late events, so
late-time coverage is conditional on estimability. Within each fit, CIs are evaluated at 2 covariate profiles
(x1 = 0 and x1 = 1.5, HR ≈ 2.9) × 5 time points: the fitting-grid endpoints
closest to where true S(t|x) crosses 0.9, 0.75, 0.5, 0.25, 0.1 (labels
S90...S10). This probes boundary proximity (S near 1 or 0) without extra
scenarios. Evaluation times are interval endpoints by construction
(`make_newdata()` snaps to endpoints; `add_surv_prob()` integrates over all
supplied intervals — the code passes *all* endpoints and filters afterwards).

Metrics per (scenario, method, profile, time label), aggregated over
replications: coverage (with Wilson CI), miss direction (above/below),
mean/median width, paired width ratio vs delta, bias and RMSE of the point
estimate, fraction of delta CIs leaving [0,1] (plus coverage of the clipped
delta CI), fit failure and non-convergence rates.

Out of scope (deliberate): model misspecification; `Vp` vs `Vc`
(`unconditional = TRUE`); smoothing-parameter selection method (REML fixed);
interval censoring — the observation behind issue #285 came from an
interval-censored example, and a coarsened-observation-time factor is the most
interesting extension, but PAMMs handle right censoring natively and the study
stays in that setting.

## Reproducibility

Master seed 285; one pre-generated L'Ecuyer-CMRG stream per (scenario, rep)
(`RESERVED_REPS = 1000` streams reserved per scenario), so results are
independent of `--cores`/scheduling, individually reproducible, and extendable
to 1000 reps per scenario without changing existing replications. Runs are
checkpointed per scenario and resumable. Scenario preparation (censoring-rate
calibration is analytic; quantile cut points come from a fixed seeded reference
simulation with n = 20000) is cached in `results/prep.rds`.

## How to run

```sh
# 0. install the package, e.g. R CMD INSTALL --no-docs .  (from the repo root)

# 1. sanity/oracle checks (do this once after any code change)
Rscript sanity-checks.R          # ~ a few minutes; --quick skips the slowest check

# 2. pilot (e.g. this container: 4 cores, ~25 reps -> MC SE for coverage ~4pp)
Rscript run-sim.R --reps 25 --cores 4
Rscript eval-sim.R --tag -pilot

# 3. full study (LRZ CoolMUC-4): see lrz/setup-lrz.md, then
#    sbatch lrz/job-full.slurm    # 500 reps, 16 cores, ~1-2 h
```

Outputs: `results/raw/scen-XX.rds` (per-replication rows, gitignored),
`output/` (committed): `coverage*.csv`, `width-ratio*.csv`, `failures*.csv`,
`summary*.txt`, and figures `fig-coverage-x*.png`, `fig-width-ratio*.png`.

## Interpretation guide

If the issue-#285 claim is correct, `default` should show coverage well above
0.95 and the largest widths, with the gap growing with time (more cumulated
intervals), with J = 40 vs J = 10, and with baseline wiggliness. If the
coauthors are right, `default` coverage should sit near (or below) 0.95 in most
cells. Delta is expected to undercover near the boundaries (S90/S10) at small
n; sim100 vs sim500 shows whether the shipped `nsim` default is adequate. With
25 pilot reps the Monte-Carlo SE of a coverage estimate is ~4–5pp — pilot
results indicate ordering, not precise levels; the 500-rep run has MC SE ~1pp.
