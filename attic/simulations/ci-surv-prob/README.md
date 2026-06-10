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
censoring/n=100 cell). Residual fit failures (13.2% and 8.6% in the two
sparsest cells — constant/nonmonotone baseline, n = 100, 60% censoring,
equidistant grid — and < 1% everywhere else in the 500-rep run) are recorded
and reported; they coincide with datasets lacking late events, so late-time
coverage is conditional on estimability. Within each fit, CIs are evaluated at 2 covariate profiles
(x1 = 0 and x1 = 1.5, HR ≈ 2.9) × 5 time points: the fitting-grid endpoints
closest to where true S(t|x) crosses 0.9, 0.75, 0.5, 0.25, 0.1 (labels
S90...S10). This probes boundary proximity (S near 1 or 0) without extra
scenarios. Evaluation times are interval endpoints by construction
(`make_newdata()` snaps to endpoints; `add_surv_prob()` integrates over all
supplied intervals — the code passes *all* endpoints and filters afterwards).

Two consequences of the endpoint-snapping design are handled explicitly in
`eval-sim.R`:

- **Endpoint availability.** Even a successful fit yields no CI at an eval
  endpoint when no subject is at risk that late (the endpoint is then absent
  from the fitted ped). Coverage therefore conditions on availability;
  `avail_frac` in `coverage*.csv` records how selective this is (worst cells:
  51% and 79% of successful fits contribute, constant baseline / 60%
  censoring / equidistant at S10/x1 = 0 for n = 100 and 250; ≥ 95%
  everywhere else). This
  conditioning is informative — only datasets with late at-risk subjects
  contribute — and inflates late-time coverage in exactly the
  heavy-censoring cells where it bites.
- **Collapsed labels.** On the coarse J = 10 grid, several S-labels can snap
  to the *same* endpoint (24 of 74 scenario × profile blocks, mostly at
  x1 = 1.5, where e.g. nominal S90 sits at a realized true S of ~0.47).
  Duplicated endpoints are deduplicated before aggregation (keeping the label
  closest to the realized true S), and marginal-by-time summaries bin by
  *realized* true S (`S_bin`), not by nominal label.

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
results indicate ordering, not precise levels; the 500-rep run has MC SE ~1pp
per cell (marginal averages share replications across cells, so their
uncertainty is *not* 1pp/sqrt(#cells)).

## Findings (full run, 500 reps; see `output/*-full*`)

Headline (share of deduplicated cells with coverage < 0.93 / 0.93–0.97 /
> 0.97; nominal 0.95):

| method | < 0.93 | 0.93–0.97 | > 0.97 | min | paired width ratio vs delta, median (q10–q90) |
|---|---|---|---|---|---|
| `default` | 2% | 62% | 36% | 0.918 | 1.07 (1.00–1.35) |
| `delta` | 36% | 64% | 0.3% | 0.723 | 1 |
| `sim100` | 57% | 43% | 0% | 0.836 | 0.95 (0.93–0.97) |
| `sim500` | 19% | 81% | 0% | 0.850 | 0.99 (0.96–1.00) |

1. **`default` over-covers, mostly mildly.** Mean coverage 0.965, rising
   toward the boundary (0.976 in the realized-S10 bin vs 0.948 at S90) and
   with baseline wiggliness (0.973 nonmonotone vs 0.961 constant). The median
   width cost vs delta is only ~7%, but right-skewed (q90 = 1.35, largest
   late + wiggly + heavy censoring). No support for a grid-resolution effect
   (0.965 both grid types — but J and k are confounded with grid type by
   design, so grid effects are not separately identified). `default` is the
   only method that never violates [0,1], has worst-case coverage 0.918, and
   holds up in the cells where all other methods fail.
2. **`delta` is anti-conservative in the tails**: 0.909 mean in the S10 bin,
   worst cell 0.723 (n = 100, 60% censoring); up to ~90% of its intervals
   leave [0,1] there, and clipping does not repair coverage (misses are
   upper bounds *below* the truth, not bounds outside [0,1]).
3. **Tail undercoverage is only partly a delta problem.** At n = 1000 with
   heavy censoring, where [0,1] violations are rare, delta still covers 0.834
   and even sim500 only 0.858 at S10: the fitted late hazard is biased
   (downward-biased S, ~0.5 SE) under heavy censoring, which caps *any*
   correctly-sized symmetric/equal-tailed interval near ~92%. `default`
   escapes via much wider, asymmetric intervals.
4. **`sim` with the shipped `nsim = 100` undercovers everywhere** (mean
   0.924; 57% of cells below 0.93). Much of this is the type-7
   `stats::quantile()` estimator, whose 2.5%/97.5% sample quantiles at
   B = 100 have expected central mass ≈ 0.93, i.e. intervals are
   systematically too narrow (sim100 is 5% narrower than delta paired);
   `nsim = 500` shrinks the deficit (0.99 ratio) but remains below nominal
   in the tails for the bias reason above.
5. **Package implications (suggested, not implemented here):** document
   `default`'s conservatism (it is the safest worst-case choice) rather than
   silently relying on it; warn when delta bounds leave [0,1]; fix the
   sim-CI quantile estimator (e.g. `type = 6`) and/or raise the default
   `nsim` — raising `nsim` alone is brute force and does not reach nominal.
   Any change of the *default* `ci_type` should wait for an unconfounded
   J/k arm and interval-censored scenarios.

Caveats: pointwise (not simultaneous) coverage; right censoring only;
correctly specified model; `Vp`-based (conditional-on-smoothing-parameters)
uncertainty throughout; late-time coverage conditional on endpoint
availability (see above).
