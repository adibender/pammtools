# Verification rerun for issue #288 (type-6 sim-CI quantiles)

Reruns scenarios 1, 30, 33 (worst mean sim100 coverage in the full type-7
run) x 100 reps with the patched pammtools (branch `issue-288`, commit
4624263: `quantile(..., type = 6)` at all sim-CI sites). Because RNG streams
are deterministic per (scenario, rep), the rerun is paired with reps 1-100 of
the stored full run in `../results/raw` (type 7, LRZ job 5256116): identical
data, identical fits, identical posterior draws up to BLAS noise — bounds
differ only by quantile type.

Produced via (from the study directory, 2026-06-10):

```sh
R CMD INSTALL -l /tmp/pammlib-288 <worktree:issue-288>
Rscript -e '.libPaths(c("/tmp/pammlib-288", .libPaths())); source("run-sim.R")' \
  --scenarios 1,30,33 --reps 100 --cores 6 --out results-type6/raw
Rscript results-type6/compare.R
```

## Results (compare.R output, pooled over 2700 paired cells per method)

- Pairing check: point estimates and delta/default bounds agree to <= 1e-4
  between runs; only the sim bounds moved.
- sim100: paired coverage 0.907 -> 0.928 (+2.1pp; theory predicts ~+1.9pp
  from expected central mass 0.931 -> 0.950). 59 cells flip miss->cover,
  2 cover->miss (residual BLAS draw noise). Median width ratio 1.081.
- sim500: 0.919 -> 0.925 (+0.6pp; theory ~+0.4pp). Width ratio 1.016.
- Remaining sub-nominal coverage (esp. scenario 30: 0.919 at nsim = 100)
  is the late-hazard fit bias under heavy censoring discussed in the #285
  study README — out of scope for #288, affects delta and sim alike.
