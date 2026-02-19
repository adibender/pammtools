# Robustness Failure Log

Generated on 2026-02-19 from robustness contexts:

- `test-ped-join-semantics.R`
- `test-penalized-lag-lead-robustness.R`
- `test-plot-data-consistency.R`
- `test-transition-prob-robustness.R`

Test run summary for these contexts: `PASS 35`, `FAIL 8`, `WARN 0`, `SKIP 0`.

## Failure Records

### 1) `right_join.ped` drops unmatched right-hand rows

- Test: `tests/testthat/test-ped-join-semantics.R` (`right_join.ped keeps unmatched right-hand keys`)
- Observed failure:
  - `nrow(out_right) == 2` while `nrow(ref_right) == 3`
  - unmatched key `id = 999` is missing in `out_right`
- Code path:
  - `R/tidyverse-methods.R` (`right_join.ped`)
- Why failure is correct under current code:
  - Implementation calls `inner_join(...)` internally instead of `right_join(...)`, which necessarily removes unmatched right-side keys.
- Classification: `Code bug likely`
- Pending PR verification:
  - Open PR `#262` (`Fix ped tidyverse semantics and attribute sanitization`) updates `right_join.ped` to use `right_join(...)`.
  - Diff evidence: `R/tidyverse-methods.R` replacement from `inner_join` to `right_join`.
- Resolution status: `covered/fixed by pending PR #262`

### 2) `fdl` smooth basis not discoverable by `mgcv::gam`

- Test: `tests/testthat/test-penalized-lag-lead-robustness.R` (`fdl basis is discoverable by mgcv smooth constructor`)
- Observed failure:
  - `no applicable method for 'smooth.construct' applied to an object of class "fdl.smooth.spec"`
- Code path:
  - Constructor exists in `R/penalized-lag-lead.R` as `smooth.construct.fdl.smooth.spec`.
- Why failure is correct under current code:
  - The method is not registered as an S3 method in `NAMESPACE`, so `mgcv` dispatch cannot find it when `bs = "fdl"` is requested.
- Classification: `Code bug likely`
- Pending PR verification:
  - No currently open PR found that touches `R/penalized-lag-lead.R` or registers this method.
- Resolution status: `open issue needed`

### 3) `gg_re()` mutates global ggplot theme

- Test: `tests/testthat/test-plot-data-consistency.R` (`gg_re does not mutate global ggplot theme`)
- Observed failure:
  - `ggplot2::theme_get()` differs after calling `gg_re()`.
- Code path:
  - `R/convenience-plots.R` (`gg_re`)
- Why failure is correct under current code:
  - `gg_re()` currently calls `theme_set(theme_bw())`, which mutates global theme state as a side effect.
- Classification: `Code bug likely`
- Pending PR verification:
  - Open PR `#263` (`Fix mgcv convenience plot side effects and CI scaling`) replaces `theme_set(theme_bw())` with `theme_bw()`.
- Resolution status: `covered/fixed by pending PR #263`

### 4) `tidy_smooth()` and `tidy_smooth2d()` CI scaling too narrow

- Test: `tests/testthat/test-plot-data-consistency.R` (`tidy smooth confidence intervals use standard normal scaling`)
- Observed failure:
  - CI width equals `1 * se`, not `qnorm(0.975) * se`.
- Code path:
  - `R/tidiers.R` (`tidy_smooth`, `tidy_smooth2d`)
- Why failure is correct under current code:
  - Current implementation computes `ci_lower = fit - se` and `ci_upper = fit + se`.
  - Robustness expectation uses 95% normal approximation (`fit ± 1.96 * se`) unless otherwise specified.
- Classification: `Code bug likely`
- Pending PR verification:
  - Open PR `#263` adds `conf_level` and scales CI with `z = qnorm((1 + conf_level)/2)`.
- Resolution status: `covered/fixed by pending PR #263`

### 5) `add_trans_prob()` depends on input row order

- Test: `tests/testthat/test-transition-prob-robustness.R` (`transition probabilities are invariant to row order within groups`)
- Observed failure:
  - All rows mismatch between sorted and shuffled inputs (`380/380` mismatches).
- Code path:
  - `R/add-functions.R` (`add_trans_prob`, `get_trans_prob`)
- Why failure is correct under current code:
  - Computation relies on cumulative hazard differences and matrix products without enforcing stable ordering in all paths.
  - Shuffling rows changes inferred increments and therefore transition probabilities.
- Classification: `Code bug likely`
- Pending PR verification:
  - Open PR `#261` (`Fix trans prob`) introduces explicit ordering by grouping/transition/time and rewrites transition matrix accumulation logic.
- Resolution status: `covered/fixed by pending PR #261`

## Notes

- Failing robustness tests are intentionally retained as high-signal regression checks.
- For CI gating, these failures can be temporarily marked as expected only if linked to tracked issues/PRs above.
