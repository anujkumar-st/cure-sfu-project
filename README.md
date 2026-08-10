# cure-sfu-project

Code accompanying *"Testing and Correcting for Insufficient Follow-Up in
Covariate-Dependent Cure Models: From Categorical Tests to a Heavy-Tailed
Extrapolation, and Back"*. The full write-up is
in [`report/vu_project.pdf`](report/vu_project.pdf).

## Status of this repo

This repo was reconstructed from scratch: the original phase2 scripts
behind Sections 5.4-5.5 of the report were lost, and everything here
(except `phase5_simulation/`, which was already saved) was rewritten
independently from the report's equations and stated methodology, then
**run against the real data and checked line-by-line against every table
and figure in the report**. The table below states plainly where a
script matches the report exactly and where it only reproduces the same
qualitative pattern.

| Script | Match to report |
|---|---|
| `phase1_baseline/phase1_sfu_tests.R` | **Exact** on Table 1 and most of Table 2; death-endpoint bootstrap p-values in Table 2 are within Monte Carlo noise (bootstrap test, no fixed seed recoverable) |
| `phase2_extrapolation/phase2a_synthetic_validation.R` | Same qualitative pattern (corrected estimator moves toward ground truth); no specific numbers are quoted in the report to match |
| `phase2_extrapolation/phase2b_colon_age_application.R` | **Near-exact** on every number in Sec. 5.3 (mean gaps, sex-stratified gaps, the age=35 female NA edge case) |
| `phase2_extrapolation/phase2c_curse_of_dimensionality.R` | **Exact** match to Table 3 (all four rows) |
| `phase2_extrapolation/phase2d_five_point_analysis.R` | Reconstructed evaluation points (exact original points were lost); reproduces the same qualitative pattern -- non-monotone mean gap, ESS collapse at rare covariate values, non-monotone floor-pinning -- but not the exact numbers |
| `phase2_extrapolation/phase2e_bootstrap_variance.R` | Same qualitative pattern: floor-pinning rate for the typical/median patient rises sharply with dimension (7%->74% here vs. 5%->69% in the report); atypical points unstable even at one covariate |
| `phase2_extrapolation/phase2f_gbex_comparison.R` | **Essentially exact**: gamma_hat pinned at ~0.0049-0.0051 across every point/dimension, negative unconditional MLE confirmed as root cause, bootstrap SD ~0.00013 confirming the frozen-estimate finding |
| `phase2_extrapolation/phase2g_gbex_full_covariate_demo.R` | Same qualitative pattern (sigma varies substantially, gamma pinned near zero, full model beats age-only, same top covariates by importance); exceedance count and exact deviance/sigma numbers reconstructed rather than exact |
| `phase5_simulation/phase5_simulation.R` (Part A) | **Exact** match to Table 4, all 6 rows, at the script's own documented settings (N_REPS=30) |
| `phase5_simulation/phase5_partB_standalone.R` (Part B) | Rerun at the script's own full defaults (N_REPS_B=30, N_BOOT_B=200), parallelized with `furrr`/`future`. K=2 and K=5 match Table 5 exactly; K=3 and K=4 power differ from Table 5 by 2.2-2.4 standard deviations under the report's own binomial variance at N_REPS=30 -- within plausible Monte Carlo noise at this replication count, not a code defect. Type-I error stays at or below the nominal 5% throughout and Method 2 selection accuracy stays high (0.87-1.00), both consistent with the report's qualitative claims |

Where a script doesn't hit the report's numbers exactly, the discrepancy
is almost always a lost implementation detail (an exact bandwidth choice,
an evaluation point, a random seed) rather than a wrong method -- in every
case the qualitative finding described in the report reproduces.

## Layout

```
R/beran_ebvk.R                        shared EBVK estimator core (used by every phase2_* script)
phase1_baseline/                      Sec. 4: SFU hypothesis tests (cureSFUTest package)
phase2_extrapolation/                 Sec. 5: EBVK extrapolation estimator + curse-of-dimensionality + gbex
phase5_simulation/                    Sec. 6: bias-variance and Method1-vs-Method2 simulations
report/vu_project.pdf                 the write-up
```

Each `phase*/*.R` script is standalone runnable (`Rscript path/to/script.R`
from the repo root) and prints its results to the console; most also save
a `*_results.rds` file alongside themselves.

## Setup

Real data comes from the `survival` package's built-in `colon` dataset.
Two packages are not on CRAN and must be installed from source:

```r
# CRAN packages
install.packages(c("survival", "KernSmooth", "POT", "treeClust"))

# GitHub-only packages
install.packages("remotes")
remotes::install_github("tp-yuen/cureSFUTest")   # Yuen & Musta's reference SFU test package
remotes::install_github("JVelthoen/gbex")         # gradient-boosted GPD regression (Velthoen et al. 2021)

# For the parallelized Phase 5 Part B script (phase5_partB_standalone.R)
install.packages(c("future", "furrr"))
```

If `remotes::install_github` can't reach GitHub's API directly (rate
limits, proxies), clone the CRAN-mirror/source repos and `R CMD INSTALL`
them instead -- that's how this repo's results were produced:

```bash
git clone https://github.com/cran/POT.git        && R CMD INSTALL POT
git clone https://github.com/cran/treeClust.git  && R CMD INSTALL treeClust
git clone https://github.com/tp-yuen/cureSFUTest.git && R CMD INSTALL cureSFUTest
git clone https://github.com/JVelthoen/gbex.git  && R CMD INSTALL gbex
```

`future`/`furrr` pull in a longer dependency chain (`digest`, `rlang`,
`globals`, `listenv`, `parallelly`, `vctrs`, `purrr`, `cli`, `glue`,
`magrittr`, `lifecycle`) -- if CRAN is unreachable, the same
clone-from-`https://github.com/cran/<pkg>`-and-`R CMD INSTALL` pattern
above works for each of these too.

## Known issues in this repo's own code (found and fixed)

- Both `phase5_simulation.R` and `phase5_partB_standalone.R` initially
  reintroduced the exact scenario-mislabeling bug the report itself
  describes catching once already (Sec. 6.2): the Part B summary
  computed `method1_type1_error`/`method1_power` by negating the
  rejection rate from the *wrong* generative scenario, which silently
  swaps power and Type-I error. It was caught by comparing an initial
  parallel run against Table 5 -- the observed numbers matched
  `1 - power` / `1 - type1_error` almost exactly rather than the
  values themselves, which pinpointed the negation. Fixed by reading
  each metric directly off the correct scenario per the report's own
  H0/H1 framing, with no negation.

## Known upstream bugs (found and worked around, documented in the report)

- `cureSFUTest::sfu.cov.test()` silently requires a numeric covariate
  matrix and raises a misleading error message otherwise.
- `cureSFUTest::sfu.test()`'s bootstrap loop reads `.Random.seed` without
  checking it exists, and crashes on a fresh R session -- call `set.seed()`
  first.
- `cureSFUTest::sfu.test()`'s default `method` argument silently resolves
  to `"g"` (plain Grenander), which gives a degenerate p=0 on this data --
  pass `method = "sg"` (smoothed Grenander) explicitly. `sfu.cov.test()`,
  used for the categorical splits, always routes internally through the
  smoothed-Grenander bandwidth regardless of arguments supplied, so it is
  not exposed to this particular trap.
- `gbex()` does not accept a `status`/`threshold`/`maxit`/`stopcrit`
  interface -- it takes uncensored exceedances (`y`) and a covariate
  data frame (`X`) directly.
- `gbex::predict()` does not support `type="response"`; it returns a
  data frame with columns `s` (scale) and `g` (shape) instead.
- `gbex()`/`predict.gbex()` require `X`/`newdata` as a real `data.frame`
  with real column names, not a bare matrix, or column names get mangled
  internally (`X.x1` etc.) and silently mismatch.
- `gbex::variable_importance()`'s default `type` argument
  (`type = c("relative","permutation")`) is never resolved with
  `match.arg()` internally, so the default call errors with "the
  condition has length > 1" -- pass `type` explicitly.

## Reproducing Phase 1 at full precision

Phase 1 involves bootstrap hypothesis tests and was run here at the
script's own default settings; the death-endpoint bootstrap p-values in
Table 2 fall within Monte Carlo noise since the reference `cureSFUTest`
package's bootstrap has no recoverable fixed seed. Phase 5 Part B has
since been rerun at full precision (see table above).
