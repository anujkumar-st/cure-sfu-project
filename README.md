# cure-sfu-project

Code accompanying *"Testing and Correcting for Insufficient Follow-Up in
Covariate-Dependent Cure Models: From Categorical Tests to a Heavy-Tailed
Extrapolation, and Back"* (Anuj Kumar, Aug 2026). The full write-up is
in [`report/vu_project.pdf`](report/vu_project.pdf).

## Project status

This repository was reconstructed from the report methodology and the available
saved analyses. The numerical artifacts in the repository are the regenerated
results used for the current version of the report.

The **report retains its original phase numbering** (Phase 1, Phase 2, Phase 5)
because those labels refer to the study/report chronology. The repository itself
uses a cleaner three-stage layout (`01`, `02`, `03`) so that a reader does not
encounter unexplained gaps such as a missing Phase 3 or Phase 4 directory.

| Repository stage | Report phase | Purpose | Status |
|---|---|---|---|
| `01_sfu_tests/` | Phase 1 | Sufficient-follow-up hypothesis tests using `cureSFUTest` | Regenerated in R; verified results stored under `01_sfu_tests/results/` |
| `02_ebvk_extrapolation/` | Phase 2 | From-scratch EBVK extrapolation, dimensionality diagnostics, and `gbex` comparison | Core estimator and exploratory extensions audited |
| `03_simulation_study/` | Phase 5 | Controlled simulation studies for finite-sample behavior and categorical methods | Part A and Part B regenerated |

### Analysis scripts

| Script | What it does |
|---|---|
| `01_sfu_tests/01_sfu_tests.R` | Baseline SFU analysis on `survival::colon`; explicit 0/1 event coding and numeric covariate-matrix input |
| `02_ebvk_extrapolation/01_synthetic_validation.R` | Synthetic EBVK sanity check; 20-rep check is exploratory rather than a full Monte Carlo validation |
| `02_ebvk_extrapolation/02_colon_age_application.R` | EBVK application to colon recurrence/death data over an age grid |
| `02_ebvk_extrapolation/03_curse_of_dimensionality.R` | Multivariate product-kernel dimensionality demonstration; exploratory extension |
| `02_ebvk_extrapolation/04_five_point_analysis.R` | Reconstructed five-point covariate-space analysis |
| `02_ebvk_extrapolation/05_bootstrap_variance.R` | Bootstrap diagnostic of tail-index floor-pinning and local instability |
| `02_ebvk_extrapolation/06_gbex_comparison.R` | Small-tail `gbex` comparison and initialization-sensitivity diagnostic |
| `02_ebvk_extrapolation/07_gbex_full_covariate_demo.R` | Ten-covariate `gbex` fit; reports in-sample deviance comparison, not predictive superiority |
| `03_simulation_study/01_bias_vs_tauc.R` | Part A: bias/MSE versus censoring horizon |
| `03_simulation_study/02_method1_vs_method2.R` | Part B: parallel Method 1/Method 2 simulation study |

### Regenerated results

The current Part A results use the seven-point censoring-horizon sweep
`4, 6, 8, 12, 20, 40, 80`; the simulated uncured 99% quantile at the target
covariate value is approximately 31.62.

The current Part B study uses `N_REPS = 30`, `N_BOOT = 200`, `K = 2:5`, and
`N_TOTAL = 400`. Empirical power is 0.90--0.97, empirical Type-I error is
0--0.067, and Method 2 selection accuracy is 0.967--1.000 in this study.
These are empirical simulation results; increasing the outer Monte Carlo count
later would give more precise operating-characteristic estimates.

## Layout

```text
R/
└── beran_ebvk.R                         shared EBVK estimator core

01_sfu_tests/
├── 01_sfu_tests.R
└── results/
    └── phase1_results.rds

02_ebvk_extrapolation/
├── 01_synthetic_validation.R
├── 02_colon_age_application.R
├── 03_curse_of_dimensionality.R
├── 04_five_point_analysis.R
├── 05_bootstrap_variance.R
├── 06_gbex_comparison.R
├── 07_gbex_full_covariate_demo.R
└── results/
    ├── 01_synthetic_validation_results.rds
    ├── 02_colon_age_application_results.rds
    ├── 03_curse_of_dimensionality_results.rds
    ├── 04_five_point_results.rds
    ├── 05_bootstrap_variance_results.rds
    ├── 06_gbex_comparison_results.rds
    └── 07_gbex_full_covariate_results.rds

03_simulation_study/
├── 01_bias_vs_tauc.R
├── 02_method1_vs_method2.R
└── results/
    ├── partA_results.csv
    ├── partA_results.rds
    ├── partA_bias_vs_tauc.png
    ├── partB_results.csv
    ├── partB_results.rds
    └── partB_level_power_vs_K.png

report/
├── vu_project.tex
└── vu_project.pdf
```

The analysis scripts are intended to be run with the repository root as the
working directory, so shared paths such as `R/beran_ebvk.R` resolve consistently.

## Setup

Real data comes from the `survival` package's built-in `colon` dataset.
Required packages include:

```r
# CRAN packages
install.packages(c("survival", "KernSmooth", "POT", "treeClust"))

# GitHub packages
install.packages(c("remotes", "future", "furrr"))
remotes::install_github("tp-yuen/cureSFUTest")
remotes::install_github("JVelthoen/gbex")
```

For the parallel simulation, `future`, `furrr`, and their dependencies must be
available. The Part B script automatically uses multiple workers when possible.

## Reproducibility notes

Phase 1 uses bootstrap hypothesis tests. The checked result is stored in
`01_sfu_tests/results/phase1_results.rds`.

Phase 2's multivariate categorical/continuous product-kernel layer is an
exploratory extension of the continuous-covariate EBVK construction; the report
does not claim that the original EBVK asymptotic theory automatically covers
this mixed-data implementation.

Phase 5 Part B records the study-level empirical power, Type-I error, and Method
2 category-selection accuracy. The present study uses 30 outer replications;
that is adequate for this exploratory study, while a larger `N_REPS` can be used
later when tighter Monte Carlo precision is desired.

## Known upstream/package behaviors worked around

- `cureSFUTest::sfu.cov.test()` requires a numeric covariate matrix despite its
  potentially misleading error message.
- `cureSFUTest::sfu.test()` uses `.Random.seed` in its bootstrap path, so the
  analysis scripts set a seed explicitly before bootstrap-based runs.
- `cureSFUTest::sfu.test()` defaults to the plain Grenander estimator (`"g"`).
  This project explicitly passes `method = "sg"` because the study uses the
  smoothed-Grenander test; this is package behavior, not an upstream defect.
- `gbex()` takes uncensored exceedances (`y`) and a covariate data frame (`X`)
  rather than a `status`/`threshold`/`maxit`/`stopcrit` interface.
- `gbex::predict()` returns scale (`s`) and shape (`g`) rather than a
  `type = "response"` prediction.
- `gbex()`/`predict.gbex()` require real data frames with column names; bare
  matrices can lead to internal name mangling.
- `gbex::variable_importance()` requires an explicit `type` argument because
  its default vector is not resolved with `match.arg()` internally.

## Report and repository conventions

The report's Phase 1 / Phase 2 / Phase 5 labels are retained where they are
part of the written study narrative. Repository directories instead use the
ordered workflow names `01_sfu_tests`, `02_ebvk_extrapolation`, and
`03_simulation_study`. This keeps the public repository navigable while
preserving the report's original structure and cross-references.
