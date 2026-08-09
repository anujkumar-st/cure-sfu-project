## ============================================================================
## Phase 1: baseline replication of the sufficient-follow-up (SFU) tests
## on survival::colon, using the real reference package cureSFUTest
## (Yuen & Musta 2024/2026; Yuen, Musta & Van Keilegom 2025).
##
## Reproduces report Section 4: Table 1 (unconditional test, both endpoints)
## and Table 2 (categorical test by sex and by age, Method 1 = corrected
## intersection-union aggregation).
## ============================================================================

need <- c("survival")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("cureSFUTest", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("tp-yuen/cureSFUTest")
}
suppressMessages({
  library(survival)
  library(cureSFUTest)
})

data(colon)

## Bug fix #1 (documented in report Sec. 4.2): the package's bootstrap loop
## reads .Random.seed without checking it exists, which crashes on a fresh
## R session. Calling set.seed() first avoids this.
set.seed(2026)

## ---- data setup: two endpoints, time in YEARS ----
rec   <- colon[colon$etype == 1, ]   # recurrence
death <- colon[colon$etype == 2, ]   # death (all-cause)

rec$time_yr   <- rec$time   / 365.25
death$time_yr <- death$time / 365.25

## tau_G = max observed follow-up (sfu.test's own default), same for both
## endpoints since they come from the same trial
tau_G <- max(rec$time_yr)
cat(sprintf("tau_G (max observed follow-up) = %.3f years\n", tau_G))

## Clinically justified tau choices (report Sec. 3):
##  - recurrence: tau_G + 2 years, justified by Pastorino et al. (2025,
##    JAMA Oncology): pooled n=35,213 across 15 trials shows recurrence
##    risk falls below 0.5% by 6 years post-surgery
##  - death: tau_G + 10 years (no disease-specific cure horizon exists for
##    all-cause mortality, so a more generous horizon is used)
tau_rec   <- tau_G + 2
tau_death <- tau_G + 10
cat(sprintf("tau (recurrence) = %.3f years, tau (death) = %.3f years\n",
            tau_rec, tau_death))

## ============================================================================
## Table 1: unconditional SFU test, both endpoints
## ============================================================================
cat("\n=== Table 1: unconditional sufficient-follow-up test ===\n")

## method = "sg" (smoothed-Grenander) is Yuen & Musta's actual test
## statistic; the package's default match.arg order would silently pick
## "g" (plain Grenander) instead, which gives degenerate p=0 for both
## endpoints on this data -- method="sg" must be passed explicitly.
fit_rec   <- sfu.test(rec$time_yr,   rec$status,   tau = tau_rec,   method = "sg")
fit_death <- sfu.test(death$time_yr, death$status, tau = tau_death, method = "sg")

tbl1 <- data.frame(
  endpoint = c("Recurrence", "Death (all-cause)"),
  p_value  = c(fit_rec$p.value, fit_death$p.value),
  conclusion = ifelse(c(fit_rec$p.value, fit_death$p.value) < 0.05,
                       "Sufficient follow-up", "Insufficient follow-up")
)
print(tbl1, row.names = FALSE)

## ============================================================================
## Table 2: categorical SFU test (Method 1 = intersection-union test),
## by sex and by age (binned at the median)
## ============================================================================
cat("\n=== Table 2: categorical sufficient-follow-up test (Method 1) ===\n")

## Bug fix #2 (documented in report Sec. 4.2): Method 1's per-subgroup
## p-values are read directly off the SAME sfu.cov.test() call used for
## Method 2, which internally always uses the correct undersmoothed
## bandwidth bw = min(tau_n^(-7/30), tau/2). An earlier hand-written
## version used the wrong bandwidth AND aggregated with min(p_x) instead
## of the correct max(p_x) implied by intersection-union logic
## (H0: insufficient for SOME x; H1: sufficient for ALL x).
run_categorical <- function(dat, cov, tau, n.boot = 1000) {
  res <- sfu.cov.test(dat$time_yr, dat$status, dat[[cov]], tau = tau, n.boot = n.boot)
  ## Method 1 (IUT): H0 = insufficient for SOME x, H1 = sufficient for ALL x.
  ## "Sufficient for all" is only declared if EVERY subgroup individually
  ## clears significance, i.e. max(p) < 0.05 -- getting this direction
  ## backwards would flip every conclusion in Table 2 (see report Sec. 4.2).
  list(p_by_subgroup = res$p.value, x.unique = res$x.unique,
       method1_p = max(res$p.value),   # intersection-union aggregation
       method1_conclusion = ifelse(max(res$p.value) < 0.05,
                                    "Sufficient for all subgroups",
                                    "Not sufficient for all subgroups"),
       full = res)
}

rec$age_bin <- as.numeric(rec$age > median(rec$age))
death$age_bin <- as.numeric(death$age > median(death$age))

rec_sex   <- run_categorical(rec,   "sex",     tau_rec)
rec_age   <- run_categorical(rec,   "age_bin", tau_rec)
death_sex <- run_categorical(death, "sex",     tau_death)
death_age <- run_categorical(death, "age_bin", tau_death)

fmt_p <- function(p) paste(sprintf("%.3f", p), collapse = ", ")

tbl2 <- data.frame(
  split = c("Recurrence x sex", "Recurrence x age", "Death x sex", "Death x age"),
  subgroup_p_values = c(fmt_p(rec_sex$p_by_subgroup),  fmt_p(rec_age$p_by_subgroup),
                         fmt_p(death_sex$p_by_subgroup), fmt_p(death_age$p_by_subgroup)),
  method1_conclusion = c(rec_sex$method1_conclusion,  rec_age$method1_conclusion,
                          death_sex$method1_conclusion, death_age$method1_conclusion)
)
print(tbl2, row.names = FALSE)

## ---- Method 2: bootstrap worst-category selection instability check ----
## (report Sec. 4.3: reproduces the source paper's own flagged instability
## where x.ast.hat's bootstrap selection can disagree with a naive
## largest-p-value ranking; checked here on death x age and death x sex,
## across two separate (re-bootstrapped) runs. Run 1 reuses the death_sex/
## death_age fits already computed above; run 2 re-runs with fresh
## bootstrap draws.)
cat("\n=== Method 2 selection-instability check (death endpoint) ===\n")
report_instability <- function(label, fit) {
  cat(sprintf("  %s: x.ast.hat = %s, naive largest-p subgroup = %s\n",
              label, fit$full$x.ast.hat, which.max(fit$p_by_subgroup)))
}
cat("-- run 1 --\n")
report_instability("death x sex", death_sex)
report_instability("death x age", death_age)
cat("-- run 2 --\n")
death_sex_r2 <- run_categorical(death, "sex", tau_death)
death_age_r2 <- run_categorical(death, "age_bin", tau_death)
report_instability("death x sex", death_sex_r2)
report_instability("death x age", death_age_r2)

## ---- save results ----
out_dir <- "phase1_baseline"
if (!dir.exists(out_dir)) out_dir <- "."
saveRDS(list(table1 = tbl1, table2 = tbl2), file.path(out_dir, "phase1_results.rds"))
cat("\nSaved phase1_results.rds\n")
