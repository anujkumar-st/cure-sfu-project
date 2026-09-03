## ============================================================================
## Repository stage 01: baseline replication of the sufficient-follow-up (SFU) tests
## on survival::colon, using the real reference package cureSFUTest.
## ============================================================================

options(stringsAsFactors = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("cureSFUTest", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("tp-yuen/cureSFUTest")
}

suppressMessages({
  library(survival)
  library(cureSFUTest)
})

colon_dat <- survival::colon
set.seed(2026)

rec   <- colon_dat[colon_dat$etype == 1, , drop = FALSE]
death <- colon_dat[colon_dat$etype == 2, , drop = FALSE]

## survival::colon uses status = 0 for censored and status = 1 for event.
## cureSFUTest expects delta = 0 for censored and delta = 1 for event.
rec$event   <- as.numeric(rec$status == 1)
death$event <- as.numeric(death$status == 1)

if (sum(rec$event) == 0 || sum(death$event) == 0) {
  stop("Event coding check failed: zero events in one endpoint.")
}

rec$time_yr   <- rec$time / 365.25
death$time_yr <- death$time / 365.25

tau_G <- max(rec$time_yr, na.rm = TRUE)
tau_rec <- tau_G + 2
tau_death <- tau_G + 10

cat(sprintf("tau_G = %.6f years\n", tau_G))
cat(sprintf("tau_rec = %.6f years; tau_death = %.6f years\n", tau_rec, tau_death))
cat(sprintf("Events: recurrence=%d, death=%d\n\n", sum(rec$event), sum(death$event)))

## ---------------------------------------------------------------------------
## Table 1: unconditional tests
## ---------------------------------------------------------------------------

fit_rec <- sfu.test(rec$time_yr, rec$event, tau = tau_rec, method = "sg")
fit_death <- sfu.test(death$time_yr, death$event, tau = tau_death, method = "sg")

tbl1 <- data.frame(
  endpoint = c("Recurrence", "Death (all-cause)"),
  p_value = c(fit_rec$p.value, fit_death$p.value),
  conclusion = ifelse(
    c(fit_rec$p.value, fit_death$p.value) < 0.05,
    "Sufficient follow-up",
    "Insufficient follow-up"
  )
)

cat("=== Table 1 ===\n")
print(tbl1, row.names = FALSE)

## ---------------------------------------------------------------------------
## Table 2: categorical tests
## ---------------------------------------------------------------------------

rec$age_bin <- as.numeric(rec$age > median(rec$age, na.rm = TRUE))
death$age_bin <- as.numeric(death$age > median(death$age, na.rm = TRUE))

run_categorical <- function(dat, covariate, tau, n.boot = 1000) {
  X <- matrix(as.numeric(dat[[covariate]]), ncol = 1)
  colnames(X) <- covariate

  result <- sfu.cov.test(
    y = dat$time_yr,
    delta = dat$event,
    X = X,
    tau = tau,
    n.boot = n.boot
  )

  pvals <- as.numeric(result$p.value)
  method1_p <- max(pvals, na.rm = TRUE)

  list(
    p_by_subgroup = pvals,
    x.unique = result$x.unique,
    method1_p = method1_p,
    method1_conclusion = ifelse(
      method1_p < 0.05,
      "Sufficient for all subgroups",
      "Not sufficient for all subgroups"
    ),
    full = result
  )
}

rec_sex <- run_categorical(rec, "sex", tau_rec)
rec_age <- run_categorical(rec, "age_bin", tau_rec)
death_sex <- run_categorical(death, "sex", tau_death)
death_age <- run_categorical(death, "age_bin", tau_death)

fmt_p <- function(p) paste(sprintf("%.3f", p), collapse = ", ")

tbl2 <- data.frame(
  split = c("Recurrence x sex", "Recurrence x age", "Death x sex", "Death x age"),
  subgroup_p_values = c(
    fmt_p(rec_sex$p_by_subgroup),
    fmt_p(rec_age$p_by_subgroup),
    fmt_p(death_sex$p_by_subgroup),
    fmt_p(death_age$p_by_subgroup)
  ),
  method1_p = c(
    rec_sex$method1_p,
    rec_age$method1_p,
    death_sex$method1_p,
    death_age$method1_p
  ),
  method1_conclusion = c(
    rec_sex$method1_conclusion,
    rec_age$method1_conclusion,
    death_sex$method1_conclusion,
    death_age$method1_conclusion
  )
)

cat("\n=== Table 2 ===\n")
print(tbl2, row.names = FALSE)

## ---------------------------------------------------------------------------
## Method 2 repeatability check (descriptive; not a claim of instability)
## ---------------------------------------------------------------------------

set.seed(2026)
death_sex_r2 <- run_categorical(death, "sex", tau_death)
death_age_r2 <- run_categorical(death, "age_bin", tau_death)

method2_repeat <- data.frame(
  analysis = c("Death x sex", "Death x age"),
  run1 = c(as.integer(death_sex$full$x.ast.hat), as.integer(death_age$full$x.ast.hat)),
  run2 = c(as.integer(death_sex_r2$full$x.ast.hat), as.integer(death_age_r2$full$x.ast.hat))
)

cat("\n=== Method 2 repeatability check ===\n")
print(method2_repeat, row.names = FALSE)

## ---------------------------------------------------------------------------
## Save complete verified results
## ---------------------------------------------------------------------------

phase1_results <- list(
  tau_G = tau_G,
  tau_rec = tau_rec,
  tau_death = tau_death,
  table1 = tbl1,
  table2 = tbl2,
  method2_repeat = method2_repeat,
  fit_rec = fit_rec,
  fit_death = fit_death,
  recurrence_sex = rec_sex,
  recurrence_age = rec_age,
  death_sex = death_sex,
  death_age = death_age,
  death_sex_run2 = death_sex_r2,
  death_age_run2 = death_age_r2
)

saveRDS(phase1_results, file.path("01_sfu_tests", "results", "phase1_results.rds"))
cat("\nSaved 01_sfu_tests/results/phase1_results.rds\n")
