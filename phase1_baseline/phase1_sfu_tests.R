## ============================================================================
## Phase 1: baseline replication of the sufficient-follow-up (SFU) tests
## on survival::colon, using the real reference package cureSFUTest
## (Yuen & Musta 2024/2026; Yuen, Musta & Van Keilegom 2025).

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

set.seed(2026)


rec   <- colon[colon$etype == 1, ]  
death <- colon[colon$etype == 2, ]   

rec$time_yr   <- rec$time   / 365.25
death$time_yr <- death$time / 365.25

tau_G <- max(rec$time_yr)
cat(sprintf("tau_G (max observed follow-up) = %.3f years\n", tau_G))


tau_rec   <- tau_G + 2
tau_death <- tau_G + 10
cat(sprintf("tau (recurrence) = %.3f years, tau (death) = %.3f years\n",
            tau_rec, tau_death))

## Table 1: unconditional SFU test, both endpoints
cat("\n=== Table 1: unconditional sufficient-follow-up test ===\n")

fit_rec   <- sfu.test(rec$time_yr,   rec$status,   tau = tau_rec,   method = "sg")
fit_death <- sfu.test(death$time_yr, death$status, tau = tau_death, method = "sg")

tbl1 <- data.frame(
  endpoint = c("Recurrence", "Death (all-cause)"),
  p_value  = c(fit_rec$p.value, fit_death$p.value),
  conclusion = ifelse(c(fit_rec$p.value, fit_death$p.value) < 0.05,
                       "Sufficient follow-up", "Insufficient follow-up")
)
print(tbl1, row.names = FALSE)

cat("\n=== Table 2: categorical sufficient-follow-up test (Method 1) ===\n")

## (H0: insufficient for SOME x; H1: sufficient for ALL x).
run_categorical <- function(dat, cov, tau, n.boot = 1000) {
  res <- sfu.cov.test(dat$time_yr, dat$status, dat[[cov]], tau = tau, n.boot = n.boot)
  
  list(p_by_subgroup = res$p.value, x.unique = res$x.unique,
       method1_p = max(res$p.value),   
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

out_dir <- "phase1_baseline"
if (!dir.exists(out_dir)) out_dir <- "."
saveRDS(list(table1 = tbl1, table2 = tbl2), file.path(out_dir, "phase1_results.rds"))
cat("\nSaved phase1_results.rds\n")
