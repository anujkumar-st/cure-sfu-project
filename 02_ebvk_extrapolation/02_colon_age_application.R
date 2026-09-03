## Repository stage 02.2: apply the EBVK extrapolation estimator to real
## survival::colon data, age as the (single) covariate, 9-point grid over
## ages 35-75.

source(if (file.exists("R/beran_ebvk.R")) "R/beran_ebvk.R" else "../R/beran_ebvk.R")
need <- c("survival", "KernSmooth")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(survival); library(KernSmooth) })

colon <- survival::colon
rec   <- colon[colon$etype == 1, ]
death <- colon[colon$etype == 2, ]
rec$time_yr   <- rec$time   / 365.25
death$time_yr <- death$time / 365.25

age_grid <- seq(35, 75, length.out = 9)


run_grid <- function(dat, grid) {
  h <- dpik(dat$age)
  out <- t(sapply(grid, function(x0) {
    f <- tryCatch(fit_one(dat$time_yr, dat$status, dat$age, x0, h),
                   error = function(e) list(beran = NA, corrected = NA, gamma = NA))
    c(beran = f$beran, corrected = f$corrected, gamma = f$gamma)
  }))
  data.frame(age = grid, out, gap = out[, "corrected"] - out[, "beran"])
}

cat("=== Recurrence: naive vs corrected cure-rate curve over age ===\n")
rec_curve <- run_grid(rec, age_grid)
print(rec_curve, row.names = FALSE)
cat(sprintf("Mean correction gap (recurrence): %.4f\n", mean(rec_curve$gap, na.rm = TRUE)))

cat("\n=== Death: naive vs corrected cure-rate curve over age ===\n")
death_curve <- run_grid(death, age_grid)
print(death_curve, row.names = FALSE)
cat(sprintf("Mean correction gap (death): %.4f\n", mean(death_curve$gap, na.rm = TRUE)))

cat(sprintf("\nCorrected >= naive everywhere (recurrence): %s\n",
            all(rec_curve$gap >= -1e-8, na.rm = TRUE)))
cat(sprintf("Corrected >= naive everywhere (death):      %s\n",
            all(death_curve$gap >= -1e-8, na.rm = TRUE)))


## Sex-stratified check against Phase 1's recurrence x sex finding
cat("\n=== Sex-stratified recurrence age-grid (check against Phase 1) ===\n")
rec_male   <- rec[rec$sex == 1, ]
rec_female <- rec[rec$sex == 0, ]
cat(sprintf("male n=%d, female n=%d\n", nrow(rec_male), nrow(rec_female)))

male_curve   <- run_grid(rec_male, age_grid)
female_curve <- run_grid(rec_female, age_grid)

cat("male:\n");   print(male_curve, row.names = FALSE)
cat("female:\n"); print(female_curve, row.names = FALSE)

cat(sprintf("\nMean correction gap: male=%.4f, female=%.4f\n",
            mean(male_curve$gap, na.rm = TRUE), mean(female_curve$gap, na.rm = TRUE)))
cat("(Report finding: this does NOT replicate Phase 1's direction -- female\n")
cat(" gap was larger than male's, opposite of the Phase 1 SFU-test result;\n")
cat(" attributed to sample thinning shrinking the local effective sample.)\n")

out_dir <- if (dir.exists("02_ebvk_extrapolation/results")) "02_ebvk_extrapolation/results" else "."
saveRDS(list(rec_curve = rec_curve, death_curve = death_curve,
             male_curve = male_curve, female_curve = female_curve),
        file.path(out_dir, "02_colon_age_application_results.rds"))
cat("\nSaved 02_colon_age_application_results.rds\n")
