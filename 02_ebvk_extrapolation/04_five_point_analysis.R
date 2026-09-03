## Phase 2d: companion curse-of-dimensionality analysis on the DEATH endpoint (report Sec. 5.4, paragraph after Table 3), with a different
## covariate order (age, +nodes, +sex, +differ) and checked across FIVE evaluation points instead of just the median patient, to test whether
## the single-point "gap shrinks monotonically with dimension" story generalizes.


source(if (file.exists("R/beran_ebvk.R")) "R/beran_ebvk.R" else "../R/beran_ebvk.R")
need <- c("survival", "KernSmooth")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(survival); library(KernSmooth) })

colon <- survival::colon
death <- colon[colon$etype == 2, ]
death <- death[complete.cases(death[, c("nodes", "differ")]), ]  
death$time_yr <- death$time / 365.25
tau_n <- max(death$time_yr)
cat(sprintf("n = %d\n", nrow(death)))


silverman <- function(x) 1.06 * sd(x) * length(x)^(-1 / 5)
hvec <- c(age    = dpik(death$age),
          nodes  = silverman(death$nodes),
          sex    = 0.15,
          differ = silverman(death$differ))
cat("Fixed bandwidths:\n"); print(hvec)


med_sex <- median(death$sex); med_differ <- median(death$differ)
points <- list(
  median_patient  = list(age = median(death$age), nodes = median(death$nodes),
                          sex = med_sex, differ = med_differ),
  young_low_nodes  = list(age = 55, nodes = 2, sex = med_sex, differ = med_differ),
  young_high_nodes = list(age = 55, nodes = 8, sex = med_sex, differ = med_differ),
  old_low_nodes    = list(age = 68, nodes = 2, sex = med_sex, differ = med_differ),
  old_high_nodes   = list(age = 68, nodes = 8, sex = med_sex, differ = med_differ)
)

combos <- list(c("age"), c("age", "nodes"), c("age", "nodes", "sex"),
               c("age", "nodes", "sex", "differ"))

results <- do.call(rbind, lapply(names(points), function(pname) {
  x0 <- points[[pname]]
  do.call(rbind, lapply(seq_along(combos), function(d) {
    combo <- combos[[d]]
    Fn <- tryCatch(
      make_beran_md(death$time_yr, death$status, death[, combo, drop = FALSE],
                     x0[combo], hvec[combo],
                     categorical = intersect(combo, "sex")),
      error = function(e) NULL)
    if (is.null(Fn)) {
      return(data.frame(point = pname, dim = d, covariates = paste(combo, collapse = "+"),
                         ess = NA, naive_p = NA, corrected_p = NA, gamma_hat = NA, gap = NA))
    }
    fit <- fit_at_x0(Fn, tau_n)
    data.frame(point = pname, dim = d, covariates = paste(combo, collapse = "+"),
               ess = attr(Fn, "ess"), naive_p = fit$beran, corrected_p = fit$corrected,
               gamma_hat = fit$gamma, gap = fit$corrected - fit$beran)
  }))
}))

cat("\n=== 5-point x 4-dimension results ===\n")
print(results, row.names = FALSE, digits = 4)

cat("\n=== Mean correction gap by dimension (across all 5 points) ===\n")
print(tapply(results$gap, results$dim, mean, na.rm = TRUE))

cat("\n=== gamma_hat floor-pinning count by dimension (out of 5 points) ===\n")
print(tapply(results$gamma_hat, results$dim, function(g) sum(abs(g - GAMMA_FLOOR) < 1e-6, na.rm = TRUE)))

cat("\n=== ESS for the two rare/atypical high-nodes points ===\n")
print(results[results$point %in% c("young_high_nodes", "old_high_nodes"), c("point", "dim", "ess")])

out_dir <- if (dir.exists("02_ebvk_extrapolation/results")) "02_ebvk_extrapolation/results" else "."
saveRDS(results, file.path(out_dir, "04_five_point_results.rds"))
cat("\nSaved 04_five_point_results.rds\n")
