## ============================================================================
## Phase 2g: gbex applied to the RECURRENCE endpoint with all ten available
## real covariates simultaneously (report Sec. 5.5, final two paragraphs).
##
## Explicitly a demonstration of covariate-scaling for the SCALE parameter
## alone (recurrence has adequate follow-up per Phase 1, so this is not a
## claim of solving the censoring/cure problem -- that remains the open
## methodological gap named in Future Work).
## ============================================================================

need <- c("survival")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("gbex", quietly = TRUE)) {
  stop("gbex not installed -- install from https://github.com/JVelthoen/gbex")
}
suppressMessages({ library(survival); library(gbex) })

data(colon)
rec <- colon[colon$etype == 1, ]
rec$time_yr <- rec$time / 365.25
covars <- c("age", "sex", "obstruct", "perfor", "adhere", "differ",
            "extent", "surg", "nodes", "node4")
rec <- rec[complete.cases(rec[, covars]), ]
cat(sprintf("n complete cases = %d\n", nrow(rec)))

## Low threshold (~1.1 years) chosen to land close to a usable uncensored
## exceedance count of roughly 220-230 (report: 223) while keeping enough
## signal for a 10-covariate boosted fit; exact threshold reconstructed
## from the target count, since the original threshold value wasn't
## recorded.
u <- 1.10
exc <- rec[rec$status == 1 & rec$time_yr > u, ]
excess <- exc$time_yr - u
cat(sprintf("threshold u = %.2f years, n uncensored exceedances = %d\n", u, nrow(exc)))

X10 <- as.data.frame(exc[, covars])
X_age_only <- as.data.frame(exc[, "age", drop = FALSE])

set.seed(2026)
fit_full <- gbex(y = excess, X = X10, silent = TRUE)
fit_age  <- gbex(y = excess, X = X_age_only, silent = TRUE)

pr_full <- predict(fit_full, newdata = X10)
pr_age  <- predict(fit_age,  newdata = X_age_only)

cat(sprintf("\nFull 10-covariate model: sigma range [%.3f, %.3f]\n",
            min(pr_full$s), max(pr_full$s)))
cat(sprintf("Shape (gamma) range: [%.4f, %.4f] (near zero throughout: %s)\n",
            min(pr_full$g), max(pr_full$g), all(abs(pr_full$g) < 0.05)))

## ---- deviance comparison: full 10-covariate vs age-only ----
gp_negloglik <- function(y, sigma, gamma) {
  z <- 1 + gamma * y / sigma
  if (any(z <= 0)) return(Inf)
  sum(log(sigma) + (1 + 1 / gamma) * log(z))
}
dev_full <- 2 * gp_negloglik(excess, pr_full$s, pr_full$g) / length(excess)
dev_age  <- 2 * gp_negloglik(excess, pr_age$s,  pr_age$g)  / length(excess)
cat(sprintf("\nMean deviance: age-only = %.3f, full 10-covariate = %.3f\n", dev_age, dev_full))
cat(sprintf("Full model improves on age-only: %s\n", dev_full < dev_age))

## ---- variable importance ----
## Note: gbex::variable_importance()'s default `type` argument
## (type = c("relative","permutation")) is never resolved with
## match.arg() internally, so the default call errors with
## "the condition has length > 1" -- a real bug in the published
## package, worked around here by passing type explicitly.
cat("\n=== Variable importance (full model) ===\n")
vi <- tryCatch(variable_importance(fit_full, type = "relative"), error = function(e) NULL)
if (!is.null(vi)) print(vi)

out_dir <- if (dir.exists("phase2_extrapolation")) "phase2_extrapolation" else "."
saveRDS(list(sigma_range = range(pr_full$s), gamma_range = range(pr_full$g),
             dev_full = dev_full, dev_age = dev_age, vi = vi, n = nrow(exc)),
        file.path(out_dir, "phase2g_gbex_full_covariate_results.rds"))
cat("\nSaved phase2g_gbex_full_covariate_results.rds\n")
