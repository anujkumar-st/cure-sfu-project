## ============================================================================
## Phase 2f: gbex comparison under matched conditions (report Sec. 5.5).
##
## gbex (Velthoen et al. 2021) has no notion of censoring, so it can only
## be fit on uncensored exceedances above a threshold. At the death
## endpoint's 50th-percentile threshold there are exactly 33 uncensored
## exceedances -- matched to the same evaluation points / covariate
## dimensions (d1-d4: age, +nodes, +sex, +differ) used for the Beran
## comparison in Phase 2d, to test whether gbex resolves the
## dimensionality-driven failure mode found there.
## ============================================================================

need <- c("survival")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("gbex", quietly = TRUE)) {
  stop("gbex not installed -- install from https://github.com/JVelthoen/gbex ",
       "(source install; also requires POT and treeClust)")
}
suppressMessages({ library(survival); library(gbex) })

data(colon)
death <- colon[colon$etype == 2, ]
death <- death[complete.cases(death[, c("nodes", "differ")]), ]
death$time_yr <- death$time / 365.25

u <- quantile(death$time_yr, 0.5)
exc <- death[death$status == 1 & death$time_yr > u, ]
cat(sprintf("Threshold u (50th pct) = %.3f years, n uncensored exceedances = %d\n",
            u, nrow(exc)))
excess <- exc$time_yr - u   # gbex has no built-in peaks-over-threshold step

med_sex <- median(death$sex); med_differ <- median(death$differ)
points <- list(
  median_patient   = list(age = median(death$age), nodes = median(death$nodes),
                           sex = med_sex, differ = med_differ),
  young_high_nodes = list(age = 55, nodes = 8, sex = med_sex, differ = med_differ),
  old_high_nodes   = list(age = 68, nodes = 8, sex = med_sex, differ = med_differ)
)
combos <- list(c("age"), c("age", "nodes"), c("age", "nodes", "sex"),
               c("age", "nodes", "sex", "differ"))

fit_gbex_at <- function(combo, x0) {
  Xdf <- as.data.frame(exc[, combo, drop = FALSE])   # real data.frame, real column names
  fit <- gbex(y = excess, X = Xdf, silent = TRUE)
  newdata <- as.data.frame(x0[combo])
  names(newdata) <- combo
  pr <- predict(fit, newdata = newdata)
  list(sigma = pr$s, gamma = pr$g, fit = fit)
}

cat("\n=== gbex point estimates across points x dimensions ===\n")
results <- do.call(rbind, lapply(names(points), function(pname) {
  do.call(rbind, lapply(seq_along(combos), function(d) {
    combo <- combos[[d]]
    r <- tryCatch(fit_gbex_at(combo, points[[pname]]),
                  error = function(e) list(sigma = NA, gamma = NA))
    data.frame(point = pname, dim = d, covariates = paste(combo, collapse = "+"),
               sigma = r$sigma, gamma = r$gamma)
  }))
}))
print(results, row.names = FALSE, digits = 4)

## ---- root cause: check the unconditional first-guess gamma ----
cat("\n=== Root cause: gbex's internal unconditional MLE first-guess ===\n")
if (!requireNamespace("POT", quietly = TRUE)) install.packages("POT")
uncond_fit <- POT::fitgpd(excess, threshold = 0)
cat(sprintf("Unconditional GPD MLE on the 33 exceedances: gamma_hat = %.5f\n",
            uncond_fit$fitted.values["shape"]))
cat("If negative, gbex::first_guess() silently falls back to a fixed\n")
cat("default initial gamma of 0.01, from which boosting shrinkage is tiny.\n")

## ---- 30-rep bootstrap: is the flat estimate genuine stability or a frozen init? ----
cat("\n=== 30-rep bootstrap of gamma_hat at age (d=1), evaluated at median age ===\n")
set.seed(2026)
B <- 30
gam_boot <- rep(NA_real_, B)
n_exc <- nrow(exc)
for (b in seq_len(B)) {
  idx <- sample.int(n_exc, n_exc, replace = TRUE)
  eb <- excess[idx]; Xb <- as.data.frame(exc$age[idx]); names(Xb) <- "age"
  fit_b <- tryCatch(gbex(y = eb, X = Xb, silent = TRUE), error = function(e) NULL)
  if (is.null(fit_b)) next
  pr <- predict(fit_b, newdata = data.frame(age = median(death$age)))
  gam_boot[b] <- pr$g
}
cat(sprintf("mean=%.5f, sd=%.5f, range=[%.5f, %.5f] (n_valid=%d/%d)\n",
            mean(gam_boot, na.rm = TRUE), sd(gam_boot, na.rm = TRUE),
            min(gam_boot, na.rm = TRUE), max(gam_boot, na.rm = TRUE),
            sum(!is.na(gam_boot)), B))
cat("A tiny SD here indicates the estimate is FROZEN under resampling\n")
cat("(an initialization artifact), not genuinely stable/converged.\n")

out_dir <- if (dir.exists("phase2_extrapolation")) "phase2_extrapolation" else "."
saveRDS(list(results = results, uncond_gamma = uncond_fit$fitted.values["shape"],
             bootstrap = gam_boot),
        file.path(out_dir, "phase2f_gbex_comparison_results.rds"))
cat("\nSaved phase2f_gbex_comparison_results.rds\n")
