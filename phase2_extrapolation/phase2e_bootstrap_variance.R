## Phase 2e: bootstrap variance check (report Sec. 5.4, final paragraph).
## B bootstrap resamples of gamma_hat/p_hat at 3 representative points (median patient, young/high-nodes, old/high-nodes) x 4 covariate
## dimensions (age, +nodes, +sex, +differ; death endpoint, n=888), to distinguish genuine dimensionality-driven instability (typical patient) from pre-existing small-sample fragility (atypical/rare covariate
## profiles) -- rather than trusting a single-sample point estimate.


source(if (file.exists("R/beran_ebvk.R")) "R/beran_ebvk.R" else "../R/beran_ebvk.R")
need <- c("survival", "KernSmooth")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(survival); library(KernSmooth) })

data(colon)
death <- colon[colon$etype == 2, ]
death <- death[complete.cases(death[, c("nodes", "differ")]), ]
death$time_yr <- death$time / 365.25
n <- nrow(death)

silverman <- function(x) 1.06 * sd(x) * length(x)^(-1 / 5)
hvec <- c(age = dpik(death$age), nodes = silverman(death$nodes),
          sex = 0.15, differ = silverman(death$differ))

med_sex <- median(death$sex); med_differ <- median(death$differ)
points <- list(
  median_patient   = list(age = median(death$age), nodes = median(death$nodes),
                           sex = med_sex, differ = med_differ),
  young_high_nodes = list(age = 55, nodes = 8, sex = med_sex, differ = med_differ),
  old_high_nodes   = list(age = 68, nodes = 8, sex = med_sex, differ = med_differ)
)
combos <- list(c("age"), c("age", "nodes"), c("age", "nodes", "sex"),
               c("age", "nodes", "sex", "differ"))

B <- 150
set.seed(2026)

run_bootstrap <- function(x0) {
  out <- vector("list", length(combos))
  for (d in seq_along(combos)) {
    combo <- combos[[d]]
    gam_boot <- rep(NA_real_, B)
    p_boot   <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      dboot <- death[idx, ]
      Fn <- tryCatch(
        make_beran_md(dboot$time_yr, dboot$status, dboot[, combo, drop = FALSE],
                       x0[combo], hvec[combo], categorical = intersect(combo, "sex")),
        error = function(e) NULL)
      if (is.null(Fn)) next
      fit <- fit_at_x0(Fn, max(dboot$time_yr))
      gam_boot[b] <- fit$gamma; p_boot[b] <- fit$corrected
    }
    out[[d]] <- data.frame(dim = d, covariates = paste(combo, collapse = "+"),
                            mean_gamma = mean(gam_boot, na.rm = TRUE),
                            sd_gamma = sd(gam_boot, na.rm = TRUE),
                            floor_pin_rate = mean(abs(gam_boot - GAMMA_FLOOR) < 1e-6, na.rm = TRUE),
                            n_valid = sum(!is.na(gam_boot)))
  }
  do.call(rbind, out)
}

results <- do.call(rbind, lapply(names(points), function(pname) {
  cat(sprintf("Bootstrapping at %s (B=%d) ...\n", pname, B))
  r <- run_bootstrap(points[[pname]])
  r$point <- pname
  r
}))
results <- results[, c("point", "dim", "covariates", "mean_gamma", "sd_gamma",
                        "floor_pin_rate", "n_valid")]

cat("\n=== Bootstrap variance results ===\n")
print(results, row.names = FALSE, digits = 4)

cat("\nFloor-pinning rate by point and dimension (fraction of B resamples):\n")
print(tapply(results$floor_pin_rate, list(results$point, results$dim), identity))

out_dir <- if (dir.exists("phase2_extrapolation")) "phase2_extrapolation" else "."
saveRDS(results, file.path(out_dir, "phase2e_bootstrap_variance_results.rds"))
cat("\nSaved phase2e_bootstrap_variance_results.rds\n")
