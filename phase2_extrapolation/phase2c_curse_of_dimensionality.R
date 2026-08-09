## ============================================================================
## Phase 2c: curse-of-dimensionality demonstration (report Sec. 5.4, Table 3).
##
## Extends the Beran estimator to a multivariate product-Epanechnikov
## kernel and adds real colon-recurrence covariates one at a time
## (age, then +nodes, +differ, +extent), evaluated at the componentwise-
## median "typical patient", with a FIXED per-dimension bandwidth reused
## unchanged as dimensions are added -- so degradation reflects
## dimensionality alone, not bandwidth retuning.
##
## Bandwidth rule: KernSmooth::dpik() for covariates it can handle (age,
## nodes); a Silverman rule-of-thumb fallback (1.06 * sd * n^(-1/5)) for
## near-degenerate ordinal covariates where dpik fails/degenerates
## (differ has 3 levels, extent has 4).
## ============================================================================

source(if (file.exists("R/beran_ebvk.R")) "R/beran_ebvk.R" else "../R/beran_ebvk.R")
need <- c("survival", "KernSmooth")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(survival); library(KernSmooth) })

data(colon)
rec <- colon[colon$etype == 1, ]
rec <- rec[complete.cases(rec[, c("nodes", "differ")]), ]   # n=888
rec$time_yr <- rec$time / 365.25
tau_n <- max(rec$time_yr)
cat(sprintf("n = %d (dropped %d rows with NA nodes/differ)\n",
            nrow(rec), sum(colon$etype == 1) - nrow(rec)))

## ---- componentwise-median "typical patient" evaluation point ----
x0 <- list(age = median(rec$age), nodes = median(rec$nodes),
           differ = median(rec$differ), extent = median(rec$extent))
cat("Evaluation point (median patient):\n"); str(x0)

## ---- fixed per-dimension bandwidths, computed once, reused unchanged ----
silverman <- function(x) 1.06 * sd(x) * length(x)^(-1 / 5)
h_age    <- dpik(rec$age)
h_nodes  <- dpik(rec$nodes)
h_differ <- silverman(rec$differ)   # dpik degenerates on 3-level ordinal data
h_extent <- silverman(rec$extent)   # dpik degenerates on 4-level ordinal data
hvec <- c(age = h_age, nodes = h_nodes, differ = h_differ, extent = h_extent)
cat("\nFixed bandwidths:\n"); print(hvec)

## ---- add covariates one at a time, same bandwidths throughout ----
combos <- list(c("age"), c("age", "nodes"), c("age", "nodes", "differ"),
               c("age", "nodes", "differ", "extent"))

table3 <- do.call(rbind, lapply(combos, function(combo) {
  Fn <- make_beran_md(rec$time_yr, rec$status, rec[, combo, drop = FALSE],
                       x0[combo], hvec[combo])
  fit <- fit_at_x0(Fn, tau_n)
  data.frame(covariates = paste(combo, collapse = " + "),
             ess = attr(Fn, "ess"), naive_p = fit$beran,
             corrected_p = fit$corrected, gamma_hat = fit$gamma)
}))

cat("\n=== Table 3: curse-of-dimensionality demonstration ===\n")
print(table3, row.names = FALSE, digits = 4)

cat(sprintf("\ngamma_hat hits the %.1f truncation floor at the final step: %s\n",
            GAMMA_FLOOR, isTRUE(all.equal(table3$gamma_hat[nrow(table3)], GAMMA_FLOOR))))

out_dir <- if (dir.exists("phase2_extrapolation")) "phase2_extrapolation" else "."
saveRDS(table3, file.path(out_dir, "phase2c_curse_of_dimensionality_results.rds"))
cat("\nSaved phase2c_curse_of_dimensionality_results.rds\n")
