## Phase 2a: synthetic sanity check of the from-scratch EBVK extrapolation
## estimator against known ground truth (report Sec. 5.2).

source(if (file.exists("R/beran_ebvk.R")) "R/beran_ebvk.R" else "../R/beran_ebvk.R")
if (!requireNamespace("KernSmooth", quietly = TRUE)) install.packages("KernSmooth")
suppressMessages(library(KernSmooth))

set.seed(2026)


gamma_fun <- function(x) (x + 1) / 2
p_fun <- function(x, b1 = 0.4, b2 = 2) {
  lp <- b1 + b2 * (2 * x - 1)
  exp(lp) / (1 + exp(lp))
}

simulate_one <- function(n, tau_c, eps_fixed_frac = 0.1) {
  X <- runif(n, 0, 1)
  gam_i <- gamma_fun(X); p_i <- p_fun(X)
  cured <- rbinom(n, 1, 1 - p_i) == 1
  U <- runif(n)
  Y <- ifelse(cured, Inf, U^(-gam_i))   # Pareto-tailed uncured survival
  C <- ifelse(runif(n) < 1 - eps_fixed_frac, runif(n, 0, tau_c), tau_c)
  Tobs <- pmin(Y, C)
  delta <- as.numeric(Y <= C)
  list(Tobs = Tobs, delta = delta, X = X)
}

## Single-run sanity check at x0 = 0.5
x0 <- 0.5
true_gamma <- gamma_fun(x0)
true_p <- p_fun(x0)
cat(sprintf("Ground truth at x0=%.1f: gamma=%.4f, p=%.4f\n", x0, true_gamma, true_p))

n <- 2000
tau_c_short <- 6 
dat <- simulate_one(n, tau_c_short)
h <- dpik(dat$X)
fit <- fit_one(dat$Tobs, dat$delta, dat$X, x0, h)

cat(sprintf("Naive Beran p_hat    = %.4f (underestimates, as predicted)\n", fit$beran))
cat(sprintf("EBVK corrected p_hat = %.4f (gamma_hat = %.4f)\n", fit$corrected, fit$gamma))
cat(sprintf("Corrected moves toward truth: %s\n",
            abs(fit$corrected - true_p) < abs(fit$beran - true_p)))

## Repeated-sample check: bias/MSE of naive vs corrected over R reps
R <- 20
reps <- t(sapply(seq_len(R), function(r) {
  d <- simulate_one(n, tau_c_short)
  hh <- dpik(d$X)
  f <- fit_one(d$Tobs, d$delta, d$X, x0, hh)
  c(beran = f$beran, corrected = f$corrected)
}))

cat(sprintf("\n%d-rep check at x0=%.1f (true p=%.4f):\n", R, x0, true_p))
cat(sprintf("  naive Beran:     mean p_hat=%.4f (bias %+.4f)\n",
            mean(reps[, "beran"], na.rm = TRUE), mean(reps[, "beran"], na.rm = TRUE) - true_p))
cat(sprintf("  EBVK corrected:  mean p_hat=%.4f (bias %+.4f)\n",
            mean(reps[, "corrected"], na.rm = TRUE), mean(reps[, "corrected"], na.rm = TRUE) - true_p))
cat(sprintf("  MSE naive=%.4f, MSE corrected=%.4f\n",
            mean((reps[, "beran"] - true_p)^2, na.rm = TRUE),
            mean((reps[, "corrected"] - true_p)^2, na.rm = TRUE)))

out_dir <- if (dir.exists("phase2_extrapolation")) "phase2_extrapolation" else "."
saveRDS(list(single_run = fit, reps = reps, true_p = true_p, true_gamma = true_gamma),
        file.path(out_dir, "phase2a_synthetic_validation_results.rds"))
cat("\nSaved phase2a_synthetic_validation_results.rds\n")
