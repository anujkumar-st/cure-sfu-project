## Repository stage 03 simulation study: Part A
## Bias/MSE versus censoring horizon for the from-scratch EBVK estimator.

source("R/beran_ebvk.R")
if (!requireNamespace("KernSmooth", quietly = TRUE)) install.packages("KernSmooth")
suppressMessages(library(KernSmooth))

set.seed(20260807)

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
  Y <- ifelse(cured, Inf, U^(-gam_i))
  C <- ifelse(runif(n) < 1 - eps_fixed_frac, runif(n, 0, tau_c), tau_c)
  Tobs <- pmin(Y, C)
  delta <- as.numeric(Y <= C)
  list(Tobs = Tobs, delta = delta, X = X)
}

N_REPS_A <- 30
N_A <- 2000
X0_A <- 0.5
EPS <- 0.01
TRUE_P_A <- p_fun(X0_A)
TRUE_GAMMA_A <- gamma_fun(X0_A)
Q_99_A <- (1 / EPS)^TRUE_GAMMA_A
TAU_C_GRID <- c(4, 6, 8, 12, 20, 40, 80)

cat(sprintf("Part A: true p(x0=%.1f) = %.6f; gamma=%.3f; uncured 99%% quantile=%.6f\n", X0_A, TRUE_P_A, TRUE_GAMMA_A, Q_99_A))

partA_results <- do.call(rbind, lapply(TAU_C_GRID, function(tau_c) {
  cat(sprintf("  tau_c = %g ...\n", tau_c))
  reps <- t(sapply(seq_len(N_REPS_A), function(r) {
    dat <- simulate_one(N_A, tau_c)
    h <- dpik(dat$X)
    fit_one(dat$Tobs, dat$delta, dat$X, X0_A, h)
  }))
  data.frame(
    tau_c = tau_c,
    status = ifelse(tau_c < Q_99_A, "insufficient", "at_or_beyond_99pct"),
    mean_beran = mean(reps[, "beran"], na.rm = TRUE),
    mean_corrected = mean(reps[, "corrected"], na.rm = TRUE),
    bias_beran = mean(reps[, "beran"], na.rm = TRUE) - TRUE_P_A,
    bias_corrected = mean(reps[, "corrected"], na.rm = TRUE) - TRUE_P_A,
    mse_beran = mean((reps[, "beran"] - TRUE_P_A)^2, na.rm = TRUE),
    mse_corrected = mean((reps[, "corrected"] - TRUE_P_A)^2, na.rm = TRUE),
    n_valid_corrected = sum(!is.na(reps[, "corrected"])),
    correction_valid_rate = mean(!is.na(reps[, "corrected"]))
  )
}))

print(partA_results, row.names = FALSE)
write.csv(partA_results, "03_simulation_study/results/partA_results.csv", row.names = FALSE)
saveRDS(partA_results, "03_simulation_study/results/partA_results.rds")

png("03_simulation_study/results/partA_bias_vs_tauc.png", width = 900, height = 650)
plot(partA_results$tau_c, partA_results$bias_beran, type = "b",
     ylim = range(c(partA_results$bias_beran, partA_results$bias_corrected), na.rm = TRUE),
     xlab = expression(tau[c]), ylab = "Bias (estimate - true p(x0))",
     main = "Naive Beran vs EBVK-corrected bias")
lines(partA_results$tau_c, partA_results$bias_corrected, type = "b")
abline(h = 0, lty = 2)
abline(v = Q_99_A, lty = 3)
legend("topright", legend = c("Naive Beran", "EBVK corrected", "uncured 99% quantile"),
       lty = c(1, 1, 3), pch = c(1, 1, NA), bty = "n")
dev.off()
cat("Saved Part A outputs.\n")

cat("\nPart B is run by 02_method1_vs_method2.R because it is parallelized.\n")
