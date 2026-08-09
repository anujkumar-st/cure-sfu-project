## Phase 5 simulation study
## Part A: sweep the censoring horizon tau_c from severely insufficient to close-to-sufficient follow-up, and show the naive Beran bias vs. the  EBVK-corrected estimator bias as a function of insufficiency. (Uses the from-scratch estimator validated in Phase 2 -- all
##         functions are reproduced here so this script is self-contained.)

## Part B: Method 1 (intersection-union test) vs Method 2 (bootstrap worst-category selection) from Yuen/Musta/Van Keilegom's categorical-covariate SFU test, as the number of covariate
##         categories K grows with a FIXED total sample size (so per-category n shrinks as K grows—mirroring the real-data sex-stratification noise we saw in Phase 2).
##           - Level (Type-I error) of Method 1 under a scenario where every category truly has sufficient follow-up.
##           - Power of Method 1 to detect that at least one category has insufficient follow-up.
##           - Selection accuracy of Method 2's x.ast.hat: how often it correctly identifies the true worst category.

need <- c("KernSmooth", "survival")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("cureSFUTest", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("tp-yuen/cureSFUTest")
}
suppressMessages({
  library(KernSmooth)
  library(survival)
  library(cureSFUTest)
})

set.seed(20260807)


## PART A: tau_c sweep for the from-scratch EBVK extrapolation estimator

epan <- function(u) ifelse(abs(u) < 1, 0.75 * (1 - u^2), 0)

make_beran <- function(T, delta, X, x0, h) {
  w_raw <- epan((x0 - X) / h)
  W <- w_raw / sum(w_raw)
  ord <- order(T)
  Ts <- T[ord]; ds <- delta[ord]; Ws <- W[ord]
  cumWprev <- c(0, cumsum(Ws)[-length(Ws)])
  denom <- 1 - cumWprev
  denom[denom <= 1e-12] <- 1e-12
  factor <- 1 - Ws / denom
  factor[factor < 0] <- 0
  logfac <- ifelse(ds == 1, log(factor), 0)
  cumlog <- cumsum(logfac)
  function(t) {
    idx <- which(Ts <= t)
    if (length(idx) == 0) return(0)
    1 - exp(cumlog[max(idx)])
  }
}

gamma_hat_fn <- function(Fn, tau_n, y2) {
  num <- Fn(y2^2 * tau_n) - Fn(y2 * tau_n)
  den <- Fn(y2 * tau_n) - Fn(tau_n)
  ratio <- num / den
  if (!is.finite(ratio) || ratio <= 0) return(NA_real_)
  -1 / (log(ratio) / log(y2))
}

p_hat_fn <- function(Fn, tau_n, y1, gam) {
  if (is.na(gam)) return(NA_real_)
  Fn_tau <- Fn(tau_n)
  Fn_tau + (Fn_tau - Fn(y1 * tau_n)) / (y1^(-1/gam) - 1)
}

G_GRID <- seq(0.25, 0.89, by = 0.02)

fit_one <- function(Tobs, delta, X, x0, h) {
  Fn <- make_beran(Tobs, delta, X, x0, h)
  tau_n <- max(Tobs)
  beran_p <- Fn(tau_n)
  gam_grid <- pmax(sapply(G_GRID, function(y2) gamma_hat_fn(Fn, tau_n, y2)), 0.1)
  pmat <- matrix(NA_real_, length(G_GRID), length(G_GRID))
  for (j in seq_along(G_GRID)) for (i in seq_along(G_GRID))
    pmat[i, j] <- p_hat_fn(Fn, tau_n, G_GRID[i], gam_grid[j])
  valid <- is.finite(pmat) & pmat >= beran_p - 1e-8 & pmat <= 1
  if (!any(valid)) return(c(beran = beran_p, corrected = NA_real_))
  mean_p <- mean(pmat[valid])
  d <- abs(pmat - mean_p); d[!valid] <- Inf
  best <- which(d == min(d), arr.ind = TRUE)[1, ]
  c(beran = beran_p, corrected = pmat[best[1], best[2]])
}


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
TAU_C_GRID <- c(4, 6, 8, 12, 20, 40)
N_A <- 2000
X0_A <- 0.5
TRUE_P_A <- p_fun(X0_A)

cat(sprintf("Part A: true p(x0=%.1f) = %.4f\n", X0_A, TRUE_P_A))

partA_results <- do.call(rbind, lapply(TAU_C_GRID, function(tau_c) {
  cat(sprintf("  tau_c = %g ...\n", tau_c))
  reps <- t(sapply(seq_len(N_REPS_A), function(r) {
    dat <- simulate_one(N_A, tau_c)
    h <- dpik(dat$X)
    fit_one(dat$Tobs, dat$delta, dat$X, X0_A, h)
  }))
  data.frame(
    tau_c = tau_c,
    mean_beran = mean(reps[, "beran"], na.rm = TRUE),
    mean_corrected = mean(reps[, "corrected"], na.rm = TRUE),
    bias_beran = mean(reps[, "beran"], na.rm = TRUE) - TRUE_P_A,
    bias_corrected = mean(reps[, "corrected"], na.rm = TRUE) - TRUE_P_A,
    mse_beran = mean((reps[, "beran"] - TRUE_P_A)^2, na.rm = TRUE),
    mse_corrected = mean((reps[, "corrected"] - TRUE_P_A)^2, na.rm = TRUE),
    n_valid_corrected = sum(!is.na(reps[, "corrected"]))
  )
}))

cat("\n=== Part A results: bias/MSE vs tau_c (severity of insufficient follow-up) ===\n")
print(partA_results, row.names = FALSE)

png("/mnt/user-data/outputs/partA_bias_vs_tauc.png", width = 800, height = 600)
plot(partA_results$tau_c, partA_results$bias_beran, type = "b", col = "firebrick",
     ylim = range(c(partA_results$bias_beran, partA_results$bias_corrected), na.rm = TRUE),
     xlab = "tau_c (censoring horizon)", ylab = "bias (mean estimate - true p(x0))",
     main = "Naive Beran vs EBVK-corrected bias as follow-up lengthens")
lines(partA_results$tau_c, partA_results$bias_corrected, type = "b", col = "steelblue")
abline(h = 0, lty = 2, col = "gray50")
legend("topright", legend = c("naive Beran", "EBVK corrected"),
       col = c("firebrick", "steelblue"), lty = 1, pch = 1)
dev.off()
cat("saved plot: partA_bias_vs_tauc.png\n")

saveRDS(partA_results, "/mnt/user-data/outputs/partA_results.rds")

## PART B: Method 1 (IUT) vs Method 2 (bootstrap x* selection) as K grows

simulate_categorical <- function(K, n_total, hard_tau_c, easy_tau_c,
                                  all_sufficient = FALSE) {
  n_per <- floor(n_total / K)
  gam <- 0.9         
  p_true <- 0.55      
  Y_list <- list(); C_list <- list(); X_list <- list()
  for (k in seq_len(K)) {
    tau_c_k <- if (all_sufficient || k != 1) easy_tau_c else hard_tau_c
    cured <- rbinom(n_per, 1, 1 - p_true) == 1
    U <- runif(n_per)
    Y <- ifelse(cured, Inf, U^(-gam))
    C <- ifelse(runif(n_per) < 0.9, runif(n_per, 0, tau_c_k), tau_c_k)
    Y_list[[k]] <- Y; C_list[[k]] <- C; X_list[[k]] <- rep(k, n_per)
  }
  Y <- unlist(Y_list); C <- unlist(C_list); X <- unlist(X_list)
  Tobs <- pmin(Y, C)
  delta <- as.numeric(Y <= C)
  list(Tobs = Tobs, delta = delta, X = X)
}

N_REPS_B <- 30          
N_BOOT_B <- 200         
K_GRID <- 2:5
N_TOTAL_B <- 400         
HARD_TAU_C <- 3          
EASY_TAU_C <- 30        
TAU_TEST_B <- 60       

run_scenario <- function(K, all_sufficient) {
  method1_reject <- logical(N_REPS_B)   
  method2_correct <- logical(N_REPS_B) 
  for (r in seq_len(N_REPS_B)) {
    set.seed(1000 * K + r + (if (all_sufficient) 500000 else 0))
    dat <- simulate_categorical(K, N_TOTAL_B, HARD_TAU_C, EASY_TAU_C, all_sufficient)
    res <- tryCatch(
      cureSFUTest::sfu.cov.test(dat$Tobs, dat$delta, dat$X, tau = TAU_TEST_B,
                                 n.boot = N_BOOT_B),
      error = function(e) NULL)
    if (is.null(res)) { method1_reject[r] <- NA; method2_correct[r] <- NA; next }
    agg_p <- max(res$p.value)
    method1_reject[r] <- agg_p < 0.05  
    if (!all_sufficient) {
    
      method2_correct[r] <- (res$x.ast.hat == which(res$x.unique == 1))
    }
  }
  list(method1_reject = method1_reject, method2_correct = method2_correct)
}

partB_results <- do.call(rbind, lapply(K_GRID, function(K) {
  cat(sprintf("K = %d ...\n", K))
  cat("  scenario: all categories sufficient (ground truth H1 -> power)\n")
  scen_null <- run_scenario(K, all_sufficient = TRUE)
  cat("  scenario: category 1 insufficient (ground truth H0 -> Type-I error / selection accuracy)\n")
  scen_alt <- run_scenario(K, all_sufficient = FALSE)

  data.frame(
    K = K,
    n_per_category = floor(N_TOTAL_B / K),
    
    method1_power = mean(scen_null$method1_reject, na.rm = TRUE),
    
    method1_type1_error = mean(scen_alt$method1_reject, na.rm = TRUE),
   
    method2_selection_accuracy = mean(scen_alt$method2_correct, na.rm = TRUE)
  )
}))

cat("\n=== Part B results: Method 1 / Method 2 behavior as K grows ===\n")
print(partB_results, row.names = FALSE)

png("/mnt/user-data/outputs/partB_level_power_vs_K.png", width = 800, height = 600)
plot(partB_results$K, partB_results$method1_type1_error, type = "b", col = "firebrick",
     ylim = c(0, 1), xlab = "number of covariate categories (K)", ylab = "probability",
     main = "Method 1 level & power, Method 2 selection accuracy vs K")
lines(partB_results$K, partB_results$method1_power, type = "b", col = "steelblue")
lines(partB_results$K, partB_results$method2_selection_accuracy, type = "b", col = "darkgreen")
abline(h = 0.05, lty = 2, col = "gray50")
legend("bottomleft", legend = c("Method 1 Type-I error (nominal alpha=0.05 dashed)",
                                 "Method 1 power", "Method 2 selection accuracy"),
       col = c("firebrick", "steelblue", "darkgreen"), lty = 1, pch = 1, cex = 0.8)
dev.off()
cat("saved plot: partB_level_power_vs_K.png\n")

saveRDS(partB_results, "/mnt/user-data/outputs/partB_results.rds")

cat("\nDone. Results and plots saved to /mnt/user-data/outputs/\n")
cat("(if running outside the sandboxed environment, change the two output\n")
cat(" paths above to a local directory of your choice)\n")
