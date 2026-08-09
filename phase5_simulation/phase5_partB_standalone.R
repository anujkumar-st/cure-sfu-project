## Phase 5, Part B (standalone, parallelized): Method 1 (IUT) vs Method 2 (bootstrap x* selection) as the number of covariate categories K grows,
## with a FIXED total sample size split evenly across categories (so per-category n shrinks as K grows).

required_pkgs <- c("survival", "remotes", "future", "furrr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
if (!requireNamespace("cureSFUTest", quietly = TRUE)) {
  remotes::install_github("tp-yuen/cureSFUTest")
}
suppressMessages({
  library(survival)
  library(cureSFUTest)
  library(future)
  library(furrr)
})

plan(multisession, workers = availableCores())
cat(sprintf("Parallel plan initialized with %d workers.\n", availableCores()))

OUTPUT_DIR <- "phase5_outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE)


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


run_one_task <- function(K, scenario, rep) {
  all_suff <- identical(scenario, "null")
  set.seed(1000 * K + rep + (if (all_suff) 500000 else 0))
  dat <- simulate_categorical(K, N_TOTAL_B, HARD_TAU_C, EASY_TAU_C, all_suff)
  res <- tryCatch(
    cureSFUTest::sfu.cov.test(dat$Tobs, dat$delta, dat$X, tau = TAU_TEST_B,
                               n.boot = N_BOOT_B),
    error = function(e) NULL)
  if (is.null(res)) {
    return(data.frame(K = K, scenario = scenario, rep = rep,
                       method1_reject = NA, method2_correct = NA))
  }
  agg_p <- max(res$p.value)              
  m1 <- agg_p < 0.05                    
  m2 <- NA
  if (!all_suff) {
    
    m2 <- (res$x.ast.hat == which(res$x.unique == 1))
  }
  data.frame(K = K, scenario = scenario, rep = rep,
             method1_reject = m1, method2_correct = m2)
}


task_grid <- expand.grid(K = K_GRID, scenario = c("null", "alt"),
                          rep = seq_len(N_REPS_B), stringsAsFactors = FALSE)
cat(sprintf("Running %d tasks (%d K values x 2 scenarios x %d reps) across %d workers...\n",
            nrow(task_grid), length(K_GRID), N_REPS_B, availableCores()))

t0 <- Sys.time()
results_list <- furrr::future_pmap(
  task_grid,
  run_one_task,
  .options = furrr::furrr_options(seed = NULL),   
  .progress = TRUE
)
results_flat <- do.call(rbind, results_list)
cat(sprintf("\ndone in %.1f minutes\n", as.numeric(Sys.time() - t0, units = "mins")))

plan(sequential) 


partB_results <- do.call(rbind, lapply(K_GRID, function(K) {
  scen_null <- results_flat[results_flat$K == K & results_flat$scenario == "null", ]
  scen_alt  <- results_flat[results_flat$K == K & results_flat$scenario == "alt", ]
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


png(file.path(OUTPUT_DIR, "partB_level_power_vs_K.png"), width = 800, height = 600)
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
cat("saved plot:", file.path(OUTPUT_DIR, "partB_level_power_vs_K.png"), "\n")

saveRDS(partB_results, file.path(OUTPUT_DIR, "partB_results.rds"))
cat("saved results:", file.path(OUTPUT_DIR, "partB_results.rds"), "\n")
