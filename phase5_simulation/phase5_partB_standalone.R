## ============================================================================
## Phase 5, Part B (standalone, parallelized): Method 1 (IUT) vs Method 2
## (bootstrap x* selection) as the number of covariate categories K grows,
## with a FIXED total sample size split evenly across categories (so
## per-category n shrinks as K grows).
##
## This is self-contained -- run it on its own, no dependency on Part A.
## Parallelized with furrr/future across all (K, scenario, rep) combinations
## at once, so it scales with however many cores are available.
## ============================================================================

## ---- setup ------------------------------------------------------------
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

## ---- synthetic categorical-covariate cure data generator ----
## K categories, fixed total n split evenly. Category 1 is deliberately the
## "hard" one (short censoring horizon -> insufficient follow-up); categories
## 2..K are comfortably sufficient. all_sufficient=TRUE overrides category 1
## to also be easy, for the Type-I error / level check.
simulate_categorical <- function(K, n_total, hard_tau_c, easy_tau_c,
                                  all_sufficient = FALSE) {
  n_per <- floor(n_total / K)
  gam <- 0.9          # fixed tail shape, same across categories (keep it simple)
  p_true <- 0.55       # fixed cure probability across categories
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

## ---- parameters (tune these first if it's too slow) ----
N_REPS_B <- 30           # reduce to ~10 for a quick smoke test
N_BOOT_B <- 200          # reduce to ~100 for speed, increase for precision
K_GRID <- 2:5
N_TOTAL_B <- 400         # fixed total n, split evenly across K categories
HARD_TAU_C <- 3          # short horizon -> category 1 insufficient
EASY_TAU_C <- 30         # long horizon -> comfortably sufficient
TAU_TEST_B <- 60         # test tau, must exceed max(Tobs) always

## ---- one (K, scenario, rep) task ----
## scenario "null" = all categories sufficient (ground truth H1 -> power)
## scenario "alt"  = category 1 insufficient   (ground truth H0 -> Type-I error / selection accuracy)
## Each task sets its own explicit seed (same formula as the original
## sequential loop), so results are bit-for-bit reproducible regardless of
## how many workers run it or in what order -- no reliance on furrr's
## automatic parallel-safe RNG streams.
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
  agg_p <- max(res$p.value)               # intersection-union aggregate p-value
  m1 <- agg_p < 0.05                      # reject H0 -> conclude sufficient for ALL x
  m2 <- NA
  if (!all_suff) {
    # true worst category is category 1 (x.unique is sorted, so its index
    # is where the X value == 1)
    m2 <- (res$x.ast.hat == which(res$x.unique == 1))
  }
  data.frame(K = K, scenario = scenario, rep = rep,
             method1_reject = m1, method2_correct = m2)
}

## ---- build the full task grid and run it in parallel ----
task_grid <- expand.grid(K = K_GRID, scenario = c("null", "alt"),
                          rep = seq_len(N_REPS_B), stringsAsFactors = FALSE)
cat(sprintf("Running %d tasks (%d K values x 2 scenarios x %d reps) across %d workers...\n",
            nrow(task_grid), length(K_GRID), N_REPS_B, availableCores()))

t0 <- Sys.time()
results_list <- furrr::future_pmap(
  task_grid,
  run_one_task,
  .options = furrr::furrr_options(seed = NULL),   # we seed manually per task, above --
                                                    # NULL disables furrr's own RNG-safety
                                                    # check/warning since we handle it ourselves
  .progress = TRUE
)
results_flat <- do.call(rbind, results_list)
cat(sprintf("\ndone in %.1f minutes\n", as.numeric(Sys.time() - t0, units = "mins")))

plan(sequential)  # clean up workers

## ---- aggregate per K (same corrected formulas as the sequential version) ----
partB_results <- do.call(rbind, lapply(K_GRID, function(K) {
  scen_null <- results_flat[results_flat$K == K & results_flat$scenario == "null", ]
  scen_alt  <- results_flat[results_flat$K == K & results_flat$scenario == "alt", ]
  data.frame(
    K = K,
    n_per_category = floor(N_TOTAL_B / K),
    # Power to confirm true full sufficiency: scen_null is the "all
    # categories sufficient" generative scenario, where the ground truth is
    # H1. Correctly detecting it means method1_reject == TRUE (IUT rejects
    # H0 and concludes "sufficient for all"), so power is the *direct*
    # rejection rate here, not its complement. Expect this to fall as K
    # grows, since all K per-category tests must clear 5% simultaneously.
    method1_power = mean(scen_null$method1_reject, na.rm = TRUE),
    # Type I error (false "all clear"): scen_alt is the "one category
    # insufficient" generative scenario, where the ground truth is H0.
    # Wrongly concluding "sufficient for all" means method1_reject == TRUE,
    # so this is also the *direct* rejection rate, not its complement.
    # Expect this to stay at or near 0 (IUTs are conservative).
    method1_type1_error = mean(scen_alt$method1_reject, na.rm = TRUE),
    # Method 2 selection accuracy: how often x.ast.hat correctly flags the
    # true worst category
    method2_selection_accuracy = mean(scen_alt$method2_correct, na.rm = TRUE)
  )
}))

cat("\n=== Part B results: Method 1 / Method 2 behavior as K grows ===\n")
print(partB_results, row.names = FALSE)

## ---- plot + save ----
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
