## Uses purrr + furrr and can be run from ANY working directory.
##
## IMPORTANT HYPOTHESES / INTERPRETATION
##
## The categorical cureSFUTest procedure is an intersection-union test:
##
##   H0: q_(1-eps)(x) >= tau_G(x) for SOME category x
##       ("insufficient follow-up for at least one category")
##
##   H1: q_(1-eps)(x) <  tau_G(x) for ALL categories x
##       ("sufficient follow-up for all categories")
##
## Therefore:
##   all_sufficient scenario:
##       H1 is true -> rejection of H0 = POWER
##
##   one_insufficient scenario:
##       H0 is true -> rejection of H0 = false "all clear" = TYPE-I ERROR
##
## This mirrors the interpretation in the Project report.


options(stringsAsFactors = FALSE)

required <- c("survival", "purrr", "furrr", "future")

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("cureSFUTest", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("tp-yuen/cureSFUTest")
}

suppressMessages({
  library(survival)
  library(cureSFUTest)
  library(purrr)
  library(furrr)
  library(future)
})



simulate_categorical <- function(K, n_total, hard_tau_c, easy_tau_c,
                                  all_sufficient = FALSE) {

  n_per <- floor(n_total / K)

  gamma <- 0.9
  p_true <- 0.55

  Y_list <- vector("list", K)
  C_list <- vector("list", K)
  X_list <- vector("list", K)

  for (k in seq_len(K)) {

    tau_c_k <- if (all_sufficient || k != 1) {
      easy_tau_c
    } else {
      hard_tau_c
    }

    cured <- rbinom(
      n_per,
      size = 1,
      prob = 1 - p_true
    ) == 1

    U <- runif(n_per)

    Y <- ifelse(
      cured,
      Inf,
      U^(-gamma)
    )

    C <- ifelse(
      runif(n_per) < 0.9,
      runif(n_per, 0, tau_c_k),
      tau_c_k
    )

    Y_list[[k]] <- Y
    C_list[[k]] <- C
    X_list[[k]] <- rep(k, n_per)
  }

  Y <- unlist(Y_list)
  C <- unlist(C_list)
  X <- unlist(X_list)

  Tobs <- pmin(Y, C)
  delta <- as.numeric(Y <= C)

  list(
    Tobs = Tobs,
    delta = delta,
    X = X
  )
}


N_REPS_B <- 30
N_BOOT_B <- 200

K_GRID <- 2:5
N_TOTAL_B <- 400

HARD_TAU_C <- 3
EASY_TAU_C <- 80
TAU_TEST_B <- 100

EPS <- 0.01
GAMMA <- 0.9

Q_99 <- (1 / EPS)^GAMMA

N_WORKERS <- max(
  1L,
  future::availableCores() - 1L
)

plan(
  multisession,
  workers = N_WORKERS
)

cat("============================================================\n")
cat("PHASE 5B: PARALLEL METHOD 1 vs METHOD 2\n")
cat("============================================================\n\n")

cat(sprintf(
  "gamma = %.3f\n",
  GAMMA
))

cat(sprintf(
  "uncured 99%% quantile = %.6f\n",
  Q_99
))

cat(sprintf(
  "hard tau_c = %.3f\n",
  HARD_TAU_C
))

cat(sprintf(
  "easy tau_c = %.3f\n",
  EASY_TAU_C
))

cat(sprintf(
  "test tau = %.3f\n",
  TAU_TEST_B
))

cat(sprintf(
  "N_TOTAL = %d\n",
  N_TOTAL_B
))

cat(sprintf(
  "replications = %d\n",
  N_REPS_B
))

cat(sprintf(
  "bootstrap resamples = %d\n",
  N_BOOT_B
))

cat(sprintf(
  "parallel workers = %d\n\n",
  N_WORKERS
))


run_one_task <- function(K, scenario, rep) {

  all_sufficient <- identical(
    scenario,
    "all_sufficient"
  )

  set.seed(
    1000 * K +
      rep +
      if (all_sufficient) 500000 else 0
  )

  dat <- simulate_categorical(
    K = K,
    n_total = N_TOTAL_B,
    hard_tau_c = HARD_TAU_C,
    easy_tau_c = EASY_TAU_C,
    all_sufficient = all_sufficient
  )

 
  fit <- tryCatch(
    cureSFUTest::sfu.cov.test(
      dat$Tobs,
      dat$delta,
      dat$X,
      tau = TAU_TEST_B,
      n.boot = N_BOOT_B
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {

    return(
      data.frame(
        K = K,
        scenario = scenario,
        rep = rep,
        method1_reject = NA,
        method2_correct = NA,
        fit_failed = TRUE,
        error_message = conditionMessage(fit)
      )
    )
  }

 ## max subgroup p-value < alpha means ALL subgroup tests reject H0, hence the global H0 ("insufficient for some category") is rejected.
  agg_p <- max(
    as.numeric(fit$p.value),
    na.rm = TRUE
  )

  method1_reject <- agg_p < 0.05

  ## Method 2 is meaningful in the one-insufficient scenario: category 1 is the known worst category.
  method2_correct <- NA_real_

  if (!all_sufficient) {

    worst_index <- which(
      as.numeric(fit$x.unique) == 1
    )[1]

    method2_correct <- as.numeric(
      as.integer(fit$x.ast.hat) ==
        as.integer(worst_index)
    )
  }

  data.frame(
    K = K,
    scenario = scenario,
    rep = rep,
    method1_reject = method1_reject,
    method2_correct = method2_correct,
    fit_failed = FALSE,
    error_message = NA_character_
  )
}

task_grid <- expand.grid(
  K = K_GRID,
  scenario = c(
    "all_sufficient",
    "one_insufficient"
  ),
  rep = seq_len(N_REPS_B),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

cat(sprintf(
  "Total tasks: %d\n\n",
  nrow(task_grid)
))

t0 <- Sys.time()

results_list <- furrr::future_pmap(
  task_grid,
  run_one_task,
  .options = furrr::furrr_options(
    seed = NULL,
    scheduling = 1
  ),
  .progress = TRUE
)

elapsed <- Sys.time() - t0

results_flat <- do.call(
  rbind,
  results_list
)

cat(sprintf(
  "\nElapsed time: %.2f minutes\n\n",
  as.numeric(elapsed, units = "mins")
))


partB_results <- do.call(
  rbind,
  lapply(
    K_GRID,
    function(K) {

      scen_all_sufficient <- results_flat[
        results_flat$K == K &
          results_flat$scenario == "all_sufficient",
        ,
        drop = FALSE
      ]

      scen_one_insufficient <- results_flat[
        results_flat$K == K &
          results_flat$scenario == "one_insufficient",
        ,
        drop = FALSE
      ]

      data.frame(
        K = K,

        n_per_category = floor(
          N_TOTAL_B / K
        ),

        
        method1_power = mean(
          scen_all_sufficient$method1_reject,
          na.rm = TRUE
        ),

       
        method1_type1_error = mean(
          scen_one_insufficient$method1_reject,
          na.rm = TRUE
        ),

       
        method2_selection_accuracy = mean(
          scen_one_insufficient$method2_correct,
          na.rm = TRUE
        ),

        all_sufficient_failures = sum(
          scen_all_sufficient$fit_failed,
          na.rm = TRUE
        ),

        one_insufficient_failures = sum(
          scen_one_insufficient$fit_failed,
          na.rm = TRUE
        )
      )
    }
  )
)

cat("============================================================\n")
cat("PART B RESULTS\n")
cat("============================================================\n\n")

print(
  partB_results,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("INTERPRETATION\n")
cat("============================================================\n\n")

cat(
  "all_sufficient scenario: H1 true -> Method 1 rejection = power.\n"
)

cat(
  "one_insufficient scenario: H0 true -> Method 1 rejection = Type-I error.\n"
)

cat(
  "Method 2 accuracy = probability of correctly identifying category 1\n"
)

cat(
  "as the true insufficient-follow-up category.\n\n"
)


OUTPUT_DIR <- if (dir.exists("03_simulation_study/results")) {
  "03_simulation_study/results"
} else {
  file.path(getwd(), "phase5B_standalone_outputs")
}

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  partB_results,
  file.path(
    OUTPUT_DIR,
    "partB_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  results_flat,
  file.path(
    OUTPUT_DIR,
    "partB_replication_level_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    summary = partB_results,
    replication_level = results_flat,

    settings = list(
      gamma = GAMMA,
      epsilon = EPS,
      q_99 = Q_99,
      N_REPS_B = N_REPS_B,
      N_BOOT_B = N_BOOT_B,
      K_GRID = K_GRID,
      N_TOTAL_B = N_TOTAL_B,
      HARD_TAU_C = HARD_TAU_C,
      EASY_TAU_C = EASY_TAU_C,
      TAU_TEST_B = TAU_TEST_B,
      workers = N_WORKERS
    )
  ),
  file.path(
    OUTPUT_DIR,
    "partB_results.rds"
  )
)
##   Type-I error should remain around/below nominal 0.05.
##   Power is confirmatory power under the all-sufficient scenario.
##   Method 2 accuracy is localization accuracy under one-insufficient.

png(
  file.path(
    OUTPUT_DIR,
    "partB_level_power_vs_K.png"
  ),
  width = 900,
  height = 650
)

plot(
  partB_results$K,
  partB_results$method1_type1_error,
  type = "b",
  ylim = c(0, 1),
  xlab = "Number of covariate categories (K)",
  ylab = "Probability",
  main = "Method 1 power / Type-I error and Method 2 accuracy"
)

lines(
  partB_results$K,
  partB_results$method1_power,
  type = "b"
)

lines(
  partB_results$K,
  partB_results$method2_selection_accuracy,
  type = "b"
)

abline(
  h = 0.05,
  lty = 2
)

legend(
  "bottomleft",
  legend = c(
    "Method 1 Type-I error",
    "Method 1 power",
    "Method 2 selection accuracy",
    "nominal alpha = 0.05"
  ),
  lty = c(1, 1, 1, 2),
  pch = c(1, 1, 1, NA),
  bty = "n"
)

dev.off()

plan(sequential)

cat("\nSaved outputs to:\n")
cat(
  normalizePath(
    OUTPUT_DIR
  ),
  "\n"
)

cat("\nDone.\n")
