#######################################################################
# scenarios.R
#
# High-level scenario wrapper functions:
#
#   - run_scenario_adm()
#   - run_scenario_random()
#
# These functions:
#   * run N simulation replicates (using run_simulation_* from simulation_core.R)
#   * optionally use parallel computation via foreach/doMC
#   * compute mean rejection probabilities
#   * save results to data/results
#
#######################################################################
# -------------------------------------------------------------------
# Required packages
# -------------------------------------------------------------------
if (!requireNamespace("foreach", quietly = TRUE)) {
  stop(
    "The package 'foreach' is required but not installed.\n",
    "Please run: install.packages('foreach')"
  )
}
library(foreach)

# Optional: doMC for parallel backend (macOS / Linux only)
if (requireNamespace("doMC", quietly = TRUE)) {
  library(doMC)
  detected_cores <- parallel::detectCores()
  cores_to_use <- if (is.na(detected_cores)) 2 else max(1, detected_cores - 1)
  doMC::registerDoMC(cores = cores_to_use)
  message("Parallel backend 'doMC' loaded. Using ", cores_to_use, " cores.")
} else {
  message("Package 'doMC' not installed — running in serial mode.")
}


##############################################
# Scenario: Type I Censoring
##############################################

#' @title Run Full Simulation Scenario (Administrative Censoring)
#' @export
run_scenario_adm <- function(
    scen_name,
    n1, n2,
    alpha1_true,
    alpha2_true,
    threshold_null,
    alpha  = 0.05,
    B      = 500,
    N      = 1000,
    T_max  = 90,
    start  = rep(0.001, 3)
) {
  
  message("------------------------------------------------------------")
  message(" Running scenario (administrative censoring): ", scen_name)
  message("------------------------------------------------------------")
 
  message("Starting ", N, " simulation replicates...")
  
  results <- foreach::foreach(i = 1:N) %dopar% {
    run_simulation_adm(
      i,
      n1, n2,
      alpha1_true, alpha2_true,
      alpha, B, T_max,
      start,
      threshold_null
    )
  }
  
  results <- unlist(results)
  mean_rej <- mean(results)
  
  out <- list(
    rej       = results,
    meanval   = mean_rej,
    threshold = threshold_null,
    scen_name = scen_name,
    params    = list(
      n1            = n1,
      n2            = n2,
      alpha1_true   = alpha1_true,
      alpha2_true   = alpha2_true,
      threshold_null = threshold_null,
      alpha         = alpha,
      B             = B,
      N             = N,
      T_max         = T_max
    )
  )
  
  save_dir <- file.path("data", "results")
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Save RData ----
  save_path <- file.path(save_dir, paste0(scen_name, ".RData"))
  save(out, file = save_path)
  
  # ---- Save CSV summary ----
  csv_path <- file.path(save_dir, paste0(scen_name, "_summary.csv"))
  
  summary_df <- data.frame(
    scen_name        = scen_name,
    mean_rejection   = mean_rej,
    threshold_null   = threshold_null,
    n1               = n1,
    n2               = n2,
    B                = B,
    N                = N,
    T_max            = T_max
  )
  
  write.csv(summary_df, csv_path, row.names = FALSE)
  
  message("Scenario '", scen_name, "' completed.")
  message("Saved files:")
  message("  - ", save_path)
  message("  - ", csv_path)
  
  invisible(out)
}


##############################################
# Scenario: Random right censoring
##############################################

#' @title Run Full Simulation Scenario (Random Right Censoring)
#' @export
run_scenario_random <- function(
    scen_name,
    n1, n2,
    alpha1_true,
    alpha2_true,
    censoring_rate,
    threshold_null,
    alpha  = 0.05,
    B      = 500,
    N      = 1000,
    T_max  = 90,
    start  = rep(0.001, 4)
) {
  
  message("------------------------------------------------------------")
  message(" Running scenario (random right censoring): ", scen_name)
  message("------------------------------------------------------------")
  

  message("Starting ", N, " simulation replicates...")
  
  results <- foreach::foreach(i = 1:N) %dopar% {
    run_simulation_random(
      i,
      n1, n2,
      alpha1_true, alpha2_true,
      censoring_rate,
      alpha, B, T_max,
      start,
      threshold_null
    )
  }
  
  results <- unlist(results)
  mean_rej <- mean(results)
  
  out <- list(
    rej        = results,
    meanval    = mean_rej,
    threshold  = threshold_null,
    scen_name  = scen_name,
    params     = list(
      n1             = n1,
      n2             = n2,
      alpha1_true    = alpha1_true,
      alpha2_true    = alpha2_true,
      censoring_rate = censoring_rate,
      threshold_null  = threshold_null,
      alpha          = alpha,
      B              = B,
      N              = N,
      T_max          = T_max
    )
  )
  
  save_dir <- file.path("data", "results")
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Save RData ----
  save_path <- file.path(save_dir, paste0(scen_name, ".RData"))
  save(out, file = save_path)
  
  
  # ---- Save CSV ----
  csv_path <- file.path(save_dir, paste0(scen_name, "_summary.csv"))
  
  summary_df <- data.frame(
    scen_name        = scen_name,
    mean_rejection   = mean_rej,
    threshold_null   = threshold_null,
    n1               = n1,
    n2               = n2,
    censoring_rate   = censoring_rate,
    B                = B,
    N                = N,
    T_max            = T_max
  )
  
  write.csv(summary_df, csv_path, row.names = FALSE)
  
  message("Scenario '", scen_name, "' completed.")
  message("Saved files:")
  message("  - ", save_path)
  message("  - ", csv_path)
  
  invisible(out)
}

#######################################################################
# End of scenarios.R
#######################################################################
