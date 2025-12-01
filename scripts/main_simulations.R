###############################################################
# Author: Maryam Farhadizadeh
# main_simulations.R
#
# Entry point to reproduce the example simulation scenarios from:
# "Testing similarity of competing risks models by comparing transition probabilities"
#
# Running this script:
#   1. Loads all necessary functions
#   2. Runs two example scenarios (admin + random)
#   3. Saves results in data/results/
###############################################################

message("Running example simulation scenarios...")

# Load project files
source("R/utilities.R")
source("R/simulation_core.R")
source("R/scenarios.R")

message("Functions loaded.")

#############
# Example 1 #
#############
message("Running Example 1: Administrative censoring (Margin case)")

admin_margin_example <- run_scenario_adm(
  scen_name     = "admin_margin_example",
  n1            = 100,
  n2            = 100,
  alpha1_true   = c(0.0009, 0.0011, 0.0004),
  alpha2_true   = c(0.0008, 0.0012, 0.0007),
  threshold_null = 0.11805,
  N             = 1000,
  B             = 500,
  T_max         = 90
)

message("Example 1 (Margin) — Mean rejection probability: ",
        admin_margin_example$meanval)
message("Saved to: data/results/admin_margin_example.*")
#############
# Example 2 #
#############
message("Running Example 2: Random censoring (Alternative 4, ψ = 0.005)")

random_alternative4_exp005 <- run_scenario_random(
  scen_name      = "random_alternative4_exp005",
  n1             = 50,
  n2             = 50,
  alpha1_true    = c(0.0023, 0.0011, 0.0004),
  alpha2_true    = c(0.0008, 0.0026, 0.0019),
  censoring_rate = 0.005,
  threshold_null = 0.0960,
  N              = 1000,
  B              = 500,
  T_max          = 90
)


message("Example 2 (Alternative 4) — Mean rejection probability: ",
        random_alternative4_exp005$meanval)
message("Saved to: data/results/random_alternative4_exp005.*")

message("All example scenarios completed.")
