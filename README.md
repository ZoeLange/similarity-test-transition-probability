# Testing Similarity of Competing Risks Models
### Reproducible Code for the Paper  
Lange, Z.K., Farhadizadeh, M., Dette, H., & Binder, N. (2025). Testing similarity of competing risks models by comparing transition probabilities. arXiv. https://arxiv.org/abs/2512.00583

---

## Overview

This repository provides the full simulation framework developed and presented in our paper **“Testing similarity of competing risks models by comparing transition probabilities”**.

The repository includes:

-  Data generation for competing risks with constant intensities
-  Likelihood estimation under administrative and random censoring  
-  Transition probability distance computation  
-  Constrained parametric bootstrap test  
-  Replication of two exemplary simulation scenarios  
-  Results saved automatically for each scenario

All code is written in **R**, modularized, documented, and easy to run.

---
## Running the Example Simulations

From the project root, run:

```r
source("scripts/main_simulations.R")



This script:

1. Loads all required functions

2. Runs two examplary simulation scenarios:

   - Administrative censoring Alternative 4

   - Random right censoring (Margin, ψ = 0.005)

3. Saves results in data/results/ as:

  - csv and RData files containing the full simulation results
