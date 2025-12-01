# Testing Similarity in Competing Risks Models
### Reproducible Code for the Paper  
**“Testing similarity of competing risks models by comparing transition probabilities”**  
Lange, Farhadizadeh, Dette & Binder (2025)

---

## Overview

This repository contains the full simulation framework used in our paper on  
**testing similarity of competing risks models** based on a  
**transition-probability–based distance**.

The repository includes:

-  Data generation for competing risks with constant hazards  
- Likelihood estimation under administrative and random censoring  
-  Transition probability distance computation  
-  Constrained parametric bootstrap test  
-  Full replication of the simulation scenarios  
-  Results saved automatically for each scenario

All code is written in **R**, modularized, documented, and easy to run.

---
## Running the Example Simulations

From the project root, run:

```r
source("scripts/main_simulations.R")



This script will:

1. Load all required functions

2. Run the example scenarios:

   - Administrative censoring

   - Random right censoring

Save results in data/results/ as:

- A file containing the full simulation results

