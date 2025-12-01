#######################################################################
# utilities.R
#
# This file contains all lower-level utility functions required for:
#   - Data generation (simulate_competing_risks)
#   - Likelihood construction (loglik_adm, loglik_random)
#   - Transition probability distance computation
#
# These functions are called by simulation_core.R and scenarios.R.
#
#######################################################################

#############################
# 1) Data generation
#############################
#' @author Maryam Farhadizadeh
#' @title Simulate Competing Risks Data (Constant Hazards)
#'
#' @description
#' Generates event times, cause indicators, and observed times
#' for a 3-state competing risks model with constant hazards.
#'
#' Supports administrative, exponential, uniform, or no censoring.
#'
#' @export
simulate_competing_risks <- function(
    n,
    hazards,
    beta       = c(0, 0, 0),
    cov_prob   = 0.5,
    censoring  = list(model = "unif", param = NULL)
) {
  
  stopifnot(length(hazards) == 3)
  
  # Covariate
  covariate <- rbinom(n, size = 1, prob = cov_prob)
  
  # Individual hazards
  csh <- matrix(NA, nrow = n, ncol = 3)
  for (j in 1:3) {
    csh[, j] <- hazards[j] * exp(beta[j] * covariate)
  }
  total_haz <- rowSums(csh)
  
  # Event times
  event_time <- rexp(n, rate = total_haz)
  
  # Event cause (multinomial)
  cause <- numeric(n)
  for (i in seq_len(n)) {
    p <- csh[i, ] / total_haz[i]
    p[p < 0 | p > 1 | is.nan(p)] <- 1e-12
    cause[i] <- which(rmultinom(1, 1, p) == 1)
  }
  
  # Censoring
  if (is.null(censoring)) {
    censor_time <- event_time
  } else {
    model <- censoring$model
    param <- censoring$param
    censor_time <- switch(model,
                          "unif" = runif(n, 0, param),
                          "exp"  = rexp(n, rate = param),
                          "adm"  = rep(param, n),
                          stop("Unknown censoring type: ", model))
  }
  
  # Observed outcomes
  obs_time  <- pmin(event_time, censor_time)
  obs_event <- as.numeric(event_time <= censor_time) * cause
  
  data.frame(
    id        = seq_len(n),
    obs_time  = obs_time,
    obs_event = obs_event,
    covariate = covariate
  )
}

#############################
# 2) Likelihood functions
#############################

#' @title Negative Log-Likelihood (Administrative Censoring)
#' @export
loglik_adm <- function(data) {
  
  eps <- 1e-10
  safe_log <- function(x) ifelse(x <= 0, log(eps), log(x))
  
  # Split observations by event type
  state1 <- data[data$obs_event == 1, ]
  state2 <- data[data$obs_event == 2, ]
  state3 <- data[data$obs_event == 3, ]
  
  f <- function(theta) {
    theta <- pmax(theta, eps)
    total_time <- sum(data$obs_time)
    
    val <- -(
      nrow(state1) * safe_log(theta[1]) - theta[1] * total_time +
        nrow(state2) * safe_log(theta[2]) - theta[2] * total_time +
        nrow(state3) * safe_log(theta[3]) - theta[3] * total_time
    )
    if (!is.finite(val)) 1e10 else val
  }
  
  f
}

#' @title Negative Log-Likelihood (Random Right Censoring)
#' @export
loglik_random <- function(data) {
  
  eps <- 1e-10
  safe_log <- function(x) ifelse(x <= 0, log(eps), log(x))
  
  state1 <- data[data$obs_event == 1, ]
  state2 <- data[data$obs_event == 2, ]
  state3 <- data[data$obs_event == 3, ]
  cens   <- data[data$obs_event == 0, ]
  
  f <- function(theta) {
    theta <- pmax(theta, eps)
    
    haz1 <- theta[1]
    haz2 <- theta[2]
    haz3 <- theta[3]
    cens_rate <- theta[4]
    
    total_time <- sum(data$obs_time)
    
    val <- -(
      nrow(state1) * safe_log(haz1) - haz1 * total_time +
        nrow(state2) * safe_log(haz2) - haz2 * total_time +
        nrow(state3) * safe_log(haz3) - haz3 * total_time +
        nrow(cens)   * safe_log(cens_rate) - cens_rate * total_time
    )
    
    if (!is.finite(val)) 1e10 else val
  }
  
  f
}

#############################
# 3) Transition probability distances
#############################

#' @title Max TP Distance (Administrative)
#' @export
compute_tp_distance_adm <- function(alpha1, alpha2, T_max) {
  
  stopifnot(length(alpha1) == length(alpha2))
  k <- length(alpha1)
  
  lambda1 <- sum(alpha1)
  lambda2 <- sum(alpha2)
  
  compute_j <- function(j) {
    a1 <- alpha1[j]
    a2 <- alpha2[j]
    
    D <- function(t) {
      P1 <- (a1 / lambda1) * (1 - exp(-lambda1 * t))
      P2 <- (a2 / lambda2) * (1 - exp(-lambda2 * t))
      abs(P1 - P2)
    }
    
    if (lambda1 != lambda2) {
      t_star <- log(a2 / a1) / (lambda2 - lambda1)
    } else {
      t_star <- NA
    }
    
    D_T <- D(T_max)
    
    if (!is.na(t_star) && t_star > 0 && t_star <= T_max) {
      D_star <- D(t_star)
    } else {
      D_star <- -Inf
    }
    
    max(D_T, D_star)
  }
  
  max_dist <- sapply(1:k, compute_j)
  list(
    max_distances = max_dist,
    global_max_distance = max(max_dist)
  )
}

#' @title Max TP Distance (Random Censoring)
#' @export
compute_tp_distance_random <- function(alpha1, alpha2, cens1, cens2, T_max) {
  
  stopifnot(length(alpha1) == length(alpha2))
  k <- length(alpha1)
  
  lambda1 <- sum(alpha1) + cens1
  lambda2 <- sum(alpha2) + cens2
  
  compute_j <- function(j) {
    a1 <- alpha1[j]
    a2 <- alpha2[j]
    
    D <- function(t) {
      P1 <- (a1 / lambda1) * (1 - exp(-lambda1 * t))
      P2 <- (a2 / lambda2) * (1 - exp(-lambda2 * t))
      abs(P1 - P2)
    }
    
    if (lambda1 != lambda2) {
      t_star <- (log(a2) - log(a1)) / (lambda2 - lambda1)
    } else {
      t_star <- NA
    }
    
    D_T <- D(T_max)
    D_star <- if (!is.na(t_star) && t_star > 0 && t_star <= T_max) D(t_star) else -Inf
    
    max(D_T, D_star)
  }
  
  max_dist <- sapply(1:k, compute_j)
  list(
    max_distances = max_dist,
    global_max_distance = max(max_dist)
  )
}

#######################################################################
# End of utilities.R
#######################################################################
