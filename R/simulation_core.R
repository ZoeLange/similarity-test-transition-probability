#######################################################################
# simulation_core.R
#
# Contains the two main simulation engines:
#
#   - run_simulation_adm()    (administrative censoring)
#   - run_simulation_random() (random exponential censoring)
#
# Each performs one replicate of the constrained parametric bootstrap
# similarity test described in the paper:
#
#   "Testing similarity of competing risks models by comparing
#    transition probabilities"
#
#######################################################################

##############################################
# 1) Simulation under administrative censoring
##############################################
#' @author Maryam Farhadizadeh
#' @title One Simulation Replicate (Administrative Censoring)
#'
#' @description
#' Performs ONE run of the similarity test under administrative censoring.
#'
#' Steps:
#'  1. Generate data for both groups
#'  2. Estimate hazards via MLE
#'  3. Compute test statistic (TP distance)
#'  4. Perform constrained estimation if needed
#'  5. Run bootstrap under constrained null
#'  6. Return rejection indicator (0/1)
#' @export
run_simulation_adm <- function(
    i,
    n1, n2,
    alpha1_true, alpha2_true,
    alpha, B, T_max,
    start,
    threshold_null
) {
  
  # Reproducibility
  set.seed(123 * i)
  
  boot_vals <- numeric(B)
  
  # -------------------------------------------------------------
  # Step 1: Generate datasets under administrative censoring
  # -------------------------------------------------------------
  dat1 <- simulate_competing_risks(
    n         = n1,
    hazards   = alpha1_true,
    beta      = c(0, 0, 0),
    cov_prob  = 0.5,
    censoring = list(model = "adm", param = T_max)
  )
  
  dat2 <- simulate_competing_risks(
    n         = n2,
    hazards   = alpha2_true,
    beta      = c(0, 0, 0),
    cov_prob  = 0.5,
    censoring = list(model = "adm", param = T_max)
  )
  
  # -------------------------------------------------------------
  # Step 2: Estimate transition intensities via MLE
  # -------------------------------------------------------------
  alpha1_hat <- optim(par = start, fn = loglik_adm(dat1))$par
  alpha2_hat <- optim(par = start, fn = loglik_adm(dat2))$par
  
  # -------------------------------------------------------------
  # Step 3: Compute test statistic (max TP distance)
  # -------------------------------------------------------------
  t_stat <- compute_tp_distance_adm(
    alpha1 = alpha1_hat,
    alpha2 = alpha2_hat,
    T_max  = T_max
  )$global_max_distance
  
  # -------------------------------------------------------------
  # Step 4: Constrained maximization if t_stat < threshold_null
  # -------------------------------------------------------------
  
  # Joint negative log-likelihood
  joint_loglik <- function(theta) {
    theta1 <- theta[1:3]
    theta2 <- theta[4:6]
    loglik_adm(dat1)(theta1) + loglik_adm(dat2)(theta2)
  }
  
  # Equality constraint: distance == threshold_null
  const_distance <- function(theta) {
    dist <- compute_tp_distance_adm(
      alpha1 = theta[1:3],
      alpha2 = theta[4:6],
      T_max  = T_max
    )$global_max_distance
    abs(dist - threshold_null)
  }
  
  # Positivity constraints
  const_positive <- function(theta) theta - 0
  
  # Choose constrained or unconstrained model
  if (t_stat >= threshold_null) {
    theta_cons <- c(alpha1_hat, alpha2_hat)
  } else {
    res <- alabama::auglag(
      par    = c(alpha1_hat, alpha2_hat),
      fn     = joint_loglik,
      hin    = const_positive,
      heq    = const_distance,
      control.outer = list(trace = FALSE)
    )
    theta_cons <- res$par
  }
  
  # Constrained hazards
  alpha1_cons <- theta_cons[1:3]
  alpha2_cons <- theta_cons[4:6]
  
  # -------------------------------------------------------------
  # Step 5: Bootstrap under constrained model
  # -------------------------------------------------------------
  for (b in seq_len(B)) {
    
    dat1b <- simulate_competing_risks(
      n         = n1,
      hazards   = alpha1_cons,
      beta      = c(0, 0, 0),
      cov_prob  = 0.5,
      censoring = list(model = "adm", param = T_max)
    )
    
    dat2b <- simulate_competing_risks(
      n         = n2,
      hazards   = alpha2_cons,
      beta      = c(0, 0, 0),
      cov_prob  = 0.5,
      censoring = list(model = "adm", param = T_max)
    )
    
    alpha1b <- optim(par = alpha1_cons, fn = loglik_adm(dat1b))$par
    alpha2b <- optim(par = alpha2_cons, fn = loglik_adm(dat2b))$par
    
    boot_vals[b] <- compute_tp_distance_adm(
      alpha1 = alpha1b,
      alpha2 = alpha2b,
      T_max  = T_max
    )$global_max_distance
  }
  
  # -------------------------------------------------------------
  # Step 6: Rejection decision
  # -------------------------------------------------------------
  crit_val <- quantile(boot_vals, alpha)
  reject   <- ifelse(t_stat < crit_val, 1, 0)
  
  return(reject)
}



##############################################
# 2) Simulation under exponential random censoring
##############################################

#' @title One Simulation Replicate (Random Right Censoring)
#'
#' @description
#' Performs ONE run of the similarity test under exponential right censoring.
#'
#' Equivalent mathematically to your original runSimu_random().
#'
#' @export
run_simulation_random <- function(
    i,
    n1, n2,
    alpha1_true, alpha2_true,
    censoring_rate,
    alpha, B, T_max,
    start,
    threshold_null
) {
  
  set.seed(123 * i)
  boot_vals <- numeric(B)
  
  # -------------------------------------------------------------
  # Step 1: Generate data with exponential censoring
  # -------------------------------------------------------------
  dat1 <- simulate_competing_risks(
    n         = n1,
    hazards   = alpha1_true,
    beta      = c(0, 0, 0),
    cov_prob  = 0.5,
    censoring = list(model = "exp", param = censoring_rate)
  )
  
  dat2 <- simulate_competing_risks(
    n         = n2,
    hazards   = alpha2_true,
    beta      = c(0, 0, 0),
    cov_prob  = 0.5,
    censoring = list(model = "exp", param = censoring_rate)
  )
  
  # -------------------------------------------------------------
  # Step 2: MLE for hazards + censoring
  # -------------------------------------------------------------
  alpha1_hat <- optim(par = start, fn = loglik_random(dat1))$par
  alpha2_hat <- optim(par = start, fn = loglik_random(dat2))$par
  
  # Test statistic
  t_stat <- compute_tp_distance_random(
    alpha1 = alpha1_hat[1:3],
    alpha2 = alpha2_hat[1:3],
    cens1  = alpha1_hat[4],
    cens2  = alpha2_hat[4],
    T_max  = T_max
  )$global_max_distance
  
  # -------------------------------------------------------------
  # Step 3: Constrained optimization
  # -------------------------------------------------------------
  joint_loglik <- function(theta) {
    theta1 <- theta[1:4]
    theta2 <- theta[5:8]
    loglik_random(dat1)(theta1) + loglik_random(dat2)(theta2)
  }
  
  const_distance <- function(theta) {
    dist <- compute_tp_distance_random(
      alpha1 = theta[1:3],
      alpha2 = theta[5:7],
      cens1  = theta[4],
      cens2  = theta[8],
      T_max  = T_max
    )$global_max_distance
    abs(dist - threshold_null)
  }
  
  const_positive <- function(theta) theta - 0
  
  if (t_stat >= threshold_null) {
    theta_cons <- c(alpha1_hat, alpha2_hat)
  } else {
    res <- alabama::auglag(
      par    = c(alpha1_hat, alpha2_hat),
      fn     = joint_loglik,
      hin    = const_positive,
      heq    = const_distance,
      control.outer = list(trace = FALSE)
    )
    theta_cons <- res$par
  }
  
  alpha1_cons <- theta_cons[1:3]
  alpha2_cons <- theta_cons[5:7]
  cens1_cons  <- theta_cons[4]
  cens2_cons  <- theta_cons[8]
  
  # -------------------------------------------------------------
  # Step 4: Bootstrap
  # -------------------------------------------------------------
  for (b in seq_len(B)) {
    
    dat1b <- simulate_competing_risks(
      n         = n1,
      hazards   = alpha1_cons,
      beta      = c(0, 0, 0),
      cov_prob  = 0.5,
      censoring = list(model = "exp", param = cens1_cons)
    )
    
    dat2b <- simulate_competing_risks(
      n         = n2,
      hazards   = alpha2_cons,
      beta      = c(0, 0, 0),
      cov_prob  = 0.5,
      censoring = list(model = "exp", param = cens2_cons)
    )
    
    alpha1b <- optim(par = c(alpha1_cons, cens1_cons), fn = loglik_random(dat1b))$par
    alpha2b <- optim(par = c(alpha2_cons, cens2_cons), fn = loglik_random(dat2b))$par
    
    boot_vals[b] <- compute_tp_distance_random(
      alpha1 = alpha1b[1:3],
      alpha2 = alpha2b[1:3],
      cens1  = alpha1b[4],
      cens2  = alpha2b[4],
      T_max  = T_max
    )$global_max_distance
  }
  
  # -------------------------------------------------------------
  # Step 5: Decision
  # -------------------------------------------------------------
  crit_val <- quantile(boot_vals, alpha)
  reject   <- ifelse(t_stat < crit_val, 1, 0)
  return(reject)
}


#######################################################################
# End of simulation_core.R
#######################################################################
