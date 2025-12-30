#' Helper Functions for BOSE Method
#' 
#' Utility functions for computing performance metrics and preparing data.
#' 
#' @author Wenqisi (Lydia) Pan

#===============================================================================
#                         PERFORMANCE METRICS
#===============================================================================

#' Squared Error
#'
#' @param estimate Estimated value
#' @param true_value True parameter value
#' @return Squared error
#' @export
SError <- function(estimate, true_value) {
  (estimate - true_value)^2
}

#' Absolute Relative Error
#'
#' @param estimate Estimated value
#' @param true_value True parameter value
#' @return Absolute relative error
#' @export
ARE <- function(estimate, true_value) {
  (estimate - true_value) / true_value
}

#===============================================================================
#                      DATA GENERATION FUNCTIONS
#===============================================================================

#' Generate Order Statistics from Normal Distribution
#'
#' @param n Sample size
#' @param mu True mean
#' @param sigma True standard deviation
#' @param S Scenario type (1, 2, or 3)
#' @return Vector of order statistics
#' @export
generate_order_stats_normal <- function(n, mu, sigma, S) {
  # Generate full sample
  sample_data <- rnorm(n, mean = mu, sd = sigma)
  sorted_sample <- sort(sample_data)
  
  # Extract order statistics based on Type 1 indices
  idx <- get_type1_indices(n, S)
  
  if (S == 1) {
    X <- c(sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3])
  } else if (S == 2) {
    X <- c(sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3])
  } else if (S == 3) {
    X <- c(sorted_sample[idx$k0],
           sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3],
           sorted_sample[idx$k4])
  }
  
  return(X)
}

#' Generate Order Statistics from t Distribution
#'
#' @param n Sample size
#' @param df Degrees of freedom
#' @param S Scenario type (1, 2, or 3)
#' @return Vector of order statistics
#' @export
generate_order_stats_t <- function(n, df, S) {
  sample_data <- rt(n, df = df)
  sorted_sample <- sort(sample_data)
  idx <- get_type1_indices(n, S)
  
  if (S == 1) {
    X <- c(sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3])
  } else if (S == 2) {
    X <- c(sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3])
  } else if (S == 3) {
    X <- c(sorted_sample[idx$k0],
           sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3],
           sorted_sample[idx$k4])
  }
  
  return(X)
}

#' Generate Order Statistics from Log-Normal Distribution
#'
#' @param n Sample size
#' @param meanlog Mean of log-normal distribution
#' @param sdlog Standard deviation of log-normal distribution
#' @param S Scenario type (1, 2, or 3)
#' @return Vector of order statistics
#' @export
generate_order_stats_lognormal <- function(n, meanlog, sdlog, S) {
  sample_data <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  sorted_sample <- sort(sample_data)
  idx <- get_type1_indices(n, S)
  
  if (S == 1) {
    X <- c(sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3])
  } else if (S == 2) {
    X <- c(sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3])
  } else if (S == 3) {
    X <- c(sorted_sample[idx$k0],
           sorted_sample[idx$k1], 
           sorted_sample[idx$k2], 
           sorted_sample[idx$k3],
           sorted_sample[idx$k4])
  }
  
  return(X)
}

#===============================================================================
#                      COMPARISON METHODS (WRAPPERS)
#===============================================================================

#' Get Comparison Method Estimates
#'
#' Wrapper function to call all comparison methods (Luo, Wan, BC, QE, MLN, BLUE).
#' Requires estmeansd and metaBLUE packages.
#'
#' @param X Vector of order statistics
#' @param n Sample size
#' @param S Scenario type
#' @return List with all method estimates
#' @export
get_comparison_estimates <- function(X, n, S) {
  
  # Initialize result list
  result <- list(
    mu_luo = NA,
    sigma_wan = NA,
    sigma_shi = NA,
    mu_bc = NA,
    sigma_bc = NA,
    mu_qe = NA,
    sigma_qe = NA,
    mu_mln = NA,
    sigma_mln = NA,
    mu_blue = NA,
    sigma_blue = NA
  )
  
  # Luo & Wan methods (available for all scenarios)
  if (S == 1) {
    result$sigma_wan <- Wan.std(X, n, type = "S1")$sigmahat
    result$mu_luo <- Luo.mean(X, n, type = "S1")$muhat
  } else if (S == 2) {
    result$sigma_wan <- Wan.std(X, n, type = "S2")$sigmahat
    result$mu_luo <- Luo.mean(X, n, type = "S2")$muhat
  } else if (S == 3) {
    result$sigma_wan <- Wan.std(X, n, type = "S3")$sigmahat
    result$mu_luo <- Luo.mean(X, n, type = "S3")$muhat
  }
  
  # Shi's method (only for S=3, five-number summary)
  if (S == 3) {
    phi_inv_1 <- qnorm((n - 0.375) / (n + 0.25))
    phi_inv_2 <- qnorm((0.75 * n - 0.125) / (n + 0.25))
    result$sigma_shi <- ((X[5] - X[1]) / ((2 + 0.14 * n^0.6) * phi_inv_1)) + 
                        ((X[4] - X[2]) / ((2 + 2 / (0.07 * n^0.6)) * phi_inv_2))
  }
  
  # BC, QE, MLN methods (S=2 and S=3)
  if (S == 2) {
    # Three quartiles: Q1, median, Q3
    tryCatch({
      bc <- bc.mean.sd(q1.val = X[1], med.val = X[2], q3.val = X[3], n = n, 
                       preserve.tail = FALSE, avoid.mc = FALSE)
      result$mu_bc <- bc$est.mean
      result$sigma_bc <- bc$est.sd
    }, error = function(e) {})
    
    tryCatch({
      qe <- qe.mean.sd(q1.val = X[1], med.val = X[2], q3.val = X[3], n = n)
      result$mu_qe <- qe$est.mean
      result$sigma_qe <- qe$est.sd
    }, error = function(e) {})
    
    tryCatch({
      mln <- mln.mean.sd(q1.val = X[1], med.val = X[2], q3.val = X[3], n = n)
      result$mu_mln <- mln$est.mean
      result$sigma_mln <- mln$est.sd
    }, error = function(e) {})
  }
  
  if (S == 3) {
    # Five-number summary
    tryCatch({
      bc <- bc.mean.sd(min.val = X[1], q1.val = X[2], med.val = X[3], 
                       q3.val = X[4], max.val = X[5], n = n, 
                       preserve.tail = FALSE, avoid.mc = FALSE)
      result$mu_bc <- bc$est.mean
      result$sigma_bc <- bc$est.sd
    }, error = function(e) {})
    
    tryCatch({
      qe <- qe.mean.sd(min.val = X[1], q1.val = X[2], med.val = X[3], 
                       q3.val = X[4], max.val = X[5], n = n)
      result$mu_qe <- qe$est.mean
      result$sigma_qe <- qe$est.sd
    }, error = function(e) {})
    
    tryCatch({
      mln <- mln.mean.sd(min.val = X[1], q1.val = X[2], med.val = X[3], 
                         q3.val = X[4], max.val = X[5], n = n)
      result$mu_mln <- mln$est.mean
      result$sigma_mln <- mln$est.sd
    }, error = function(e) {})
  }
  
  # BLUE method
  if (S == 1) {
    blue <- BLUE_s(X, n, type = "S1")
    result$mu_blue <- blue$muhat
    result$sigma_blue <- blue$sigmahat
  } else if (S == 2) {
    blue <- BLUE_s(X, n, type = "S2")
    result$mu_blue <- blue$muhat
    result$sigma_blue <- blue$sigmahat
  } else if (S == 3) {
    blue <- BLUE_s(X, n, type = "S3")
    result$mu_blue <- blue$muhat
    result$sigma_blue <- blue$sigmahat
  }
  
  return(result)
}
