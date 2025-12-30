#===============================================================================
#' BOSE Method: Type 1 Simulation Study
#' 
#' Compares BOSE (Bayesian Order Statistics Estimator) with existing methods
#' (Luo, Wan, BC, QE, MLN, BLUE) across different sample sizes and distributions.
#' 
#' Key finding: Type 1 quantiles maintain proper 95% credible interval coverage
#' while Type 7 quantiles exhibit systematic under-coverage.
#' 
#' @author Wenqisi (Lydia) Pan
#===============================================================================

# Load required packages
library(future.apply)
library(stringr)
library(estmeansd)
library(metaBLUE)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

# Source BOSE core functions
# source("R/bose_core.R")
# source("R/helper_functions.R")

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)
set.seed(42)

#===============================================================================
#                         SIMULATION PARAMETERS
#===============================================================================

# Prior parameters for inverse-gamma on sigma^2
alpha <- 0.01
beta <- 0.01
z <- 5  # Factor for setting adaptive grid bounds

# Sample sizes to test (no restriction on n = 4q+1)
n_values <- seq(5, 200, by = 1)
n_q_pairs <- data.frame(n = n_values)

# Number of Monte Carlo replications per sample size
n_reps <- 1000

# Grid sizes for posterior approximation
coarse_size <- 128
fine_size <- 256
n_post_samples <- 1000

# Output directory (modify as needed)
outpath <- "./simulation_results"  # Change this to your desired path

#===============================================================================
#                         DISTRIBUTION PARAMETERS
#===============================================================================

# Define parameter sets for each distribution
param_sets <- list(
  normal = list(
    list(mu_true = 5, sigma_true = 1),
    list(mu_true = 50, sigma_true = 17)
  )
  # Uncomment to add t and lognormal distributions
  # lognormal = list(
  #   list(lg_mu = 5, lg_sig = 0.1),
  #   list(lg_mu = 5, lg_sig = 0.2),
  #   list(lg_mu = 5, lg_sig = 0.3)
  # ),
  # t = list(
  #   list(df = 3),
  #   list(df = 5),
  #   list(df = 10),
  #   list(df = 25)
  # )
)

#===============================================================================
#                    MAIN COMPARISON FUNCTION
#===============================================================================

#' Run Single Comparison for Given n, S, and Distribution
#'
#' @param S Scenario type (1, 2, or 3)
#' @param dist Distribution name ("normal", "t", "lognormal")
#' @param n Sample size
#' @param params List of distribution parameters
#' @return Data frame with RMSE and coverage results
compute_comparison <- function(S, dist, n, params) {
  
  # Initialize storage for point estimates
  mu_bayes <- sigma_bayes <- rep(NA, n_reps)
  sigma_wan <- mu_luo <- sigma_shi <- rep(NA, n_reps)
  mu_bc <- sigma_bc <- rep(NA, n_reps)
  mu_qe <- sigma_qe <- rep(NA, n_reps)
  mu_mln <- sigma_mln <- rep(NA, n_reps)
  mu_blue <- sigma_blue <- rep(NA, n_reps)
  sample_mu <- sample_sigma <- rep(NA, n_reps)
  
  # Initialize storage for MSE
  mu_mse_bayes <- sigma_mse_bayes <- rep(NA, n_reps)
  sigma_mse_wan <- mu_mse_luo <- sigma_mse_shi <- rep(NA, n_reps)
  mu_mse_bc <- sigma_mse_bc <- rep(NA, n_reps)
  mu_mse_qe <- sigma_mse_qe <- rep(NA, n_reps)
  mu_mse_mln <- sigma_mse_mln <- rep(NA, n_reps)
  mu_mse_blue <- sigma_mse_blue <- rep(NA, n_reps)
  mu_mse_sample <- sigma_mse_sample <- rep(NA, n_reps)
  
  # Initialize coverage counters
  coverage_mu <- 0
  coverage_sigma <- 0
  
  # Extract distribution parameters
  if (dist == "normal") {
    mu_true <- params$mu_true
    sigma_true <- params$sigma_true
  } else if (dist == "t") {
    df <- params$df
  } else if (dist == "lognormal") {
    lg_mu <- params$lg_mu
    lg_sig <- params$lg_sig
  }
  
  # ========== Monte Carlo Loop ==========
  for (rep in 1:n_reps) {
    
    # Generate order statistics based on distribution
    if (dist == "normal") {
      X <- generate_order_stats_normal(n, mu_true, sigma_true, S)
      sample_data <- rnorm(n, mean = mu_true, sd = sigma_true)
    } else if (dist == "t") {
      X <- generate_order_stats_t(n, df, S)
      sample_data <- rt(n, df = df)
    } else if (dist == "lognormal") {
      X <- generate_order_stats_lognormal(n, lg_mu, lg_sig, S)
      sample_data <- rlnorm(n, meanlog = lg_mu, sdlog = lg_sig)
    }
    
    # Sample mean and SD (benchmark)
    sample_mu[rep] <- mean(sample_data)
    sample_sigma[rep] <- sd(sample_data)
    
    # ========== Get estimates from all comparison methods ==========
    comp_est <- get_comparison_estimates(X, n, S)
    
    mu_luo[rep] <- comp_est$mu_luo
    sigma_wan[rep] <- comp_est$sigma_wan
    sigma_shi[rep] <- comp_est$sigma_shi
    mu_bc[rep] <- comp_est$mu_bc
    sigma_bc[rep] <- comp_est$sigma_bc
    mu_qe[rep] <- comp_est$mu_qe
    sigma_qe[rep] <- comp_est$sigma_qe
    mu_mln[rep] <- comp_est$mu_mln
    sigma_mln[rep] <- comp_est$sigma_mln
    mu_blue[rep] <- comp_est$mu_blue
    sigma_blue[rep] <- comp_est$sigma_blue
    
    # ========== BOSE Method: Compute posterior without interpolation ==========
    L <- mu_luo[rep] - z * sigma_wan[rep]
    U <- mu_luo[rep] + z * sigma_wan[rep]
    L_sig <- max(sigma_wan[rep] / z, 1e-8)
    U_sig <- max(z * sigma_wan[rep], L_sig * 1.0001)
    
    samp <- get_posterior_samples_adaptive(
      X = X, n = n, S = S,  
      L = L, U = U, L_sig = L_sig, U_sig = U_sig,
      coarse = coarse_size, fine = fine_size, 
      mass = 0.99, n_samples = n_post_samples
    )
    
    post.mu <- samp$mu
    post.sig <- samp$sig
    
    # Posterior estimates
    mu_bayes[rep] <- mean(post.mu)
    sigma_bayes[rep] <- median(post.sig)
    
    # ========== Compute MSE and coverage ==========
    if (dist == "normal") {
      ci_mu <- quantile(post.mu, c(0.025, 0.975))
      ci_sigma <- quantile(post.sig, c(0.025, 0.975))
      coverage_mu <- coverage_mu + ifelse(ci_mu[1] <= mu_true & ci_mu[2] >= mu_true, 1, 0)
      coverage_sigma <- coverage_sigma + ifelse(ci_sigma[1] <= sigma_true & ci_sigma[2] >= sigma_true, 1, 0)
      mean_true <- mu_true
      sd_true <- sigma_true
    } else if (dist == "t") {
      mean_true <- 0
      sd_true <- sqrt(df / (df - 2))
    } else if (dist == "lognormal") {
      mean_true <- exp(lg_mu + (lg_sig^2) / 2)
      sd_true <- sqrt((exp(lg_sig^2) - 1) * exp(2 * lg_mu + lg_sig^2))
    }
    
    # MSE for all methods
    mu_mse_bayes[rep] <- SError(mu_bayes[rep], mean_true)
    sigma_mse_bayes[rep] <- SError(sigma_bayes[rep], sd_true)
    sigma_mse_wan[rep] <- SError(sigma_wan[rep], sd_true)
    mu_mse_luo[rep] <- SError(mu_luo[rep], mean_true)
    sigma_mse_shi[rep] <- SError(sigma_shi[rep], sd_true)
    mu_mse_bc[rep] <- SError(mu_bc[rep], mean_true)
    sigma_mse_bc[rep] <- SError(sigma_bc[rep], sd_true)
    mu_mse_qe[rep] <- SError(mu_qe[rep], mean_true)
    sigma_mse_qe[rep] <- SError(sigma_qe[rep], sd_true)
    mu_mse_mln[rep] <- SError(mu_mln[rep], mean_true)
    sigma_mse_mln[rep] <- SError(sigma_mln[rep], sd_true)
    mu_mse_blue[rep] <- SError(mu_blue[rep], mean_true)
    sigma_mse_blue[rep] <- SError(sigma_blue[rep], sd_true)
    mu_mse_sample[rep] <- SError(sample_mu[rep], mean_true)
    sigma_mse_sample[rep] <- SError(sample_sigma[rep], sd_true)
    
  } # end of replication loop
  
  # ========== Aggregate results across replications ==========
  mu_coverage <- coverage_mu / n_reps
  sigma_coverage <- coverage_sigma / n_reps
  
  mu_mse_Bayes <- mean(mu_mse_bayes)
  sigma_mse_Bayes <- mean(sigma_mse_bayes)
  mu_mse_Luo <- mean(mu_mse_luo)
  sigma_mse_Wan <- mean(sigma_mse_wan)
  sigma_mse_Shi <- mean(sigma_mse_shi)
  mu_mse_BC <- mean(mu_mse_bc, na.rm = TRUE)
  sigma_mse_BC <- mean(sigma_mse_bc, na.rm = TRUE)
  mu_mse_QE <- mean(mu_mse_qe, na.rm = TRUE)
  sigma_mse_QE <- mean(sigma_mse_qe, na.rm = TRUE)
  mu_mse_MLN <- mean(mu_mse_mln, na.rm = TRUE)
  sigma_mse_MLN <- mean(sigma_mse_mln, na.rm = TRUE)
  mu_mse_Blue <- mean(mu_mse_blue)
  sigma_mse_Blue <- mean(sigma_mse_blue)
  mu_mse_Sample <- mean(mu_mse_sample)
  sigma_mse_Sample <- mean(sigma_mse_sample)
  
  # Relative MSE (RMSE ratio to sample mean/SD)
  mu_rmse_Bayes <- mu_mse_Bayes / mu_mse_Sample
  sigma_rmse_Bayes <- sigma_mse_Bayes / sigma_mse_Sample
  mu_rmse_Luo <- mu_mse_Luo / mu_mse_Sample
  sigma_rmse_Wan <- sigma_mse_Wan / sigma_mse_Sample
  sigma_rmse_Shi <- sigma_mse_Shi / sigma_mse_Sample
  mu_rmse_BC <- mu_mse_BC / mu_mse_Sample
  sigma_rmse_BC <- sigma_mse_BC / sigma_mse_Sample
  mu_rmse_QE <- mu_mse_QE / mu_mse_Sample
  sigma_rmse_QE <- sigma_mse_QE / sigma_mse_Sample
  mu_rmse_MLN <- mu_mse_MLN / mu_mse_Sample
  sigma_rmse_MLN <- sigma_mse_MLN / sigma_mse_Sample
  mu_rmse_Blue <- mu_mse_Blue / mu_mse_Sample
  sigma_rmse_Blue <- sigma_mse_Blue / sigma_mse_Sample
  
  return(data.frame(
    n = n,
    S = S,
    mu_coverage = mu_coverage,
    sigma_coverage = sigma_coverage,
    mu_rmse_Bayes = mu_rmse_Bayes,
    sigma_rmse_Bayes = sigma_rmse_Bayes,
    sigma_rmse_Wan = sigma_rmse_Wan,
    mu_rmse_Luo = mu_rmse_Luo,
    sigma_rmse_Shi = sigma_rmse_Shi,
    mu_rmse_BC = mu_rmse_BC,
    sigma_rmse_BC = sigma_rmse_BC,
    mu_rmse_QE = mu_rmse_QE,
    sigma_rmse_QE = sigma_rmse_QE,
    mu_rmse_MLN = mu_rmse_MLN,
    sigma_rmse_MLN = sigma_rmse_MLN,
    mu_rmse_Blue = mu_rmse_Blue,
    sigma_rmse_Blue = sigma_rmse_Blue
  ))
}

#===============================================================================
#                    RUN SIMULATIONS AND SAVE RESULTS
#===============================================================================

#' Generate and Save Simulation Results
#'
#' @param n_q_pairs Data frame with sample sizes
#' @param param_sets List of parameter sets for each distribution
#' @param S_values Vector of scenario types to test (default: c(1, 2, 3))
#' @return List of result data frames
generate_and_save_results <- function(n_q_pairs, param_sets, S_values = c(1, 2, 3)) {
  results_list <- list()
  
  for (dist in names(param_sets)) {
    for (params in param_sets[[dist]]) {
      param_str <- paste(names(params), params, sep = "_", collapse = "_")
      
      for (S in S_values) {
        message("🚀 Running for S = ", S, ", Distribution = ", dist, 
                ", Params = ", param_str, " ...")
        start_time <- Sys.time()
        
        # Parallel computation across sample sizes
        results_list_i <- future_lapply(1:nrow(n_q_pairs), function(i) {
          n <- n_q_pairs$n[i]
          compute_comparison(S = S, dist = dist, n = n, params = params)
        }, future.seed = TRUE)
        
        results <- dplyr::bind_rows(results_list_i)
        end_time <- Sys.time()
        execution_time <- end_time - start_time
        
        var_name <- paste0("results_s", S, "_", dist, "_", param_str)
        results_list[[var_name]] <- results
        message("✅ Done: ", var_name, " | ⏳ Time taken: ", execution_time)
      }
    }
  }
  return(results_list)
}

#' Save All Simulation Results
#'
#' @param results_list List of result data frames
#' @param output_dir Output directory path
#' @param save_csv Save as CSV (default: TRUE)
#' @param save_rds Save as RDS (default: TRUE)
save_all_results <- function(results_list,
                              output_dir = outpath,
                              save_csv = TRUE,
                              save_rds = TRUE) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (name_prefix in names(results_list)) {
    result <- results_list[[name_prefix]]
    
    if (save_csv) {
      csv_path <- file.path(output_dir, paste0(name_prefix, ".csv"))
      write.csv(result, csv_path, row.names = FALSE)
      message("✅ CSV saved: ", csv_path)
    }
    
    if (save_rds) {
      rds_path <- file.path(output_dir, paste0(name_prefix, ".rds"))
      saveRDS(result, rds_path)
      message("✅ RDS saved: ", rds_path)
    }
  }
}

#===============================================================================
#                              MAIN EXECUTION
#===============================================================================

# Run simulations
cat("\n", strrep("=", 80), "\n")
cat("     BOSE METHOD: TYPE 1 SIMULATION STUDY\n")
cat(strrep("=", 80), "\n\n")

results_all <- generate_and_save_results(n_q_pairs, param_sets)

# Save all results
save_all_results(results_all)

cat("\n✅ All simulations complete!\n")
cat(strrep("=", 80), "\n\n")
