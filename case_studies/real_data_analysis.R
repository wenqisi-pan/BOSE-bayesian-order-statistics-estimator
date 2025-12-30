#' Real Data Analysis for BOSE Method
#' 
#' Main script for applying BOSE to real cardiovascular studies.
#' Estimates mean and SD from reported order statistics and computes
#' mean differences (MD) and standardized mean differences (SMD).
#' 
#' @author Wenqisi (Lydia) Pan

# Required packages
library(tidyverse)
library(parallel)
library(future.apply)
library(stringr)
library(estmeansd)
library(metaBLUE)

# Source core BOSE functions
# source("R/bose_core.R")
# source("R/helper_functions.R")
# source("case_studies/skewness_assessment.R")

# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)
set.seed(100)

#===============================================================================
#                         GLOBAL PARAMETERS
#===============================================================================

alpha <- 0.01
beta <- 0.01
z <- 5

n_posterior_samples <- 1000  # Number of posterior samples for inference
grid_size_default <- 1000    # Grid size for real data analysis

#===============================================================================
#                         EXAMPLE STUDY DATA
#===============================================================================

#' Example cardiovascular studies data
#' 
#' Three case-control studies with three-number summaries {min, median, max}.
#' Studies 4-5 have incomplete data (S=0) for demonstration purposes.
study_inputs <- tibble::tibble(
  Study = paste0("Study ", rep(1:5, each = 2), c(" Case", " Control")),
  S = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0),
  X_list = list(
    c(2.25, 16.0, 74.25),   # Study 1 Case
    c(9.0, 27.25, 132.5),   # Study 1 Control
    c(43.75, 65.75, 130.5), # Study 2 Case
    c(48.5, 69.5, 125.0),   # Study 2 Control
    c(16.75, 39.75, 89.25), # Study 3 Case
    c(26.25, 65.5, 114.75), # Study 3 Control
    NA, NA, NA, NA          # Studies 4-5 (incomplete data)
  ),
  n = c(40, 40, 40, 38, 15, 15, 51, 51, 22, 23)
)

#===============================================================================
#                    BAYESIAN ESTIMATION FOR REAL DATA
#===============================================================================

#' Estimate Parameters for Real Data Using BOSE
#'
#' Wrapper function that applies BOSE estimation to a single study
#' with optional contour plot data generation.
#'
#' @param S Scenario type (1, 2, or 3)
#' @param X Vector of observed order statistics
#' @param n Sample size
#' @param grid_size Grid size for posterior approximation (default: 1000)
#' @param return_contour_data Logical; if TRUE, return grid data for plotting
#' @return List with parameter estimates and optional contour data
#'
#' @export
estimate_bayes_real_data <- function(S, X, n, 
                                     grid_size = 1000,
                                     return_contour_data = FALSE) {
  
  # Get initial estimates using frequentist methods
  if (S == 1) {
    sigma_wan <- Wan.std(X, n, type = "S1")$sigmahat
    mu_luo <- Luo.mean(X, n, type = "S1")$muhat
  } else if (S == 2) {
    sigma_wan <- Wan.std(X, n, type = "S2")$sigmahat
    mu_luo <- Luo.mean(X, n, type = "S2")$muhat
  } else if (S == 3) {
    sigma_wan <- Wan.std(X, n, type = "S3")$sigmahat
    mu_luo <- Luo.mean(X, n, type = "S3")$muhat
  }
  
  # Set adaptive grid bounds
  L <- mu_luo - z * sigma_wan
  U <- mu_luo + z * sigma_wan
  L_sig <- max(sigma_wan / z, 1e-8)
  U_sig <- max(z * sigma_wan, L_sig * 1.0001)
  
  # Draw posterior samples
  samp <- get_posterior_samples_adaptive(
    X = X, n = n, S = S,
    L = L, U = U, L_sig = L_sig, U_sig = U_sig,
    coarse = 128, fine = 256,
    mass = 0.99, n_samples = n_posterior_samples
  )
  
  post_mu <- samp$mu
  post_sig <- samp$sig
  
  # Point estimates
  mu_mean <- mean(post_mu)
  sigma_mean <- median(post_sig)
  
  # Standard errors (posterior SD)
  mu_se <- sd(post_mu)
  sigma_se <- sd(post_sig)
  
  # 95% credible intervals
  mu_CI <- quantile(post_mu, c(0.025, 0.975))
  sigma_CI <- quantile(post_sig, c(0.025, 0.975))
  
  # Get Type 1 indices for reporting
  idx <- get_type1_indices(n, S)
  
  result <- list(
    mu_mean = mu_mean,
    sigma_mean = sigma_mean,
    mu_se = mu_se,
    sigma_se = sigma_se,
    mu_CI = as.numeric(mu_CI),
    sigma_CI = as.numeric(sigma_CI),
    mu_samples = post_mu,
    sigma_samples = post_sig,
    type1_indices = idx,
    initial_estimates = list(mu_luo = mu_luo, sigma_wan = sigma_wan),
    X_observed = X,
    study_info = paste0("n=", n, ", S=", S)
  )
  
  # Optionally return grid data for contour plots
  if (return_contour_data) {
    # Create finer grid for visualization
    mu_grid <- seq(L, U, length.out = grid_size)
    sigma_grid <- seq(L_sig, U_sig, length.out = grid_size)
    
    logpost_grid <- compute_logpost_grid_fast(mu_grid, sigma_grid, X, n, S, L, U)
    posterior_density <- exp(logpost_grid - max(logpost_grid, na.rm = TRUE))
    posterior_density <- posterior_density / sum(posterior_density)
    
    result$mu_grid <- mu_grid
    result$sigma_grid <- sigma_grid
    result$posterior_density <- posterior_density
  }
  
  return(result)
}

#===============================================================================
#                    MAIN ANALYSIS FUNCTION
#===============================================================================

#' Run Real Data Analysis on Multiple Studies
#'
#' Applies BOSE estimation to multiple studies and computes
#' mean differences (MD) and standardized mean differences (SMD).
#'
#' @param study_data Data frame with study information (default: study_inputs)
#' @return List with results data frame and formatted summary table
#'
#' @details For case-control studies, computes:
#' - MD: Mean difference (Control - Case)
#' - SMD: Standardized mean difference (Cohen's d)
#' - CV: Coefficient of variation (SD/mean)
#'
#' @examples
#' \dontrun{
#' results <- run_real_data_analysis()
#' print(results$final_results)
#' }
#'
#' @export
run_real_data_analysis <- function(study_data = study_inputs) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("              Type 1 Real Data Analysis\n")
  cat(strrep("=", 60), "\n\n")
  
  start_time <- Sys.time()
  
  # Apply BOSE to each study
  results <- study_data %>%
    rowwise() %>%
    mutate(result = list(
      if (S > 0 && !is.na(X_list[[1]])) {
        estimate_bayes_real_data(S = S, X = X_list, n = n, 
                                 grid_size = grid_size_default,
                                 return_contour_data = FALSE)
      } else {
        NULL
      }
    )) %>%
    ungroup() %>%
    mutate(
      # Extract point estimates
      mu_mean     = sapply(result, function(res) if (!is.null(res)) res$mu_mean else NA_real_),
      mu_se       = sapply(result, function(res) if (!is.null(res)) res$mu_se else NA_real_),
      mu_lower    = sapply(result, function(res) if (!is.null(res)) res$mu_CI[1] else NA_real_),
      mu_upper    = sapply(result, function(res) if (!is.null(res)) res$mu_CI[2] else NA_real_),
      sigma_mean  = sapply(result, function(res) if (!is.null(res)) res$sigma_mean else NA_real_),
      sigma_se    = sapply(result, function(res) if (!is.null(res)) res$sigma_se else NA_real_),
      sigma_lower = sapply(result, function(res) if (!is.null(res)) res$sigma_CI[1] else NA_real_),
      sigma_upper = sapply(result, function(res) if (!is.null(res)) res$sigma_CI[2] else NA_real_),
      
      # Type 1 index information
      k_median    = sapply(result, function(res) if (!is.null(res)) res$type1_indices$k2 else NA_integer_),
      
      # Coefficient of variation (CV)
      cv_draws    = lapply(result, function(res) if (!is.null(res)) res$sigma_samples / res$mu_samples else NULL),
      cv_mean     = sapply(cv_draws, function(cv) if (!is.null(cv)) mean(cv, na.rm = TRUE) else NA_real_),
      cv_lower    = sapply(cv_draws, function(cv) if (!is.null(cv)) quantile(cv, 0.025) else NA_real_),
      cv_upper    = sapply(cv_draws, function(cv) if (!is.null(cv)) quantile(cv, 0.975) else NA_real_)
    ) %>%
    mutate(
      # Formatted output strings
      mu_fmt    = ifelse(!is.na(mu_mean), sprintf("%.3f (%.3f, %.3f)", mu_mean, mu_lower, mu_upper), "N/A"),
      sigma_fmt = ifelse(!is.na(sigma_mean), sprintf("%.3f (%.3f, %.3f)", sigma_mean, sigma_lower, sigma_upper), "N/A"),
      cv_fmt    = ifelse(!is.na(cv_mean), sprintf("%.1f%% (%.1f%%, %.1f%%)", cv_mean*100, cv_lower*100, cv_upper*100), "N/A")
    )
  
  end_time <- Sys.time()
  
  # Print Type 1 indices information
  cat("---------- Type 1 Indices (ceiling) ----------\n")
  results %>%
    filter(!is.na(k_median)) %>%
    select(Study, n, k_median) %>%
    mutate(
      k_median_formula = paste0("ceiling(0.5 * ", n, ") = ", k_median)
    ) %>%
    print()
  
  # Calculate MD and SMD for case-control pairs
  case_control <- results %>%
    mutate(Group = ifelse(grepl("Case", Study), "Case", "Control"),
           StudyID = gsub(" (Case|Control)", "", Study)) %>%
    select(StudyID, Group, mu_mean, sigma_mean, n)
  
  wide_data <- case_control %>%
    pivot_wider(names_from = Group, values_from = c(mu_mean, sigma_mean, n)) %>%
    mutate(
      MD = mu_mean_Control - mu_mean_Case,
      pooled_sd = sqrt(((n_Case - 1) * sigma_mean_Case^2 + (n_Control - 1) * sigma_mean_Control^2) /
                         (n_Case + n_Control - 2)),
      SMD = MD / pooled_sd
    )
  
  # Merge final results
  final_results <- results %>%
    mutate(StudyID = gsub(" (Case|Control)", "", Study)) %>%
    left_join(wide_data %>% select(StudyID, MD, SMD), by = "StudyID") %>%
    mutate(MD = sprintf("%.3f", MD), SMD = sprintf("%.3f", SMD)) %>%
    select(Study, n, k_median, mu_fmt, mu_se, sigma_fmt, sigma_se, cv_fmt, MD, SMD)
  
  cat("\n---------- Final Results ----------\n")
  print(final_results)
  
  cat("\n⏱️ Time taken: ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds\n")
  
  return(list(results = results, final_results = final_results))
}

#===============================================================================
#                              MAIN EXECUTION
#===============================================================================

# Uncomment to run analysis
# cat("\n", strrep("=", 80), "\n")
# cat("     TYPE 1 REAL DATA ANALYSIS WITH SKEWNESS ASSESSMENT\n")
# cat(strrep("=", 80), "\n")
# 
# # Step 1: Skewness Analysis
# cat("\n📊 Step 1: Skewness Analysis\n")
# skewness_analysis <- run_skewness_analysis(study_inputs)
# 
# # Step 2: Bayesian Parameter Estimation
# cat("\n📊 Step 2: Bayesian Parameter Estimation\n")
# analysis_results <- run_real_data_analysis(study_inputs)
# 
# cat("\n✅ Analysis complete!\n")
# cat(strrep("=", 80), "\n\n")
