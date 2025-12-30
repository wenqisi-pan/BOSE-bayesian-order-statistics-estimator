#' BOSE Full Workflow Example
#' 
#' Complete workflow demonstrating:
#' 1. Data preparation
#' 2. Skewness assessment
#' 3. Bayesian estimation
#' 4. Result visualization
#' 5. Meta-analysis application
#' 
#' @author Wenqisi (Lydia) Pan

#===============================================================================
#                              SETUP
#===============================================================================

# Load all required functions
source("R/bose_core.R")
source("R/helper_functions.R")
source("case_studies/skewness_assessment.R")
source("case_studies/real_data_analysis.R")
source("case_studies/contour_plots.R")

# Load required packages
library(tidyverse)
library(estmeansd)
library(metaBLUE)

# Set random seed
set.seed(2024)

cat("\n")
cat(strrep("=", 80), "\n")
cat("                    BOSE METHOD: FULL WORKFLOW EXAMPLE\n")
cat(strrep("=", 80), "\n\n")

#===============================================================================
#                       STEP 1: PREPARE DATA
#===============================================================================

cat("STEP 1: Prepare Study Data\n")
cat(strrep("-", 80), "\n\n")

# Example: Three case-control studies reporting three-number summaries
# (Real cardiovascular disease data from literature)

study_data <- tibble(
  Study = c("Study 1 Case", "Study 1 Control",
            "Study 2 Case", "Study 2 Control",
            "Study 3 Case", "Study 3 Control"),
  S = c(1, 1, 1, 1, 1, 1),  # All use three-number summary
  X_list = list(
    c(2.25, 16.0, 74.25),    # Study 1 Case:    min, median, max
    c(9.0, 27.25, 132.5),    # Study 1 Control
    c(43.75, 65.75, 130.5),  # Study 2 Case
    c(48.5, 69.5, 125.0),    # Study 2 Control
    c(16.75, 39.75, 89.25),  # Study 3 Case
    c(26.25, 65.5, 114.75)   # Study 3 Control
  ),
  n = c(40, 40, 40, 38, 15, 15)
)

print(study_data)
cat("\n")

#===============================================================================
#                    STEP 2: ASSESS SKEWNESS
#===============================================================================

cat("STEP 2: Assess Distributional Skewness\n")
cat(strrep("-", 80), "\n\n")

cat("Checking if normal distribution assumption is reasonable...\n\n")

skewness_results <- run_skewness_analysis(study_data)

cat("\nInterpretation:\n")
cat("  - T1 skewness uses min, median, max\n")
cat("  - |T1| < 0.2: Approximately symmetric (normal assumption OK)\n")
cat("  - |T1| > 0.5: Consider robust methods or transformations\n\n")

#===============================================================================
#                  STEP 3: BAYESIAN ESTIMATION
#===============================================================================

cat("STEP 3: Apply BOSE Method\n")
cat(strrep("-", 80), "\n\n")

cat("Estimating population mean and SD for each study...\n\n")

# Apply BOSE to all studies
results <- study_data %>%
  rowwise() %>%
  mutate(
    bose_result = list(estimate_bose(X = X_list, n = n, S = S))
  ) %>%
  ungroup() %>%
  mutate(
    mu_est = sapply(bose_result, function(r) r$mu_est),
    mu_lower = sapply(bose_result, function(r) r$mu_CI[1]),
    mu_upper = sapply(bose_result, function(r) r$mu_CI[2]),
    sigma_est = sapply(bose_result, function(r) r$sigma_est),
    sigma_lower = sapply(bose_result, function(r) r$sigma_CI[1]),
    sigma_upper = sapply(bose_result, function(r) r$sigma_CI[2])
  )

# Display results
results_table <- results %>%
  mutate(
    mu_fmt = sprintf("%.2f (%.2f, %.2f)", mu_est, mu_lower, mu_upper),
    sigma_fmt = sprintf("%.2f (%.2f, %.2f)", sigma_est, sigma_lower, sigma_upper)
  ) %>%
  select(Study, n, mu_fmt, sigma_fmt)

print(results_table)
cat("\n")

#===============================================================================
#                  STEP 4: COMPUTE EFFECT SIZES
#===============================================================================

cat("STEP 4: Calculate Mean Differences (MD) and Standardized Mean Differences (SMD)\n")
cat(strrep("-", 80), "\n\n")

# Reshape data for case-control comparison
effects <- results %>%
  mutate(
    Group = ifelse(grepl("Case", Study), "Case", "Control"),
    StudyID = gsub(" (Case|Control)", "", Study, fixed = FALSE)
  ) %>%
  select(StudyID, Group, mu_est, sigma_est, n) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(mu_est, sigma_est, n)
  ) %>%
  mutate(
    # Mean Difference (Control - Case)
    MD = mu_est_Control - mu_est_Case,
    
    # Pooled SD
    pooled_sd = sqrt(
      ((n_Case - 1) * sigma_est_Case^2 + (n_Control - 1) * sigma_est_Control^2) /
      (n_Case + n_Control - 2)
    ),
    
    # Standardized Mean Difference (Cohen's d)
    SMD = MD / pooled_sd
  )

print(effects %>% select(StudyID, MD, SMD, pooled_sd))
cat("\n")

cat("Interpretation:\n")
cat("  - MD > 0: Control group has higher mean than Case group\n")
cat("  - SMD (Cohen's d): |SMD| < 0.2 (small), 0.5 (medium), 0.8 (large)\n\n")

#===============================================================================
#                  STEP 5: COMPARE WITH OTHER METHODS
#===============================================================================

cat("STEP 5: Compare BOSE with Existing Methods\n")
cat(strrep("-", 80), "\n\n")

# Take Study 1 Case as example
example_study <- study_data[1, ]
X <- example_study$X_list[[1]]
n <- example_study$n

cat("Example: Study 1 Case\n")
cat("  Data: min =", X[1], ", median =", X[2], ", max =", X[3], "\n")
cat("  Sample size n =", n, "\n\n")

# BOSE
bose <- estimate_bose(X, n, S = 1)

# Luo-Wan
luo <- Luo.mean(X, n, type = "S1")
wan <- Wan.std(X, n, type = "S1")

# BLUE
blue <- blue.mean.sd(min.val = X[1], med.val = X[2], max.val = X[3], n = n)

# Comparison table
comparison <- tibble(
  Method = c("BOSE", "Luo-Wan", "BLUE"),
  Mean = c(bose$mu_est, luo$muhat, blue$est.mean),
  SD = c(bose$sigma_est, wan$sigmahat, blue$est.sd),
  `95% CI Available` = c("Yes", "No", "No")
)

print(comparison)
cat("\n")

cat("Key Advantage of BOSE:\n")
cat("  ✓ Full Bayesian uncertainty quantification\n")
cat("  ✓ Proper 95% credible intervals (maintains coverage)\n")
cat("  ✓ Type 1 quantiles avoid interpolation issues\n\n")

#===============================================================================
#                  STEP 6: VISUALIZATION (OPTIONAL)
#===============================================================================

cat("STEP 6: Create Visualizations (Optional)\n")
cat(strrep("-", 80), "\n\n")

cat("To create contour plots of posterior distributions:\n\n")

cat("  # Uncomment and run:\n")
cat("  # contour_results <- run_contour_analysis(study_data, n_studies = 3)\n")
cat("  # print(contour_results$final_grid)\n\n")

cat("To visualize simulation results (RMSE and coverage):\n\n")

cat("  # Load simulation results\n")
cat("  # results <- read.csv('simulations/results/simulation_results.csv')\n\n")

cat("  # Create RMSE plots\n")
cat("  # source('visualization/plot_rmse.R')\n")
cat("  # plots <- plot_rmse_by_scenario(results, dist = 'normal')\n")
cat("  # print(plots$S1)\n\n")

cat("  # Create coverage plots\n")
cat("  # source('visualization/plot_coverage.R')\n")
cat("  # p <- plot_coverage_by_scenario(results, dist = 'normal', case_num = 1)\n")
cat("  # print(p)\n\n")

#===============================================================================
#                          STEP 7: EXPORT RESULTS
#===============================================================================

cat("STEP 7: Export Results for Meta-Analysis\n")
cat(strrep("-", 80), "\n\n")

# Prepare data for standard meta-analysis packages (metafor, meta, etc.)
meta_data <- effects %>%
  mutate(
    # Standard error of MD (for fixed-effect meta-analysis)
    se_MD = sqrt(
      (sigma_est_Case^2 / n_Case) + (sigma_est_Control^2 / n_Control)
    ),
    
    # Standard error of SMD (approximate)
    se_SMD = sqrt(
      ((n_Case + n_Control) / (n_Case * n_Control)) +
      (SMD^2 / (2 * (n_Case + n_Control)))
    )
  ) %>%
  select(StudyID, MD, se_MD, SMD, se_SMD, n_Case, n_Control)

print(meta_data)
cat("\n")

cat("This data is now ready for:\n")
cat("  - Random-effects meta-analysis (metafor::rma())\n")
cat("  - Forest plots\n")
cat("  - Heterogeneity assessment (I², τ²)\n")
cat("  - Publication bias tests\n\n")

# Optional: Save to CSV
# write.csv(meta_data, "bose_meta_analysis_results.csv", row.names = FALSE)

#===============================================================================
#                            SUMMARY
#===============================================================================

cat(strrep("=", 80), "\n")
cat("                              WORKFLOW COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("What we did:\n")
cat("  ✓ Prepared case-control study data\n")
cat("  ✓ Assessed distributional assumptions (skewness)\n")
cat("  ✓ Applied BOSE method to estimate parameters\n")
cat("  ✓ Calculated effect sizes (MD and SMD)\n")
cat("  ✓ Compared with existing methods\n")
cat("  ✓ Prepared data for meta-analysis\n\n")

cat("Key findings:\n")
cat("  - BOSE provides full uncertainty quantification\n")
cat("  - Type 1 quantiles maintain proper 95% CI coverage\n")
cat("  - Ready for downstream meta-analysis\n\n")

cat("For more information:\n")
cat("  - See README.md for detailed documentation\n")
cat("  - Check case_studies/ for real data examples\n")
cat("  - Run simulations/ to validate the method\n\n")

cat(strrep("=", 80), "\n\n")
