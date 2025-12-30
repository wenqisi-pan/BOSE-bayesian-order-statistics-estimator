#' BOSE Quick Start Example
#' 
#' This is a 5-minute introduction to the BOSE method.
#' Shows the absolute minimum code needed to estimate mean and SD
#' from order statistics (quantiles).
#' 
#' @author Wenqisi (Lydia) Pan

#===============================================================================
#                         SETUP
#===============================================================================

# Load core BOSE functions
source("R/bose_core.R")
source("R/helper_functions.R")

# Set random seed for reproducibility
set.seed(123)

#===============================================================================
#                    EXAMPLE 1: THREE-NUMBER SUMMARY
#===============================================================================

cat("\n", strrep("=", 70), "\n")
cat("  EXAMPLE 1: Estimate Mean and SD from Three-Number Summary\n")
cat(strrep("=", 70), "\n\n")

# Suppose a study reports: min = 45, median = 52, max = 68, n = 50
# We want to estimate the population mean (mu) and SD (sigma)

X <- c(45, 52, 68)  # Observed quantiles: {min, median, max}
n <- 50              # Sample size
S <- 1               # Scenario 1: three-number summary

# Run BOSE estimation
result <- estimate_bose(X, n, S)

# Display results
cat("Observed data:\n")
cat("  Min = ", X[1], ", Median = ", X[2], ", Max = ", X[3], "\n")
cat("  Sample size n = ", n, "\n\n")

cat("BOSE Estimates:\n")
cat("  Mean (mu):       ", sprintf("%.2f", result$mu_est), "\n")
cat("  95% CI for mu:   (", sprintf("%.2f", result$mu_CI[1]), ", ", 
    sprintf("%.2f", result$mu_CI[2]), ")\n\n", sep = "")

cat("  SD (sigma):      ", sprintf("%.2f", result$sigma_est), "\n")
cat("  95% CI for sigma: (", sprintf("%.2f", result$sigma_CI[1]), ", ", 
    sprintf("%.2f", result$sigma_CI[2]), ")\n\n", sep = "")

#===============================================================================
#                    EXAMPLE 2: FIVE-NUMBER SUMMARY
#===============================================================================

cat(strrep("=", 70), "\n")
cat("  EXAMPLE 2: Estimate from Five-Number Summary\n")
cat(strrep("=", 70), "\n\n")

# Study reports: min = 10, Q1 = 15, median = 20, Q3 = 28, max = 45, n = 100
X <- c(10, 15, 20, 28, 45)
n <- 100
S <- 3  # Scenario 3: five-number summary

result <- estimate_bose(X, n, S)

cat("Observed data:\n")
cat("  Five-number summary: [", paste(X, collapse = ", "), "]\n")
cat("  Sample size n = ", n, "\n\n")

cat("BOSE Estimates:\n")
cat("  Mean (mu):   ", sprintf("%.2f (%.2f, %.2f)", 
    result$mu_est, result$mu_CI[1], result$mu_CI[2]), "\n")
cat("  SD (sigma):  ", sprintf("%.2f (%.2f, %.2f)", 
    result$sigma_est, result$sigma_CI[1], result$sigma_CI[2]), "\n\n")

#===============================================================================
#                    EXAMPLE 3: COMPARE WITH OTHER METHODS
#===============================================================================

cat(strrep("=", 70), "\n")
cat("  EXAMPLE 3: Compare BOSE with Existing Methods\n")
cat(strrep("=", 70), "\n\n")

# Install required packages if needed
# install.packages(c("estmeansd", "metaBLUE"))

library(estmeansd)
library(metaBLUE)

# Use the three-number summary from Example 1
X <- c(45, 52, 68)
n <- 50
S <- 1

# BOSE estimate
bose_result <- estimate_bose(X, n, S)

# Luo's method (2018) for mean
luo_result <- Luo.mean(X, n, type = "S1")

# Wan's method (2014) for SD
wan_result <- Wan.std(X, n, type = "S1")

# BLUE method
blue_result <- blue.mean.sd(min.val = X[1], med.val = X[2], 
                             max.val = X[3], n = n)

# Compare results
cat("Method Comparison:\n")
cat(sprintf("  %-15s  Mean: %6.2f   SD: %6.2f\n", 
            "BOSE", bose_result$mu_est, bose_result$sigma_est))
cat(sprintf("  %-15s  Mean: %6.2f   SD: %6.2f\n", 
            "Luo-Wan", luo_result$muhat, wan_result$sigmahat))
cat(sprintf("  %-15s  Mean: %6.2f   SD: %6.2f\n", 
            "BLUE", blue_result$est.mean, blue_result$est.sd))
cat("\n")

cat("Key Difference:\n")
cat("  BOSE provides full uncertainty quantification with 95% credible intervals,\n")
cat("  while other methods only provide point estimates.\n\n")

#===============================================================================
#                         WHAT'S NEXT?
#===============================================================================

cat(strrep("=", 70), "\n")
cat("  Next Steps:\n")
cat(strrep("=", 70), "\n\n")

cat("1. Run simulation study:\n")
cat("   source('simulations/run_simulation_type1.R')\n\n")

cat("2. Analyze real data:\n")
cat("   source('case_studies/real_data_analysis.R')\n\n")

cat("3. Create visualization:\n")
cat("   source('visualization/plot_rmse.R')\n")
cat("   # Then load your simulation results and plot\n\n")

cat("4. Assess skewness:\n")
cat("   source('case_studies/skewness_assessment.R')\n\n")

cat("See inst/examples/full_workflow.R for a complete analysis example.\n\n")

cat(strrep("=", 70), "\n")
cat("  Quick Start Complete!\n")
cat(strrep("=", 70), "\n\n")
