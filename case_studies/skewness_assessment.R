#' Skewness Assessment Functions for BOSE Method
#' 
#' Functions to assess distributional skewness from order statistics.
#' Implements three skewness measures (T1, T2, T3) suitable for meta-analysis.
#' 
#' @author Wenqisi (Lydia) Pan
#' @references Shi, J., Luo, D., et al. (2020). Research Synthesis Methods.

# Required packages
# library(tidyverse)

#===============================================================================
#                        SKEWNESS CALCULATION FUNCTIONS
#===============================================================================

#' Calculate T1 Skewness (Three-Number Summary)
#'
#' Computes skewness from min, median, and max using the formula:
#' T1 = (min + max - 2*median) / (max - min)
#'
#' @param min_val Minimum value
#' @param median_val Median value
#' @param max_val Maximum value
#' @return Numeric skewness value (negative = left-skewed, positive = right-skewed)
#'
#' @details 
#' T1 measures asymmetry of the distribution. 
#' - T1 = 0: Symmetric distribution
#' - T1 > 0: Right-skewed (long tail on right)
#' - T1 < 0: Left-skewed (long tail on left)
#'
#' @examples
#' calculate_T1_skewness(10, 15, 25)  # Right-skewed: 0.333
#' calculate_T1_skewness(10, 17.5, 25)  # Symmetric: 0
#'
#' @export
calculate_T1_skewness <- function(min_val, median_val, max_val) {
  if (is.na(min_val) || is.na(median_val) || is.na(max_val)) {
    return(NA_real_)
  }
  if (max_val == min_val) {
    return(0)  # No variation, perfectly symmetric
  }
  T1 <- (max_val + min_val - 2 * median_val) / (max_val - min_val)
  return(T1)
}

#' Calculate T2 Skewness (Quartile-Based)
#'
#' Computes skewness from quartiles using the formula:
#' T2 = (Q1 + Q3 - 2*median) / (Q3 - Q1)
#'
#' @param q1 First quartile (25th percentile)
#' @param median_val Median (50th percentile)
#' @param q3 Third quartile (75th percentile)
#' @return Numeric skewness value
#'
#' @details 
#' T2 is less sensitive to extreme values than T1, focusing on the
#' middle 50% of the distribution (interquartile range).
#'
#' @examples
#' calculate_T2_skewness(12, 15, 20)  # Right-skewed: 0.125
#'
#' @export
calculate_T2_skewness <- function(q1, median_val, q3) {
  if (is.na(q1) || is.na(median_val) || is.na(q3)) {
    return(NA_real_)
  }
  if (q3 == q1) {
    return(0)  # IQR = 0, no variation in middle 50%
  }
  T2 <- (q3 + q1 - 2 * median_val) / (q3 - q1)
  return(T2)
}

#' Calculate T3 Skewness (Sample-Size Adjusted)
#'
#' Comprehensive skewness measure combining T1 and T2 with sample size adjustment:
#' T3 = max{[2.65 × ln(0.6n) / √n] × |T1|, |T2|}
#'
#' @param min_val Minimum value
#' @param max_val Maximum value
#' @param q1 First quartile
#' @param q3 Third quartile
#' @param median_val Median value
#' @param n Sample size
#' @return Numeric skewness value
#'
#' @details 
#' T3 accounts for sampling variability by adjusting T1 based on sample size.
#' For large samples, T1 is weighted more heavily; for small samples, T2 is
#' preferred due to its robustness.
#'
#' @references Shi, J., Luo, D., et al. (2020). Optimally estimating the sample
#'   standard deviation from the five-number summary.
#'
#' @export
calculate_T3_skewness <- function(min_val, max_val, q1, q3, median_val, n) {
  if (any(is.na(c(min_val, max_val, q1, q3, median_val, n)))) {
    return(NA_real_)
  }
  
  # Calculate T1 and T2
  T1 <- calculate_T1_skewness(min_val, median_val, max_val)
  T2 <- calculate_T2_skewness(q1, median_val, q3)
  
  if (is.na(T1) || is.na(T2)) {
    return(NA_real_)
  }
  
  # Sample size adjustment factor
  adjustment <- 2.65 * log(0.6 * n) / sqrt(n)
  
  # T3 = max{adjustment × |T1|, |T2|}
  T3 <- max(adjustment * abs(T1), abs(T2))
  
  return(T3)
}

#===============================================================================
#                        SKEWNESS CLASSIFICATION
#===============================================================================

#' Classify Skewness Magnitude
#'
#' Categorizes skewness into interpretable levels based on absolute value.
#'
#' @param skew_value Numeric skewness value
#' @return Character string classification
#'
#' @details Classification thresholds:
#' - |skewness| < 0.2: Approximately symmetric
#' - 0.2 ≤ |skewness| < 0.5: Slightly skewed
#' - 0.5 ≤ |skewness| < 0.8: Moderately skewed
#' - |skewness| ≥ 0.8: Highly skewed
#'
#' @examples
#' classify_skewness(0.1)   # "Approximately symmetric"
#' classify_skewness(0.6)   # "Moderately skewed"
#'
#' @export
classify_skewness <- function(skew_value) {
  if (is.na(skew_value)) {
    return("Data unavailable")
  }
  abs_skew <- abs(skew_value)
  
  if (abs_skew < 0.2) {
    return("Approximately symmetric")
  } else if (abs_skew < 0.5) {
    return("Slightly skewed")
  } else if (abs_skew < 0.8) {
    return("Moderately skewed")
  } else {
    return("Highly skewed")
  }
}

#' Get Skewness Direction
#'
#' Determines the direction of skewness (left vs right).
#'
#' @param skew_value Numeric skewness value
#' @return Character string: "Left-skewed", "Right-skewed", or "Symmetric"
#'
#' @export
get_skewness_direction <- function(skew_value) {
  if (is.na(skew_value)) {
    return("N/A")
  }
  if (skew_value > 0) {
    return("Right-skewed")
  } else if (skew_value < 0) {
    return("Left-skewed")
  } else {
    return("Symmetric")
  }
}

#===============================================================================
#                        UTILITY FUNCTIONS
#===============================================================================

#' Extract Quantiles from Order Statistics
#'
#' Maps observed order statistics to five-number summary based on scenario type.
#'
#' @param X Vector of observed order statistics
#' @param S Scenario type (1, 2, or 3)
#' @return List with elements: min_val, q1, median_val, q3, max_val
#'
#' @details Scenario mapping:
#' - S=1: X = {min, median, max}
#' - S=2: X = {Q1, median, Q3}
#' - S=3: X = {min, Q1, median, Q3, max}
#'
#' @export
extract_quantiles_from_X <- function(X, S) {
  quantiles <- list(
    min_val = NA_real_,
    q1 = NA_real_,
    median_val = NA_real_,
    q3 = NA_real_,
    max_val = NA_real_
  )
  
  if (S == 1) {
    # Three-number summary: {min, median, max}
    quantiles$min_val <- X[1]
    quantiles$median_val <- X[2]
    quantiles$max_val <- X[3]
  } else if (S == 2) {
    # Quartile summary: {Q1, median, Q3}
    quantiles$q1 <- X[1]
    quantiles$median_val <- X[2]
    quantiles$q3 <- X[3]
  } else if (S == 3) {
    # Five-number summary: {min, Q1, median, Q3, max}
    quantiles$min_val <- X[1]
    quantiles$q1 <- X[2]
    quantiles$median_val <- X[3]
    quantiles$q3 <- X[4]
    quantiles$max_val <- X[5]
  }
  
  return(quantiles)
}

#' Calculate All Skewness Metrics
#'
#' Comprehensive function that computes all available skewness measures
#' based on the observed order statistics.
#'
#' @param X Vector of observed order statistics
#' @param S Scenario type (1, 2, or 3)
#' @param n Sample size
#' @return Data frame with skewness metrics and classifications
#'
#' @details Returns:
#' - T1_skewness, T2_skewness, T3_skewness: Numeric skewness values
#' - T1_category, T2_category, T3_category: Classification strings
#' - T1_direction, T2_direction: Direction strings
#'
#' @examples
#' # Three-number summary
#' calculate_all_skewness(c(10, 15, 25), S = 1, n = 30)
#'
#' # Five-number summary
#' calculate_all_skewness(c(10, 12, 15, 18, 25), S = 3, n = 30)
#'
#' @export
calculate_all_skewness <- function(X, S, n) {
  # Extract quantiles based on scenario
  quants <- extract_quantiles_from_X(X, S)
  
  # Initialize skewness metrics
  T1 <- NA_real_
  T2 <- NA_real_
  T3 <- NA_real_
  
  # T1: Requires min, median, max
  if (!is.na(quants$min_val) && !is.na(quants$median_val) && !is.na(quants$max_val)) {
    T1 <- calculate_T1_skewness(quants$min_val, quants$median_val, quants$max_val)
  }
  
  # T2: Requires Q1, median, Q3
  if (!is.na(quants$q1) && !is.na(quants$median_val) && !is.na(quants$q3)) {
    T2 <- calculate_T2_skewness(quants$q1, quants$median_val, quants$q3)
  }
  
  # T3: Requires complete five-number summary
  if (!is.na(quants$min_val) && !is.na(quants$max_val) && 
      !is.na(quants$q1) && !is.na(quants$q3) && !is.na(quants$median_val)) {
    T3 <- calculate_T3_skewness(quants$min_val, quants$max_val, 
                                quants$q1, quants$q3, quants$median_val, n)
  }
  
  return(data.frame(
    T1_skewness = T1,
    T2_skewness = T2,
    T3_skewness = T3,
    T1_category = classify_skewness(T1),
    T2_category = classify_skewness(T2),
    T3_category = classify_skewness(T3),
    T1_direction = get_skewness_direction(T1),
    T2_direction = get_skewness_direction(T2)
  ))
}

#===============================================================================
#                        BATCH ANALYSIS FUNCTIONS
#===============================================================================

#' Run Skewness Analysis on Multiple Studies
#'
#' Analyzes skewness for a dataset of multiple studies with case-control design.
#'
#' @param study_data Data frame with columns: Study, S, X_list, n
#' @return List with results, summary statistics, and paired comparisons
#'
#' @details Input data frame should contain:
#' - Study: Study identifier (character)
#' - S: Scenario type (1, 2, or 3)
#' - X_list: List column of observed order statistics
#' - n: Sample size (integer)
#'
#' @examples
#' \dontrun{
#' study_data <- tibble(
#'   Study = c("Study 1 Case", "Study 1 Control"),
#'   S = c(1, 1),
#'   X_list = list(c(10, 15, 25), c(12, 18, 28)),
#'   n = c(30, 35)
#' )
#' results <- run_skewness_analysis(study_data)
#' }
#'
#' @export
run_skewness_analysis <- function(study_data) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("                    Skewness Analysis\n")
  cat(strrep("=", 70), "\n\n")
  
  # Calculate skewness for each study
  skewness_results <- study_data %>%
    filter(!is.na(X_list) & S > 0) %>%
    rowwise() %>%
    mutate(
      skewness_metrics = list(calculate_all_skewness(X_list, S, n))
    ) %>%
    unnest(skewness_metrics) %>%
    mutate(
      group_type = ifelse(str_detect(Study, "Case"), "Case", "Control"),
      study_num = str_extract(Study, "\\d+")
    ) %>%
    ungroup()
  
  # Print detailed results
  cat("=== Detailed Skewness Results ===\n\n")
  
  detailed_table <- skewness_results %>%
    select(Study, n, S, 
           T1_skewness, T1_category, T1_direction,
           T2_skewness, T2_category, T2_direction,
           T3_skewness, T3_category) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  
  print(detailed_table, n = Inf)
  
  # Summary by group
  cat("\n=== Summary by Group ===\n\n")
  
  summary_by_group <- skewness_results %>%
    group_by(group_type) %>%
    summarise(
      n_studies = n(),
      across(c(T1_skewness, T2_skewness, T3_skewness),
             list(
               mean = ~round(mean(., na.rm = TRUE), 4),
               sd = ~round(sd(., na.rm = TRUE), 4),
               min = ~round(min(., na.rm = TRUE), 4),
               max = ~round(max(., na.rm = TRUE), 4)
             ),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  print(summary_by_group)
  
  # Paired comparison
  cat("\n=== Paired Comparison by Study ===\n\n")
  
  paired_comparison <- skewness_results %>%
    group_by(study_num) %>%
    summarise(
      n_case = n[group_type == "Case"],
      n_control = n[group_type == "Control"],
      case_T1 = round(T1_skewness[group_type == "Case"], 4),
      control_T1 = round(T1_skewness[group_type == "Control"], 4),
      diff_T1 = round(case_T1 - control_T1, 4),
      .groups = "drop"
    )
  
  print(paired_comparison)
  
  cat("\n", strrep("=", 70), "\n\n")
  
  return(list(
    results = skewness_results,
    summary = summary_by_group,
    paired = paired_comparison
  ))
}
