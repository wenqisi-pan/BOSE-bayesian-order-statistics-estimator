#' Coverage Probability Visualization Functions for BOSE Method
#' 
#' Functions to create boxplots showing the distribution of coverage
#' probabilities across different sample sizes. Used to assess whether
#' 95% credible intervals maintain proper coverage.
#' 
#' @author Wenqisi (Lydia) Pan

# Required packages
library(ggplot2)
library(dplyr)
library(tidyr)

#===============================================================================
#                         UTILITY FUNCTIONS
#===============================================================================

#' Get Case Title by Distribution
#'
#' Creates formatted title strings for different parameter combinations.
#'
#' @param dist Distribution type: "normal", "t", or "lognormal"
#' @param case_num Case number (integer)
#' @return Character string for plot title
#' @export
get_case_title <- function(dist = "normal", case_num = 1) {
  if (dist == "normal") {
    case_titles <- c("N(5, 1²)", "N(50, 17²)")
  } else if (dist == "t") {
    case_titles <- c("df = 3", "df = 5", "df = 10", "df = 25")
  } else if (dist == "lognormal") {
    case_titles <- c("Log-Normal(5, 0.1²)", "Log-Normal(5, 0.2²)", "Log-Normal(5, 0.3²)")
  } else {
    stop("dist must be 'normal', 't', or 'lognormal'")
  }
  
  return(case_titles[case_num])
}

#===============================================================================
#                    COVERAGE PLOTTING FUNCTIONS
#===============================================================================

#' Plot Coverage Probability Boxplots by Scenario
#'
#' Creates boxplots showing the distribution of coverage probabilities
#' for mu and sigma across different scenarios (S1, S2, S3).
#'
#' @param df Data frame with simulation results containing:
#'   - S: Scenario type (1, 2, or 3)
#'   - case: Parameter combination index
#'   - n: Sample size
#'   - mu_coverage: Coverage probability for mu
#'   - sigma_coverage: Coverage probability for sigma
#' @param dist Distribution type: "normal", "t", or "lognormal"
#' @param case_num Which parameter combination to plot (integer)
#' @param n_type Optional filter for sample size type:
#'   - NULL: All sample sizes (default)
#'   - 1: n = 4q+1
#'   - 2: n = 4q+2
#'   - 3: n = 4q+3
#'   - 4: n = 4q (divisible by 4)
#' @return ggplot object with coverage boxplots
#'
#' @details 
#' Plot features:
#' - Black dashed line at 0.95 (target coverage)
#' - Red dashed lines showing mean coverage for each scenario
#' - Boxplots showing distribution across different sample sizes
#' - Separate panels for mu and sigma coverage
#'
#' @examples
#' \dontrun{
#' # Load simulation results
#' results <- read.csv("simulation_results.csv")
#' 
#' # Plot coverage for Normal(5,1), all sample sizes
#' p <- plot_coverage_by_scenario(results, dist = "normal", case_num = 1)
#' print(p)
#' 
#' # Plot coverage for n = 4q+1 only
#' p_4q1 <- plot_coverage_by_scenario(results, dist = "normal", 
#'                                     case_num = 1, n_type = 1)
#' print(p_4q1)
#' 
#' # Save plot
#' ggsave("coverage_boxplot.png", p, width = 14, height = 12, dpi = 600)
#' }
#'
#' @export
plot_coverage_by_scenario <- function(df, dist = "normal", case_num = 1, 
                                       n_type = NULL) {
  
  # Check required columns
  required_cols <- c("S", "case", "mu_coverage", "sigma_coverage", "n")
  if (!all(required_cols %in% names(df))) {
    stop(paste("Data frame must contain columns:", 
               paste(required_cols, collapse = ", ")))
  }
  
  # Filter by case first, then compute n_type
  df <- df %>%
    filter(case == case_num) %>%
    mutate(n_type = ifelse(n %% 4 == 0, 4, n %% 4))
  
  # Optionally filter by n_type
  if (!is.null(n_type)) {
    if (!n_type %in% c(1, 2, 3, 4)) {
      stop("n_type must be NULL (all data), 1, 2, 3, or 4")
    }
    
    df <- df %>% filter(n_type == !!n_type)
    
    # Check if any data remains after filtering
    if (nrow(df) == 0) {
      warning(paste0("No data found for case=", case_num, ", n=4q+", n_type))
      return(NULL)
    }
  }
  
  # Build plot title
  case_title <- get_case_title(dist, case_num)
  
  if (!is.null(n_type)) {
    n_type_label <- paste0("n = 4q+", n_type)
    case_title <- paste0(case_title, " | ", n_type_label)
  }
  
  # Prepare data for plotting
  data_long <- df %>%
    mutate(Scenario = paste0("S", S)) %>%
    select(Scenario, mu_coverage, sigma_coverage) %>%
    pivot_longer(
      cols = c(mu_coverage, sigma_coverage),
      names_to = "Parameter",
      values_to = "Coverage"
    ) %>%
    mutate(
      Parameter = case_when(
        Parameter == "mu_coverage" ~ "Coverage~(mu)",
        Parameter == "sigma_coverage" ~ "Coverage~(sigma)"
      ),
      Parameter = factor(Parameter, levels = c("Coverage~(mu)", "Coverage~(sigma)"))
    )
  
  # Calculate mean coverage for each group
  means_data <- data_long %>%
    group_by(Scenario, Parameter) %>%
    summarise(mean_coverage = mean(Coverage, na.rm = TRUE), .groups = "drop")
  
  # Create plot
  p <- ggplot(data_long, aes(x = Scenario, y = Coverage)) +
    # Boxplots
    geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0.5, 
                 fill = "grey80") +
    # Target coverage line (95%)
    geom_hline(yintercept = 0.95, linetype = "dashed", 
               color = "black", linewidth = 0.3) +
    # Mean coverage lines (red)
    geom_segment(data = means_data,
                 aes(x = as.numeric(factor(Scenario)) - 0.3,
                     xend = as.numeric(factor(Scenario)) + 0.3,
                     y = mean_coverage, yend = mean_coverage),
                 color = "red", linetype = "dashed", linewidth = 0.8) +
    # Faceting
    facet_wrap(~ Parameter,
               nrow = 2, scales = "free_y",
               strip.position = "left",
               labeller = label_parsed) +
    # Labels
    labs(x = "", y = "", title = case_title) +
    # Theming
    theme_minimal() +
    theme(
      strip.text = element_text(size = 28, face = "bold", angle = 90),
      axis.title = element_text(size = 28, color = "black"),      
      axis.text.x = element_text(size = 28, color = "black"),
      axis.text.y = element_text(size = 22, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      panel.spacing = unit(2, "lines"),
      strip.placement = "outside",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      axis.text.y.left = element_text(margin = margin(r = 5)),
      panel.spacing.y = unit(1, "lines")
    ) +
    # Y-axis settings
    coord_cartesian(ylim = c(0.7, 1)) +
    scale_y_continuous(breaks = seq(0.7, 1, by = 0.05))
  
  return(p)
}

#' Plot Coverage Probability Boxplots by n Type
#'
#' Creates boxplots showing coverage probabilities grouped by sample size type
#' (n = 4q+1, 4q+2, 4q+3, 4q) for a specific scenario.
#'
#' @param df Data frame with simulation results
#' @param dist Distribution type: "normal", "t", or "lognormal"
#' @param case_num Which parameter combination to plot
#' @param scenario_num Which scenario to plot (1, 2, or 3)
#' @return ggplot object
#'
#' @details 
#' Useful for examining Type 1 vs Type 7 quantile differences.
#' Type 7 quantiles show under-coverage for certain n types, while
#' Type 1 quantiles maintain proper coverage across all n types.
#'
#' @examples
#' \dontrun{
#' # Compare coverage across n types for S1
#' p <- plot_coverage_by_ntype(results, dist = "normal", 
#'                              case_num = 1, scenario_num = 1)
#' print(p)
#' }
#'
#' @export
plot_coverage_by_ntype <- function(df, dist = "normal", case_num = 1, 
                                    scenario_num = 1) {
  
  # Check required columns
  required_cols <- c("S", "case", "mu_coverage", "sigma_coverage", "n")
  if (!all(required_cols %in% names(df))) {
    stop(paste("Data frame must contain columns:", 
               paste(required_cols, collapse = ", ")))
  }
  
  # Filter by case and scenario
  df <- df %>%
    filter(case == case_num, S == scenario_num) %>%
    mutate(n_type = ifelse(n %% 4 == 0, 4, n %% 4))
  
  # Check if data exists
  if (nrow(df) == 0) {
    warning(paste0("No data found for case=", case_num, ", S=", scenario_num))
    return(NULL)
  }
  
  # Build plot title
  case_title <- get_case_title(dist, case_num)
  scenario_label <- paste0("S", scenario_num)
  plot_title <- paste0(case_title, " | ", scenario_label)
  
  # Prepare data for plotting
  data_long <- df %>%
    mutate(n_type_label = paste0("n = 4q+", n_type)) %>%
    select(n_type_label, mu_coverage, sigma_coverage) %>%
    pivot_longer(
      cols = c(mu_coverage, sigma_coverage),
      names_to = "Parameter",
      values_to = "Coverage"
    ) %>%
    mutate(
      Parameter = case_when(
        Parameter == "mu_coverage" ~ "Coverage~(mu)",
        Parameter == "sigma_coverage" ~ "Coverage~(sigma)"
      ),
      Parameter = factor(Parameter, levels = c("Coverage~(mu)", "Coverage~(sigma)")),
      n_type_label = factor(n_type_label, 
                            levels = c("n = 4q+1", "n = 4q+2", "n = 4q+3", "n = 4q+4"))
    )
  
  # Calculate mean coverage
  means_data <- data_long %>%
    group_by(n_type_label, Parameter) %>%
    summarise(mean_coverage = mean(Coverage, na.rm = TRUE), .groups = "drop")
  
  # Create plot
  p <- ggplot(data_long, aes(x = n_type_label, y = Coverage)) +
    # Boxplots
    geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0.5, 
                 fill = "grey80") +
    # Target coverage line
    geom_hline(yintercept = 0.95, linetype = "dashed", 
               color = "black", linewidth = 0.3) +
    # Mean coverage lines
    geom_segment(data = means_data,
                 aes(x = as.numeric(factor(n_type_label)) - 0.3,
                     xend = as.numeric(factor(n_type_label)) + 0.3,
                     y = mean_coverage, yend = mean_coverage),
                 color = "red", linetype = "dashed", linewidth = 0.8) +
    # Faceting
    facet_wrap(~ Parameter,
               nrow = 2, scales = "free_y",
               strip.position = "left",
               labeller = label_parsed) +
    # Labels
    labs(x = "", y = "", title = plot_title) +
    # Theming
    theme_minimal() +
    theme(
      strip.text = element_text(size = 28, face = "bold", angle = 90),
      axis.title = element_text(size = 28, color = "black"),      
      axis.text.x = element_text(size = 28, color = "black"),
      axis.text.y = element_text(size = 22, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      panel.spacing = unit(2, "lines"),
      strip.placement = "outside",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      axis.text.y.left = element_text(margin = margin(r = 5)),
      panel.spacing.y = unit(1, "lines")
    ) +
    # Y-axis settings
    coord_cartesian(ylim = c(0.7, 1)) +
    scale_y_continuous(breaks = seq(0.7, 1, by = 0.05))
  
  return(p)
}

#===============================================================================
#                         BATCH PLOTTING FUNCTION
#===============================================================================

#' Create All Coverage Plots for a Dataset
#'
#' Convenience function that generates all standard coverage plots for
#' a simulation study.
#'
#' @param df Data frame with simulation results
#' @param dist Distribution type
#' @param output_dir Directory to save plots (optional)
#' @return List of plots
#'
#' @details Creates:
#' - Coverage by scenario plots for each case
#' - Coverage by n_type plots for each case and scenario combination
#'
#' @examples
#' \dontrun{
#' results <- read.csv("simulation_results.csv")
#' all_plots <- create_all_coverage_plots(results, dist = "normal",
#'                                         output_dir = "./plots")
#' }
#'
#' @export
create_all_coverage_plots <- function(df, dist = "normal", output_dir = NULL) {
  
  all_plots <- list()
  
  # Determine number of cases based on distribution
  n_cases <- if (dist == "normal") 2 else if (dist == "t") 4 else 3
  
  # Coverage by scenario (for each case)
  for (case_num in 1:n_cases) {
    plot_name <- paste0("coverage_scenario_case", case_num)
    all_plots[[plot_name]] <- plot_coverage_by_scenario(df, dist, case_num)
    
    # Optionally save
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      filename <- file.path(output_dir, paste0(plot_name, ".png"))
      ggsave(filename, all_plots[[plot_name]], 
             width = 14, height = 12, dpi = 600)
      cat("Saved:", filename, "\n")
    }
  }
  
  # Coverage by n_type (for each case and scenario)
  for (case_num in 1:n_cases) {
    for (scenario_num in 1:3) {
      plot_name <- paste0("coverage_ntype_case", case_num, "_S", scenario_num)
      all_plots[[plot_name]] <- plot_coverage_by_ntype(df, dist, 
                                                        case_num, scenario_num)
      
      # Optionally save
      if (!is.null(output_dir)) {
        filename <- file.path(output_dir, paste0(plot_name, ".png"))
        ggsave(filename, all_plots[[plot_name]], 
               width = 14, height = 12, dpi = 600)
        cat("Saved:", filename, "\n")
      }
    }
  }
  
  return(all_plots)
}
