#' RMSE Visualization Functions for BOSE Method
#' 
#' Functions to create publication-quality RMSE comparison plots for
#' simulation study results. Compares BOSE with existing methods (Luo, Wan,
#' BC, QE, MLN, BLUE, Shi).
#' 
#' @author Wenqisi (Lydia) Pan

# Required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)

#===============================================================================
#                         CONFIGURATION & STYLING
#===============================================================================

#' Get Method Colors
#'
#' Standard color palette for method comparison plots.
#'
#' @return Named vector of colors
#' @export
get_method_colors <- function() {
  c(
    "LW" = "red",
    "BC" = "orange",
    "QE" = "green",
    "MLN" = "darkgrey",
    "BLUE" = "blue",
    "BOSE" = "black",
    "Shi" = "pink"
  )
}

#' Get Method Shapes
#'
#' Standard shape palette for method comparison plots.
#'
#' @return Named vector of shape codes
#' @export
get_method_shapes <- function() {
  c(
    "LW" = 16,     # Filled circle
    "BC" = 17,     # Filled triangle
    "QE" = 15,     # Filled square
    "MLN" = 18,    # Filled diamond
    "BLUE" = 4,    # Cross
    "BOSE" = 19,   # Large filled circle
    "Shi" = 8      # Star
  )
}

#' Get RMSE Columns by Distribution
#'
#' Determines which RMSE columns to plot based on distribution type.
#'
#' @param dist Distribution type: "normal", "t", or "lognormal"
#' @return Character vector of column names
#' @export
get_rmse_columns <- function(dist = "normal") {
  if (dist == "normal") {
    c("mu_rmse_Blue", "mu_rmse_Luo", "mu_rmse_Bayes",
      "sigma_rmse_Bayes", "sigma_rmse_Wan", "sigma_rmse_Blue", "sigma_rmse_Shi")
  } else if (dist == "t") {
    c("mu_rmse_Blue", "mu_rmse_Luo", "mu_rmse_Bayes",
      "sigma_rmse_Bayes", "sigma_rmse_Wan", "sigma_rmse_Blue", "sigma_rmse_Shi")
  } else if (dist == "lognormal") {
    c("mu_rmse_BC", "mu_rmse_MLN", "mu_rmse_Blue", "mu_rmse_Luo", "mu_rmse_Bayes",
      "sigma_rmse_Bayes", "sigma_rmse_Wan", "sigma_rmse_BC", "sigma_rmse_MLN",
      "sigma_rmse_Blue", "sigma_rmse_Shi")
  } else {
    stop("dist must be 'normal', 't', or 'lognormal'")
  }
}

#' Get Case Labels by Distribution
#'
#' Creates formatted labels for different parameter combinations.
#'
#' @param dist Distribution type
#' @return Character vector of plot-ready labels (with expression syntax)
#' @export
get_case_labels <- function(dist = "normal") {
  if (dist == "normal") {
    c("paste('N(5, 1'^2, ')')", 
      "paste('N(50, 17'^2, ')')")
  } else if (dist == "t") {
    c("paste(italic(df), '= 3')",
      "paste(italic(df), '= 5')",
      "paste(italic(df), '= 10')",
      "paste(italic(df), '= 25')")
  } else if (dist == "lognormal") {
    c("paste('Log-Normal(5, 0.1'^2, ')')",
      "paste('Log-Normal(5, 0.2'^2, ')')",
      "paste('Log-Normal(5, 0.3'^2, ')')")
  }
}

#===============================================================================
#                    DATA PREPARATION FUNCTIONS
#===============================================================================

#' Prepare RMSE Data for Plotting
#'
#' Reshapes simulation results into long format and adds method labels.
#'
#' @param df Data frame with simulation results
#' @param dist Distribution type
#' @return Data frame in long format ready for plotting
#' @export
prepare_rmse_data <- function(df, dist = "normal") {
  
  rmse_cols <- get_rmse_columns(dist)
  case_levels <- get_case_labels(dist)
  
  long_data <- df %>%
    mutate(Scenario = paste0("S", S)) %>%
    pivot_longer(cols = all_of(rmse_cols), 
                 names_to = "Estimator", 
                 values_to = "RMSE_Value") %>%
    mutate(
      # Determine parameter type (mu vs sigma)
      Parameter = ifelse(grepl("mu_", Estimator), 
                         "RMSE~(hat(mu))", 
                         "RMSE~(hat(sigma))"),
      
      # Map estimator to method name
      Method = case_when(
        Estimator %in% c("mu_rmse_Luo", "sigma_rmse_Wan") ~ "LW",
        Estimator %in% c("mu_rmse_BC", "sigma_rmse_BC") ~ "BC",
        Estimator %in% c("mu_rmse_QE", "sigma_rmse_QE") ~ "QE",
        Estimator %in% c("mu_rmse_MLN", "sigma_rmse_MLN") ~ "MLN",
        Estimator %in% c("mu_rmse_Blue", "sigma_rmse_Blue") ~ "BLUE",
        Estimator %in% c("mu_rmse_Bayes", "sigma_rmse_Bayes") ~ "BOSE",
        Estimator %in% c("sigma_rmse_Shi") ~ "Shi"
      )
    ) %>%
    # Shi's method only applies to S=3
    filter(!(Method == "Shi" & S != 3)) %>%
    # Make BOSE the reference level (appears first in legend)
    mutate(Method = fct_relevel(Method, "BOSE"))
  
  # Add formatted case labels
  long_data <- long_data %>%
    mutate(case_label = case_levels[case]) %>%
    mutate(case_label = factor(case_label,
                               levels = case_levels,
                               ordered = TRUE))
  
  return(long_data)
}

#===============================================================================
#                    MAIN PLOTTING FUNCTIONS
#===============================================================================

#' Plot RMSE Comparison Split by Scenario
#'
#' Creates stacked plots (mu on top, sigma on bottom) for each scenario (S1, S2, S3).
#' Returns a list of plots, one for each scenario.
#'
#' @param df Data frame with simulation results containing:
#'   - S: Scenario type (1, 2, or 3)
#'   - case: Parameter combination index
#'   - n: Sample size
#'   - mu_rmse_*, sigma_rmse_*: RMSE values for different methods
#' @param dist Distribution type: "normal", "t", or "lognormal"
#' @return List of ggplot objects, one per scenario (named "S1", "S2", "S3")
#'
#' @details Creates publication-quality plots with:
#' - Top panel: RMSE for mu estimates
#' - Bottom panel: RMSE for sigma estimates
#' - Horizontal reference line at RMSE = 1 (sample mean/SD performance)
#' - BOSE highlighted with thicker lines
#'
#' @examples
#' \dontrun{
#' # Load simulation results
#' results <- read.csv("simulation_results.csv")
#' 
#' # Create plots for Normal distribution
#' plots <- plot_rmse_by_scenario(results, dist = "normal")
#' 
#' # Display S1 plot
#' print(plots$S1)
#' 
#' # Save all plots
#' ggsave("S1_rmse.png", plots$S1, width = 24, height = 12, dpi = 600)
#' ggsave("S2_rmse.png", plots$S2, width = 24, height = 12, dpi = 600)
#' ggsave("S3_rmse.png", plots$S3, width = 24, height = 12, dpi = 600)
#' }
#'
#' @export
plot_rmse_by_scenario <- function(df, dist = "normal") {
  
  # Prepare data
  long_data <- prepare_rmse_data(df, dist)
  
  # Get styling
  method_colors <- get_method_colors()
  method_shapes <- get_method_shapes()
  
  # Create one plot per scenario
  scenario_list <- unique(long_data$Scenario)
  plots <- list()
  
  for (sc in scenario_list) {
    data_s <- filter(long_data, Scenario == sc)
    
    mu_data <- filter(data_s, Parameter == "RMSE~(hat(mu))")
    sigma_data <- filter(data_s, Parameter == "RMSE~(hat(sigma))")
    
    # Top panel: mu RMSE
    p_mu <- ggplot(mu_data, aes(x = n, y = RMSE_Value, 
                                 color = Method, shape = Method, 
                                 group = Estimator)) +
      geom_line(linewidth = 0.3) +
      geom_point(size = 2.5) +
      # Emphasize BOSE with thicker line
      geom_line(data = filter(mu_data, Method == "BOSE"), linewidth = 0.4) +
      geom_point(data = filter(mu_data, Method == "BOSE"), size = 2.5) +
      # Reference line at 1 (sample mean performance)
      geom_hline(yintercept = 1, linetype = "dashed", 
                 color = "black", linewidth = 0.3) +
      labs(x = NULL, 
           y = expression(RMSE(hat(mu))), 
           color = "Methods",
           shape = "Methods") +
      facet_grid(. ~ case_label, labeller = label_parsed) +
      scale_color_manual(values = method_colors) +
      scale_shape_manual(values = method_shapes) +
      theme_bw() +
      theme(
        strip.text.x = element_text(size = 42, face = "bold"),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, size = 28),
        axis.title.x = element_text(size = 28),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 2.0),
        legend.position = "none"
      )
    
    # Bottom panel: sigma RMSE
    p_sigma <- ggplot(sigma_data, aes(x = n, y = RMSE_Value, 
                                       color = Method, shape = Method, 
                                       group = Estimator)) +
      geom_line(linewidth = 0.3) +
      geom_point(size = 2.5) +
      # Emphasize BOSE
      geom_line(data = filter(sigma_data, Method == "BOSE"), linewidth = 0.4) +
      geom_point(data = filter(sigma_data, Method == "BOSE"), size = 2.5) +
      # Reference line
      geom_hline(yintercept = 1, linetype = "dashed", 
                 color = "black", linewidth = 0.3) +
      labs(x = "Sample Size (n)", 
           y = expression(RMSE(hat(sigma))), 
           color = "Methods",
           shape = "Methods") +
      facet_grid(. ~ case_label, labeller = label_parsed) +
      scale_color_manual(values = method_colors) +
      scale_shape_manual(values = method_shapes) +
      theme_bw() +
      theme(
        strip.text = element_blank(),  # Hide facet labels on second row
        strip.background = element_rect(fill = "lightgray"),
        axis.title.y = element_text(angle = 90, size = 28),
        axis.title.x = element_text(size = 28),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 2.0),
        legend.position = "right",
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26)
      )
    
    # Combine plots (mu on top, sigma on bottom)
    plots[[sc]] <- (p_mu / p_sigma) + 
      plot_layout(guides = "collect") +
      plot_annotation(theme = theme(plot.title = element_text(size = 38, 
                                                               face = "bold", 
                                                               hjust = 0.5)))
  }
  
  return(plots)
}

#' Plot RMSE Comparison Split by Case
#'
#' Creates stacked plots for each parameter combination (case).
#' Useful when you want to compare scenarios within each parameter setting.
#'
#' @param df Data frame with simulation results
#' @param dist Distribution type: "normal", "t", or "lognormal"
#' @return List of ggplot objects, one per case (named "Case1", "Case2", etc.)
#'
#' @details Similar to plot_rmse_by_scenario() but organizes plots by
#'   parameter combination instead of scenario type.
#'
#' @export
plot_rmse_by_case <- function(df, dist = "normal") {
  
  # Prepare data
  long_data <- prepare_rmse_data(df, dist)
  
  # Get styling
  method_colors <- get_method_colors()
  method_shapes <- get_method_shapes()
  
  # Create one plot per case
  case_list <- unique(long_data$case)
  plots <- list()
  
  for (cs in case_list) {
    data_c <- filter(long_data, case == cs)
    
    mu_data <- filter(data_c, Parameter == "RMSE~(hat(mu))")
    sigma_data <- filter(data_c, Parameter == "RMSE~(hat(sigma))")
    
    # Get case label for title
    case_label_text <- unique(data_c$case_label)[1]
    
    # Top panel: mu RMSE
    p_mu <- ggplot(mu_data, aes(x = n, y = RMSE_Value, 
                                 color = Method, shape = Method, 
                                 group = Estimator)) +
      geom_line(linewidth = 0.3) +
      geom_point(size = 2.5) +
      geom_line(data = filter(mu_data, Method == "BOSE"), linewidth = 0.4) +
      geom_point(data = filter(mu_data, Method == "BOSE"), size = 2.5) +
      geom_hline(yintercept = 1, linetype = "dashed", 
                 color = "black", linewidth = 0.3) +
      labs(x = NULL, 
           y = expression(RMSE(hat(mu))), 
           color = "Methods",
           shape = "Methods") +
      facet_grid(. ~ Scenario) +
      scale_color_manual(values = method_colors) +
      scale_shape_manual(values = method_shapes) +
      theme_bw() +
      theme(
        strip.text.x = element_text(size = 42, face = "bold"),
        strip.background = element_blank(),
        axis.title.y = element_text(angle = 90, size = 28),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 2.0),
        legend.position = "none"
      )
    
    # Bottom panel: sigma RMSE
    p_sigma <- ggplot(sigma_data, aes(x = n, y = RMSE_Value, 
                                       color = Method, shape = Method, 
                                       group = Estimator)) +
      geom_line(linewidth = 0.3) +
      geom_point(size = 2.5) +
      geom_line(data = filter(sigma_data, Method == "BOSE"), linewidth = 0.4) +
      geom_point(data = filter(sigma_data, Method == "BOSE"), size = 2.5) +
      geom_hline(yintercept = 1, linetype = "dashed", 
                 color = "black", linewidth = 0.3) +
      labs(x = "Sample Size (n)", 
           y = expression(RMSE(hat(sigma))), 
           color = "Methods",
           shape = "Methods") +
      facet_grid(. ~ Scenario) +
      scale_color_manual(values = method_colors) +
      scale_shape_manual(values = method_shapes) +
      theme_bw() +
      theme(
        strip.text = element_blank(),
        strip.background = element_rect(fill = "lightgray"),
        axis.title.y = element_text(angle = 90, size = 28),
        axis.title.x = element_text(size = 28),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 2.0),
        legend.position = "right",
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26)
      )
    
    # Combine plots
    plots[[paste0("Case", cs)]] <- (p_mu / p_sigma) + 
      plot_layout(guides = "collect")
  }
  
  return(plots)
}
