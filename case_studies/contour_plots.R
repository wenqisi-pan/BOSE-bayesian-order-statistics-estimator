#' Contour Plot Visualization for BOSE Method
#' 
#' Functions to create publication-quality contour plots showing the
#' joint posterior distribution of (mu, sigma) for real data applications.
#' 
#' @author Wenqisi (Lydia) Pan

# Required packages
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)

# Source core functions if needed
# source("R/bose_core.R")
# source("case_studies/real_data_analysis.R")

#===============================================================================
#                         GLOBAL PARAMETERS
#===============================================================================

n_show <- 50               # Number of posterior samples to display on plot
grid_size_contour <- 200   # Grid size for contour plots (smaller for speed)

#===============================================================================
#                         CONTOUR PLOT CREATION
#===============================================================================

#' Create Contour Plot of Joint Posterior
#'
#' Generates a contour plot showing the posterior distribution of (mu, sigma)
#' with credible intervals, point estimates, and sample points.
#'
#' @param result_list List returned from estimate_bayes_real_data() with
#'   return_contour_data = TRUE
#' @param study_name Character string for plot title (e.g., "Study 1")
#' @param group_type Character string for group (e.g., "Case", "Control")
#' @param xlim Numeric vector of length 2 for x-axis limits (optional)
#' @param ylim Numeric vector of length 2 for y-axis limits (optional)
#' @return ggplot object
#'
#' @details Plot elements:
#' - Filled contours: Posterior density levels (10%, 25%, 50%, 75%, 90%)
#' - Black × points: Random posterior samples (n=50)
#' - Red × point: Posterior mean (mu, sigma)
#' - Blue triangle: Initial estimates (Luo, Wan)
#' - Red dashed lines: 95% credible intervals
#'
#' @examples
#' \dontrun{
#' result <- estimate_bayes_real_data(S=1, X=c(10,15,25), n=30, 
#'                                    return_contour_data=TRUE)
#' p <- create_contour_plot(result, "Study 1", "Case")
#' print(p)
#' }
#'
#' @export
create_contour_plot <- function(result_list, study_name = "", group_type = "", 
                                xlim = NULL, ylim = NULL) {
  if (is.null(result_list)) {
    return(NULL)
  }
  
  # Prepare contour data
  contour_data <- expand.grid(
    mu = result_list$mu_grid,
    sigma = result_list$sigma_grid
  )
  contour_data$density <- as.vector(result_list$posterior_density)
  
  # Define contour levels
  max_density <- max(contour_data$density)
  breaks <- c(0.1, 0.25, 0.5, 0.75, 0.9) * max_density
  
  # Sample posterior points for display
  sample_points <- data.frame(
    mu = result_list$mu_samples,
    sigma = result_list$sigma_samples
  ) %>% dplyr::slice_sample(n = n_show)
  
  # Build subtitle with Type 1 index information
  n_val <- as.integer(str_extract(result_list$study_info, "\\d+"))
  k_med <- ceiling(0.5 * n_val)
  
  subtitle_text <- paste0(
    "Data: [", paste(round(result_list$X_observed, 1), collapse = ", "), 
    "] | n = ", n_val,
    " | k_med = ", k_med, " (Type 1)"
  )
  
  # Create plot
  p <- ggplot(contour_data, aes(x = mu, y = sigma)) +
    # Filled contours (posterior density)
    geom_contour_filled(aes(z = density), 
                        breaks = breaks,
                        alpha = 0.6) +
    # Contour lines (for clarity)
    geom_contour(aes(z = density), 
                 breaks = breaks,
                 color = "black", 
                 linewidth = 0.6,
                 alpha = 0.8) +
    # Posterior samples (black × points)
    geom_point(data = sample_points, aes(x = mu, y = sigma),
               color = "black", size = 1.5, shape = 4, alpha = 0.8) +
    # Posterior mean (red × point)
    geom_point(aes(x = result_list$mu_mean, y = result_list$sigma_mean),
               color = "#D32F2F", size = 4, shape = 4, stroke = 2.5) +
    # Initial estimates (blue triangle)
    geom_point(aes(x = result_list$initial_estimates$mu_luo, 
                   y = result_list$initial_estimates$sigma_wan),
               color = "#1976D2", size = 3, shape = 17) +
    # 95% credible intervals (red dashed lines)
    geom_vline(xintercept = result_list$mu_CI, 
               color = "#D32F2F", linetype = "dashed", 
               linewidth = 1, alpha = 0.8) +
    geom_hline(yintercept = result_list$sigma_CI, 
               color = "#D32F2F", linetype = "dashed", 
               linewidth = 1, alpha = 0.8) +
    # Styling
    scale_fill_grey(start = 0.9, end = 0.3, name = "Density") +
    labs(
      title = paste(study_name, group_type),
      subtitle = subtitle_text,
      x = expression(mu),
      y = expression(sigma)
    ) +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
    )
  
  # Apply axis limits if provided
  if (!is.null(xlim)) {
    p <- p + scale_x_continuous(limits = xlim, expand = c(0.02, 0))
  }
  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim, expand = c(0.02, 0))
  }
  
  return(p)
}

#===============================================================================
#                    BATCH CONTOUR PLOT GENERATION
#===============================================================================

#' Run Contour Plot Analysis for Multiple Studies
#'
#' Generates contour plots for all studies in a dataset with matched axis
#' ranges for case-control pairs.
#'
#' @param study_data Data frame with study information (default: study_inputs)
#' @param n_studies Number of study pairs to plot (default: 3)
#' @return List with results, individual plots, study ranges, and final grid
#'
#' @details 
#' For case-control studies, this function:
#' 1. Computes posterior densities for all studies
#' 2. Determines common axis ranges for each case-control pair
#' 3. Creates matched contour plots with consistent scales
#' 4. Arranges plots in a grid (cases on top row, controls on bottom)
#'
#' @examples
#' \dontrun{
#' contour_results <- run_contour_analysis()
#' print(contour_results$final_grid)
#' }
#'
#' @export
run_contour_analysis <- function(study_data = study_inputs, n_studies = 3) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("              Type 1 Contour Plot Analysis\n")
  cat(strrep("=", 60), "\n\n")
  
  start_time <- Sys.time()
  
  results_list <- list()
  plots_list <- list()
  
  # Step 1: Collect results for all studies
  for (i in 1:(2 * n_studies)) {
    study_row <- study_data[i, ]
    
    study_parts <- str_split(study_row$Study, " ")[[1]]
    study_name <- paste(study_parts[1], study_parts[2])
    group_type <- study_parts[3]
    
    if (study_row$S > 0 && !is.na(study_row$X_list[[1]][1])) {
      result <- estimate_bayes_real_data(
        S = study_row$S, 
        X = study_row$X_list[[1]], 
        n = study_row$n,
        grid_size = grid_size_contour,
        return_contour_data = TRUE
      )
      
      if (!is.null(result)) {
        key <- paste0(gsub(" ", "", study_name), "_", group_type)
        results_list[[key]] <- result
        
        idx <- result$type1_indices
        cat(sprintf("%s: n=%d, k_median=%d (Type 1: ceiling)\n", 
                    key, study_row$n, idx$k2))
      }
    }
  }
  
  # Step 2: Calculate common axis ranges for each study pair
  study_ranges <- list()
  for (study_num in 1:n_studies) {
    case_key <- paste0("Study", study_num, "_Case")
    control_key <- paste0("Study", study_num, "_Control")
    
    case_result <- results_list[[case_key]]
    control_result <- results_list[[control_key]]
    
    # Helper function to find effective range of posterior density
    get_effective_range <- function(result) {
      threshold <- max(result$posterior_density) * 0.01
      
      # Find mu range where density > threshold
      mu_ranges <- apply(result$posterior_density, 2, function(col) {
        valid_indices <- which(col > threshold)
        if (length(valid_indices) > 0) {
          c(result$mu_grid[min(valid_indices)], result$mu_grid[max(valid_indices)])
        } else {
          c(NA, NA)
        }
      })
      
      # Find sigma range where density > threshold
      sigma_ranges <- apply(result$posterior_density, 1, function(row) {
        valid_indices <- which(row > threshold)
        if (length(valid_indices) > 0) {
          c(result$sigma_grid[min(valid_indices)], result$sigma_grid[max(valid_indices)])
        } else {
          c(NA, NA)
        }
      })
      
      # Compute overall ranges
      mu_dist_range <- c(min(mu_ranges[1,], na.rm = TRUE), max(mu_ranges[2,], na.rm = TRUE))
      sigma_dist_range <- c(min(sigma_ranges[1,], na.rm = TRUE), max(sigma_ranges[2,], na.rm = TRUE))
      
      # Include credible intervals in range
      mu_with_ci <- range(c(mu_dist_range, result$mu_CI))
      sigma_with_ci <- range(c(sigma_dist_range, result$sigma_CI))
      
      list(mu = mu_with_ci, sigma = sigma_with_ci)
    }
    
    case_range <- get_effective_range(case_result)
    control_range <- get_effective_range(control_result)
    
    # Combine case and control ranges
    study_mu_range <- range(c(case_range$mu, control_range$mu))
    study_sigma_range <- range(c(case_range$sigma, control_range$sigma))
    
    # Add padding
    mu_padding <- diff(study_mu_range) * 0.05 / 2
    sigma_padding <- diff(study_sigma_range) * 0.1 / 2
    
    study_mu_range <- study_mu_range + c(-mu_padding, mu_padding)
    study_sigma_range <- study_sigma_range + c(-sigma_padding, sigma_padding)
    
    study_ranges[[paste0("Study", study_num)]] <- list(
      mu = study_mu_range,
      sigma = study_sigma_range
    )
  }
  
  # Step 3: Generate plots with matched axis ranges
  for (key in names(results_list)) {
    study_parts <- str_split(key, "_")[[1]]
    study_name <- study_parts[1]
    study_num <- str_extract(study_name, "\\d+")
    study_name <- paste("Study", study_num)
    group_type <- study_parts[2]
    
    ranges <- study_ranges[[paste0("Study", study_num)]]
    
    plot <- create_contour_plot(
      results_list[[key]], 
      study_name, 
      group_type,
      xlim = ranges$mu,
      ylim = ranges$sigma
    )
    plots_list[[key]] <- plot
  }
  
  # Step 4: Create final grid arrangement
  final_grid <- NULL
  if (length(plots_list) >= 2 * n_studies) {
    # Order: cases on top row, controls on bottom row
    plot_order <- c(
      paste0("Study", 1:n_studies, "_Case"),
      paste0("Study", 1:n_studies, "_Control")
    )
    
    ordered_plots <- plots_list[plot_order]
    
    final_grid <- grid.arrange(
      grobs = ordered_plots,
      ncol = n_studies, 
      nrow = 2,
      widths = unit(rep(1, n_studies), "null"),
      heights = unit(rep(1, 2), "null"),
      top = NULL,
      respect = TRUE
    )
  }
  
  end_time <- Sys.time()
  cat("\n⏱️ Time taken: ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds\n")
  
  return(list(
    results = results_list, 
    plots = plots_list,
    study_ranges = study_ranges,
    final_grid = final_grid
  ))
}

#===============================================================================
#                              USAGE EXAMPLE
#===============================================================================

# Uncomment to run contour plot analysis
# cat("\n📊 Step 3: Contour Plot Generation\n")
# contour_results <- run_contour_analysis()
# 
# # Display the final grid
# if (!is.null(contour_results$final_grid)) {
#   print(contour_results$final_grid)
# }
# 
# # Optional: Save to file
# # ggsave("type1_contour_plots.png", contour_results$final_grid, 
# #        width = 10, height = 8, dpi = 600, bg = "white")
