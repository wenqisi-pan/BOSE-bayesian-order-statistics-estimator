library(shiny)
library(bslib)

source("00_bose_core.R")

weighted_quantile <- function(values, weights, probs) {
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord]
  cum_w <- cumsum(weights) / sum(weights)
  v_out <- sapply(probs, function(p) {
    if (p <= cum_w[1]) return(values[1])
    idx <- which(cum_w >= p)[1]
    return(values[idx])
  })
  return(v_out)
}

ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  
  titlePanel("Recovering the sample mean and standard deviation"),
  
  div(style = "margin-bottom: 20px;",
      p("This webpage applies the Bayesian Order Statistics-Based Estimator (BOSE) method, along with the Luo/Wan/Shi (LWS) approaches, the Best Linear Unbiased Estimator (BLUE), the Quantile Estimation (QE) method, the Box-Cox (BC) method, and the Method for Unknown Non-Normal Distributions (MLN) approach to estimate the sample mean and standard deviation from a study that presents one of the following sets of summary statistics:"),
      tags$ul(
        tags$li("S1: median, minimum and maximum values, and sample size"),
        tags$li("S2: median, first and third quartiles, and sample size"),
        tags$li("S3: median, minimum and maximum values, first and third quartiles, and sample size")
      ),
      p("The BOSE method additionally provides 95% credible intervals for the estimated mean, standard deviation, and coefficient of variation (CV)."),
      p(em("If the summary statistics are not in the appropriate form (i.e., S1, S2, or S3), values of NA will be given for the estimated mean and standard deviation."))
  ),
  
  hr(),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("min_val", "Minimum value (if available)", value = NA),
      numericInput("q1_val", "First quartile value (if available)", value = NA),
      numericInput("med_val", "Median value (required)", value = NA),
      numericInput("q3_val", "Third quartile value (if available)", value = NA),
      numericInput("max_val", "Maximum value (if available)", value = NA),
      numericInput("n_val", "Sample size (required)", value = NA, min = 1),
      hr(),
      radioButtons("method", "Method",
                   choices = c(
                     "BOSE" = "bose",
                     "LWS (Luo / Wan / Shi)" = "lws",
                     "BLUE" = "blue",
                     "Box-Cox (BC)" = "bc",
                     "Quantile Estimation (QE)" = "qe",
                     "Method for Unknown Non-Normal Distributions (MLN)" = "mln"
                   ),
                   selected = "bose")
    ),
    
    mainPanel(
      br(),
      uiOutput("results_ui")
    )
  ),
  
  hr(),
  
  div(style = "margin-top: 30px; margin-bottom: 30px; color: #555;",
      h4("References"),
      tags$ul(
        tags$li(strong("BOSE:"), " A Bayesian Order Statistics-Based Estimator for Recovering the Sample Mean and Standard Deviation (Working Paper)."),
        
        tags$li(strong("Wan:"), " Wan, X., Wang, W., Liu, J., & Tong, T. (2014). ", 
                tags$a(href="https://doi.org/10.1186/1471-2288-14-135", target="_blank", 
                       "Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range."), 
                em("BMC Medical Research Methodology"), ", 14(1), 135."),
        
        tags$li(strong("Luo:"), " Luo, D., Wan, X., Liu, J., & Tong, T. (2018). ", 
                tags$a(href="https://doi.org/10.1177/0962280216669183", target="_blank", 
                       "Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range."), 
                em("Statistical Methods in Medical Research"), ", 27(6), 1785-1805."),
        
        tags$li(strong("Shi:"), " Shi, J., Luo, D., Weng, H., Zeng, X. T., Lin, L., Chu, H., & Tong, T. (2020). ", 
                tags$a(href="https://doi.org/10.1002/jrsm.1429", target="_blank", 
                       "Optimally estimating the sample standard deviation from the five-number summary."), 
                em("Research Synthesis Methods"), ", 11(5), 641-654."),
        
        tags$li(strong("BLUE:"), " Yang, X., Hutson, A. D., & Wang, D. (2022). ", 
                tags$a(href="https://doi.org/10.1080/02664763.2021.1967890", target="_blank", 
                       "A generalized BLUE approach for combining location and scale information in a meta-analysis."), 
                em("Journal of Applied Statistics"), ", 49(15), 3846-3867."),
        
        tags$li(strong("QE & BC:"), " McGrath, S., Zhao, X., Steele, R., Thombs, B. D., Benedetti, A., & Levis, B. (2020). ", 
                tags$a(href="https://doi.org/10.1177/0962280219889080", target="_blank", 
                       "Estimating the sample mean and standard deviation from commonly reported quantiles in meta-analysis."), 
                em("Statistical Methods in Medical Research"), ", 29(9), 2520-2537."),
        
        tags$li(strong("MLN:"), " Cai, S., Zhou, J., & Pan, J. (2021). ", 
                tags$a(href="https://doi.org/10.1177/09622802211047348", target="_blank", 
                       "Estimating the sample mean and standard deviation from order statistics and sample size in meta-analysis."), 
                em("Statistical Methods in Medical Research"), ", 30(12), 2701-2719.")
      )
  )
)

server <- function(input, output, session) {
  
  parsed_data <- reactive({
    req(input$n_val, input$med_val)
    min_v <- input$min_val; q1_v <- input$q1_val; med_v <- input$med_val
    q3_v <- input$q3_val; max_v <- input$max_val
    
    has_min_max <- !is.na(min_v) && !is.na(max_v)
    has_q1_q3   <- !is.na(q1_v) && !is.na(q3_v)
    
    if (has_min_max && has_q1_q3) {
      S <- 3; X <- c(min_v, q1_v, med_v, q3_v, max_v)
    } else if (has_q1_q3 && !has_min_max) {
      S <- 2; X <- c(q1_v, med_v, q3_v)
    } else if (has_min_max && !has_q1_q3) {
      S <- 1; X <- c(min_v, med_v, max_v)
    } else { S <- NA; X <- NULL }
    
    list(S = S, X = X, n = input$n_val)
  })
  
  output$results_ui <- renderUI({
    dat <- parsed_data()
    
    if (is.na(dat$S)) {
      return(HTML("<span style='color:#e74c3c;'><strong>Please enter valid summary statistics (S1, S2, or S3).</strong></span>"))
    }
    
    raw_vals <- c(input$min_val, input$q1_val, input$med_val, input$q3_val, input$max_val)
    valid_vals <- raw_vals[!is.na(raw_vals)]
    if (is.unsorted(valid_vals)) {
      return(HTML("<span style='color:#e74c3c;'><strong>Error: The input values must be in ascending order (Min &le; Q1 &le; Median &le; Q3 &le; Max).</strong></span>"))
    }
    
    if ((dat$S == 1 || dat$S == 2) && dat$n < 3) {
      return(HTML("<span style='color:#e74c3c;'><strong>Error: For a 3-number summary (S1 or S2), the sample size must be at least 3.</strong></span>"))
    }
    if (dat$S == 3 && dat$n < 5) {
      return(HTML("<span style='color:#e74c3c;'><strong>Error: For a 5-number summary (S3), the sample size must be at least 5.</strong></span>"))
    }
    
    comp <- run_competing_methods(dat$X, dat$n, dat$S)
    est_mu <- NA; est_sig <- NA; extra_info <- ""
    
    if (input$method == "bose") {
      z <- 5
      L_b <- comp$mu_luo - z * comp$sigma_wan; U_b <- comp$mu_luo + z * comp$sigma_wan
      L_s <- max(comp$sigma_wan / z, 1e-8); U_s <- max(z * comp$sigma_wan, L_s * 1.0001)
      
      res <- get_posterior_weights_adaptive(dat$X, dat$n, dat$S, L_b, U_b, L_s, U_s)
      stats <- compute_posterior_statistics_cpp(res$mu_grid, res$sigma_grid, res$w_f)
      
      est_mu <- stats$mu_mean
      est_sig <- stats$sigma_median
      
      is_cv_valid <- TRUE
      if (any(dat$X <= 0, na.rm = TRUE) || abs(est_mu) < 1e-4) {
        is_cv_valid <- FALSE
      }
      
      if (is_cv_valid) {
        mu_mat <- matrix(res$mu_grid, nrow = length(res$mu_grid), ncol = length(res$sigma_grid))
        sig_mat <- matrix(res$sigma_grid, nrow = length(res$mu_grid), ncol = length(res$sigma_grid), byrow = TRUE)
        cv_vals <- as.vector(sig_mat / mu_mat)
        w_flat <- as.vector(res$w_f)
        
        cv_mean <- sum(cv_vals * w_flat)
        cv_ci <- weighted_quantile(cv_vals, w_flat, c(0.025, 0.975))
        
        cv_warning_html <- ""
        if (!is.na(cv_ci[2]) && cv_ci[2] > 1) {
          cv_warning_html <- "<br><p style='color: #d35400; background-color: #fcf3cf; padding: 10px; border-left: 4px solid #d35400;'><strong>⚠️ Warning:</strong> The 95% Confidence Interval for the CV is extremely large (upper bound > 1), indicating very high uncertainty in the parameter estimates relative to the mean. Please interpret the results with caution.</p>"
        }
        
        cv_output_html <- sprintf("
          <p><strong>Estimated CV:</strong> %.4f</p>
          <p><strong>95%% CI (CV):</strong> [%.4f, %.4f]</p>
          %s", cv_mean, cv_ci[1], cv_ci[2], cv_warning_html)
        
      } else {
        cv_output_html <- "<br><p style='color: #7f8c8d; background-color: #f2f4f4; padding: 10px; border-left: 4px solid #7f8c8d;'><strong>Note:</strong> The Coefficient of Variation (CV) is not calculated. CV is only interpretable for strictly positive data on a ratio scale. Your input data contains zero/negative values, or the estimated mean is too close to zero.</p>"
      }
      
      extra_info <- sprintf("
        <hr>
        <p><strong>95%% Credible Interval (Mean):</strong> [%.4f, %.4f]</p>
        <p><strong>95%% Credible Interval (SD):</strong> [%.4f, %.4f]</p>
        %s",
                            stats$mu_ci_95[1], stats$mu_ci_95[2],
                            stats$sigma_ci_95[1], stats$sigma_ci_95[2],
                            cv_output_html)
      
    } else if (input$method == "lws") {
      est_mu <- comp$mu_luo
      est_sig <- if (dat$S == 3 && !is.na(comp$sigma_shi)) comp$sigma_shi else comp$sigma_wan
      
    } else if (input$method == "blue") {
      est_mu <- comp$mu_blue
      est_sig <- comp$sigma_blue
      extra_info <- sprintf("
        <hr>
        <p><strong>Estimated standard error (Mean):</strong> %.4f</p>
        <p><strong>Estimated standard error (SD):</strong> %.4f</p>",
                            sqrt(comp$blue_Var_mu), sqrt(comp$blue_Var_sigma))
      
    } else if (input$method == "bc") {
      est_mu <- comp$mu_bc; est_sig <- comp$sigma_bc
    } else if (input$method == "qe") {
      est_mu <- comp$mu_qe; est_sig <- comp$sigma_qe
    } else if (input$method == "mln") {
      est_mu <- comp$mu_mln; est_sig <- comp$sigma_mln
    }
    
    HTML(sprintf("
      <h4><strong>Estimated mean:</strong> %s</h4>
      <h4><strong>Estimated standard deviation:</strong> %s</h4>
      %s
    ",
                 ifelse(is.na(est_mu), "NA", formatC(est_mu, format = "f", digits = 4)),
                 ifelse(is.na(est_sig), "NA", formatC(est_sig, format = "f", digits = 4)),
                 extra_info))
  })
}

shinyApp(ui = ui, server = server)
