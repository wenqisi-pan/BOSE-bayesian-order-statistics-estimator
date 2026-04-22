library(Rcpp)
library(estmeansd)
library(metaBLUE)
library(tidyverse)

# ============================================================================
# Rcpp Core Functions
# ============================================================================

rcpp_code <- '
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

const double TINY = 1e-300;
const double LOG_TINY = std::log(TINY);
const double SQRT_2 = 1.4142135623730951;
const double LOG_SQRT_2PI = 0.9189385332046727;

inline double dnorm_log(double x, double mu, double sigma) {
    double z = (x - mu) / sigma;
    return -0.5 * z * z - std::log(sigma) - LOG_SQRT_2PI;
}
inline double pnorm_fast(double x) {
    return 0.5 * (1.0 + std::erf(x / SQRT_2));
}
inline double safe_log(double x) {
    return (x > TINY) ? std::log(x) : LOG_TINY;
}

// Log-posterior on (mu, sigma) grid.
// Type 1 quantile indices: ceiling-based.
// Priors: mu ~ Uniform(L,U), sigma^2 ~ Inv-Gamma(alpha, beta).
// [[Rcpp::export]]
NumericMatrix compute_logpost_grid_cpp(
    NumericVector mu_grid, NumericVector sigma_grid, NumericVector X,
    int n, int S, double L, double U,
    double alpha_prior, double beta_prior
) {
    int G1 = mu_grid.size(), G2 = sigma_grid.size(), nX = X.size();
    NumericMatrix log_post(G1, G2);
    double log_prior_mu = -std::log(U - L);

    int k0 = 1, k1, k2, k3, k4;
    if (S == 1) { k1 = 1; k2 = (int)std::ceil(0.5*n); k3 = n; }
    else if (S == 2) { k1=(int)std::ceil(0.25*n); k2=(int)std::ceil(0.5*n); k3=(int)std::ceil(0.75*n); }
    else { k0=1; k1=(int)std::ceil(0.25*n); k2=(int)std::ceil(0.5*n); k3=(int)std::ceil(0.75*n); k4=n; }

    for (int j = 0; j < G2; j++) {
        double sig = sigma_grid[j], sig2 = sig*sig, inv_sig = 1.0/sig;
        double log_prior_sigma2 = alpha_prior*std::log(beta_prior) - std::lgamma(alpha_prior)
            - (alpha_prior+1.0)*std::log(sig2) - beta_prior/sig2;
        for (int i = 0; i < G1; i++) {
            double mu = mu_grid[i], log_lik = 0.0;
            for (int k = 0; k < nX; k++) log_lik += dnorm_log(X[k], mu, sig);
            std::vector<double> p(nX);
            for (int k = 0; k < nX; k++) p[k] = pnorm_fast((X[k]-mu)*inv_sig);
            double cdf_terms = 0.0;
            if (S == 1) {
                cdf_terms = (k2-k1-1)*safe_log(p[1]-p[0]) + (k3-k2-1)*safe_log(p[2]-p[1]);
            } else if (S == 2) {
                cdf_terms = (k1-1)*safe_log(p[0]) + (k2-k1-1)*safe_log(p[1]-p[0])
                    + (k3-k2-1)*safe_log(p[2]-p[1]) + (n-k3)*safe_log(1.0-p[2]);
            } else {
                cdf_terms = (k1-k0-1)*safe_log(p[1]-p[0]) + (k2-k1-1)*safe_log(p[2]-p[1])
                    + (k3-k2-1)*safe_log(p[3]-p[2]) + (k4-k3-1)*safe_log(p[4]-p[3]);
            }
            log_post(i,j) = log_lik + cdf_terms + log_prior_mu + log_prior_sigma2;
        }
    }
    return log_post;
}

// Quantile from discrete marginal posterior
// [[Rcpp::export]]
double quantile_from_grid_cpp(NumericVector grid_values,
                              NumericVector probabilities,
                              double quantile_level) {
    int n = grid_values.size();
    double sum_prob = 0.0;
    for (int i = 0; i < n; i++) sum_prob += probabilities[i];
    std::vector<double> cdf(n); double cum_sum = 0.0;
    for (int i = 0; i < n; i++) { cum_sum += probabilities[i]/sum_prob; cdf[i] = cum_sum; }
    if (quantile_level <= cdf[0]) return grid_values[0];
    if (quantile_level >= cdf[n-1]) return grid_values[n-1];
    int idx_before = 0;
    for (int i = 0; i < n; i++) { if (cdf[i] < quantile_level) idx_before = i; else break; }
    int idx_after = idx_before + 1;
    return grid_values[idx_before] + (quantile_level - cdf[idx_before]) /
        (cdf[idx_after] - cdf[idx_before]) * (grid_values[idx_after] - grid_values[idx_before]);
}

// Posterior summary: point estimates + credible intervals (95/90/80%).
// [[Rcpp::export]]
List compute_posterior_statistics_cpp(NumericVector mu_grid,
                                      NumericVector sigma_grid,
                                      NumericMatrix w_f) {
    int G1 = mu_grid.size(), G2 = sigma_grid.size();
    double w_sum = 0.0;
    for (int i = 0; i < G1; i++)
        for (int j = 0; j < G2; j++) w_sum += w_f(i,j);
    NumericVector post_mu(G1), post_sigma(G2);
    for (int i = 0; i < G1; i++) {
        double s=0; for (int j=0;j<G2;j++) s+=w_f(i,j); post_mu[i]=s/w_sum;
    }
    for (int j = 0; j < G2; j++) {
        double s=0; for (int i=0;i<G1;i++) s+=w_f(i,j); post_sigma[j]=s/w_sum;
    }
    double mu_mean=0, sigma_mean=0;
    for (int i=0;i<G1;i++) mu_mean += mu_grid[i]*post_mu[i];
    for (int j=0;j<G2;j++) sigma_mean += sigma_grid[j]*post_sigma[j];
    return List::create(
        Named("mu_mean")=mu_mean,
        Named("mu_median")=quantile_from_grid_cpp(mu_grid, post_mu, 0.5),
        Named("mu_ci_95")=NumericVector::create(
            quantile_from_grid_cpp(mu_grid,post_mu,0.025), quantile_from_grid_cpp(mu_grid,post_mu,0.975)),
        Named("mu_ci_90")=NumericVector::create(
            quantile_from_grid_cpp(mu_grid,post_mu,0.05), quantile_from_grid_cpp(mu_grid,post_mu,0.95)),
        Named("mu_ci_80")=NumericVector::create(
            quantile_from_grid_cpp(mu_grid,post_mu,0.10), quantile_from_grid_cpp(mu_grid,post_mu,0.90)),
        Named("sigma_mean")=sigma_mean,
        Named("sigma_median")=quantile_from_grid_cpp(sigma_grid, post_sigma, 0.5),
        Named("sigma_ci_95")=NumericVector::create(
            quantile_from_grid_cpp(sigma_grid,post_sigma,0.025), quantile_from_grid_cpp(sigma_grid,post_sigma,0.975)),
        Named("sigma_ci_90")=NumericVector::create(
            quantile_from_grid_cpp(sigma_grid,post_sigma,0.05), quantile_from_grid_cpp(sigma_grid,post_sigma,0.95)),
        Named("sigma_ci_80")=NumericVector::create(
            quantile_from_grid_cpp(sigma_grid,post_sigma,0.10), quantile_from_grid_cpp(sigma_grid,post_sigma,0.90))
    );
}

// Mixture normal random generator
// [[Rcpp::export]]
NumericVector rmixnorm_cpp(int n, double mu1, double sigma1,
                           double mu2, double sigma2, double p1) {
    NumericVector result(n);
    for (int i = 0; i < n; i++) {
        if (R::runif(0.0, 1.0) < p1) result[i] = R::rnorm(mu1, sigma1);
        else result[i] = R::rnorm(mu2, sigma2);
    }
    return result;
}
'

message("Compiling Rcpp code...")
sourceCpp(code = rcpp_code)
message("Rcpp compiled successfully.")

# ============================================================================
# Two-Stage Adaptive Grid Posterior
# ============================================================================

get_posterior_weights_adaptive <- function(X, n, S, L, U, L_sig, U_sig,
                                          coarse = 128, fine = 256, mass = 0.99,
                                          alpha_prior = 0.01, beta_prior = 0.01) {
  mu_c  <- seq(L, U, length.out = coarse)
  sig_c <- seq(L_sig, U_sig, length.out = coarse)
  lp_c  <- compute_logpost_grid_cpp(mu_c, sig_c, X, n, S, L, U, alpha_prior, beta_prior)

  max_lp <- max(lp_c, na.rm = TRUE)
  w_c <- exp(lp_c - max_lp); w_sum <- sum(w_c)
  if (!is.finite(w_sum) || w_sum <= 0) {
    w_f <- matrix(1, fine, fine)
    return(list(mu_grid = seq(L, U, length.out = fine),
                sigma_grid = seq(L_sig, U_sig, length.out = fine),
                w_f = w_f / sum(w_f)))
  }
  w_c <- w_c / w_sum

  flat_w <- sort(as.vector(w_c), decreasing = TRUE)
  thr <- flat_w[which(cumsum(flat_w) >= mass)[1]]
  idx_mat <- which(w_c >= thr, arr.ind = TRUE)
  mu_range  <- range(mu_c[unique(idx_mat[, 1])])
  sig_range <- range(sig_c[unique(idx_mat[, 2])])

  mu_pad  <- max(1e-8, 0.1 * diff(mu_range))
  sig_pad <- max(1e-8, 0.1 * diff(sig_range))
  mu_f  <- seq(max(L, mu_range[1] - mu_pad), min(U, mu_range[2] + mu_pad), length.out = fine)
  sig_f <- seq(max(L_sig, sig_range[1] - sig_pad), min(U_sig, sig_range[2] + sig_pad), length.out = fine)

  lp_f <- compute_logpost_grid_cpp(mu_f, sig_f, X, n, S, L, U, alpha_prior, beta_prior)
  w_f <- exp(lp_f - max(lp_f, na.rm = TRUE)); w_sum_f <- sum(w_f)
  if (!is.finite(w_sum_f) || w_sum_f <= 0) return(list(mu_grid = mu_c, sigma_grid = sig_c, w_f = w_c))
  list(mu_grid = mu_f, sigma_grid = sig_f, w_f = w_f / w_sum_f)
}

# ============================================================================
# Competing Methods Wrapper
# ============================================================================
# Returns: list(mu_luo, sigma_wan, sigma_shi, mu_blue, sigma_blue,
#               blue_Var_mu, blue_Var_sigma, mu_bc, sigma_bc,
#               mu_qe, sigma_qe, mu_mln, sigma_mln)

run_competing_methods <- function(X, n, S) {
  if (S == 1) {
    sigma_wan <- Wan.std(X, n, type = "S1")$sigmahat
    mu_luo    <- Luo.mean(X, n, type = "S1")$muhat
    blue      <- BLUE_s(X, n, type = "S1")
    sigma_shi <- NA
    bc <- tryCatch({ b <- bc.mean.sd(min.val=X[1],med.val=X[2],max.val=X[3],n=n,preserve.tail=FALSE,avoid.mc=FALSE); list(mu=b$est.mean,sigma=b$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
    qe <- tryCatch({ q <- qe.mean.sd(min.val=X[1],med.val=X[2],max.val=X[3],n=n); list(mu=q$est.mean,sigma=q$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
    ml <- tryCatch({ m <- mln.mean.sd(min.val=X[1],med.val=X[2],max.val=X[3],n=n); list(mu=m$est.mean,sigma=m$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
  } else if (S == 2) {
    sigma_wan <- Wan.std(X, n, type = "S2")$sigmahat
    mu_luo    <- Luo.mean(X, n, type = "S2")$muhat
    blue      <- BLUE_s(X, n, type = "S2")
    sigma_shi <- NA
    bc <- tryCatch({ b <- bc.mean.sd(q1.val=X[1],med.val=X[2],q3.val=X[3],n=n,preserve.tail=FALSE,avoid.mc=FALSE); list(mu=b$est.mean,sigma=b$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
    qe <- tryCatch({ q <- qe.mean.sd(q1.val=X[1],med.val=X[2],q3.val=X[3],n=n); list(mu=q$est.mean,sigma=q$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
    ml <- tryCatch({ m <- mln.mean.sd(q1.val=X[1],med.val=X[2],q3.val=X[3],n=n); list(mu=m$est.mean,sigma=m$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
  } else {
    phi1 <- qnorm((n - 0.375) / (n + 0.25))
    phi2 <- qnorm((0.75 * n - 0.125) / (n + 0.25))
    sigma_wan <- Wan.std(X, n, type = "S3")$sigmahat
    mu_luo    <- Luo.mean(X, n, type = "S3")$muhat
    sigma_shi <- (X[5]-X[1])/((2+0.14*n^0.6)*phi1) + (X[4]-X[2])/((2+2/(0.07*n^0.6))*phi2)
    blue      <- BLUE_s(X, n, type = "S3")
    bc <- tryCatch({ b <- bc.mean.sd(min.val=X[1],q1.val=X[2],med.val=X[3],q3.val=X[4],max.val=X[5],n=n,preserve.tail=FALSE,avoid.mc=FALSE); list(mu=b$est.mean,sigma=b$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
    qe <- tryCatch({ q <- qe.mean.sd(min.val=X[1],q1.val=X[2],med.val=X[3],q3.val=X[4],max.val=X[5],n=n); list(mu=q$est.mean,sigma=q$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
    ml <- tryCatch({ m <- mln.mean.sd(min.val=X[1],q1.val=X[2],med.val=X[3],q3.val=X[4],max.val=X[5],n=n); list(mu=m$est.mean,sigma=m$est.sd) }, error=function(e) list(mu=NA,sigma=NA))
  }
  list(mu_luo = mu_luo, sigma_wan = sigma_wan, sigma_shi = sigma_shi,
       mu_blue = as.numeric(blue$muhat), sigma_blue = as.numeric(blue$sigmahat),
       blue_Var_mu = as.numeric(blue$Var_mu), blue_Var_sigma = as.numeric(blue$Var_sigma),
       mu_bc = bc$mu, sigma_bc = bc$sigma,
       mu_qe = qe$mu, sigma_qe = qe$sigma,
       mu_mln = ml$mu, sigma_mln = ml$sigma)
}

# ============================================================================
# True Parameter Computation
# ============================================================================

get_true_params <- function(dist, params) {
  if (dist == "normal") {
    list(mean = params$mu_true, sd = params$sigma_true)
  } else if (dist == "lognormal") {
    list(mean = exp(params$meanlog + params$sdlog^2 / 2),
         sd = sqrt((exp(params$sdlog^2) - 1) * exp(2*params$meanlog + params$sdlog^2)))
  } else if (dist == "gamma_dist") {
    list(mean = params$shape / params$rate, sd = sqrt(params$shape) / params$rate)
  } else if (dist == "t_dist") {
    list(mean = 0, sd = if (params$df > 2) sqrt(params$df / (params$df - 2)) else Inf)
  } else if (dist %in% c("mixture_normal", "mixture_mild", "mixture_asym")) {
    m <- params$p1 * params$mu1 + (1 - params$p1) * params$mu2
    v <- params$p1 * (params$sigma1^2 + params$mu1^2) +
         (1 - params$p1) * (params$sigma2^2 + params$mu2^2) - m^2
    list(mean = m, sd = sqrt(v))
  } else stop(paste("Unknown distribution:", dist))
}

# ============================================================================
# Data Generation
# ============================================================================

generate_data <- function(n, dist, params) {
  switch(dist,
    "normal"     = rnorm(n, params$mu_true, params$sigma_true),
    "lognormal"  = rlnorm(n, params$meanlog, params$sdlog),
    "gamma_dist" = rgamma(n, params$shape, params$rate),
    "t_dist"     = rt(n, params$df),
    "mixture_normal" =, "mixture_mild" =, "mixture_asym" =
      rmixnorm_cpp(n, params$mu1, params$sigma1, params$mu2, params$sigma2, params$p1),
    stop(paste("Unknown distribution:", dist))
  )
}

# ============================================================================
# Extract Summary Statistics (Type 1 quantile by default)
# ============================================================================

extract_summary_stats <- function(data, S, qtype = 1L) {
  switch(S,
    `1` = as.numeric(quantile(data, c(0, 0.5, 1), type = qtype)),
    `2` = as.numeric(quantile(data, c(0.25, 0.5, 0.75), type = qtype)),
    `3` = as.numeric(quantile(data, seq(0, 1, length.out = 5), type = qtype)),
    stop("Invalid S value")
  )
}

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

# ============================================================================
# Distribution Parameter Helpers
# ============================================================================

solve_sigma_lnorm <- function(target_skew) {
  uniroot(function(s) { u <- exp(s^2); (u+2)*sqrt(u-1) - target_skew },
          interval = c(1e-6, 10))$root
}

# ============================================================================
# Logging Helper
# ============================================================================

make_logger <- function(progress_file) {
  function(msg) {
    m <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
    message(m)
    cat(m, "\n", file = progress_file, append = TRUE)
  }
}

# ============================================================================
# Main API Function
# ============================================================================

estimate_bose <- function(X, n, S, z = 5) {
  # 1. Get initial estimates to set grid bounds
  comp <- run_competing_methods(X, n, S)
  
  # 2. Determine adaptive grid range
  L_b <- comp$mu_luo - z * comp$sigma_wan
  U_b <- comp$mu_luo + z * comp$sigma_wan
  L_s <- max(comp$sigma_wan / z, 1e-8)
  U_s <- max(z * comp$sigma_wan, L_s * 1.0001)
  
  # 3. Compute posterior weights
  res <- get_posterior_weights_adaptive(X, n, S, L_b, U_b, L_s, U_s)
  
  # 4. Extract Mean and SD statistics
  stats <- compute_posterior_statistics_cpp(res$mu_grid, res$sigma_grid, res$w_f)
  
  # 5. Compute CV and its 95% CI
  mu_mat <- matrix(res$mu_grid, nrow = length(res$mu_grid), ncol = length(res$sigma_grid))
  sig_mat <- matrix(res$sigma_grid, nrow = length(res$mu_grid), ncol = length(res$sigma_grid), byrow = TRUE)
  cv_vals <- as.vector(sig_mat / mu_mat)
  w_flat <- as.vector(res$w_f)
  
  cv_mean <- sum(cv_vals * w_flat)
  cv_ci <- weighted_quantile(cv_vals, w_flat, c(0.025, 0.975))
  
  # 6. Return structured list (invisible for clean console, but usable for objects)
  results <- list(
    mean = list(est = stats$mu_mean, ci = stats$mu_ci_95),
    sd   = list(est = stats$sigma_median, ci = stats$sigma_ci_95),
    cv   = list(est = cv_mean, ci = cv_ci)
  )
  
  # Optional: Print clean summary to console
  message("\n--- BOSE Estimation Summary ---")
  cat(sprintf("Mean: %.4f (95%% CI: [%.4f, %.4f])\n", results$mean$est, results$mean$ci[1], results$mean$ci[2]))
  cat(sprintf("SD:   %.4f (95%% CI: [%.4f, %.4f])\n", results$sd$est,   results$sd$ci[1],   results$sd$ci[2]))
  cat(sprintf("CV:   %.4f (95%% CI: [%.4f, %.4f])\n", results$cv$est,   results$cv$ci[1],   results$cv$ci[2]))
  
  return(invisible(results))
}
