#' BOSE: Bayesian Order Statistics Estimator - Core Functions
#' 
#' This file contains the main functions for Bayesian estimation of mean and
#' standard deviation from order statistics using Type 1 quantiles.
#' 
#' @author Wenqisi (Lydia) Pan
#' @references Pan, W. & Wang, X. (2025). BOSE: A Bayesian Order Statistics
#'   Estimator for Recovering the Sample Mean and Standard Deviation

# Required packages
# library(estmeansd)
# library(metaBLUE)

#===============================================================================
#                         TYPE 1 QUANTILE INDICES
#===============================================================================

#' Calculate Type 1 Quantile Indices
#'
#' Computes the order statistic indices for Type 1 quantiles (ceiling method).
#' This is a key methodological contribution: Type 1 quantiles provide proper
#' credible interval coverage while Type 7 exhibits systematic under-coverage.
#'
#' @param n Sample size (integer)
#' @param S Scenario type (integer: 1, 2, or 3)
#'   - S=1: three-number summary {min, median, max}
#'   - S=2: quartile summary {Q1, median, Q3}
#'   - S=3: five-number summary {min, Q1, median, Q3, max}
#'
#' @return List of indices (k1, k2, k3, or k0-k4 depending on S)
#'
#' @details Type 1 quantiles use ceiling(p*n) for the pth quantile index.
#'   This avoids interpolation issues that affect uncertainty quantification.
#'
#' @examples
#' get_type1_indices(30, S = 1)  # Returns: k1=1, k2=15, k3=30
#' get_type1_indices(30, S = 2)  # Returns: k1=8, k2=15, k3=23
#'
#' @export
get_type1_indices <- function(n, S) {
  # Type 1: ceiling method (round up)
  if (S == 1) {
    # {min, median, max}
    list(
      k1 = 1,
      k2 = ceiling(0.5 * n),
      k3 = n
    )
  } else if (S == 2) {
    # {Q1, median, Q3}
    list(
      k1 = ceiling(0.25 * n),
      k2 = ceiling(0.5 * n),
      k3 = ceiling(0.75 * n)
    )
  } else if (S == 3) {
    # {min, Q1, median, Q3, max}
    list(
      k0 = 1,
      k1 = ceiling(0.25 * n),
      k2 = ceiling(0.5 * n),
      k3 = ceiling(0.75 * n),
      k4 = n
    )
  } else {
    stop("Invalid S value. Must be 1, 2, or 3.")
  }
}

#===============================================================================
#                         LOG POSTERIOR COMPUTATION
#===============================================================================

#' Compute Log Posterior on Grid (Vectorized, Type 1)
#'
#' Efficiently computes the log posterior density over a 2D grid of (mu, sigma)
#' using vectorized operations. Uses Type 1 quantile indices.
#'
#' @param mu_grid Vector of mu values
#' @param sigma_grid Vector of sigma values  
#' @param X Vector of observed order statistics
#' @param n Sample size (integer)
#' @param S Scenario type (1, 2, or 3)
#' @param L Lower bound for mu prior (uniform on [L, U])
#' @param U Upper bound for mu prior
#' @param alpha Shape parameter for inverse-gamma prior on sigma^2 (default: 0.01)
#' @param beta Scale parameter for inverse-gamma prior on sigma^2 (default: 0.01)
#'
#' @return Matrix (G1 x G2) of log posterior values
#'
#' @details The posterior is proportional to:
#'   likelihood × prior(mu) × prior(sigma^2)
#'   where prior(mu) ~ Uniform(L, U) and prior(sigma^2) ~ InvGamma(alpha, beta)
#'
#' @examples
#' mu_grid <- seq(4, 6, length.out = 100)
#' sigma_grid <- seq(0.5, 2, length.out = 100)
#' X <- c(3.2, 5.0, 6.8)  # min, median, max
#' logpost <- compute_logpost_grid_fast(mu_grid, sigma_grid, X, n=30, S=1, L=0, U=10)
#'
#' @export
compute_logpost_grid_fast <- function(mu_grid, sigma_grid, X, n, S, L, U,
                                      alpha = 0.01, beta = 0.01) {
  
  G1 <- length(mu_grid)
  G2 <- length(sigma_grid)
  
  # Build meshgrid: MU is G1 x G2, SIG is G1 x G2
  MU <- matrix(rep(mu_grid, times = G2), nrow = G1, ncol = G2)
  SIG <- matrix(rep(sigma_grid, each = G1), nrow = G1, ncol = G2)
  SIG2 <- SIG^2
  
  # Get Type 1 indices
  idx <- get_type1_indices(n, S)
  
  # ========== PDF Terms: Sum of log likelihoods for observed order statistics ==========
  log_lik <- matrix(0, nrow = G1, ncol = G2)
  for (x in X) {
    log_lik <- log_lik + dnorm(x, mean = MU, sd = SIG, log = TRUE)
  }
  
  # ========== CDF Terms: Based on order statistic likelihood ==========
  tiny <- 1e-300  # Prevent log(0)
  
  if (S == 1) {
    # X = {a, m, b} corresponding to {X_(k1), X_(k2), X_(k3)}
    # Exponents: (k2 - k1 - 1) and (k3 - k2 - 1)
    
    p1 <- pnorm((X[1] - MU) / SIG)  # F(a)
    p2 <- pnorm((X[2] - MU) / SIG)  # F(m)
    p3 <- pnorm((X[3] - MU) / SIG)  # F(b)
    
    exp1 <- idx$k2 - idx$k1 - 1
    exp2 <- idx$k3 - idx$k2 - 1
    
    cdf_terms <- exp1 * log(pmax(p2 - p1, tiny)) + 
                 exp2 * log(pmax(p3 - p2, tiny))
    
  } else if (S == 2) {
    # X = {q1, m, q3} corresponding to {X_(k1), X_(k2), X_(k3)}
    # Exponents: (k1-1), (k2 - k1 - 1), (k3 - k2 - 1), (n - k3)
    
    p1 <- pnorm((X[1] - MU) / SIG)  # F(q1)
    p2 <- pnorm((X[2] - MU) / SIG)  # F(m)
    p3 <- pnorm((X[3] - MU) / SIG)  # F(q3)
    
    exp0 <- idx$k1 - 1
    exp1 <- idx$k2 - idx$k1 - 1
    exp2 <- idx$k3 - idx$k2 - 1
    exp3 <- n - idx$k3
    
    cdf_terms <- exp0 * log(pmax(p1, tiny)) +
                 exp1 * log(pmax(p2 - p1, tiny)) + 
                 exp2 * log(pmax(p3 - p2, tiny)) +
                 exp3 * log(pmax(1 - p3, tiny))
    
  } else if (S == 3) {
    # X = {a, q1, m, q3, b} corresponding to {X_(k0), X_(k1), X_(k2), X_(k3), X_(k4)}
    # Exponents: (k1-k0-1), (k2-k1-1), (k3-k2-1), (k4-k3-1)
    
    p0 <- pnorm((X[1] - MU) / SIG)  # F(a)
    p1 <- pnorm((X[2] - MU) / SIG)  # F(q1)
    p2 <- pnorm((X[3] - MU) / SIG)  # F(m)
    p3 <- pnorm((X[4] - MU) / SIG)  # F(q3)
    p4 <- pnorm((X[5] - MU) / SIG)  # F(b)
    
    exp0 <- idx$k1 - idx$k0 - 1
    exp1 <- idx$k2 - idx$k1 - 1
    exp2 <- idx$k3 - idx$k2 - 1
    exp3 <- idx$k4 - idx$k3 - 1
    
    cdf_terms <- exp0 * log(pmax(p1 - p0, tiny)) +
                 exp1 * log(pmax(p2 - p1, tiny)) + 
                 exp2 * log(pmax(p3 - p2, tiny)) +
                 exp3 * log(pmax(p4 - p3, tiny))
  }
  
  # ========== Priors ==========
  # Uniform prior on mu: p(mu) = 1/(U-L)
  if (L >= U) stop("Lower bound should be smaller than upper bound!")
  log_prior_mu <- -log(U - L)
  
  # Inverse-gamma prior on sigma^2: p(sigma^2) ~ IG(alpha, beta)
  if (alpha <= 0 || beta <= 0) stop("Alpha and beta should both be positive!")
  log_prior_sigma2 <- alpha * log(beta) - lgamma(alpha) - 
                      (alpha + 1) * log(SIG2) - beta / SIG2
  
  # Return: log posterior = log likelihood + log prior
  log_lik + cdf_terms + log_prior_mu + log_prior_sigma2
}

#===============================================================================
#                    ADAPTIVE POSTERIOR SAMPLING
#===============================================================================

#' Adaptive Grid-Based Posterior Sampler
#'
#' Two-stage adaptive sampling: (1) coarse grid to find region of interest,
#' (2) fine grid within ROI for accurate posterior approximation.
#'
#' @param X Vector of observed order statistics
#' @param n Sample size
#' @param S Scenario type (1, 2, or 3)
#' @param L Lower bound for mu
#' @param U Upper bound for mu
#' @param L_sig Lower bound for sigma
#' @param U_sig Upper bound for sigma
#' @param coarse Coarse grid size (default: 128)
#' @param fine Fine grid size (default: 256)
#' @param mass Probability mass to capture in ROI (default: 0.99)
#' @param n_samples Number of posterior samples to draw (default: 1000)
#' @param alpha Inverse-gamma shape parameter (default: 0.01)
#' @param beta Inverse-gamma scale parameter (default: 0.01)
#'
#' @return List with elements:
#'   - mu: Vector of posterior samples for mu
#'   - sig: Vector of posterior samples for sigma
#'
#' @details The adaptive approach:
#'   1. Coarse grid (e.g., 128x128) to identify high-density region
#'   2. Refine grid (e.g., 256x256) within region containing 99% of mass
#'   3. Draw samples proportional to posterior density
#'
#' @examples
#' X <- c(3.2, 5.0, 6.8)
#' samp <- get_posterior_samples_adaptive(X, n=30, S=1, L=0, U=10, 
#'                                        L_sig=0.1, U_sig=5)
#' mean(samp$mu)      # Posterior mean of mu
#' median(samp$sig)   # Posterior median of sigma
#'
#' @export
get_posterior_samples_adaptive <- function(X, n, S, L, U, L_sig, U_sig,
                                           coarse = 128, fine = 256,
                                           mass = 0.99, n_samples = 1000,
                                           alpha = 0.01, beta = 0.01) {
  # ========== Stage 1: Coarse Grid ==========
  mu_c <- seq(L, U, length.out = coarse)
  sig_c <- seq(L_sig, U_sig, length.out = coarse)
  lp_c <- compute_logpost_grid_fast(mu_c, sig_c, X, n, S, L, U, alpha, beta)
  w_c <- exp(lp_c - max(lp_c, na.rm = TRUE))
  w_sum <- sum(w_c)
  
  # Fallback: if posterior is degenerate, sample uniformly
  if (!is.finite(w_sum) || w_sum <= 0) {
    mu_smpl <- runif(n_samples, min = L, max = U)
    sig_smpl <- runif(n_samples, min = L_sig, max = U_sig)
    return(list(mu = mu_smpl, sig = sig_smpl))
  }
  w_c <- w_c / w_sum
  
  # Determine region of interest (ROI) containing 'mass' of probability
  flat_order <- order(as.vector(w_c), decreasing = TRUE)
  flat_w <- as.vector(w_c)[flat_order]
  csum <- cumsum(flat_w)
  cut_idx <- which(csum >= mass)[1]
  thr <- flat_w[cut_idx]
  idx_mat <- which(w_c >= thr, arr.ind = TRUE)
  mu_idx <- unique(idx_mat[, 1])
  sig_idx <- unique(idx_mat[, 2])
  mu_min <- min(mu_c[mu_idx]); mu_max <- max(mu_c[mu_idx])
  sig_min <- min(sig_c[sig_idx]); sig_max <- max(sig_c[sig_idx])
  
  # Add padding around ROI
  mu_pad <- max(1e-8, 0.1 * (mu_max - mu_min))
  sig_pad <- max(1e-8, 0.1 * (sig_max - sig_min))
  mu_lo <- max(L, mu_min - mu_pad); mu_hi <- min(U, mu_max + mu_pad)
  sig_lo <- max(L_sig, sig_min - sig_pad); sig_hi <- min(U_sig, sig_max + sig_pad)
  
  # ========== Stage 2: Fine Grid within ROI ==========
  mu_f <- seq(mu_lo, mu_hi, length.out = fine)
  sig_f <- seq(sig_lo, sig_hi, length.out = fine)
  lp_f <- compute_logpost_grid_fast(mu_f, sig_f, X, n, S, L, U, alpha, beta)
  w_f <- exp(lp_f - max(lp_f, na.rm = TRUE))
  w_sum_f <- sum(w_f)
  
  # Fallback to coarse grid if fine grid fails
  if (!is.finite(w_sum_f) || w_sum_f <= 0) {
    post_mu_c <- rowSums(w_c)
    mu_idx_s <- sample.int(length(mu_c), size = n_samples, replace = TRUE, prob = post_mu_c)
    sig_idx_s <- vapply(mu_idx_s, function(i) {
      row <- w_c[i, ]
      if (sum(row) <= 0 || !all(is.finite(row))) sample.int(length(sig_c), 1) 
      else sample.int(length(sig_c), 1, prob = row)
    }, integer(1))
    return(list(mu = mu_c[mu_idx_s], sig = sig_c[sig_idx_s]))
  }
  w_f <- w_f / w_sum_f
  
  # ========== Draw Samples from Fine Grid ==========
  # Marginal distribution of mu
  post_mu <- rowSums(w_f)
  if (any(!is.finite(post_mu)) || sum(post_mu) <= 0) {
    post_mu[] <- 1 / length(post_mu)
  }
  
  # Sample mu indices
  mu_idx_s <- sample.int(length(mu_f), size = n_samples, replace = TRUE, prob = post_mu)
  
  # For each mu, sample sigma from conditional distribution
  sig_idx_s <- vapply(mu_idx_s, function(i) {
    row <- w_f[i, ]
    s <- sum(row)
    if (!is.finite(s) || s <= 0) {
      sample.int(length(sig_f), 1)
    } else {
      sample.int(length(sig_f), 1, prob = row)
    }
  }, integer(1))
  
  list(mu = mu_f[mu_idx_s], sig = sig_f[sig_idx_s])
}

#===============================================================================
#                         MAIN BOSE ESTIMATION
#===============================================================================

#' BOSE: Main Estimation Function
#'
#' Estimates mean and standard deviation from order statistics using Bayesian
#' methods with Type 1 quantiles.
#'
#' @param X Vector of observed order statistics (length 3 or 5)
#' @param n Sample size
#' @param S Scenario type (1, 2, or 3)
#' @param coarse_size Coarse grid size (default: 128)
#' @param fine_size Fine grid size (default: 256)
#' @param n_post_samples Number of posterior samples (default: 1000)
#' @param z Factor for initial bounds (default: 5)
#'
#' @return List with:
#'   - mu_est: Posterior mean of mu
#'   - sigma_est: Posterior median of sigma
#'   - mu_CI: 95% credible interval for mu
#'   - sigma_CI: 95% credible interval for sigma
#'   - post_samples: List of posterior samples (mu, sig)
#'
#' @details Uses Luo's method for initial mu estimate and Wan's method for
#'   initial sigma estimate to set adaptive grid bounds.
#'
#' @examples
#' # Three-number summary: {min, median, max}
#' X <- c(3.2, 5.0, 6.8)
#' result <- estimate_bose(X, n = 30, S = 1)
#' result$mu_est
#' result$sigma_est
#'
#' @export
estimate_bose <- function(X, n, S, 
                          coarse_size = 128, 
                          fine_size = 256,
                          n_post_samples = 1000,
                          z = 5) {
  
  # Get initial estimates using frequentist methods
  # Note: Requires estmeansd package
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
  
  # Set adaptive grid bounds based on initial estimates
  L <- mu_luo - z * sigma_wan
  U <- mu_luo + z * sigma_wan
  L_sig <- max(sigma_wan / z, 1e-8)
  U_sig <- max(z * sigma_wan, L_sig * 1.0001)
  
  # Draw posterior samples
  samp <- get_posterior_samples_adaptive(
    X = X, n = n, S = S,  
    L = L, U = U, L_sig = L_sig, U_sig = U_sig,
    coarse = coarse_size, fine = fine_size, 
    mass = 0.99, n_samples = n_post_samples
  )
  
  post.mu <- samp$mu
  post.sig <- samp$sig
  
  # Point estimates
  mu_est <- mean(post.mu)
  sigma_est <- median(post.sig)
  
  # 95% credible intervals
  ci_mu <- quantile(post.mu, c(0.025, 0.975))
  ci_sigma <- quantile(post.sig, c(0.025, 0.975))
  
  return(list(
    mu_est = mu_est,
    sigma_est = sigma_est,
    mu_CI = as.numeric(ci_mu),
    sigma_CI = as.numeric(ci_sigma),
    post_samples = list(mu = post.mu, sig = post.sig),
    initial_estimates = list(mu_luo = mu_luo, sigma_wan = sigma_wan)
  ))
}
