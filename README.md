# BOSE: Bayesian Order Statistics Based Estimator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> A Bayesian methodology for meta-analysis that recovers sample mean and standard deviation from order statistics (quantiles) using Type 1 quantile definitions.

## Overview

BOSE (Bayesian Order Statistics Estimator) provides a rigorous Bayesian framework for estimating population parameters from summary statistics commonly reported in medical and scientific literature. This is particularly valuable for meta-analyses where full datasets are unavailable.

### Key Features

- ✅ **Type 1 Quantile Support**: Maintains proper 95% credible interval coverage
- ✅ **Multiple Summary Formats**: Handles 3-number, 5-number, and quartile summaries
- ✅ **Adaptive Grid Sampling**: Efficient posterior approximation without MCMC
- ✅ **Interactive Shiny App**: User-friendly web interface for estimation
- ✅ **Comprehensive Validation**: Extensive simulation studies comparing with existing methods

### Critical Finding

**Type 1 quantiles provide superior uncertainty quantification compared to Type 7 quantiles.** Our research demonstrates that Type 1 quantiles (using `ceiling(p*n)`) maintain proper 95% credible interval coverage, while Type 7 quantiles (R's default) exhibit systematic under-coverage for certain sample sizes.

## Installation

```r
# Install from GitHub (recommended)
# install.packages("devtools")
devtools::install_github("yourusername/BOSE-method")

# Or download and install locally
# install.packages("path/to/BOSE-method", repos = NULL, type = "source")
```

### Required Dependencies

```r
install.packages(c(
  "estmeansd",   # For comparison methods (Luo, Wan)
  "metaBLUE",    # For BLUE estimator
  "tidyverse",   # Data manipulation and visualization
  "shiny",       # Interactive web application
  "future.apply" # Parallel computing for simulations
))
```

## Quick Start

### Basic Example: Three-Number Summary

```r
library(BOSE)

# You have a study reporting: min = 3.2, median = 5.0, max = 6.8, n = 30
X <- c(3.2, 5.0, 6.8)  # Order statistics
n <- 30                # Sample size
S <- 1                 # Scenario: three-number summary

# Estimate mean and SD
result <- estimate_bose(X, n, S)

# View results
result$mu_est        # Estimated mean: ~4.98
result$sigma_est     # Estimated SD: ~1.02
result$mu_CI         # 95% CI for mean: [4.5, 5.5]
result$sigma_CI      # 95% CI for SD: [0.8, 1.3]
```

### Scenario Types

BOSE supports three types of summary statistics:

| Scenario | Summary Type | Input Format | Example |
|----------|--------------|--------------|---------|
| **S = 1** | Three-number | `{min, median, max}` | `c(3.2, 5.0, 6.8)` |
| **S = 2** | Quartiles | `{Q1, median, Q3}` | `c(4.0, 5.0, 6.0)` |
| **S = 3** | Five-number | `{min, Q1, median, Q3, max}` | `c(3.2, 4.0, 5.0, 6.0, 6.8)` |

### Advanced Example: Quartile Summary

```r
# Study reports: Q1 = 45, Median = 52, Q3 = 60, n = 120
result <- estimate_bose(
  X = c(45, 52, 60),
  n = 120,
  S = 2  # Quartile summary
)

# Posterior samples are also returned
hist(result$post_samples$mu, 
     main = "Posterior Distribution of Mean",
     xlab = expression(mu))
```

## Interactive Shiny App

Launch the interactive web application:

```r
library(shiny)
runApp("inst/shiny/app.R")
# Or if package is installed:
# BOSE::launch_app()
```

The Shiny app provides:
- Interactive input for all three scenario types
- Real-time visualization of posterior distributions
- Contour plots showing joint posterior of (μ, σ)
- Downloadable results and plots

## Methodology

### Type 1 vs Type 7 Quantiles

The choice of quantile definition significantly impacts inference:

- **Type 1** (used by BOSE): `k = ceiling(p × n)`
  - No interpolation between order statistics
  - Maintains proper uncertainty quantification
  - Suitable for Bayesian credible intervals
  
- **Type 7** (R default): Linear interpolation
  - Can produce systematic under-coverage
  - May underestimate uncertainty for certain n values

### Bayesian Framework

Given order statistics **X** = {X_(k₁), X_(k₂), ..., X_(kₘ)} from a sample of size n:

**Likelihood:**
```
L(μ, σ | X) ∝ ∏ᵢ φ(Xᵢ; μ, σ) × ∏ⱼ [F(Xⱼ₊₁) - F(Xⱼ)]^(kⱼ₊₁ - kⱼ - 1)
```

**Priors:**
- μ ~ Uniform(L, U)  [weakly informative]
- σ² ~ InverseGamma(0.01, 0.01)  [non-informative]

**Posterior Approximation:**
1. Coarse grid (128 × 128) to identify high-density region
2. Fine grid (256 × 256) within region containing 99% of posterior mass
3. Draw samples proportional to posterior density

## Simulation Studies

Comprehensive simulations demonstrate BOSE's performance:

```r
# Run simulation for Normal(5, 1²) with n ranging from 5 to 200
source("simulations/run_simulation_type1.R")

# Results include:
# - RMSE comparisons with Luo, Wan, BC, QE, MLN, BLUE methods
# - Coverage probabilities for 95% credible intervals
# - Performance across three scenarios (S=1, S=2, S=3)
```

### Simulation Results Summary

For Normal(5, 1²) distribution with three-number summary (S=1):

| Method | RMSE(μ̂) | RMSE(σ̂) | Coverage(μ) | Coverage(σ) |
|--------|---------|---------|-------------|-------------|
| **BOSE** | **0.92** | **0.95** | **94.8%** | **95.2%** |
| Luo/Wan | 0.98 | 1.03 | N/A | N/A |
| BLUE | 0.94 | 0.97 | N/A | N/A |

## Real Data Analysis

Example using cardiovascular studies:

```r
source("case_studies/real_data_analysis.R")

# Analyze multiple studies
results <- run_real_data_analysis()

# Generate contour plots
contour_results <- run_contour_analysis()
```

## Visualization

### RMSE Plots

```r
source("visualization/plot_rmse.R")

# Compare methods across sample sizes
plots <- plot_RMSE_combined_split_by_S(
  df = simulation_results,
  dist = "normal",
  facet = "case"
)

# Display for scenario S=1
plots$S1
```

### Coverage Analysis

```r
source("visualization/plot_coverage.R")

# Boxplots of coverage by scenario
plot_coverage_by_scenario(
  df = simulation_results,
  case_num = 1
)
```

## Citation

If you use BOSE in your research, please cite:

```bibtex
@article{pan2025bose,
  title={BOSE: A Bayesian Order Statistics Estimator for Recovering the Sample Mean and Standard Deviation},
  author={Pan, Wenqisi and Wang, Xinlei},
  journal={In Preparation},
  year={2025}
}
```

## Repository Structure

```
BOSE-method/
├── R/                      # Core functions
│   ├── bose_core.R         # Main BOSE estimation
│   └── helper_functions.R  # Utilities
├── simulations/            # Simulation studies
│   └── run_simulation_type1.R
├── case_studies/           # Real data applications
│   └── real_data_analysis.R
├── visualization/          # Plotting functions
│   └── shiny_app.R         # Interactive app
├── data/                   # Example datasets
├── inst/examples/          # Usage examples
└── vignettes/              # Tutorials
```

## Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Wenqisi (Lydia) Pan**  
PhD Candidate, Statistics  
University of Texas at Arlington  
📧 wensqisi.pan@mavs.uta.edu  
🔗 [GitHub](https://github.com/yourusername)

**Advisor: Dr. Xinlei Wang**  
Jenkins-Garrett Professor  
Department of Mathematics  
University of Texas at Arlington

## Acknowledgments

- Dr. Xinlei Wang for invaluable guidance and mentorship
- UTA American Statistical Association Student Chapter
- Developers of `estmeansd` and `metaBLUE` packages for comparison methods

## References

1. Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. *Statistical Methods in Medical Research*, 27(6), 1785-1805.

2. Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. *BMC Medical Research Methodology*, 14(1), 135.

3. Shi, J., Luo, D., Weng, H., et al. (2020). Optimally estimating the sample standard deviation from the five-number summary. *Research Synthesis Methods*, 11(5), 641-654.

---

**Note:** This is research software. While extensively tested, please validate results for your specific application. For production use in meta-analyses, we recommend consulting with a statistician.
