# BOSE: Bayesian Order Statistics Based Estimator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Bayesian methodology for meta-analysis that recovers the sample mean and standard deviation from order statistics (quantiles).

## Overview

BOSE (Bayesian Order Statistics Estimator) provides a rigorous Bayesian framework for estimating population parameters from summary statistics commonly reported in medical and scientific literature. This is particularly valuable for meta-analyses where full datasets are unavailable.

### Key Features

- ✅ **Type 1 Quantile Support**: Maintains proper 95% credible interval coverage
- ✅ **Multiple Summary Formats**: Handles 3-number, 5-number, and quartile summaries
- ✅ **Adaptive Grid Sampling**: Efficient posterior approximation without MCMC
- ✅ **Comprehensive Validation**: Extensive simulation studies comparing with existing methods
- ✅ **Production-Ready Code**: Modular design with complete documentation
- ✅ **Example Data Included**: Ready-to-use cardiovascular study examples

---

## 🚀 Quick Start

### Installation

```r
# Clone the repository
git clone https://github.com/Lydia0205/BOSE-bayesian-order-statistics-estimator.git
cd BOSE-bayesian-order-statistics-estimator

# Install required packages
install.packages(c("tidyverse", "estmeansd", "metaBLUE", "ggplot2", "patchwork"))
```

### 5-Minute Example

```r
# Load core functions
source("R/bose_core.R")
source("R/helper_functions.R")

# Example: Estimate mean and SD from three-number summary
# Study reports: min = 45, median = 52, max = 68, n = 50
X <- c(45, 52, 68)  # Order statistics: {min, median, max}
n <- 50             # Sample size
S <- 1              # Scenario 1: three-number summary

# Run BOSE estimation
result <- estimate_bose(X, n, S)

# View results
print(paste("Mean:", round(result$mu_est, 2)))
print(paste("95% CI:", round(result$mu_CI[1], 2), "to", round(result$mu_CI[2], 2)))
print(paste("SD:", round(result$sigma_est, 2)))
print(paste("95% CI:", round(result$sigma_CI[1], 2), "to", round(result$sigma_CI[2], 2)))
```

**Output:**
```
[1] "Mean: 54.32"
[1] "95% CI: 50.18 to 58.46"
[1] "SD: 12.85"
[1] "95% CI: 10.23 to 16.12"
```

### Complete Workflow

For a full analysis example including skewness assessment, effect size calculation, and visualization:

```r
# Run the complete workflow example (7 steps)
source("inst/examples/full_workflow.R")

# Or try the quick start (3 examples)
source("inst/examples/quick_start.R")
```

### Using Your Own Data

```r
# Load example data to see the format
example_data <- read.csv("data/example_studies.csv")
head(example_data)

# Prepare your data in the same format
my_studies <- tibble(
  Study = c("My Study 1", "My Study 2"),
  S = c(1, 3),  # 1 = three-number, 2 = quartiles, 3 = five-number
  X_list = list(
    c(10, 15, 25),           # Study 1: min, median, max
    c(8, 12, 16, 22, 35)     # Study 2: min, Q1, median, Q3, max
  ),
  n = c(30, 45)
)

# Apply BOSE to each study
library(tidyverse)
results <- my_studies %>%
  rowwise() %>%
  mutate(
    estimates = list(estimate_bose(X = X_list, n = n, S = S))
  )
```

---

## 📊 Supported Summary Statistics

BOSE handles three types of order statistics commonly reported in literature:

| Scenario | Summary Type | Input Format | Example | Use Case |
|----------|--------------|--------------|---------|----------|
| **S = 1** | Three-number | `{min, median, max}` | `c(3.2, 5.0, 6.8)` | Oldest meta-analyses, clinical trials |
| **S = 2** | Quartiles | `{Q1, median, Q3}` | `c(4.0, 5.0, 6.0)` | Modern reporting, box plots |
| **S = 3** | Five-number | `{min, Q1, median, Q3, max}` | `c(3.2, 4.0, 5.0, 6.0, 6.8)` | Complete summary, Tukey's five-number |

---

## 🔬 Methodology

### Bayesian Framework

Given order statistics **X** = {X_(k₁), X_(k₂), ..., X_(kₘ)} from a sample of size n:

**Likelihood:**
```
L(μ, σ | X) ∝ ∏ᵢ φ(Xᵢ; μ, σ) × ∏ⱼ [F(Xⱼ₊₁) - F(Xⱼ)]^(kⱼ₊₁ - kⱼ - 1)
```

**Priors:**
- μ ~ Uniform(L, U) [weakly informative, adaptive bounds]
- σ² ~ InverseGamma(0.01, 0.01) [non-informative]

**Posterior Approximation (Two-Stage Adaptive Grid):**
1. **Coarse grid** to identify high-posterior-density region
2. **Fine grid** within region containing 99% of posterior mass
3. **Sample** proportional to posterior density (default: 1000 draws)

### Type 1 vs Type 7 Quantiles

The choice of quantile definition critically impacts credible interval coverage:

| Quantile Type | Formula | Properties | Coverage |
|---------------|---------|------------|----------|
| **Type 1** (BOSE) | `k = ceiling(p × n)` | No interpolation, conservative | ✅ Maintains 95% |
| **Type 7** (R default) | Linear interpolation | Smoother, can underestimate uncertainty | ❌ Under-coverage for certain n |

**Key Insight:** Type 7's interpolation between order statistics artificially reduces uncertainty, leading to systematic under-coverage for sample sizes where n mod 4 ≠ 1.

---

## 📈 Simulation Studies

### Running Simulations

```r
# Run Type 1 quantile simulations for Normal(5, 1²)
source("simulations/run_simulation_type1.R")

# Results include:
# - RMSE comparisons with Luo, Wan, BC, QE, MLN, BLUE, Shi methods
# - Coverage probabilities for 95% credible intervals
# - Performance across three scenarios (S=1, S=2, S=3)
# - Multiple distributions: Normal, Student's t, Log-Normal
```

### Performance Summary

For Normal(5, 1²) distribution with three-number summary (S=1):

| Method | RMSE(μ̂) | RMSE(σ̂) | Coverage(μ) | Coverage(σ) | CI Available |
|--------|---------|---------|-------------|-------------|--------------|
| **BOSE** | **0.92** | **0.95** | **94.8%** | **95.2%** | ✅ Yes |
| Luo/Wan | 0.98 | 1.03 | N/A | N/A | ❌ No |
| BLUE | 0.94 | 0.97 | N/A | N/A | ❌ No |
| BC | 1.05 | 1.12 | N/A | N/A | ❌ No |

**Key Advantage:** BOSE is the only method providing proper uncertainty quantification with validated credible interval coverage.

### Visualizing Results

```r
# Load simulation results
results <- read.csv("simulations/results/simulation_results.csv")

# Create RMSE comparison plots
source("visualization/plot_rmse.R")
plots <- plot_rmse_by_scenario(results, dist = "normal")
print(plots$S1)  # Display scenario 1 plot

# Create coverage assessment plots
source("visualization/plot_coverage.R")
p <- plot_coverage_by_scenario(results, dist = "normal", case_num = 1)
print(p)

# Save plots
ggsave("S1_rmse.png", plots$S1, width = 24, height = 12, dpi = 600)
ggsave("coverage.png", p, width = 14, height = 12, dpi = 600)
```

---

## 📚 Real Data Applications

### Case Studies Included

The repository includes complete analyses of cardiovascular studies:

```r
# Load case study functions
source("case_studies/skewness_assessment.R")
source("case_studies/real_data_analysis.R")
source("case_studies/contour_plots.R")

# Step 1: Assess skewness (check distributional assumptions)
skewness_results <- run_skewness_analysis(study_data)

# Step 2: Apply BOSE estimation
analysis_results <- run_real_data_analysis(study_data)

# Step 3: Visualize posterior distributions
contour_results <- run_contour_analysis(study_data, n_studies = 3)
print(contour_results$final_grid)
```

### Skewness Assessment

Before applying BOSE, assess whether normal distribution is appropriate:

**Skewness Measures:**
- **T1**: Based on min, median, max → `(min + max - 2×median) / (max - min)`
- **T2**: Based on Q1, median, Q3 → `(Q1 + Q3 - 2×median) / (Q3 - Q1)`
- **T3**: Sample-size adjusted combination of T1 and T2

**Interpretation:**
- |T1| < 0.2: Approximately symmetric (normal assumption OK)
- |T1| < 0.5: Slightly skewed (BOSE still applicable)
- |T1| ≥ 0.8: Highly skewed (consider transformations)

### Meta-Analysis Integration

BOSE output is ready for standard meta-analysis packages:

```r
# Prepare for metafor
library(metafor)

# Calculate effect sizes
meta_data <- results %>%
  mutate(
    yi = MD,  # Mean difference (from BOSE estimates)
    sei = sqrt((sigma_est_case^2/n_case) + (sigma_est_control^2/n_control))
  )

# Random-effects meta-analysis
model <- rma(yi, sei = sei, data = meta_data, method = "REML")
forest(model)
```

---

## 📖 Documentation

### Core Functions

**`estimate_bose(X, n, S)`** - Main estimation function
- `X`: Vector of observed order statistics
- `n`: Sample size (integer)
- `S`: Scenario type (1, 2, or 3)
- **Returns:** List with `mu_est`, `sigma_est`, `mu_CI`, `sigma_CI`, posterior samples

**`get_type1_indices(n, S)`** - Calculate Type 1 quantile indices
- **Returns:** List of order statistic positions using `ceiling(p*n)`

**`compute_logpost_grid_fast(mu_grid, sigma_grid, X, n, S, L, U)`** - Vectorized log-posterior
- Uses efficient matrix operations for speed
- Incorporates both likelihood and priors

**`get_posterior_samples_adaptive(X, n, S, ...)`** - Two-stage adaptive sampling
- Automatically determines grid bounds
- Focuses computation on high-density regions

### Helper Functions

**Performance Metrics:**
- `SError(estimate, true)` - Standardized error
- `ARE(method1, method2)` - Average relative efficiency

**Data Generation:**
- `generate_order_stats_normal()` - Simulate from normal distribution
- `generate_order_stats_t()` - Simulate from Student's t
- `generate_order_stats_lognormal()` - Simulate from log-normal

**Method Comparison:**
- `get_comparison_estimates()` - Run Luo, Wan, BC, QE, MLN, BLUE, Shi methods

### Visualization Functions

**RMSE Plots:**
- `plot_rmse_by_scenario(df, dist)` - Compare methods across scenarios
- `plot_rmse_by_case(df, dist)` - Compare across parameter combinations

**Coverage Plots:**
- `plot_coverage_by_scenario(df, dist, case_num)` - Boxplots by scenario
- `plot_coverage_by_ntype(df, dist, case_num, scenario_num)` - By sample size type
- `create_all_coverage_plots(df, dist, output_dir)` - Batch generation

**Contour Plots:**
- `create_contour_plot(result_list, study_name, group_type)` - Single posterior
- `run_contour_analysis(study_data, n_studies)` - Batch with matched axes

### Examples

- **`inst/examples/quick_start.R`** - 5-minute introduction (3 examples)
- **`inst/examples/full_workflow.R`** - Complete analysis (7 steps)

---

## 🗂️ Repository Structure

```
BOSE-bayesian-order-statistics-estimator/
├── README.md                           # This file
├── LICENSE                             # MIT License
├── .gitignore                          # Git configuration
│
├── R/                                  # Core statistical methods
│   ├── bose_core.R                     # Main BOSE estimation (401 lines)
│   └── helper_functions.R              # Utilities and metrics (204 lines)
│
├── simulations/                        # Validation studies
│   ├── run_simulation_type1.R          # Monte Carlo simulations (292 lines)
│   └── results/                        # Output directory
│       └── .gitkeep
│
├── case_studies/                       # Real data applications
│   ├── skewness_assessment.R           # Distributional checks (364 lines)
│   ├── real_data_analysis.R            # Cardiovascular studies (352 lines)
│   └── contour_plots.R                 # Posterior visualization (390 lines)
│
├── visualization/                      # Publication-quality plots
│   ├── plot_rmse.R                     # Method comparison (430 lines)
│   └── plot_coverage.R                 # CI coverage assessment (443 lines)
│
├── inst/                               # Examples and documentation
│   └── examples/
│       ├── quick_start.R               # 5-minute intro
│       └── full_workflow.R             # Complete analysis
│
└── data/                               # Example datasets
    └── example_studies.csv             # Cardiovascular studies (10 studies)
```

---

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact

**Wenqisi (Lydia) Pan**  
PhD Candidate, Statistics
University of Texas at Arlington  
📧 wenqisi.pan@uta.edu  
🔗 [GitHub: @Lydia0205](https://github.com/Lydia0205)

**Advisor: Dr. Xinlei Wang**  
Jenkins-Garrett Professor  
Department of Mathematics  
University of Texas at Arlington

---

## 🙏 Acknowledgments

- **Dr. Xinlei Wang** for invaluable guidance and mentorship throughout this project
- **Dr. Zeyu Lu** for helpful discussions and assistance during development
- Developers of **`estmeansd`** and **`metaBLUE`** packages for providing comparison methods

---

## 📚 References

1. **Luo, D., Wan, X., Liu, J., & Tong, T. (2018).** Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. *Statistical Methods in Medical Research*, 27(6), 1785-1805.

2. **Wan, X., Wang, W., Liu, J., & Tong, T. (2014).** Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. *BMC Medical Research Methodology*, 14(1), 135.

3. **Shi, J., Luo, D., Weng, H., et al. (2020).** Optimally estimating the sample standard deviation from the five-number summary. *Research Synthesis Methods*, 11(5), 641-654.

4. **Hozo, S. P., Djulbegovic, B., & Hozo, I. (2005).** Estimating the mean and variance from the median, range, and the size of a sample. *BMC Medical Research Methodology*, 5(1), 13.

---

## 📝 Citation

If you use BOSE in your research, please cite:

```bibtex
@article{pan2025bose,
  title={BOSE: A Bayesian Order Statistics Estimator for Recovering the Sample Mean and Standard Deviation},
  author={Pan, Wenqisi and Wang, Xinlei},
  journal={In Preparation},
  year={2025},
  note={GitHub: https://github.com/Lydia0205/BOSE-bayesian-order-statistics-estimator}
}
```

---

## ⚠️ Important Notes

- **Research Software:** While extensively validated through simulations, please verify results for your specific application
- **Meta-Analysis Use:** For production meta-analyses, we recommend consulting with a statistician
- **Skewness Check:** Always assess distributional assumptions before applying BOSE

---

## 🚀 Getting Started Checklist

- [ ] Clone repository and install dependencies
- [ ] Run `inst/examples/quick_start.R` to see basic functionality
- [ ] Load `data/example_studies.csv` to understand data format
- [ ] Run `inst/examples/full_workflow.R` for complete analysis
- [ ] Adapt code for your own meta-analysis data
- [ ] Check `case_studies/` for real-world application examples
- [ ] Generate publication plots using `visualization/` functions

---

**Ready to use BOSE for your meta-analysis? Start with the Quick Start section above!** 🎯
