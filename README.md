# BOSE: Bayesian Order Statistics-Based Estimator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Bayesian methodology for meta-analysis that recovers the sample mean and standard deviation from order statistics (quantiles).

BOSE provides a rigorous Bayesian framework for estimating population parameters from summary statistics commonly reported in medical and scientific literature. This is particularly valuable for meta-analyses where full datasets are unavailable.

---

## Repository Structure

-   `bose_core.R` — Shared core functions. Contains Rcpp compilation for fast log-posterior and grid computations, the two-stage adaptive grid posterior computation, and wrappers for competing estimation methods.
-   `Rshiny.R` — A Shiny web application for interactive estimation using the provided functions.
-   `README.md`

### Core Functions

-   **`run_competing_methods(X, n, S)`**
    Wrapper function to compute estimates using competing methods: LWS (Luo/Wan/Shi), BLUE, Box-Cox (BC), Quantile Estimation (QE), and MLN.
-   **`get_posterior_weights_adaptive(...)`**
    Computes the two-stage adaptive grid posterior for the BOSE method.
-   **`compute_posterior_statistics_cpp(...)`**
    C++ function via Rcpp that extracts posterior statistics (mean, median, and credible intervals) from the adaptive grid.

---

## Shiny App
To use the interactive web interface, simply run the Rshiny.R script in RStudio or standard R environment:

```r
shiny::runApp("Rshiny.R")
```

---

## 🚀 Quick Start

### Installation

```r
# Install required packages
install.packages(c("Rcpp", "tidyverse", "estmeansd", "metaBLUE", "shiny", "bslib"))
```

### 5-Minute Example

```r
# Load the core functions (this will also compile the Rcpp code)
source("bose_core.R")

# Example: Estimate mean and SD from three-number summary
# Study reports: min = 2.25, median = 16.0, max = 74.25, n = 40
X <- c(2.25, 16.0, 74.25)   # Order statistics: {min, median, max}
n <- 40                     # Sample size
S <- 1                      # Scenario 1: three-number summary

# Run BOSE estimation
estimate_bose(X, n, S)
```

**Output:**
```
Mean: 21.2665 (95% CI: [15.7232, 27.5442])
SD:   15.2626 (95% CI: [11.4644, 21.0530])
CV:   0.7414 (95% CI: [0.5365, 1.0648])
```

### 📊 Summary Statistics Scenarios Supported

BOSE handles three types of order statistics commonly reported in literature:

| Scenario | Summary Type | Input Format | Example | Use Case |
|----------|--------------|--------------|---------|----------|
| **S = 1** | Three-number | `{min, median, max}` | `c(3.2, 5.0, 6.8)` | Oldest meta-analyses, clinical trials |
| **S = 2** | Quartiles | `{Q1, median, Q3}` | `c(4.0, 5.0, 6.0)` | Modern reporting, box plots |
| **S = 3** | Five-number | `{min, Q1, median, Q3, max}` | `c(3.2, 4.0, 5.0, 6.0, 6.8)` | Complete summary, Tukey's five-number |

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact

**Wenqisi (Lydia) Pan**  
📧 wenqisi.pan@uta.edu <br>
PhD Candidate, Statistics <br>
Department of Mathematics <br>
The University of Texas at Arlington <br>

**Advisor: Dr. Xinlei Wang**  
📧 xinlei.wang@uta.edu  
Jenkins-Garrett Professor  
Department of Mathematics  
The University of Texas at Arlington

---

## 🙏 Acknowledgments

- **Dr. Xinlei Wang** for invaluable guidance and mentorship throughout this project
- **Dr. Zeyu Lu** for helpful discussions and assistance during development
- Developers of **`estmeansd`** and **`metaBLUE`** packages for providing comparison methods

---

## 📝 Citation

If you use BOSE in your research, please cite:

```bibtex
@article{pan2026bose,
  title={BOSE: A Bayesian Order Statistics Based Estimator for Recovering the Sample Mean and Standard Deviation},
  author={Pan, Wenqisi and Wang, Xinlei},
  journal={In Preparation},
  year={2026},
  note={GitHub: https://github.com/wenqisi-pan/BOSE-bayesian-order-statistics-estimator}
}
```

---

## ⚠️ Important Notes

- **Research Software:** While extensively validated through simulations, please verify results for your specific application
- **Meta-Analysis Use:** For production meta-analyses, we recommend consulting with a statistician

---

**Ready to use BOSE for your meta-analysis? Start with the Quick Start section above!** 🎯
