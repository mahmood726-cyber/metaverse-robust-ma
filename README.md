
# metaverse

Versatile evidence synthesis with robust estimation, high-dimensional variable
selection, and modern inference methods for meta-analysis.

## Overview

**metaverse** provides cutting-edge statistical methods for robust meta-analysis
in R:

- **Robust Estimation**: M, MM, S, and tau estimators with contamination models
  (t-mixture, slash, Box-Cox)
- **Variable Selection**: Knockoff filters (FDR-controlled), spike-and-slab
  priors (Bayesian), LASSO/elastic net/SCAD/MCP, and e-values
- **Modern Inference**: Conformal prediction, selective (post-selection)
  inference, permutation tests, Bayesian Gibbs sampling
- **Sensitivity Analysis**: Egger/Begg/Thompson-Sharp tests, trim-and-fill,
  p-curve, fragility index, E-values, influence diagnostics, Cook's distance,
  outlier detection
- **Visualization**: Forest, funnel, radial (Galbraith), L'Abbe, Baujat,
  influence, cumulative, and diagnostic plots (ggplot2-based)
- **Causal Helpers**: Heterogeneous treatment effect estimation, subgroup
  analysis, transportability analysis
- **Utilities**: Effect size conversion (SMD/OR/RR/RD/r/z/g/f/eta2), power
  analysis, dependent-effects aggregation, network data preparation

## Installation

```r
# Install development version from GitHub
# install.packages("devtools")
devtools::install_github("mahmood726-cyber/metaverse")

library(metaverse)
```

## Quick Start

```r
# Simulate meta-analysis data with 10% outlier contamination
set.seed(42)
data <- simulate_meta_data(k = 30, theta = 0.5, tau2 = 0.1,
                           contamination = 0.1)

# Robust meta-analysis with MM-estimator
result <- meta_robust(data$yi, data$vi, method = "MM")
print(result)

# With t-mixture contamination model
result_contam <- meta_robust(data$yi, data$vi, contamination = "t-mixture")

# Forest plot
forest_plot_complete(result)

# Funnel plot with contours
funnel_plot_complete(result)
```

### Variable Selection

```r
# Generate data with 5 moderators
data_mod <- simulate_meta_data(k = 50, moderators = 5)
X <- as.matrix(data_mod[, grep("mod", names(data_mod))])

# Knockoff filter with FDR control at 10%
selected <- select_moderators(X, data_mod$yi, data_mod$vi,
                             method = "knockoff", fdr = 0.1)
print(selected)

# Spike-and-slab Bayesian variable selection
selected_bayes <- select_moderators(X, data_mod$yi, data_mod$vi,
                                    method = "spike_slab")
```

### Sensitivity Analysis

```r
# Using the bundled evaluation dataset
data(metaverse_eval_data)

# Build model object
model <- new("metaverse",
             data = metaverse_eval_data[, c("yi", "vi", "study_id")],
             robust = meta_robust(metaverse_eval_data$yi,
                                  metaverse_eval_data$vi))

# Run multiple sensitivity analyses
sens <- assess_sensitivity(model,
                           analyses = c("publication_bias", "influence",
                                       "fragility", "evalues", "outliers"))
```

### Effect Size Conversion

```r
# Convert SMD to log odds ratio
convert_effect_sizes(0.5, vi = 0.04, from = "SMD", to = "logOR")

# Convert correlation to SMD
convert_effect_sizes(0.3, vi = 0.01, from = "r", to = "SMD")
```

### Power Analysis

```r
# How many studies for 80% power?
power_analysis(theta = 0.3, tau2 = 0.05, power = 0.80)

# Power with 20 studies
power_analysis(k = 20, theta = 0.3, tau2 = 0.05)
```

## Bundled Data

The package includes `metaverse_eval_data`, a simulated dataset of 60 studies
with 10% outlier contamination and 10 potential moderators (2 true, 8 noise).
This dataset is used in the manuscript and vignettes.

```r
data(metaverse_eval_data)
str(metaverse_eval_data)
```

## Key Dependencies

- **robustbase**: Robust regression (S- and MM-estimators)
- **glmnet**: Penalized regression (LASSO, elastic net)
- **lme4**: Mixed-effects models for multilevel meta-analysis
- **ggplot2**: Publication-quality plotting
- **boot**: Bootstrap confidence intervals

## Citation

If you use metaverse in your research, please cite:

```
@Manual{ahmad2025metaverse,
  title = {metaverse: Versatile Evidence Synthesis with Robust Estimation},
  author = {Mahmood Ahmad},
  year = {2025},
  note = {R package version 0.1.0},
  url = {https://github.com/mahmood726-cyber/metaverse}
}
```

## License

GPL (>= 3)
