# mpssAIC: Multi-Path Stepwise Selection with AIC

[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/R-4-Data-Science/Final_project_Group11)

An R package implementing multi-path forward selection using AIC, resampling-based stability estimation, and plausible model selection.

## Overview

Traditional forward selection follows a single greedy path, potentially missing competitive models. **mpssAIC** explores multiple near-optimal paths simultaneously, assesses variable stability through resampling, and identifies a "Rashomon set" of plausible models that are both:

- **Statistically sound** (low AIC)
- **Robust** (built from stable predictors)

## Installation

Install the development version from GitHub:
```r
# Install devtools if needed
install.packages("devtools")

# Install mpssAIC with vignettes
devtools::install_github(
  "R-4-Data-Science/Final_project_Group11",
  subdir = "mpssAIC",
  build_vignettes = TRUE
)

# Load the package
library(mpssAIC)
```

Other packages, such as `devtools`, `care`, and `corpcor` may be required by the package. Install them first if the installation process complains.

## Quick Start

### Example 1: Linear Regression (Gaussian)
```r
library(mpssAIC)

# Generate synthetic data
set.seed(42)
n <- 100
X <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n),
  x4 = rnorm(n)
)
y <- 2*X$x1 - 1.5*X$x2 + rnorm(n)

# Step 1: Multi-path forward selection
forest <- build_paths(
  X = X, 
  y = y,
  family = "gaussian",
  K = 4,           # Max model size
  delta = 2,       # AIC tolerance
  L = 30           # Max models per step
)

# Step 2: Compute variable stability
pi_scores <- stability(
  X = X,
  y = y,
  family = "gaussian",
  B = 30,          # Number of bootstrap resamples
  K = 4,
  delta = 2,
  L = 30
)

# Step 3: Select plausible models
plausible <- plausible_models(
  forest = forest,
  pi = pi_scores,
  Delta = 2,       # AIC window
  tau = 0.6        # Stability threshold
)

# View results
print(plausible)
```

**Expected Output:**
```
   AIC AIC_diff size      vars avg_stability
1  287.3    0.00    2  x1 + x2         0.950
2  288.1    0.80    3  x1 + x2 + x3    0.817
```

### Example 2: Logistic Regression (Binomial)
```r
# Generate synthetic binary data
set.seed(123)
n <- 150
X <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n)
)

# True model: logit(p) = 0.5 + 1.2*x1 - 0.8*x2
eta <- 0.5 + 1.2*X$x1 - 0.8*X$x2
prob <- 1 / (1 + exp(-eta))
y <- rbinom(n, 1, prob)

# Step 1: Multi-path forward selection
forest <- build_paths(
  X = X,
  y = y,
  family = "binomial",
  K = 3,
  delta = 2,
  L = 20
)

# Step 2: Compute stability
pi_scores <- stability(
  X = X,
  y = y,
  family = "binomial",
  B = 25,
  K = 3,
  delta = 2,
  L = 20
)

# Step 3: Select plausible models
plausible <- plausible_models(
  forest = forest,
  pi = pi_scores,
  Delta = 2,
  tau = 0.6
)

print(plausible)

# Evaluate best model with confusion matrix
best_idx <- as.integer(strsplit(plausible$key[1], ",")[[1]])
X_best <- X[, best_idx, drop = FALSE]
fit <- glm(y ~ ., data = data.frame(y = y, X_best), family = binomial())

# Compute confusion matrix
cm_result <- confusion_metrics(fit, y = y, cutoff = 0.5)
print(cm_result$confusion_matrix)
print(cm_result$metrics)
```

**Expected Output:**
```
        True
Predicted  0  1
        0 68  8
        1  7 67

  prevalence    accuracy sensitivity specificity         FDR         DOR 
      0.5000      0.9000      0.8933      0.9067      0.0946     73.5714
```

## Core Functions

### `build_paths()`
Multi-path forward selection exploring multiple model paths simultaneously.

**Key Parameters:**
- `family`: "gaussian" or "binomial"
- `K`: Maximum model size
- `delta`: AIC tolerance for keeping near-ties (default: 2)
- `eps`: Minimum improvement threshold (default: 1e-6)
- `L`: Max models to keep per step (default: 50)

### `stability()`
Resampling-based variable stability estimation.

**Key Parameters:**
- `B`: Number of bootstrap/subsample replicates (default: 50)
- `resample`: "bootstrap" or "subsample"

**Returns:** Stability scores π ∈ [0,1] for each predictor

### `plausible_models()`
Selects models combining AIC quality and stability.

**Key Parameters:**
- `Delta`: AIC tolerance window (default: 2)
- `tau`: Minimum average stability (default: 0.6)

**Returns:** Data frame of plausible models with AIC, size, variables, and stability

### `confusion_metrics()` (Optional Helper)
Computes confusion matrix and classification metrics for logistic models.

**Returns:** List with confusion matrix and metrics (accuracy, sensitivity, specificity, FDR, DOR)

## Vignettes

View detailed examples on real datasets:
```r
# Diabetes progression example (regression)
vignette("diabetes_mpssAIC", package = "mpssAIC")

# Browse all vignettes
browseVignettes("mpssAIC")
```

## Parameter Guidance

| Parameter | Typical Value | Description |
|-----------|---------------|-------------|
| K | min(p, 10) | Maximum model size |
| eps | 1e-6 | Minimum AIC improvement |
| delta (δ) | 2 | AIC tolerance for near-ties |
| L | 25-100 | Max models per step |
| B | 50-100 | Number of resamples |
| Delta (Δ) | 2 | AIC window for plausibility |
| tau (τ) | 0.6 | Stability threshold |

## Why Multi-Path Selection?

**Traditional Forward Selection Problems:**
- Follows only one greedy path
- Misses competitive alternative models
- Unstable under small data perturbations

**mpssAIC Advantages:**
- Explores multiple near-optimal paths
- Quantifies variable stability via resampling
- Returns a "Rashomon set" of trustworthy models
- More robust model selection

## Citation

If you use this package, please cite:
```
Philip, M., Ansah, P.M., & Ng, J. (2025). mpssAIC: Multi-Path Stepwise 
Selection with AIC. R package version 0.1.0.
https://github.com/R-4-Data-Science/Final_project_Group11
```

## License

MIT License

## Authors

- Mark Philip
- Prince Mensah Ansah
- Jonathan Ng

**Repository:** [R-4-Data-Science/Final_project_Group11](https://github.com/R-4-Data-Science/Final_project_Group11)

## References

Efron, B., Hastie, T., Johnstone, I., & Tibshirani, R. (2004). Least angle regression. *The Annals of Statistics*, 32(2), 407-499.
```

---

## **FIX #4: Update DESCRIPTION File**

### **Update: `mpssAIC/DESCRIPTION`**
```
Package: mpssAIC
Title: Multi-Path Stepwise Selection with AIC
Version: 0.1.0
Authors@R: c(
    person("Mark", "Philip", email = "mark@example.com", role = c("aut", "cre")),
    person("Prince Mensah", "Ansah", email = "prince@example.com", role = "aut"),
    person("Jonathan", "Ng", email = "jonathan@example.com", role = "aut")
  )
Description: Implements multi-path forward selection using the Akaike Information
    Criterion (AIC), resampling-based stability estimation, and plausible model 
    selection. Explores multiple near-optimal model paths simultaneously, avoiding 
    the greedy trap of traditional stepwise selection. Computes variable stability 
    scores via bootstrap or subsample resampling. Returns a "Rashomon set" of models 
    that are both statistically sound (low AIC) and robust (high stability).
License: MIT + file LICENSE
Encoding: UTF-8
LazyData: true
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.0
Depends: 
    R (>= 3.5.0)
Imports:
    stats,
    graphics
Suggests:
    care,
    knitr,
    rmarkdown,
    testthat (>= 3.0.0)
VignetteBuilder: knitr
URL: https://github.com/R-4-Data-Science/Final_project_Group11
BugReports: https://github.com/R-4-Data-Science/Final_project_Group11/issues
