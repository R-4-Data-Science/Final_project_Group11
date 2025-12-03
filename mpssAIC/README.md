# mpssAIC
Multi-Path Stepwise AIC Model Selection with Stability and Plausible Models

**mpssAIC** implements an expanded, robust version of stepwise model selection: Multi-path forward selection (explores many candidate paths instead of one)
AIC-based scoring of all visited models. Resampling-based variable stability via bootstrap/subsampling
Plausible model identification using AIC tolerance and stability thresholds. This package was developed as part of the final project for R Programming for Data Science (STAT 6210) at Auburn University.

__Install from GitHub__

You can install `mpssAIC` directly from GitHub using devtools:

```
# install.packages("devtools")
devtools::install_github(
  "R-4-Data-Science/Final_project_Group11",
  subdir = "mpssAIC",
  build_vignettes = FALSE
)

library(mpssAIC)
```

Below is a minimal, reproducible example demonstrating the full workflow:

```
library(mpssAIC)

set.seed(1)
n <- 100; p <- 5

X <- as.data.frame(matrix(rnorm(n*p), n, p))
names(X) <- paste0("x", 1:p)
beta <- c(2, -1, 0, 0, 0)
y <- as.numeric(as.matrix(X) %*% beta + rnorm(n))
```

# 1. Build multi-path forward selection forest

```
forest <- build_paths(
  X, y,
  family = "gaussian",
  K = 5,
  eps = 1e-6,
  delta = 2,
  L = 20
)
```

# 2. Compute stability of variables under resampling
```
pi_vec <- stability(
  X, y,
  family = "gaussian",
  B = 50,
  resample = "bootstrap",
  K = 5,
  eps = 1e-6,
  delta = 2,
  L = 20
)
```
# 3. Identify plausible, stable models
```
plaus <- plausible_models(
  forest,
  pi = pi_vec,
  Delta = 2,
  tau = 0.6
)
```
```
plaus
```
# Sample output:

| key |   AIC    | AIC_diff | size |   vars    | avg_stability |
|-----|----------|----------|------|-----------|----------------|
| 1,2 | 157.8527 |    0     |  2   | x1 + x2   |   0.7883333    |




