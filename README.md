<div style="font-family:'Times New Roman', serif;">

<h1>Final Project – Multi-Path Stepwise Selection with AIC</h1>
<h3>STAT Programming with R — Group 11</h3>

<p><strong>Authors:</strong> Mark Philip Castillo, Prince Mensah Ansah, Jonathan Ng</p>

<h2>Overview</h2>

<p>This repository contains our implementation of the Multi-Path Stepwise AIC Model Selection framework. The project includes all required components:</p>

<ul>
  <li><strong>Algorithm 1 – Multi-Path Forward Selection</strong><br>Builds several forward-selection paths and avoids the limitations of greedy search.</li>

  <li><strong>Algorithm 2 – Stability via Resampling</strong><br>Uses bootstrap resampling to estimate how consistently each predictor appears across model paths.</li>

  <li><strong>Algorithm 3 – Plausible Model Selection</strong><br>Combines AIC quality and variable stability to identify a final set of plausible and robust models.</li>

  <li><strong>Algorithm 4 – Full Pipeline</strong><br>Integrates Algorithms 1–3 into a single reproducible model-discovery workflow.</li>
</ul>

<h2>R Package: <code>mpssAIC</code></h2>
<p>Located in the <code>mpssAIC/</code> directory with full documentation and a vignette.</p>

<h2>Final Report</h2>
<p>Provided as <code>FinalProject.Rmd</code> and <code>FinalProject.html</code>.</p>

<h2>Installation</h2>

<pre>
Installation 

# Option 1: Using devtools (our primary method)
  
needed_pkgs <- c("devtools", "care", "corpcor")
for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

devtools::install_github(
  "R-4-Data-Science/Final_project_Group11",
  subdir = "mpssAIC",
  build_vignettes = TRUE
)

library(mpssAIC)
library(care)

# Option 2: Using remotes (also works, satisfies assignment requirement)
# install.packages("remotes")
# remotes::install_github(
#   "R-4-Data-Science/Final_project_Group11",
#   subdir = "mpssAIC",
#   build_vignettes = TRUE
# )

library(mpssAIC)

</pre>

<h2>Example Usage</h2>

<pre>
set.seed(1)
X <- data.frame(
  x1 = rnorm(50),
  x2 = rnorm(50),
  x3 = rnorm(50)
)
y <- 1 + 2*X$x1 - X$x2 + rnorm(50)

# Build model paths
forest <- build_paths(X, y)

# Compute stability
pi <- stability(X, y, B = 10)

# Determine plausible models
plausible_models(forest, pi)
</pre>

<h2>Vignette</h2>
<p>A full worked example using the diabetes dataset is available:</p>
<pre>vignette("diabetes_mpssAIC")</pre>

<h2>Project Summary</h2>
<p>
This project implements a complete model-selection workflow that explores multiple model paths, quantifies predictor stability, and identifies robust, statistically competitive models. The approach improves interpretability and addresses the limitations of traditional single-path stepwise selection.
</p>

</div>
