# Final Project – Multi-Path Stepwise Selection with AIC
STAT Programming with R — Group 11

Authors: Mark Philip Castillo, Prince Mensah Ansah, Jonathan Ng

Overview

This repository contains our implementation of the Multi-Path Stepwise AIC Model Selection framework. The project includes all required components:

Algorithm 1 – Multi-Path Forward Selection
Builds several forward-selection paths and avoids the limitations of greedy search.

Algorithm 2 – Stability via Resampling
Uses bootstrap resampling to estimate how consistently each predictor appears across model paths.

Algorithm 3 – Plausible Model Selection
Combines AIC quality and variable stability to identify a final set of plausible and robust models.

Algorithm 4 – Full Pipeline
Integrates Algorithms 1–3 into a single reproducible model-discovery workflow.

R Package: mpssAIC
Located in the mpssAIC/ directory with full documentation and a vignette.

Final Report
Provided as FinalProject.Rmd and FinalProject.html.



# Installation

To install the mpssAIC package directly from GitHub:

# install.packages("remotes")
remotes::install_github(
  "R-4-Data-Science/Final_project_Group11",
  subdir = "mpssAIC",
  build_vignettes = TRUE
)


Then load the package:

library(mpssAIC)

Example Usage
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

Vignette

A full worked example using the diabetes dataset is available:

vignette("diabetes_mpssAIC")

Project Summary

This project implements a complete model-selection workflow that:

explores multiple model paths,

quantifies predictor stability, and

identifies robust, statistically competitive models.

The approach improves interpretability and addresses the limitations of traditional single-path stepwise selection.
