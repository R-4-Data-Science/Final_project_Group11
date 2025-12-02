Final Project – Multi-Path Stepwise Selection with AIC
STAT Programming with R — Group 11

Authors: Mark Philip Castillo, Prince Mensah Ansah, Jonathan

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
