# =============================================================================
# RSV Serology Analysis: Serological Jump Model Fitting
# =============================================================================
# This script fits hierarchical Bayesian serological jump models to estimate
# antibody kinetics and infection timing from longitudinal RSV serological data.
# Multiple model variants are fitted to compare different hierarchical structures.

# =============================================================================
# Load Required Packages
# =============================================================================

library(furrr)        # For parallel processing
library(serojump)     # Custom package for serological jump models
library(tidyverse)    # Data manipulation and visualization
source("R/utils.R")   # Custom utility functions

# =============================================================================
# Model 1: Basic Model (No Hierarchical Structure)
# =============================================================================

# Load the basic model specification
models <- readRDS(here::here("models", "model_base.RDS"))

# Set up parallel processing for faster computation
plan(multisession, workers = 4)

# Configure MCMC settings for the basic model
# These settings are optimized for the transdimensional nature of the model
rj_settings <- list(
    numberChainRuns = 4,        # Number of MCMC chains
    iterations = 500000,        # Total iterations per chain
    burninPosterior = 250000,   # Burn-in period (50% of iterations)
    thin = 1000,               # Thinning interval to reduce autocorrelation
    runParallel = TRUE,        # Enable parallel processing
    onDebug = FALSE            # Disable debug mode for production runs
)

# Define save information for the basic model
save_info_base <- list(
    file_name = "fudan_e3",    # Output directory name
    model_name = "base"        # Model variant identifier
)

# Fit the basic serological jump model
# This model estimates antibody kinetics without age-specific effects
model_summary <- runSeroJump(models$model[[1]], rj_settings, save_info = models$save[[1]])


# =============================================================================
# Generate Diagnostic Plots and Results for Basic Model
# =============================================================================

# Alternative: Load previously fitted model if available
# model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base", "model_summary.RDS"))

# Generate MCMC convergence diagnostics
# This creates trace plots, R_hat values, and other convergence checks
plotMCMCDiagnosis(model_summary, save_info = models$save[[1]])

# Generate posterior summary plots
# This creates antibody kinetics trajectories and other result visualizations
plotPostFigs(model_summary, save_info = models$save[[1]])

# =============================================================================
# Model 2: Two-Level Hierarchical Model
# =============================================================================

# Load the two-level hierarchical model specification
# This model includes age group effects on antibody kinetics
models <- readRDS(here::here("models", "model_base_hier_2.RDS"))

# Note: Parallel processing already set up above
# plan(multisession, workers = 4)

# Configure MCMC settings for the hierarchical model
# Reduced iterations for faster testing (can be increased for final analysis)
rj_settings <- list(
    numberChainRuns = 4,        # Number of MCMC chains
    iterations = 1000,          # Total iterations per chain (reduced for testing)
    burninPosterior = 500,      # Burn-in period (50% of iterations)
    thin = 1,                   # No thinning for shorter runs
    runParallel = TRUE,         # Enable parallel processing
    onDebug = FALSE             # Disable debug mode
)

# Define save information for the hierarchical model
save_info_base_hier <- list(
    file_name = "fudan_e3",        # Output directory name
    model_name = "base_hier_2"     # Model variant identifier
)

# Fit the two-level hierarchical serological jump model
# This model estimates age-specific antibody kinetics parameters
model_summary <- runSeroJump(models$model[[1]], rj_settings, save_info = models$save[[1]])

# Alternative: Load previously fitted model if available
# model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier", "model_summary.RDS"))

# Generate diagnostic plots and results for the hierarchical model
# Note: Requires save_info to be properly configured
plotMCMCDiagnosis(model_summary, save_info = models$save[[1]])
plotPostFigs(model_summary, save_info = models$save[[1]])


# =============================================================================
# Model 3: Three-Level Hierarchical Model
# =============================================================================

# Load the three-level hierarchical model specification
# This model includes additional hierarchical structure for more complex
# age-specific effects and individual-level variation
models <- readRDS(here::here("models", "model_base_hier_3.RDS"))

# Note: Parallel processing already set up above
# plan(multisession, workers = 4)

# Configure MCMC settings for the three-level hierarchical model
# Same settings as the two-level model for consistency
rj_settings <- list(
    numberChainRuns = 4,        # Number of MCMC chains
    iterations = 1000,          # Total iterations per chain
    burninPosterior = 500,      # Burn-in period (50% of iterations)
    thin = 1,                   # No thinning for shorter runs
    runParallel = TRUE,         # Enable parallel processing
    onDebug = FALSE             # Disable debug mode
)

# Define save information for the three-level hierarchical model
save_info_base_hier <- list(
    file_name = "fudan_e3",        # Output directory name
    model_name = "base_hier_3"     # Model variant identifier
)

# Fit the three-level hierarchical serological jump model
# This model provides the most complex hierarchical structure for comparison
model_summary <- runSeroJump(models$model[[1]], rj_settings, save_info = models$save[[1]])

# Alternative: Load previously fitted model if available
# model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier", "model_summary.RDS"))

# Generate diagnostic plots and results for the three-level hierarchical model
# This completes the model comparison analysis
plotMCMCDiagnosis(model_summary, save_info = models$save[[1]])
plotPostFigs(model_summary, save_info = models$save[[1]])
