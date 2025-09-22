# =============================================================================
# RSV Serology Analysis: Fudan Study - Package Initialization
# =============================================================================
# This script loads all required packages and sets global variables for the
# RSV serological analysis. It should be run first before any other analysis
# scripts to ensure all dependencies are available.

# Load development tools and core packages
library(devtools)          # For package development and installation
library(serojump)          # Custom package for serological jump models
library(tidyverse)         # Core data manipulation and visualization suite
require(readxl)            # For reading Excel files
require(lubridate)         # For date/time manipulation

# Bayesian analysis packages
require(posterior)         # For posterior analysis and diagnostics
require(bayesplot)         # For Bayesian model visualization
require(ggdist)            # For distribution plots
require(tidybayes)         # For tidy Bayesian analysis
require(cmdstanr)          # R interface to Stan for Bayesian modeling

# Additional utility packages
library(patchwork)         # For combining plots
library(data.table)        # For fast data manipulation

# Load custom utility functions
source(here::here("R", "utils.R"))

# =============================================================================
# Global Variables
# =============================================================================
# Set the number of MCMC chains for Bayesian model fitting
# This ensures consistency across all model runs
N_chains <- 4