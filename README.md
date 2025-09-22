# RSV Serology Analysis: Fudan Study

## Overview

This repository contains the analysis of respiratory syncytial virus (RSV) serological data from the Fudan study, focusing on understanding antibody kinetics, infection dynamics, and correlates of protection (COP) across different age groups during two epidemic waves in 2023-2024.

## Research Goals

1. **Antibody Kinetics Modeling**: Characterize the temporal dynamics of RSV-specific antibody responses following infection
2. **Age-Stratified Analysis**: Understand how antibody responses vary across different age groups (≤5, 5-18, 19-59, 60-74, 75+ years)
3. **Correlates of Protection (COP)**: Develop statistical models to identify antibody titers associated with protection against RSV infection
4. **Epidemic Wave Analysis**: Compare antibody dynamics between two distinct epidemic waves (E1: early 2023, E2: late 2023-early 2024)
5. **Model Comparison**: Evaluate different hierarchical Bayesian models for their ability to capture serological dynamics

## Data Structure

### Raw Data
- **Source**: `data/fudan/data_titer_original_20241009.csv`
- **Content**: Longitudinal RSV PreF antibody titers from participants across two epidemic waves
- **Time Period**: January 2023 - June 2024
- **Measurements**: Log10-transformed antibody titers at multiple time points per individual

### Processed Data
- **E1 (First Wave)**: `data_clean/fudan/df_soro_model_e1.csv` - Early 2023 epidemic data
- **E2 (Second Wave)**: `data_clean/fudan/df_soro_model_e2.csv` - Late 2023-early 2024 epidemic data  
- **Combined**: `data_clean/fudan/df_soro_model.csv` - Full dataset with age group stratification
- **Exposure Priors**: `data_clean/fudan/exp_prior*.csv` - Epidemiological priors for infection timing

## Analysis Workflow

### 1. Data Preparation (`R/clean_fudan.R`)
- **Purpose**: Clean and prepare serological data for analysis
- **Key Steps**:
  - Transform raw titers to log10 scale
  - Create separate datasets for each epidemic wave
  - Calculate age groups from birth dates
  - Filter individuals with multiple measurements
  - Generate exposure priors from epidemiological data

### 2. SeroJump Model (`R/serojump_model_*.R`)
- **Purpose**: Fit hierarchical Bayesian models to serological data
- **Models**:
  - `base`: Basic model without age stratification
  - `base_hier_2`: Two-level hierarchical model
  - `base_hier_3`: Three-level hierarchical model
- **Key Parameters**:
  - Antibody kinetics parameters (peak, decay, baseline)
  - Infection timing probabilities
  - Age-specific effects

### 3. Correlates of Protection Analysis (`R/COP_*.R`)
- **Purpose**: Model the relationship between antibody titers and protection
- **Models**:
  - `logit_curve_fit.stan`: Basic logistic curve fitting
  - `logit_curve_fit_covar.stan`: Age-stratified logistic curve fitting
- **Key Outputs**:
  - Age-specific protection curves
  - Relative and absolute correlates of protection
  - Risk reduction estimates

### 4. Model Diagnostics (`R/*_post_conv.R`)
- **Purpose**: Assess model convergence and fit quality
- **Diagnostics**:
  - R_hat and ESS values
  - Trace plots for key parameters
  - Posterior predictive checks

### 5. Visualization (`R/*_post_*.R`)
- **Purpose**: Generate publication-ready figures
- **Key Figures**:
  - Antibody kinetics trajectories by age group
  - Correlates of protection curves
  - Model comparison plots
  - Convergence diagnostics

## Key Dependencies

- **R Packages**: `serojump`, `tidyverse`, `cmdstanr`, `posterior`, `tidybayes`, `patchwork`
- **Stan**: For Bayesian modeling
- **Parallel Processing**: `furrr` for multi-core computation

## Usage

1. **Setup**: Run `R/aaa_init.R` to load required packages
2. **Data Cleaning**: Execute `R/clean_fudan.R` to prepare datasets
3. **Model Fitting**: Run `R/serojump_model_run.R` to fit serological models
4. **COP Analysis**: Execute `R/COP_model.R` to analyze correlates of protection
5. **Visualization**: Run post-processing scripts to generate figures

## Model Specifications

### Serological Jump Model
- **Likelihood**: Normal distribution for antibody measurements
- **Priors**: Hierarchical structure with age-specific effects
- **Transdimensional**: Variable number of infection events per individual
- **MCMC**: 4 chains, 500,000 iterations with 250,000 burn-in

### Correlates of Protection Model
- **Function**: Logistic curve: L × (1 - 1/(1 + exp(-k × (x - x0))))
- **Parameters**:
  - L: Upper asymptote (maximum protection)
  - k: Steepness of the curve
  - x0: Midpoint (titer at 50% protection)
- **Hierarchical**: Age-specific parameters with global priors

## Outputs

### Main Results
- Age-stratified antibody kinetics trajectories
- Correlates of protection curves by age group
- Model comparison statistics
- Convergence diagnostics

### Key Findings
- Age-specific differences in antibody responses
- Optimal antibody titers for protection
- Relative effectiveness of different age groups
- Temporal dynamics across epidemic waves

## Contact

For questions about this analysis, please contact david_hodgson1991@hotmail.co.uk

## License

This analysis is part of ongoing research. Please cite appropriately if used.