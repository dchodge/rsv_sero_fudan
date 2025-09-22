# =============================================================================
# RSV Serology Analysis: Correlates of Protection (COP) Modeling
# =============================================================================
# This script analyzes the relationship between RSV antibody titers and protection
# against infection, fitting age-stratified logistic curves to estimate correlates
# of protection across different age groups.

# =============================================================================
# Load Serological Jump Model Results
# =============================================================================

# Load the fitted serological jump model results
# This contains the posterior estimates of infection timing and antibody kinetics
model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier_2", paste0("model_summary.RDS")))

# =============================================================================
# Extract Model Components and Prepare Data
# =============================================================================

# Extract the fitted model object and posterior samples
fitfull <- model_summary$fit      # Fitted model object
outputfull <- model_summary$post  # Posterior samples

# Get MCMC chain information
N_chains <- outputfull$n_chains   # Number of MCMC chains
n_post <- outputfull$n_post       # Number of posterior samples per chain

# Create chain identifier vector for plotting
chain_samples <- 1:N_chains %>% map(~c(rep(.x, n_post))) %>% unlist

# Extract original data used in the model
data_t <- fitfull$data_t

# =============================================================================
# Prepare Correlates of Protection Data
# =============================================================================

# Load the correlates of protection data generated from the serological model
# This contains the relationship between antibody titers and infection probability
cop_raw <- readRDS(here::here("outputs", "fits", "fudan_e3_hpc", "base_hier_2", "figs", "post", "plt_data", "cop_data.RDS"))

# Merge with age group information and filter for valid titers
# Filter out very low titers (log10 < 3) as these are below the detection limit
cop_age <- cop_raw[[1]] %>% 
    left_join(data_t$raw_sero %>% select(id, age_group) %>% unique) %>% 
    filter(titre_val > 3)  # Remove titers below detection limit


# =============================================================================
# Load Required Packages for Stan Modeling
# =============================================================================

require(cmdstanr)    # R interface to Stan
require(posterior)   # Posterior analysis tools
require(tidybayes)   # Tidy Bayesian analysis

# =============================================================================
# Prepare Data for Stan Models
# =============================================================================

# Create working copy of the COP data
cop_age_s <- cop_age 

# Prepare data list for Stan models
# This format is required by Stan for hierarchical modeling
data_list_A <- list(
    N = nrow(cop_age_s),                                    # Number of observations
    x = cop_age_s$titre_val,                               # Antibody titers (log10)
    y = cop_age_s$prop,                                    # Infection probability
    x_cov = cop_age_s$age_group %>% as.numeric,           # Age group as numeric (1-5)
    N_cov = 5                                              # Number of age groups
)

# =============================================================================
# Fit Basic Logistic Curve Model (No Age Stratification)
# =============================================================================

# Compile the basic Stan model
log_curve_fit <- cmdstan_model(here::here("src", "logit_curve_fit.stan"))

# Fit the model using MCMC sampling
# - 4 parallel chains for faster computation
# - 1000 warmup iterations to allow chains to converge
# - 1000 sampling iterations for posterior estimates
post_A <- log_curve_fit$sample(
    data = data_list_A, 
    parallel_chains = 4, 
    iter_warmup = 1000, 
    iter_sampling = 1000,
    adapt_delta = 0.99
)

# Save the fitted model for later analysis
post_A$save_object(here::here("outputs", "fits_stan", "base", "cop_data_post.RDS"))

# =============================================================================
# Fit Age-Stratified Logistic Curve Model
# =============================================================================

# Compile the hierarchical Stan model with age group effects
log_curve_fit <- cmdstan_model(here::here("src", "logit_curve_fit_covar.stan"))

# Fit the age-stratified model
# This model allows different logistic curve parameters for each age group
post_A <- log_curve_fit$sample(
    data = data_list_A, 
    parallel_chains = 4, 
    iter_warmup = 1000, 
    iter_sampling = 1000,
    adapt_delta = 0.99
)

# Save the age-stratified model results
post_A$save_object(here::here("outputs", "fits_stan", "age_cov", "cop_data_post.RDS"))