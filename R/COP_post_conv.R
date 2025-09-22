model_type <- "age_cov" # base

post_A <- readRDS(here::here("outputs", "fits_stan", model_type, "cop_data_post.RDS"))

# Load required libraries
library(posterior)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Extract diagnostics - use the correct method for parameter-level diagnostics
diagnostics <- post_A$summary()

# Check the structure of diagnostics
cat("Diagnostics structure:\n")
str(diagnostics)
cat("\nFirst few rows:\n")
print(head(diagnostics))

# Get all available parameters
all_params <- unique(diagnostics$variable)
cat("\nAll available parameters:\n")
print(all_params)

if (model_type == "age_cov") {
# Get R_hat and ESS for key parameters
  key_params <- c("L_global", "k_global", "x0_global", "sigma", 
                "L_deviation[1]", "L_deviation[2]", "L_deviation[3]", "L_deviation[4]", "L_deviation[5]",
                "k_deviation[1]", "k_deviation[2]", "k_deviation[3]", "k_deviation[4]", "k_deviation[5]",
                "x0_deviation[1]", "x0_deviation[2]", "x0_deviation[3]", "x0_deviation[4]", "x0_deviation[5]",
                "L_sigma", "k_sigma", "x0_sigma",
                "L_deviation", "k_deviation", "x0_deviation")
} else {
    key_params <- c("L", "k", "x0", "sigma")
}

# Check which key parameters are actually in the diagnostics
available_key_params <- key_params[key_params %in% all_params]
cat("\nAvailable key parameters:\n")
print(available_key_params)

# Filter diagnostics for available key parameters
if(length(available_key_params) > 0) {
  param_diagnostics <- diagnostics %>%
    filter(variable %in% available_key_params) %>%
    select(variable, rhat, ess_bulk, ess_tail)
} else {
  # If no key parameters found, use all parameters
  cat("No key parameters found, using all available parameters\n")
  param_diagnostics <- diagnostics %>%
    select(variable, rhat, ess_bulk, ess_tail)
}

# Check if the required columns exist
cat("\nAvailable columns in diagnostics:\n")
print(names(diagnostics))

# If rhat, ess_bulk, ess_tail don't exist, try alternative names
if(!"rhat" %in% names(diagnostics)) {
  # Try alternative column names
  if("Rhat" %in% names(diagnostics)) {
    param_diagnostics <- param_diagnostics %>% rename(rhat = Rhat)
  }
  if("ess_bulk" %in% names(diagnostics)) {
    # ess_bulk exists
  } else if("ess_bulk" %in% names(diagnostics)) {
    # ess_bulk exists
  } else if("ess" %in% names(diagnostics)) {
    param_diagnostics <- param_diagnostics %>% rename(ess_bulk = ess)
  }
  if("ess_tail" %in% names(diagnostics)) {
    # ess_tail exists
  } else {
    # Create a placeholder if ess_tail doesn't exist
    param_diagnostics$ess_tail <- param_diagnostics$ess_bulk
  }
}

cat("\nParameter diagnostics:\n")
print(param_diagnostics)

# Create R_hat plot
rhat_plot <- param_diagnostics %>%
  ggplot(aes(x = reorder(variable, rhat), y = rhat)) +
  geom_point(size = 3, color = "steelblue") +
  geom_hline(yintercept = 1.01, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 1.05, color = "orange", linetype = "dashed", size = 1) +
  coord_flip() +
  labs(title = "R_hat Diagnostics", 
       x = "Parameter", 
       y = "R_hat") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Create ESS plot
ess_plot <- param_diagnostics %>%
  pivot_longer(cols = c(ess_bulk, ess_tail), names_to = "ess_type", values_to = "ess_value") %>%
  ggplot(aes(x = reorder(variable, ess_value), y = ess_value, color = ess_type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 400, color = "red", linetype = "dashed", size = 1) +
  coord_flip() +
  labs(title = "Effective Sample Size (ESS)", 
       x = "Parameter", 
       y = "ESS",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

# Extract draws for trace plots
draws_df <- post_A$draws() %>%
  as_draws_df()

# Check which parameters are available in the draws
draws_params <- names(draws_df)
cat("\nAvailable parameters in draws:\n")
print(draws_params[!grepl("^\\.", draws_params)])  # Exclude .chain, .iteration, .draw

# Select available global parameters for trace plots
if (model_type == "age_cov") {
# Get R_hat and ESS for key parameters
  global_params <- c("L_global", "k_global", "x0_global", "sigma")
} else {
  global_params <- c("L", "k", "x0", "sigma")
}

available_global <- global_params[global_params %in% draws_params]

sigma_params <- c("L_sigma", "k_sigma", "x0_sigma")
available_sigma <- sigma_params[sigma_params %in% draws_params]

cat("\nAvailable global parameters for trace plots:\n")
print(available_global)
cat("\nAvailable sigma parameters for trace plots:\n")
print(available_sigma)

# Create trace plots for global parameters
if(length(available_global) > 0) {
  trace_global <- draws_df %>%
    select(.chain, .iteration, all_of(available_global)) %>%
    pivot_longer(cols = all_of(available_global), 
                 names_to = "parameter", values_to = "value") %>%
    ggplot(aes(x = .iteration, y = value, color = factor(.chain))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(title = "Trace Plots - Global Parameters",
         x = "Iteration", 
         y = "Value",
         color = "Chain") +
    theme_minimal() +
    theme(legend.position = "bottom")
} else {
  trace_global <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "No global parameters available", size = 6) +
    theme_void()
}

# Create trace plots for sigma parameters
if(length(available_sigma) > 0) {
  trace_sigma <- draws_df %>%
    select(.chain, .iteration, all_of(available_sigma)) %>%
    pivot_longer(cols = all_of(available_sigma), 
                 names_to = "parameter", values_to = "value") %>%
    ggplot(aes(x = .iteration, y = value, color = factor(.chain))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
    labs(title = "Trace Plots - Sigma Parameters",
         x = "Iteration", 
         y = "Value",
         color = "Chain") +
    theme_minimal() +
    theme(legend.position = "bottom")
} else {
  trace_sigma <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "No sigma parameters available", size = 6) +
    theme_void()
}

# Combine all plots
diagnostic_plot <- (rhat_plot | ess_plot) / 
  (trace_global / trace_sigma) +
  plot_annotation(title = "COP Model Convergence Diagnostics",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))

# Save the diagnostic plot
ggsave(here::here("outputs", "fits_stan", model_type, "figs", "cop_convergence_diagnostics.png"), 
       diagnostic_plot, width = 9, height = 12, dpi = 300)
ggsave(here::here("outputs", "fits_stan", model_type, "figs", "cop_convergence_diagnostics.pdf"), 
       diagnostic_plot, width = 9, height = 12)

# Extract sampler diagnostics (divergences, tree depth, etc.)
# Try multiple methods to extract diagnostics

# Method 1: Try to extract from draws directly
sampler_diagnostics <- post_A$sampler_diagnostics()
sampler_diagnostics_df <- sampler_diagnostics %>% as_draws_df

# Method 2: Try to access diagnostics through the fit object


# Extract diagnostic information
if (!is.null(sampler_diagnostics_df) && length(sampler_diagnostics_df) > 0) {
  cat("Sampler diagnostics structure:\n")
  str(sampler_diagnostics_df)
  
  # Handle different structures
  
  cat("Diagnostic data frame columns:", paste(names(sampler_diagnostics_df), collapse = ", "), "\n")
  
  # Extract divergence information
  if ("divergent__" %in% names(sampler_diagnostics_df)) {
    n_divergent <- sum(sampler_diagnostics_df$divergent__, na.rm = TRUE)
    total_transitions <- nrow(sampler_diagnostics_df)
    divergent_pct <- (n_divergent / total_transitions) * 100
  } else {
    n_divergent <- NA
    total_transitions <- NA
    divergent_pct <- NA
  }
} else {
  # Fallback to draws method
  cat("Falling back to draws method for diagnostics\n")
  diag_df <- draws_df
  n_divergent <- NA
  total_transitions <- NA
  divergent_pct <- NA
}

# Extract max tree depth information
if ("treedepth__" %in% names(sampler_diagnostics_df)) {
  max_treedepth <- max(sampler_diagnostics_df$treedepth__, na.rm = TRUE)
  mean_treedepth <- mean(sampler_diagnostics_df$treedepth__, na.rm = TRUE)
  # Count how many times max tree depth was hit (usually 10 or 12)
  max_depth_hits <- sum(sampler_diagnostics_df$treedepth__ >= 10, na.rm = TRUE)
  max_depth_pct <- (max_depth_hits / total_transitions) * 100
} else {
  max_treedepth <- NA
  mean_treedepth <- NA
  max_depth_hits <- NA
  max_depth_pct <- NA
}

# Extract energy information
if ("energy__" %in% names(sampler_diagnostics_df)) {
  min_energy <- min(sampler_diagnostics_df$energy__, na.rm = TRUE)
  max_energy <- max(sampler_diagnostics_df$energy__, na.rm = TRUE)
  mean_energy <- mean(sampler_diagnostics_df$energy__, na.rm = TRUE)
} else {
  min_energy <- NA
  max_energy <- NA
  mean_energy <- NA
}

# Extract step size information
if ("stepsize__" %in% names(sampler_diagnostics_df)) {
  mean_stepsize <- mean(sampler_diagnostics_df$stepsize__, na.rm = TRUE)
  final_stepsize <- sampler_diagnostics_df$stepsize__[nrow(diag_df)]
} else {
  mean_stepsize <- NA
  final_stepsize <- NA
}

# Print summary of convergence diagnostics
cat("Convergence Diagnostics Summary:\n")
cat("===============================\n")
cat("Parameters with R_hat > 1.01:", sum(param_diagnostics$rhat > 1.01, na.rm = TRUE), "\n")
cat("Parameters with R_hat > 1.05:", sum(param_diagnostics$rhat > 1.05, na.rm = TRUE), "\n")
cat("Parameters with ESS < 400:", sum(param_diagnostics$ess_bulk < 400, na.rm = TRUE), "\n")

cat("\nSampler Diagnostics:\n")
cat("====================\n")
cat("Total transitions:", ifelse(is.na(total_transitions), "N/A", total_transitions), "\n")
cat("Divergent transitions:", ifelse(is.na(n_divergent), "N/A", n_divergent), "\n")
cat("Divergent percentage:", ifelse(is.na(divergent_pct), "N/A", paste0(round(divergent_pct, 2), "%")), "\n")
cat("Max tree depth reached:", ifelse(is.na(max_treedepth), "N/A", max_treedepth), "\n")
cat("Mean tree depth:", ifelse(is.na(mean_treedepth), "N/A", round(mean_treedepth, 2)), "\n")
cat("Max depth hits (>=10):", ifelse(is.na(max_depth_hits), "N/A", max_depth_hits), "\n")
cat("Max depth hit percentage:", ifelse(is.na(max_depth_pct), "N/A", paste0(round(max_depth_pct, 2), "%")), "\n")
cat("Mean energy:", ifelse(is.na(mean_energy), "N/A", round(mean_energy, 2)), "\n")
cat("Energy range:", ifelse(is.na(min_energy) || is.na(max_energy), "N/A", 
                           paste0(round(min_energy, 2), " - ", round(max_energy, 2))), "\n")
cat("Mean step size:", ifelse(is.na(mean_stepsize), "N/A", round(mean_stepsize, 6)), "\n")
cat("Final step size:", ifelse(is.na(final_stepsize), "N/A", round(final_stepsize, 6)), "\n")


# Sentence to R markdown


# Warning messages based on diagnostics
cat("\nDiagnostic Warnings:\n")
cat("====================\n")
if (!is.na(divergent_pct) && divergent_pct > 1) {
  cat("WARNING: High divergence rate (", round(divergent_pct, 2), "%) - consider increasing adapt_delta\n")
}
if (!is.na(max_depth_pct) && max_depth_pct > 10) {
  cat("WARNING: High max tree depth hit rate (", round(max_depth_pct, 2), "%) - consider increasing max_treedepth\n")
}
if (!is.na(max_treedepth) && max_treedepth >= 10) {
  cat("WARNING: Maximum tree depth reached - consider increasing max_treedepth\n")
}

cat("\nDetailed parameter diagnostics:\n")
print(param_diagnostics)
