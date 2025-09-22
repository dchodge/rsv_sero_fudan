model_type <- "base" # base

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}


post_base <- readRDS(here::here("outputs", "fits_stan", "base", "cop_data_post.RDS"))
post_age <- readRDS(here::here("outputs", "fits_stan", "age_cov", "cop_data_post.RDS"))

# Extract diagnostics - use the correct method for parameter-level diagnostics
post_base_log <- post_base$draws("log_lik") 
post_age_log <- post_age$draws("log_lik") 

require(loo)
loo_base <- loo(post_base_log)
loo_age <- loo(post_age_log)

loo_compare <- loo_compare(loo_base, loo_age) # Proof that age-statifiaction is better


# prior predictive distributions
titre_range <- post_base$draws() %>% spread_draws(x_new[i]) %>% pull(x_new) %>% unique %>% range
x_new <- seq(titre_range[1], titre_range[2], length.out = 20)
midpoint = (titre_range[1] + titre_range[2]) / 2;

L_sample <- runif(100, 0, 1)
k_sample <- rnorm(100, 0, 5)
x0_sample <- rnorm(100, midpoint, midpoint/2)

prior_pred_base <- map(1:20, function(i) {

  data.frame(
      x_new = x_new[i],
      inf_risk = L_sample * (1 - 1 / (1 + exp(-k_sample * (x_new[i] - x0_sample)))),
      cop_prot = (1 / (1 + exp(-k_sample * (x_new[i] - x0_sample))))
    )
  }
) %>% bind_rows

p1 <- prior_pred_base %>% 
  ggplot() + stat_lineribbon(aes(x = x_new, y = inf_risk), .width = 0.95, alpha = 0.3) + theme_minimal() + 
  labs(x = "PreF titre (log10)", y = "Posterior probability of infection", fill = "95% credible interval")

p2 <- prior_pred_base %>% 
  ggplot() + stat_lineribbon(aes(x = x_new, y = cop_prot), .width = 0.95, alpha = 0.3) + theme_minimal() +
  labs(x = "PreF titre (log10)", y = "Posterior probability of protection", fill = "95% credible interval")

p1 / p2 + plot_layout(guides = "collect")
ggsave(here::here("outputs", "fits_stan", "base", "figs", "cop_prior_pred.png"), width = 9, height = 12, dpi = 300)
ggsave(here::here("outputs", "fits_stan", "base", "figs", "cop_prior_pred.pdf"),  width = 9, height = 12)



L_global_sample <- rnorm(100, 0, 1.5)
k_global_sample <- rnorm(100, 0, 5)
x0_global_sample <- rnorm(100, midpoint, midpoint/2)

L_deviation_sample <- rnorm(100, 0, 1)
k_deviation_sample <- rnorm(100, 0, 1)
x0_deviation_sample <- rnorm(100, 0, 1)

L_sigma_sample <- rexp(100, 1)
k_sigma_sample <- rexp(100, 1)
x0_sigma_sample <- rexp(100, 1)

L_sample_age <- inv_logit(L_global_sample + L_deviation_sample * L_sigma_sample)
k_sample_age <- k_global_sample + k_deviation_sample * k_sigma_sample
x0_sample_age <- x0_global_sample + x0_deviation_sample * x0_sigma_sample

prior_pred_age <- map(1:20, function(i) {

  data.frame(
      x_new = x_new[i],
      inf_risk = L_sample_age * (1 - 1 / (1 + exp(-k_sample_age * (x_new[i] - x0_sample_age)))),
      cop_prot = (1 / (1 + exp(-k_sample_age * (x_new[i] - x0_sample_age))))
    )
  }
) %>% bind_rows

p1 <- prior_pred_age %>% 
  ggplot() + stat_lineribbon(aes(x = x_new, y = inf_risk), .width = 0.95, alpha = 0.3) + theme_minimal() + 
  labs(x = "PreF titre (log10)", y = "Posterior probability of infection", fill = "95% credible interval")

p2 <- prior_pred_age %>% 
  ggplot() + stat_lineribbon(aes(x = x_new, y = cop_prot), .width = 0.95, alpha = 0.3) + theme_minimal() +
  labs(x = "PreF titre (log10)", y = "Posterior probability of protection", fill = "95% credible interval")

p1 / p2 + plot_layout(guides = "collect")
ggsave(here::here("outputs", "fits_stan", "age_cov", "figs", "cop_prior_pred.png"), width = 9, height = 12, dpi = 300)
ggsave(here::here("outputs", "fits_stan", "age_cov", "figs", "cop_prior_pred.pdf"),  width = 9, height = 12)