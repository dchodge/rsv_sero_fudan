model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3_hpc", "base_hier_2", paste0("model_summary.RDS")))

############################################################
################ INFECTION RATES ########################
############################################################

fitfull <- model_summary$fit    
outputfull <- model_summary$postc

filename <- outputfull$filename
modelname <- outputfull$modelname

n_chains <- outputfull$n_chains
n_post <- outputfull$n_post

chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

data_t <- fitfull$data_t

cop_raw <- readRDS(here::here("outputs", "fits", "fudan_e3_hpc", "base_hier_2", "figs", "post", "plt_data", "cop_data.RDS"))
cop_age <- cop_raw[[1]] %>% left_join(data_t$raw_sero %>% select(id, age_group) %>% unique) %>% filter(titre_val > 3)


require(cmdstanr)
require(posterior)
require(tidybayes)
cop_age_s <- cop_age %>% filter(age_group == "≤5")

data_list_A <- list(
    N = nrow(cop_age_s),
    x = cop_age_s$titre_val,
    y = cop_age_s$prop
)

log_curve_fit <- cmdstan_model( here::here("src", "logit_curve_fit.stan") )
post_A <- log_curve_fit$sample(data = data_list_A, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)

y_hat_new_post_full <- post_A$draws() %>% spread_draws(y_hat_rel[i], y_prot[i], y_hat_new[i], x_new[i]) 


y_hat_new_post_full %>% 
    ggplot() +
    stat_lineribbon(aes(x = x_new, y = y_hat_new), .width = 0.95,alpha = 0.3) + 
    geom_point(data = cop_age_s, aes(x = titre_val, y = prop), alpha = 0.6, size = 2) +
    facet_wrap(vars(biomarker)) + theme_minimal() +
    guides(color = "none", fill = "none") +
#scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "Log titre", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker") + 
    ggtitle("D. Fitted curves for infection risk")




require(cmdstanr)
require(posterior)
require(tidybayes)
cop_age_s <- cop_age 

data_list_A <- list(
    N = nrow(cop_age_s),
    x = cop_age_s$titre_val,
    y = cop_age_s$prop,
    x_cov = cop_age_s$age_group,
    N_cov = 5
)

log_curve_fit <- cmdstan_model( here::here("src", "logit_curve_fit_covar.stan") )
post_A <- log_curve_fit$sample(data = data_list_A, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)
