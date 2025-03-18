filename <- "hpc/nih_2024_inf/p3"
n_chains <- 4
require(readxl)
require(lubridate)
require(posterior)
require(bayesplot)
require(ggdist)
library(tidyverse)

model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3_hpc", "base_hier_2", paste0("model_summary.RDS")))

############################################################
################ INFECTION RATES ########################
############################################################


fitfull <- model_summary$fit    
outputfull <- model_summary$post

filename <- outputfull$filename
modelname <- outputfull$modelname

n_chains <- outputfull$n_chains
n_post <- outputfull$n_post

chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

data_t <- fitfull$data_t
N <- data_t$N
T_max <- (data_t$endTitreTime %>% max) + 10

post <- fitfull$post
par_tab <- fitfull$par_tab
model_outline <- fitfull$model
post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))



# testing hierarchical model
logit_inverse <- function(x) {
    exp(x) / (1 + exp(x))
}

# Calculate extracted parameters list based on hierarchy flag
calculate_param_list <- function(params, hier_flag, model, id_key) {
  if (isTRUE(hier_flag)) {
  #  hier_data <- model$abkineticsModel[[id_key]]$dataHier
    hier_params <- model$abkineticsModel[[id_key]]$parsHier
    base_params <- model$abkineticsModel[[id_key]]$parsBase
    num_groups <- model$abkineticsModel[[id_key]]$dataHierN
    
    param_list <- map(seq_len(num_groups), function(i) {
      map(base_params, ~ if (. %in% hier_params) c(., paste0("z_", ., "_", i), paste0("sigma_", .)) else .) %>% unlist()
    })
  } else {
    param_list <- list(params)
  }
  return(param_list)
}


# Extract parameters based on hierarchy
extract_parameters <- function(model_outline, id_key, samples, param_list, sample_idx, group_idx) {
    hier_flag <- model_outline$abkineticsModel[[id_key]]$hierFlag
    param_list_group <- param_list[[group_idx]]

  if (isTRUE(hier_flag)) {
    hier_params <- model_outline$abkineticsModel[[id_key]]$parsHier
    base_params <- model_outline$abkineticsModel[[id_key]]$parsBase
    return(apply_hierarchical_adjustment(samples, base_params, hier_params, model_outline, group_idx)[sample_idx, ])
  }
  return(samples[sample_idx, param_list[[group_idx]], drop = FALSE])
}

# Adjust hierarchical parameters
apply_hierarchical_adjustment <- function(post_fit, base_params, hier_params, model_outline, group_idx ) {
    #post_fit = mock_post_fit_1
    #param_list = result_hier[[group_idx]]
    #base_params = mock_model_1$abkineticsModel[[id_key]]$parsBase
    #hier_params =  mock_model_1$abkineticsModel[[id_key]]$parsHier
    #model_summary = mock_model_summary
  
  for (param in hier_params) {
   #param <- hier_params[1]
    boundaries <- model_outline$infoModel$logitBoundaries %>%
      filter(par_name == param)
    upper <- boundaries$ub
    lower <- boundaries$lb

    param_names <- c(param, paste0("z_", param, "_", group_idx), paste0("sigma_", param) )
    
    post_fit <- post_fit %>% mutate(
      !!sym(param_names[1]) := logit_inverse(
        !!sym(param_names[1]) +
        !!sym(param_names[2]) *
        !!sym(param_names[3])
      ) * (upper - lower) + lower
    )
  }
  return(post_fit %>% select(all_of(base_params)))
}


df_trajectory_full <- map_df(1:length(model_outline$abkineticsModel),
    function(name1) {
        pars_extract <- model_outline$abkineticsModel[[name1]]$pars
        functionalForm <- model_outline$abkineticsModel[[name1]]$funcForm
        biomarker <- model_outline$abkineticsModel[[name1]]$biomarker
        exposureType <- model_outline$abkineticsModel[[name1]]$exposureType

        # compare <- bind_rows(
        #     post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame  %>% pivot_longer(everything(), names_to = "param", values_to = "value") %>%
        #         mutate(type = "Posterior distribution") %>% filter(param %in% pars_extract)
#
        # )

        ab_function <- function(T, pars) {
            1:T %>% map( 
                ~functionalForm(0, .x, pars)
            ) %>% unlist
        }

        ## hierarchical effects
        hierFlag_value <- model_outline$abkineticsModel[[name1]]$hierFlag

        pars_extract_list <- calculate_param_list(pars_extract, hierFlag_value, model_outline, name1)

        df_trajectories <- map_df(1:length(pars_extract_list), 
            function(j) {
                df_adjusted_par <- extract_parameters(model_outline, name1, post_fit, pars_extract_list, 1:nrow(post_fit), j)# %>% mutate(covar_value = j, exposure_type = exposureType, biomarker = biomarker)
                T <- T_max
                traj_post <- 1:(1000) %>% purrr::map_df(
                    ~data.frame(
                        time = 1:T,
                        value = ab_function(T, as.numeric(df_adjusted_par[.x, ])),
                        sample = .x
                    )
                )  %>% group_by(time) %>% mutate(exposure_type = exposureType) %>% 
                            mutate(biomarker = biomarker, covar = j)
            }
        ) 
    }
)

df_trajectory_summaries <- df_trajectory_full  %>% group_by(time, exposure_type, biomarker, covar) %>% mean_qi(value)

df_trajectory_summaries_trim <- df_trajectory_summaries %>% filter(exposure_type != "none") %>% 
    mutate(age_group = recode(covar, `1` = "<=5", `2` = "5-18", `3` = "19-59", `4` = "60-74", `5` = "75+")) %>% 
    mutate(age_group = factor(age_group, levels = c("<=5", "5-18", "19-59", "60-74", "75+"))) %>% 
    filter(time <= 380)


df_trajectory_max <- df_trajectory_full %>% group_by(exposure_type, biomarker, covar, sample) %>% filter(value == max(value)) %>% filter(exposure_type != "none")

library(pracma)
df_trajectory_auc <- df_trajectory_full %>% group_by(exposure_type, biomarker, covar, sample) %>% filter(exposure_type != "none") %>% 
    summarise(auc = trapz(value))

df_trajectory_compare <- df_trajectory_max %>% left_join(df_trajectory_auc) %>% 
    mutate(age_group = recode(covar, `1` = "<=5", `2` = "5-18", `3` = "19-59", `4` = "60-74", `5` = "75+")) %>% 
    mutate(age_group = factor(age_group, levels = c("<=5", "5-18", "19-59", "60-74", "75+"))) 

sample_fitted_states <- function(outputfull) {
    require(data.table)
    fit_states_dt <- data.table::as.data.table(outputfull$fit_states)
    fit_states_dt
}

fit_states_dt <- sample_fitted_states(outputfull)


df_difference_boosting <- fit_states_dt %>% filter(inf_ind > 0) %>% group_by(id) %>% 
    summarise(inf_time_mean = mean(inf_time), PreF_inf = mean(PreF), n = n() / 2000) %>% left_join(data_t$raw_sero) %>% 
        mutate(time_after = time - inf_time_mean) %>% filter(time_after > 0) %>% 
        mutate(preF_boost = PreF - PreF_inf)  

p1 <- df_difference_boosting %>% 
    ggplot() + 
        geom_point(aes(x = time_after, y = preF_boost, fill = age_group, alpha = n), shape = 21, size = 3) +
        geom_ribbon(data = df_trajectory_summaries_trim, aes(ymin = .lower, ymax = .upper, x = time, fill = age_group), size = 3, alpha = 0.3) + 
        geom_line(data = df_trajectory_summaries_trim, aes(y = value, x = time, color = age_group), size = 2) +
        theme_minimal() +
        labs(x = "Time post infection (days)", y = "PreF fold-rise", 
            alpha = "Posterior probability of infection", color = "Age group", fill = "Age group") +
        facet_wrap(~age_group)  + theme_ft() +  theme(legend.position = "top") + 
        scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), labels = 10^c(0, 1, 2, 3, 4, 5)) 


df_trajectory_compare_sum <- df_trajectory_compare %>% group_by(age_group) %>% summarise(auc = mean(auc), value = mean(value))
p2 <- df_trajectory_compare %>% ggplot() + 
    geom_point(aes(x = value, y = auc, color = age_group), alpha = 0.5, size = 0.8) + 
    geom_point(data = df_trajectory_compare_sum, aes(x = value, y = auc, fill = age_group), shape = 21, alpha = 1, size = 8) + 
    xlim(0.5, 2) + ylim(0, 500) + theme_bw() + 
    labs(x = "PreF max titre", y = "AUC", color = "Age group", fill = "Age group") +  theme(legend.position = "bottom") 


p2 <- df_trajectory_compare %>% 
    ggplot() + stat_pointinterval(aes(x = age_group, y = auc, color = age_group)) + theme_bw() + 
    labs(x = "Age group (years)", y = "Area Under Curve\n of antibody kinetics") + theme_ft() + theme(legend.position = "none") 

p3 <- df_trajectory_compare %>% 
    ggplot() + stat_pointinterval(aes(x = age_group, y = value, color = age_group)) + theme_bw() + theme(legend.position = "none") + 
    labs(x = "Age group (years)", y = "Peak fold-rise of antibody kinetics") + theme_ft() + theme(legend.position = "none")  + 
        scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), labels = 10^c(0, 1, 2, 3, 4, 5)) 
 

require(patchwork)

p1/ (p2 + p3) + plot_annotation(tag_levels = "A")
ggsave(here::here("outputs", "manu", "main", "age_abkinetics.png"), width = 12, height = 12)
ggsave(here::here("outputs", "manu", "main", "age_abkinetics.pdf"), width = 12, height = 12)