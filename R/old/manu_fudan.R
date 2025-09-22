modelsummary_e3 <- readRDS(here::here("outputs", "fits", "fudan_e3", "e3_base", paste0("model_summary.RDS")))

############################################################
################ INFECTION RATES ########################
############################################################

model_summary <- modelsummary_e3

fitfull <- model_summary$fit    
outputfull <- model_summary$post

fit_states <- outputfull$fit_states
filename <- outputfull$filename
modelname <- outputfull$modelname

data_t <- fitfull$data_t
N <- data_t$N

fitfull <- model_summary$fit  
outputfull <- model_summary$post
fit_states <- outputfull$fit_states
        model_outline <- fitfull$model
biomarker <- model_outline$abkineticsModel[[1]]$biomarker

fit_states_exp_prob <- fit_states %>% select(id, inf_ind) %>% summarise(inf_post = mean(inf_ind), .by = "id")

dataplt <- fit_states %>% filter(inf_ind == 1) %>% select(id, inf_time) %>% group_by(id) %>% mean_qi(inf_time) %>% left_join(fit_states_exp_prob)

datasets <- fit_states %>% filter(inf_time > -1, inf_ind == 1)  %>% arrange(sample) %>% 
    split(cumsum(c(TRUE, diff(.$sample) != 0))) %>% map(~pull(.x, inf_time) )
minmax <- datasets %>% unlist() %>% range
break_vec <- ceiling(seq(0, minmax[2], length.out = 30))

hist_data <- datasets %>%
    map(~ hist(., breaks = break_vec, plot = FALSE)) %>%
    map(~ {
        data.frame(
        count = .$counts,
        bin = .$mids
        ) %>% mutate(biomarker = biomarker)
    }) 

mean_counts <- hist_data %>%
    map(~ .x %>% pull(count)) %>%  do.call(rbind, .) %>% colMeans
std_counts <- do.call(rbind, hist_data %>% map(~ .$count)) %>% apply(2, sd)

#dataplt_inf <- fit_states %>% filter(inf_time > -1, inf_ind == 1) %>% select(id, inf_time) %>% group_by(id) %>% mean_qi(inf_time)

df_his_counts <- data.frame(bin = hist_data[[1]]$bin, count = mean_counts, lower = mean_counts - std_counts,  upper = mean_counts + std_counts) %>% 
    mutate(biomarker = biomarker) %>% mutate(epidemic = if_else(bin < 350, "Wave 1", "Wave 2"))

figC <- df_his_counts %>% ggplot() +
    geom_line(aes(x = bin, y = count, color = epidemic), size = 2) +
    geom_ribbon(
                aes(x = bin, ymin = lower, ymax = upper, fill = epidemic), 
                    alpha = 0.5) +
    theme_bw()  + 
            labs(x = "Time in study", y = "Number of people infected", fill = "") + 
            ggtitle("Recovery of infection timings") + 
    theme(legend.position = "bottom") 
ggsave(here::here("outputs", "manu", "main", "ar_2016_A.pdf"))


df_mean <- data.frame(
    s = 1:400,
    wave1 = (datasets %>% map(~sum(. < 350)) %>% unlist %>% as.numeric %>% `/`(., N)) ,
    wave2 = (datasets %>% map(~sum(. >= 350)) %>% unlist %>% as.numeric  %>% `/`(., N)) 
) %>% pivot_longer(c(wave1, wave2), names_to = "wave", values_to = "ar") %>% group_by(wave) %>% mean_qi(ar)


df_mean %>% 
    ggplot() +
  geom_bar(aes(x = wave, y = ar, fill = wave), stat = "identity", position = position_dodge(width = 0.0), alpha = 1) +
  geom_errorbar(aes(x = wave, ymin = .lower, ymax = .upper), width = 0.4, size = 2) + 
  labs(title = "Known and inferred attack rates",
       x = "RSV",
       y = "Attack rate", fill = "") +
    theme_minimal() +
    theme(text = element_text(size = 18, color = "black")) 
ggsave(here::here("outputs", "manu", "main", "ar_2016_B.pdf"))



############################################################
################ COR ########################
############################################################


fitfull <- model_summary$fit    
outputfull <- model_summary$post

fit_states <- outputfull$fit_states
filename <- outputfull$filename
modelname <- outputfull$modelname

post <- fitfull$post
data_t <- fitfull$data_t

n_chains <- outputfull$n_chains
n_post <- outputfull$n_post
n_length <- n_chains * n_post
chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

model_outline <- fitfull$model

post_fit <- post$mcmc %>% lapply(as.data.frame) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(chain = as.character(chain_samples ))

n_post <- outputfull$n_post
n_length <- n_chains * n_post

cop_exp_sum_plot_all <- map_df(1:length(model_outline$observationalModel), 
    function(i) {
        i <- 1
    biomarker <- model_outline$observationalModel[[i]]$biomarker

    df_full_info <- fit_states %>% mutate(titre_type = case_when(
        inf_time == -1 ~ "No inf",
        inf_time != -1 ~ "Inf")
    ) %>%  group_by(id, titre_type) %>% add_count(name = "n_rows") %>% 
    group_by(id, titre_type, n_rows) %>% mutate(prop = n_rows / (n_length)) %>%
    group_by(id, titre_type, prop) %>%
    mean_qi(!!sym(biomarker) ) %>% 
    complete(id = 1:data_t$N, titre_type, fill = list(prop = 0)) %>%
    arrange(!!biomarker) %>% as.data.frame 

    # Get response for each individual
    df_full_info_res <- df_full_info %>% filter(titre_type == "Inf") %>% select(id, prop)

    # Get covariate for each individual
    df_full_info_x <- df_full_info %>% group_by(id) %>% mutate(!!sym(biomarker) := weighted.mean(x = !!sym(biomarker), w = prop)) %>% 
    select(id, !!sym(biomarker))

    df_data <- data.frame(
        id = 1:data_t$N,
        start_titre = data_t$initialTitreValue %>% as.data.frame %>% pull(!!biomarker),
        known_inf = data_t$knownInfsVec
    ) %>% mutate(known_inf = ifelse(known_inf == 0, "Not known", "Known"))


    df_cor_info <- df_full_info_res %>% left_join(df_full_info_x) %>% left_join(df_data)  %>%
        rename(titre_val = !!(biomarker)) %>% mutate(biomarker = !!biomarker) 

    }) %>% unique


cop_exp_sum_plot_all %>% filter(titre_val > 3) %>%
    ggplot() + 
        geom_count(aes(x = titre_val, y = prop), alpha = 0.5) + theme_bw() + 
        geom_smooth(method = "loess", span = 0.8, aes(x = titre_val, y = prop), size = 2) + 
        facet_grid(rows = vars(biomarker)) + 
        labs(x = "Titre value at exposure", y = "Probability of infection") + 
        theme(text = element_text(size = 15))
ggsave(here::here("outputs", "manu", "main", "cor_rsv.pdf"))