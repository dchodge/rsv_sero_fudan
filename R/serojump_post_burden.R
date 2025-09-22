
model_summary <- readRDS(file = "./outputs/fits/fudan_e3_hpc/base_hier_2/model_summary.RDS")

#plotMCMCDiagnosis(model_summary, save_info = save_info)
# takes ages to run these so comment out!
#plotPostFigs(model_summary, save_info = save_info) 

fitfull <- model_summary$fit    
outputfull <- model_summary$post

filename <- outputfull$filename
modelname <- outputfull$modelname

n_chains <- outputfull$n_chains
n_post <- outputfull$n_post

chain_samples <- 1:n_chains %>% map(~c(rep(.x, n_post))) %>% unlist

data_t <- fitfull$data_t

# Make Figure 4
## Get antibody kinetics
abkin_df_ind <- readRDS(file = "./outputs/fits/fudan_e3_hpc/base_hier_2/figs/post/plt_data/ab_kinetics_recov_indiv_ High.RDS")


abkin_df_ind[[1]] <- abkin_df_ind[[1]] %>% mutate(date = days(t) + ymd("2023-01-01") )

pA <- abkin_df_ind[[1]] %>% 
    ggplot() +
        geom_line(aes(x = date, y = titre_trajectory, group = sample), color = "#960319", alpha = 0.05) + facet_wrap(vars(id)) + 
        geom_point(data = data_t$raw_sero %>% filter(id %in% unique(abkin_df_ind[[1]]$id)) %>% mutate(date = days(time) + ymd("2023-01-01") ), 
            aes(x = date, y = PreF), shape = 21, size = 3, fill = "gray") + theme_bw() + 
        labs(x = "Date", y = "PreF titre value (log10)", color = "Exposure type") + ggtitle(paste0("A. Individual-level model-fits for antibody kinetics"))  + 
        theme(text = element_text(size = 12))  + 
    theme_ft()

exptime_df <- readRDS(file = "./outputs/fits/fudan_e3_hpc/base_hier_2/figs/post/plt_data/exposure_time.RDS")


inferred_wave <- map(1:length(exptime_df[[2]]), 
    ~exptime_df[[2]][[.x]] %>% mutate(wave_number = case_when(bin < 350 ~ "Wave 1", TRUE ~ "Wave 2")) %>% mutate(date = days(round(bin, 0)) + ymd("2023-01-01"))
)

pB <- exptime_df[[1]] %>% ggplot() +
     # geom_histogram(aes(x = date, y = after_stat(count), fill = "gray40"), size = 2, alpha = 0.5, color = "black") + 
    #geom_line(data = bind_rows(hist_data), aes(x = bin, y = count), alpha = 0.1, color = "blue") +
    geom_line(data = data.frame(date = inferred_wave[[1]]$date, count = exptime_df[[3]], wave_number = inferred_wave[[1]]$wave_number), 
        aes(x = date, y = count, color = wave_number), size = 1.5) +
    geom_ribbon(data = data.frame(date = inferred_wave[[1]]$date, 
                                    lower = exptime_df[[3]] - exptime_df[[4]], 
                                    upper = exptime_df[[3]] + exptime_df[[4]], 
                                    wave_number = inferred_wave[[1]]$wave_number), 
                aes(x = date, ymin = lower, ymax = upper, fill = wave_number), 
                 alpha = 0.5) +
    scale_fill_manual(
         values =c('Wave 1'='#1f78b4','Wave 2'='#e6550d')) +
    scale_color_manual(
         values =c('Wave 1'='#1f78b4','Wave 2'='#e6550d')) +
    theme_minimal() +
            labs(x = "Date", y = "Number of people infected", fill = "") + 
            ggtitle("B. Inferred epidemic wave") + 
    guides(color = "none") +
    theme(legend.position = "bottom") + 
    theme_ft()

known_delta <- nrow(exptime_df[[1]])

df_plot <- 
    map_df(1:length(inferred_wave), ~
    data.frame(
        wave_number = c("Wave 1", "Wave 2"),
        tot = (inferred_wave[[.x]] %>% summarise(sum = sum(count), .by = wave_number) %>% pull(sum)))
    ) %>%   group_by(wave_number) %>%
        summarise(
        q025 = quantile(tot, 0.025),
        q50  = quantile(tot, 0.5),
        q975 = quantile(tot, 0.975)
    )
    
  
colnames(df_plot) <- c("wave_number", "lb", "mean", "ub")
N <- model_summary$fit$data_t$N
pC <- df_plot %>% 
    ggplot() +
  geom_bar(aes(x = wave_number, y = mean / N, fill = wave_number), size = 2, alpha = 0.5, stat = "identity", position = position_dodge(width = 0.0), width = 0.5) +
  geom_errorbar( aes(x = wave_number, ymin = lb /N, ymax = ub/ N), width = 0.4, size = 2) + 
    scale_fill_manual(
         values =c('Wave 1'='#1f78b4','Wave 2'='#e6550d')) +
  labs(title = "C. Inferred attack rates",
       x = "RSV wave",
       y = "Attack rate", fill = "")  + theme_minimal()  + 
    theme_ft() + theme(legend.position = "none")

pBC <- pB + pC + plot_layout(widths = c(2, 1))

## age group
outputfull <- model_summary$post
fit_states <- outputfull$fit_states
n_chains <- outputfull$n_chains
n_post <- outputfull$n_post
n_length <- n_chains * n_post

df_sero_model_age <- data_t$raw_sero %>% select(id, age_group) %>% unique %>% 
    mutate(age_group = recode(age_group, "<=5" = "≤5", "5-18" = "5-18", "19-59" = "19-59", "60-74" = "60-74", "75+" = "75+"))

length_i <- 4 * (fit_states %>% left_join(df_sero_model_age) %>% pull(sample_c) %>% max)


fit_states_age_1 <- fit_states %>% left_join(df_sero_model_age) %>%
    summarise(inf_tot = sum(inf_ind), n = n(), .by = c("age_group")) %>% 
    mutate(inf_rate = inf_tot / sum(n)) 

fit_states_age_alt <- fit_states %>% left_join(df_sero_model_age) %>%
    summarise(inf_rate = mean(inf_ind), n = n(), .by = c("age_group")) 

pDi <- fit_states_age_1 %>% 
    ggplot() + geom_col(aes(x = age_group, y = inf_rate))  + 
    labs(x = "Age group", y = "Attack rate (both waves)") +
    ggtitle("D. Attack rate stratified by age group")   + 
    theme_ft()
    
pDii <- fit_states_age_alt %>% 
    ggplot() + geom_col(aes(x = age_group, y = inf_rate))  + 
    labs(x = "Age group", y = "Attack rate (both waves)") +
    ggtitle("E. Attack rate per age group")   + 
    theme_ft()

pA / pBC / (pDi + pDii) + plot_layout(heights = c(2, 1, 1))
ggsave(here::here("outputs", "manu", "main", "manu_fudan_wave2.png"), width = 12, height = 15)
ggsave(here::here("outputs", "manu", "main", "manu_fudan_wave2.pdf"), width = 12, height = 15)







# Extract values for paragraph 1 of results
cat("=== VALUES FOR PARAGRAPH 1 OF RESULTS ===\n")

# Total individuals included
total_individuals <- N
cat("Total individuals included:", total_individuals, "\n")

# Calculate excluded subjects (those with only single sample)
# This would need to be calculated from your data - assuming 6 based on your text
excluded_subjects <- 6  # You may need to calculate this from your actual data
cat("Excluded subjects (single sample):", excluded_subjects, "\n")

get_median_95ci <- function(x) {
    data.frame(
        median = median(x),
        lb = quantile(x, 0.025),
        ub = quantile(x, 0.975)
    )
}

# Total infected individuals across both waves
extract_infections <- map_df(1:length(inferred_wave), ~
    data.frame(
        sample = .x,
        wave_number = c("Wave 1", "Wave 2"),
        tot = (inferred_wave[[.x]] %>% summarise(sum = sum(count), .by = wave_number) %>% pull(sum)))
) 


inf_all <- extract_infections %>% group_by(sample) %>% summarise(tot = sum(tot)) %>%
    summarise(get_median_95ci(tot))

inf_wave <- extract_infections %>% group_by(wave_number) %>% summarise(get_median_95ci(tot)) 
inf_rate_wave <- extract_infections %>% group_by(wave_number) %>% summarise(get_median_95ci(tot) / 502) 

# Wave-specific attack rates
wave1_rate <- round(inf_wave$median[inf_wave$wave_number == "Wave 1"] / N * 100, 1)
wave1_rate_lb <- round(inf_wave$lb[inf_wave$wave_number == "Wave 1"] / N * 100, 1)
wave1_rate_ub <- round(inf_wave$ub[inf_wave$wave_number == "Wave 1"] / N * 100, 1)

wave2_rate <- round(inf_wave$median[inf_wave$wave_number == "Wave 2"] / N * 100, 1)
wave2_rate_lb <- round(inf_wave$lb[inf_wave$wave_number == "Wave 2"] / N * 100, 1)
wave2_rate_ub <- round(inf_wave$ub[inf_wave$wave_number == "Wave 2"] / N * 100, 1)

# Age-specific rates
age_rates <- fit_states_age_alt %>% 
    mutate(rate_pct = round(inf_rate * 100, 1)) %>%
    select(age_group, rate_pct)

cat("Age-specific attack rates (both waves combined):\n")
print(age_rates)

# Young children rate (≤5 years)
young_children_rate <- age_rates$rate_pct[age_rates$age_group == "≤5"]
cat("Young children (≤5 years) rate:", young_children_rate, "%\n")

# Adults rates (assuming 19-59 and 60-74 are adults)
adults_19_59_rate <- age_rates$rate_pct[age_rates$age_group == "19-59"]
adults_60_74_rate <- age_rates$rate_pct[age_rates$age_group == "60-74"]
adults_75plus_rate <- age_rates$rate_pct[age_rates$age_group == "75+"]

# Calculate overall adult rate (weighted average or simple average)
adults_wave1_rate <- round(mean(c(adults_19_59_rate, adults_60_74_rate, adults_75plus_rate)), 1)

cat("Adults rate (19+ years):", adults_wave1_rate, "%\n")

# Create a summary for easy copying
cat("\n=== SUMMARY FOR MARKDOWN ===\n")
inf_info <- list(
    total_individuals = total_individuals,
    excluded_subjects = excluded_subjects,
    total_infected = inf_all$median,
    total_infected_lb = inf_all$lb,
    total_infected_ub = inf_all$ub,
    wave1_rate = wave1_rate,
    wave1_rate_lb = wave1_rate_lb,
    wave1_rate_ub = wave1_rate_ub,
    wave2_rate = wave2_rate,
    wave2_rate_lb = wave2_rate_lb,
    wave2_rate_ub = wave2_rate_ub,
    young_children_rate = young_children_rate,
    adults_19_59_rate = adults_19_59_rate,
    adults_60_74_rate = adults_60_74_rate,
    adults_75plus_rate = adults_75plus_rate,
    adults_wave1_rate = adults_wave1_rate
)

saveRDS(inf_info, file = here::here("markdown", "inf_info.RDS"))