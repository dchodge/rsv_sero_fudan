# Load required libraries
library(dplyr)
library(ggplot2)
library(pROC)
library(here)
library(tidyr)
library(purrr)
library(patchwork)
library(gt)

get_auc <- function(true_labels, predicted_probs) {
        roc_obj <- roc(true_labels, predicted_probs)
        data.frame(
            auc = as.numeric(roc_obj$auc)
        )
}

get_coverage <- function(x, k, x0) {
     coverage_prob <- function(beta_i, x0) {
        c(x0 - qlogis(0.95) / beta_i, x0 - qlogis(0.05) / beta_i)
    }

    covarage_limits <- coverage_prob(k, x0)
    cov_per <- (x > covarage_limits[1]) & (x < covarage_limits[2])
    cov_amount <- cov_per %>% mean
}
### AMG 
get_improvement <- function(x, k, x0) {
    protection <- function(x, beta, x0) {
        1 - (plogis(beta * (x0 - x)))
    }
    protection_gain <- function(x, delta, beta, x0) {
        protection(x + delta, beta, x0) - protection(x, beta, x0)
    }

    improvement <- protection_gain(x, log10(4), k, x0) 
}

get_thresholds <- function(df_cop) {
    get_mean_95CI <- function(x) {
        data.frame(
            mean = median(x),
            lb = quantile(x, 0.025),
            ub = quantile(x, 0.975)
        )

    }
    get_boundaries <- bind_rows(
        df_cop %>% filter(y_prot > 0.49 & y_prot < 0.51) %>% ungroup %>% summarise(get_mean_95CI(x_new)) %>% mutate(y_prot = "50%"),
        df_cop %>% filter(y_prot > 0.74 & y_prot < 0.76) %>% ungroup %>% summarise(get_mean_95CI(x_new)) %>% mutate(y_prot = "75%"),
        df_cop %>% filter(y_prot > 0.89 & y_prot < 0.91) %>% ungroup %>% summarise(get_mean_95CI(x_new)) %>% mutate(y_prot = "90%")
    )
    get_boundaries
}


### Look at some important metrics to consider when we think about COP
# Load the correlates of protection data generated from the serological model
# This contains the relationship between antibody titers and infection probability
cop_raw <- readRDS(here::here("outputs", "fits", "fudan_e3_hpc", "base_hier_2", "figs", "post", "plt_data", "cop_data.RDS"))

# Merge with age group information and filter for valid titers
# Filter out very low titers (log10 < 3) as these are below the detection limit
cop_age <- cop_raw[[1]] %>% 
    left_join(data_t$raw_sero %>% select(id, age_group) %>% unique) %>% 
    filter(titre_val > 3)  # Remove titers below detection limit

cop_age_s <- cop_age %>% mutate(age_group = recode(age_group, "<=5" = "≤5", "5-18" = "5-18", "19-59" = "19-59", "60-74" = "60-74", "75+" = "75+"))


N = nrow(cop_age_s)                                    # Number of observations
x = cop_age_s$titre_val                               # Antibody titers (log10)
y = cop_age_s$prop                                    # Infection probability
x_cov = cop_age_s$age_group %>% as.numeric           # Age group as numeric (1-5)
N_cov = 5  

x_med <-  median(x)                                            # Number of age groups

cop_age_s_sum <- cop_age_s %>% mutate(titre_group = 
    cut(titre_val, breaks = seq(3, 5, 0.1), labels  = seq(3.05, 5, 0.1))) %>%
    mutate(titre_group = as.numeric(as.character(titre_group))) %>%
    summarise(prop = mean(prop), n = n(), .by = c("titre_group")) %>% arrange(titre_group)

post_A <- readRDS(here::here("outputs", "fits_stan", "age_cov", "cop_data_post.RDS"))

df_cop <- post_A$draws() %>% spread_draws(y_hat_new[i], y_prot[i], x_new[i]) 

df_thresholds <- get_thresholds(df_cop)

p1 <- df_cop %>% 
    ggplot() +
    geom_vline(xintercept = x_med, color = "red", linetype = "dashed", size = 1) +
    stat_lineribbon(aes(x = x_new, y = y_hat_new), .width = 0.95,alpha = 0.3) + 
    geom_point(data = cop_age_s_sum, aes(x = titre_group, y = prop, size = n), alpha = 0.8) +
    theme_minimal() +
    guides(color = "none", fill = "none") +
#scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "PreF titre (log10)", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker") + 
    ggtitle("A. Fitted curves for infection risk") + theme_ft()

p2 <- df_cop %>% 
    ggplot() +
    geom_vline(xintercept = x_med, color = "red", linetype = "dashed", size = 1) +
    stat_lineribbon(aes(x = x_new, y = y_prot), .width = 0.95 ,alpha = 0.3) + 
    theme_minimal() +
    guides(color = "none", fill = "none") +
#scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "PreF titre (log10)", y = "Posterior probability of protection", color = "Biomarker", fill = "Biomarker") + 
    ggtitle("B. Fitted curves for correlate of protection") + theme_ft()

p3 <- df_thresholds %>% 
    ggplot() +
    geom_point(aes(y = y_prot, x = mean), size = 5) +
    geom_errorbar(aes(y = y_prot, xmin = lb, xmax = ub)) +
    theme_minimal() +
    labs(y = "Protection threshold", x = "PreF titre (log10)") + 
    ggtitle("C. Thresholds for protection") + theme_ft()


p1 / p2 / p3
ggsave(here::here("outputs", "manu", "main", "base_cop.pdf"), width = 10, height = 12)
ggsave(here::here("outputs", "manu", "main", "base_cop.png"), width = 10, height = 12)


### METRICS
#### AUC
y_hat_est <- post_A$draws() %>% spread_draws(y_hat_est[i]) %>% summarise(y_hat_est = mean(y_hat_est))
log_pars <- post_A$draws() %>% spread_draws(x0_global, k_global) %>% summarise(x0 = mean(x0_global), k = mean(k_global))

auc <- get_auc(as.numeric(y  > 0.5), y_hat_est$y_hat_est) 
coverage <- get_coverage(x, log_pars$k, log_pars$x0)
improvement <- get_improvement(x_med, log_pars$k, log_pars$x0)
improvement_4 <- get_improvement(4, log_pars$k, log_pars$x0)

df_metrics_full <- data.frame(
    auc = auc,
    coverage = coverage,
    improvement = improvement,
    improvement_4 = improvement_4
)


###### ###### ###### ###### ###### 
###### AGE-SPECIFIC DATA ###### 
###### ###### ###### ###### ###### 

age_groups_i <- c("≤5", "5-18", "19-59", "60-74", "75+")


y_hat_new_post_age <- post_A$draws() %>% spread_draws(y_hat_rel_age[i, j], y_prot_age[i, j], y_hat_new_age[i, j], x_new[i]) %>% 
    mutate(age_group = recode(j, `1` = "≤5", `2` = "5-18", `3` = "19-59", `4` = "60-74", `5` = "75+")) %>% 
    mutate(age_group = factor(age_group, levels = c("≤5", "5-18", "19-59", "60-74", "75+"))) 

x_med_age <- cop_age_s %>% group_by(age_group) %>% summarise(titre_val = median(titre_val), n = n())

cop_age_s_sum <- cop_age_s %>% mutate(titre_group = 
    cut(titre_val, breaks = seq(3, 5, 0.2), labels  = seq(3.05, 5, 0.2))) %>%
    group_by(age_group, titre_group) %>%
    mutate(titre_group = as.numeric(as.character(titre_group))) %>%
    summarise(prop = mean(prop), n = n()) %>% arrange(age_group, titre_group)

df_cop_age <- y_hat_new_post_age %>% rename(y_prot = y_prot_age)
df_thresholds_age <- map_df(1:N_cov, function(i_i) {
    df_thresholds_1 <- get_thresholds(df_cop_age %>% filter(j == i_i)) %>% mutate(age_group = age_groups_i[i_i])
    df_thresholds_1
})



### METRICS
log_pars_age <- post_A$draws() %>% spread_draws(k_temp[k], x0_temp[k]) %>% summarise(x0 = mean(x0_temp), k = mean(k_temp))


age_groups_i <- c("≤5", "5-18", "19-59", "60-74", "75+")
df_metrics_age_auc <- map_df (1:N_cov, 
    function(i) {
        selector <- x_cov == i
        auc <- get_auc(as.numeric(y[selector]  > 0.5), y_hat_est$y_hat_est[selector]) 
        data.frame(
            age_group = age_groups_i[i],
            auc = auc
        )
    }
)

df_metrics_age_covarage <- map_df (1:N_cov, 
    function(i) {
        selector <- x_cov == i
        data.frame(
            age_group = age_groups_i[i],
            coverage = get_coverage(x[selector], log_pars_age[i, ]$k, log_pars_age[i, ]$x0)
        )
    }
)

df_metrics_age_AMG <- map_df (1:N_cov, 
    function(i) {
        selector <- x_cov == i
        data.frame(
            age_group = age_groups_i[i],
            AMG = get_improvement(mean(x[selector]), log_pars_age[i, ]$k, log_pars_age[i, ]$x0)
        )
    }
)


df_metrics_age_AMG_4 <- map_df (1:N_cov, 
    function(i) {
        selector <- x_cov == i
        data.frame(
            age_group = age_groups_i[i],
            AMG_4 = get_improvement(4, log_pars_age[i, ]$k, log_pars_age[i, ]$x0)
        )
    }
)


df_metrics_age_full <- left_join(df_metrics_age_auc, df_metrics_age_covarage, by = "age_group") %>% 
    left_join(df_metrics_age_AMG, by = "age_group") %>%
    left_join(df_metrics_age_AMG_4, by = "age_group") %>% 
    mutate(age_group = factor(age_group, levels = c("≤5", "5-18", "19-59", "60-74", "75+")))


require(gt)
# Create gt table with color coding

# Create gt table with color coding
create_gt_table <- function(df, columns_to_color, title = "C. Age-Specific Metrics for COP") {
  
  # Rename columns for display
  df <- df %>%
    rename(
      `Age Group` = age_group,
      AUC = auc,
      Coverage = coverage,
      `AMG (median)` = AMG,
      `AMG (log_10 4)` = AMG_4
    )
  
  # Create the gt table
  gt_table <- df %>%
    gt() %>%
    tab_header(
      title = title,
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_title("title")
    ) %>%
    # Style the table
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    # Add color coding for each specified column
    data_color(
      columns = all_of(columns_to_color),
      colors = scales::col_numeric(
        palette = c("#C8C8C8", "#DC1414"),
        domain = NULL
      ),
      alpha = 0.8
    ) %>%
    # Format numbers
    fmt_number(
      columns = all_of(columns_to_color),
      decimals = 4
    ) %>%
    # Add borders
    tab_style(
      style = cell_borders(sides = "all", color = "white", weight = px(1)),
      locations = cells_body()
    ) %>%
    # Style the header
    tab_style(
      style = cell_borders(sides = "all", color = "gray", weight = px(2)),
      locations = cells_column_labels()
    ) %>%
    # Adjust column widths
    cols_width(
      `Age Group` ~ px(100),
      everything() ~ px(140)
    )
  
  return(gt_table)
}


# Create the gt table
gt_table <- create_gt_table(df_metrics_age_full, c("AUC", "Coverage", "AMG (median)", "AMG (log_10 4)"))

# Save as various formats using gtsave
gtsave(gt_table, 
       filename = "df_metrics_age_full_table.png",
       path = here::here("outputs", "manu", "main"))

gtsave(gt_table, 
       filename = "df_metrics_age_full_table.pdf",
       path = here::here("outputs", "manu", "main"))

# Print the table to console for immediate viewing
cat("Color-coded table saved as:\n")
cat("- PNG image: outputs/manu/main/df_metrics_age_full_table.png\n")
cat("- PDF image: outputs/manu/main/df_metrics_age_full_table.pdf\n")
cat("- HTML table: outputs/manu/main/df_metrics_age_full_table.html\n")
cat("- LaTeX table: outputs/manu/main/df_metrics_age_full_table.tex\n")
cat("- CSV: outputs/manu/main/df_metrics_age_full_colored.csv\n")
cat("- Basic CSV: outputs/manu/main/df_metrics_age_full.csv\n\n")

# Display the gt table in console
print(gt_table)

# Display the raw data
df_metrics_age_full

# Extract COP information for results text
cat("=== COP VALUES FOR RESULTS ===\n")

# Overall COP metrics
cat("Overall COP metrics:\n")
cat("AUC:", round(df_metrics_full$auc, 3), "\n")
cat("Coverage:", round(df_metrics_full$coverage, 3), "\n")
cat("AMG (median):", round(df_metrics_full$improvement * 100, 0), "%\n")
cat("AMG (log10 4):", round(df_metrics_full$improvement_4 * 100, 0), "%\n")

# Protection thresholds
cat("\nProtection thresholds:\n")
print(df_thresholds)

# Extract 75% and 90% protection thresholds
threshold_75 <- df_thresholds %>% filter(y_prot == "75%")
threshold_90 <- df_thresholds %>% filter(y_prot == "90%")

threshold_75_median <- round(threshold_75$mean, 1)
threshold_75_lb <- round(threshold_75$lb, 1)
threshold_75_ub <- round(threshold_75$ub, 1)

threshold_90_median <- round(threshold_90$mean, 1)
threshold_90_lb <- round(threshold_90$lb, 1)
threshold_90_ub <- round(threshold_90$ub, 1)

cat("75% protection threshold:", threshold_75_median, "(95% CrI:", threshold_75_lb, "-", threshold_75_ub, ")\n")
cat("90% protection threshold:", threshold_90_median, "(95% CrI:", threshold_90_lb, "-", threshold_90_ub, ")\n")

# Age-specific COP metrics
cat("\nAge-specific COP metrics:\n")
print(df_metrics_age_full)

# Extract age-specific values
age_leq5_metrics <- df_metrics_age_full %>% filter(age_group == "≤5")
age_5_18_metrics <- df_metrics_age_full %>% filter(age_group == "5-18")
age_19_59_metrics <- df_metrics_age_full %>% filter(age_group == "19-59")
age_60_74_metrics <- df_metrics_age_full %>% filter(age_group == "60-74")
age_75plus_metrics <- df_metrics_age_full %>% filter(age_group == "75+")

# Create summary for markdown
cop_info <- list(
    # Overall COP
    cop_auc = round(df_metrics_full$auc, 3),
    cop_coverage = round(df_metrics_full$coverage, 3),
    amg_median = round(df_metrics_full$improvement * 100, 0),
    amg_log4 = round(df_metrics_full$improvement_4 * 100, 0),
    
    # Thresholds
    threshold_75_median = threshold_75_median,
    threshold_75_lb = threshold_75_lb,
    threshold_75_ub = threshold_75_ub,
    threshold_90_median = threshold_90_median,
    threshold_90_lb = threshold_90_lb,
    threshold_90_ub = threshold_90_ub,
    
    # Age-specific AUC
    age_leq5_auc = round(age_leq5_metrics$auc, 2),
    age_5_18_auc = round(age_5_18_metrics$auc, 2),
    age_19_59_auc = round(age_19_59_metrics$auc, 2),
    age_60_74_auc = round(age_60_74_metrics$auc, 2),
    age_75plus_auc = round(age_75plus_metrics$auc, 2),
    
    # Age-specific Coverage
    age_leq5_coverage = round(age_leq5_metrics$coverage, 2),
    age_5_18_coverage = round(age_5_18_metrics$coverage, 2),
    age_19_59_coverage = round(age_19_59_metrics$coverage, 2),
    age_60_74_coverage = round(age_60_74_metrics$coverage, 2),
    age_75plus_coverage = round(age_75plus_metrics$coverage, 2)
)

# Save the information
saveRDS(cop_info, file = here::here("markdown", "cop_info.RDS"))



p1 <- y_hat_new_post_age %>% 
    ggplot() +
    geom_vline(data = x_med_age, aes(xintercept = titre_val, color = age_group), linetype = "dashed", size = 1) +
    stat_lineribbon(aes(x = x_new, y = y_hat_new_age, fill = age_group), .width = 0.95, alpha = 0.2) + 
    geom_point(data = cop_age_s_sum, aes(x = titre_group, y = prop, color = age_group, size = n), alpha = 1) +
    facet_wrap(vars(biomarker)) + theme_minimal() +
    guides(color = "none", fill = "none", size = "none") +
    facet_wrap(vars(age_group)) +
#scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "PreF titre (log10)", y = "Posterior probability of infection", color = "Biomarker", fill = "Biomarker") + 
    ggtitle("A. Age-specific fitted curves for infection risk") + theme_ft()


p2 <- y_hat_new_post_age %>% 
    ggplot() +
    geom_vline(data = x_med_age, aes(xintercept = titre_val, color = age_group), linetype = "dashed", size = 1) +
    stat_lineribbon(aes(x = x_new, y = y_prot_age, fill = age_group), .width = 0.95, alpha = 0.2) + 
     facet_wrap(vars(age_group)) +theme_minimal() +
        guides(color = "none", fill = "none") +
    labs(x = "PreF titre (log10)", y = "Probability of protection given exposure", color = "Biomarker") + 
    ggtitle("B. Age-specific fitted curves for absolute COP") + theme_ft()

p1 / p2
ggsave(here::here("outputs", "manu", "main", "age_cop_2.pdf"), width = 10, height = 10)
ggsave(here::here("outputs", "manu", "main", "age_cop_2.png"), width = 10, height = 10)


p3 <- y_hat_new_post_age %>% 
    ggplot() +
  #  geom_vline(data = x_med_age, aes(xintercept = titre_val, color = age_group), linetype = "dashed", size = 1) +
    stat_summary(aes(x = x_new, y = y_prot_age, color = age_group), geom = "line",alpha = 1, size = 2) + 
    labs(x = "PreF titre (log10)", y = "Probability of protection given exposure", color = "Biomarker") + 
    ggtitle("B. Fitted curves for absolute COP") + theme_ft()

p4 <- df_thresholds_age %>% mutate(age_group = factor(age_group, levels = c("≤5", "5-18", "19-59", "60-74", "75+"))) %>%
    ggplot() +
    geom_point(aes(y = age_group, x = mean, group = y_prot, color = y_prot), position = position_dodge(width = 0.5), size = 5) +
    geom_linerange(aes(y = age_group, xmin = lb, xmax = ub, group = y_prot), position = position_dodge(width = 0.5)) +
    theme_minimal() +
    #    facet_wrap(vars(y_prot)) +
    labs(y = "Protection threshold", x = "PreF titre (log10)") + 
    ggtitle("C. Thresholds for protection") + theme_ft()

p3 / p4
ggsave(here::here("outputs", "manu", "main", "age_cop.pdf"), width = 10, height = 10)
ggsave(here::here("outputs", "manu", "main", "age_cop.png"), width = 10, height = 10)


post_A$draws() %>% spread_draws( x0_deviation[j]) 
post_A$draws() %>% spread_draws( x0_sigma) 

post_params <- post_A$draws() %>% spread_draws(L_temp[j], k_temp[j], x0_temp[j]) 

p1A <- post_params %>% ggplot() + 
  stat_pointinterval(aes(x = as.character(j), y = L_temp)) 

p2A <- post_params %>% ggplot() + 
  stat_pointinterval(aes(x = as.character(j), y = k_temp)) 

p3A <- post_params %>% ggplot() + 
  stat_pointinterval(aes(x = as.character(j), y = x0_temp)) 

require(patchwork)
p1A / p2A / p3A


p3 <- y_hat_new_post_age %>% 
    ggplot() +
    stat_summary(aes(x = x_new, y = y_hat_rel_age, color = age_group), geom = "line",alpha = 1, size = 2) + 
   # geom_point(data = cop_age_s, aes(x = titre_val, y = prop), alpha = 0.6, size = 2) +
   # facet_wrap(vars(j)) + theme_minimal() +

   # facet_wrap(vars(j)) +
#scale_x_continuous(breaks = seq(0, 1, length.out = 10), labels = labels_x_spike) + 
    labs(x = "PreF titre (log10)", y = "Risk reduction", color = "Biomarker") + 
    ggtitle("C. Fitted curves for relative COP") + theme_ft() 

