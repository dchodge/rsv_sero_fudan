# =============================================================================
# RSV Serology Analysis: Data Cleaning and Preparation
# =============================================================================
# This script processes raw serological data from the Fudan study, creating
# clean datasets for analysis of two epidemic waves (E1: early 2023, E2: late 2023-early 2024)
# and preparing age-stratified datasets for hierarchical modeling.

# =============================================================================
# Load and Initial Data Processing
# =============================================================================

# Load raw serological data
# Contains longitudinal RSV PreF antibody titers from participants
raw_sero <- read.csv(here::here("data", "fudan", "data_titer_original_20241009.csv"))

# Transform data from wide to long format and process dates
# - Convert titer columns (Titer_2303 to Titer_2406) to long format
# - Apply log10 transformation to antibody titers (standard practice)
# - Parse date information from column names (format: Titer_MMYY)
# - Create proper date objects for temporal analysis
df_sero <- raw_sero %>% 
    pivot_longer(Titer_2303:Titer_2406, names_to = "times", values_to = "values") %>% 
    mutate(values = log10(values)) %>%  # Log10 transform antibody titers
    separate(times, c("blank", "day")) %>%  # Extract date from column name
    mutate(day = paste0("01-", substr(day, 3, 4), "-", substr(day, 1, 2))) %>%  # Format as DD-MM-YY
    select(!blank) %>%  # Remove blank column
    mutate(time = dmy(day)) %>%  # Convert to proper date format
    group_by(No) %>%  # Group by participant ID
    mutate(r = row_number())  # Add row number for each participant

# Split data into two epidemic waves based on temporal patterns
# E1: First 4 measurements (early 2023 epidemic wave)
# E2: Measurements 4 onwards (late 2023-early 2024 epidemic wave)
# Note: Measurement 4 is included in both waves for continuity
df_sero_case1 <- df_sero %>% filter(r <= 4)  # First epidemic wave
df_sero_case2 <- df_sero %>% filter(r >= 4)  # Second epidemic wave

# =============================================================================
# E1: First Epidemic Wave Data Preparation (Early 2023)
# =============================================================================

# Convert dates to days since start of study (2023-01-01)
# This creates a continuous time scale for modeling
df_sero_model_e1_temp <- df_sero_case1 %>% 
    mutate(day = as.numeric(time - ymd("2023-01-01"))) %>% 
    group_by(No) %>% 
    filter(!is.na(values)) %>%  # Remove missing values
    mutate(r = row_number())    # Renumber after filtering

# Filter out individuals with only one measurement
# Serological jump models require multiple measurements per person to estimate
# infection timing and antibody kinetics
solo_entries <- df_sero_model_e1_temp %>% 
    filter(max(r) == 1) %>% 
    pull(No)

# Create final dataset for E1 analysis
df_sero_model_e1 <- df_sero_model_e1_temp %>% 
    filter(!No %in% solo_entries) %>%  # Remove single-measurement individuals
    mutate(id = factor(No, levels = unique(No))) %>%  # Create factor ID
    mutate(id = as.numeric(id)) %>%  # Convert to numeric for Stan
    select(id, day, value = values) %>% 
    mutate(biomarker = "PreF") %>%  # Specify biomarker type
    rename(PreF = value, time = day) %>% 
    ungroup %>%
    select(!c(No))  # Remove original ID column

# CRITICAL: Ensure there is more than one entry per person for valid modeling
# This is essential for the serological jump model to estimate infection timing

# =============================================================================
# Create Epidemiological Prior for E1
# =============================================================================

# Define incidence data for the first epidemic wave
# These values represent monthly incidence rates (per 1000 population)
# Based on epidemiological surveillance data
inci <- c(0, 0, 0, 0, 1.2, 2.65, 3.33, 3.01, 5.26, 2.06, 1.09, 0, 0, 2.5, 9.17, 4, 1.67, 0, 0.41, 1.3)
time <- (0:(length(inci) - 1) * 30) + 1  # Convert to days (30-day months)

# Create exposure prior distribution for E1
# This provides the model with information about when infections were likely to occur
library(zoo)

exp_prior_e1 <- data.frame(
    time = time,
    inci = inci
) %>% 
    complete(time = 1:600) %>%  # Fill in missing days
    mutate(inci = na.approx(inci, x = time, na.rm = FALSE)) %>%  # Interpolate missing values
    filter(time < 305) %>%  # Focus on E1 period (first ~10 months)
    mutate(prob = inci / sum(inci))  # Convert to probability distribution

exp_prior_e1 <- exp_prior_e1 %>% select(day = time, prob)

# Save processed datasets
write.csv(df_sero_model_e1, file = here::here("data_clean", "fudan", "df_soro_model_e1.csv"))
write.csv(exp_prior_e1, file = here::here("data_clean", "fudan", "exp_prior_e1.csv"))


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## E2: Get serological dynamics for the first epidemic wave
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 



# Convert the date to day
df_sero_model_e2_temp <- df_sero_case2 %>% mutate(day = as.numeric(time - ymd("2023-11-01") )) %>% 
    group_by(No) %>% filter(!is.na(values)) %>% mutate(r = row_number())

solo_entries <- df_sero_model_e2_temp %>% filter(max(r) == 1) %>% pull(No)
df_sero_model_e2 <- df_sero_model_e2_temp %>% filter(!No %in% solo_entries) %>% 
    mutate(id = factor(No, levels = unique(No))) %>% mutate(id = as.numeric(id)) %>% 
    select(id, day, value = values) %>% mutate(biomarker = "PreF") %>% rename(PreF = value, time = day) %>% 
    ungroup %>%
    select(!c(No))

exp_prior_e2 <- data.frame(
    time = time,
    inci = inci
) %>% complete(time = 1:600) %>%
 mutate(
    inci = na.approx(inci, x = time, na.rm = FALSE)
  ) %>% filter(time > 403) %>% mutate(prob = inci / sum(inci))

exp_prior_e2 <- exp_prior_e2 %>% select(day = time, prob)


write.csv(df_sero_model_e2, file = here::here("data_clean", "fudan", "df_soro_model_e2.csv"))
write.csv(exp_prior_e2, file = here::here("data_clean", "fudan", "exp_prior_e2.csv"))



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## E3: Get serological dynamics for both waves
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# Convert the date to day
df_sero_model_temp <- df_sero %>% mutate(day = as.numeric(time - ymd("2023-01-01") )) %>% 
    group_by(No) %>% filter(!is.na(values)) %>% mutate(r = row_number())

df_linker_age <- df_sero %>%  mutate(date_of_birth = as.Date(date_of_birth, format = "%Y/%m/%d")) %>%
  mutate(age = as.integer(difftime(time, date_of_birth, units = "days") / 365.25)) %>%
  mutate(age_group = case_when(
    age <= 5 ~ "≤5",
    age >= 6 & age <= 18 ~ "5-18",
    age >= 19 & age <= 59 ~ "19-59",
    age >= 60 & age <= 74 ~ "60-74",
    age >= 75 ~ "75+"
  )) %>% mutate(age_group = factor(age_group, levels = c("≤5", "5-18", "19-59", "60-74", "75+"))) %>% select(No, age_group) %>% unique %>% 
  group_by(No) %>% mutate(r = row_number()) %>% filter( r == 1) %>% select(No, age_group)

solo_entries <- df_sero_model_temp %>% filter(max(r) == 1) %>% pull(No)

df_sero_model <- df_sero_model_temp %>% filter(!No %in% solo_entries) %>% 
    mutate(id = factor(No, levels = unique(No))) %>% mutate(id = as.numeric(id)) %>% 
    select(id, day, value = values) %>% mutate(biomarker = "PreF") %>% rename(PreF = value, time = day) %>% 
    ungroup %>% left_join(df_linker_age) %>%
    select(!c(No)) 

## MUST ENSURE THERE IS MORE THAN ONE ENTRY PER PERSON!

#df_sero_model_e1[complete.cases(df_sero_model_e1)]
# full data
inci <- c(0, 0, 0, 0, 1.2, 2.65, 3.33, 3.01, 5.26, 2.06, 1.09, 0, 0, 2.5, 9.17, 4, 1.67, 0, 0.41, 1.3)
time <- (0:(length(inci) - 1) * 30) + 1

# extract data from epidemic 1
library(zoo)

exp_prior <- data.frame(
    time = time,
    inci = inci
) %>% complete(time = 1:600) %>%
 mutate(
    inci = na.approx(inci, x = time, na.rm = FALSE)
  ) %>% filter(!is.na(inci)) %>% mutate(prob = inci / sum(inci)) 

exp_prior <- exp_prior %>% select(time, prob)


write.csv(df_sero_model, file = here::here("data_clean", "fudan", "df_soro_model.csv"))
write.csv(exp_prior, file = here::here("data_clean", "fudan", "exp_prior.csv"))


boost_id_case1 <- df_sero_model_e1 %>% group_by(id) %>% mutate(change = PreF - lag(PreF)) %>% filter(change > 0.3) %>% pull(id)
p1 <- df_sero_model_e1 %>% filter(id %in% boost_id_case1) %>% 
    ggplot() + 
        geom_line(aes(x = time, y = PreF, group = id))

boost_id_case2 <- df_sero_case2 %>% mutate(change = values - lag(values)) %>% filter(change > 0.3) %>% pull(No)
p2 <- df_sero_case2 %>% filter(No %in% boost_id_case2) %>% 
    ggplot() + 
        geom_line(aes(x = time, y = values, group = No))

p1 / p2
