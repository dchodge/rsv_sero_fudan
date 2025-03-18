
# get serological data and clean
raw_sero <- read.csv(here::here("data", "fudan", "data_titer_original_20241009.csv"))

df_sero <- raw_sero %>% pivot_longer(Titer_2303:Titer_2406, names_to = "times", values_to = "values") %>% 
    mutate(values = log10(values)) %>% separate(times, c("blank", "day")) %>%
    mutate(day = paste0("01-", substr(day, 3, 4), "-", substr(day, 1, 2))) %>% select(!blank) %>%
    mutate(time = dmy(day)) %>% group_by(No) %>% mutate(r = row_number()) 
    
df_sero_case1 <- df_sero %>% filter(r <= 4)
df_sero_case2 <- df_sero %>% filter(r >= 4)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## E1: Get serological dynamics for the first epidemic wave
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# Convert the date to day
df_sero_model_e1_temp <- df_sero_case1 %>% mutate(day = as.numeric(time - ymd("2023-01-01") )) %>% 
    group_by(No) %>% filter(!is.na(values)) %>% mutate(r = row_number())

solo_entries <- df_sero_model_e1_temp %>% filter(max(r) == 1) %>% pull(No)
df_sero_model_e1 <- df_sero_model_e1_temp %>% filter(!No %in% solo_entries) %>% 
    mutate(id = factor(No, levels = unique(No))) %>% mutate(id = as.numeric(id)) %>% 
    select(id, day, value = values) %>% mutate(biomarker = "PreF") %>% rename(PreF = value, time = day) %>% 
    ungroup %>%
    select(!c(No))

## MUST ENSURE THERE IS MORE THAN ONE ENTRY PER PERSON!

#df_sero_model_e1[complete.cases(df_sero_model_e1)]
# full data
inci <- c(0, 0, 0, 0, 1.2, 2.65, 3.33, 3.01, 5.26, 2.06, 1.09, 0, 0, 2.5, 9.17, 4, 1.67, 0, 0.41, 1.3)
time <- (0:(length(inci) - 1) * 30) + 1

# extract data from epidemic 1
library(zoo)

exp_prior_e1 <- data.frame(
    time = time,
    inci = inci
) %>% complete(time = 1:600) %>%
 mutate(
    inci = na.approx(inci, x = time, na.rm = FALSE)
  ) %>% filter(time < 305) %>% mutate(prob = inci / sum(inci))

exp_prior_e1 <- exp_prior_e1 %>% select(day = time, prob)


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
    age <= 5 ~ "<=5",
    age >= 6 & age <= 18 ~ "5-18",
    age >= 19 & age <= 59 ~ "19-59",
    age >= 60 & age <= 74 ~ "60-74",
    age >= 75 ~ "75+"
  )) %>% mutate(age_group = factor(age_group, levels = c("<=5", "5-18", "19-59", "60-74", "75+"))) %>% select(No, age_group) %>% unique %>% 
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
