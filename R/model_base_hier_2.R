library(devtools)
library(serojump)
library(tidyverse)


obsLogLikelihood <- function(titre_val, titre_est, pars) {

    ll <- dnorm(titre_val, titre_est, pars[1], log = TRUE)
}

noInfSerumKinetics <- function(titre_est, timeSince, pars) {
    titre_est_new <- titre_est - (pars[1] * titre_est) * (timeSince)
    titre_est_new <- max(titre_est_new, 0)
    titre_est_new
}

infTuenisPower2016 <- function(titre_est, timeSince, pars) {

    y1 <- pars[1]
    t1 <- pars[2]
    r <- pars[3]
   # alpha <- pars[4]

    v <- 0.001
    mu <- 1 / t1 * y1

    if (timeSince < t1) {
        titre_est_boost <- exp(mu * timeSince)
    } else {
        titre_est_boost <- exp(y1) * (1 + (r - 1) * exp(y1)^{r - 1} * v * (timeSince - t1)) ^ {-1 / (r - 1)}
    }

    titre_est_log <- titre_est + log(titre_est_boost) #* max(0, 1 - titre_est * alpha)
    titre_est_log 
}
#write.csv(df_sero_model, file = here::here("data_clean", "fudan", "df_soro_model.csv"))


data_titre <- read.csv( file = here::here("data_clean", "fudan", "df_soro_model.csv")) %>%  mutate(age_group = factor(age_group, levels = c("<=5", "5-18", "19-59", "60-74", "75+"))) %>% select(!X) 
exp_prior <- read.csv(file = here::here("data_clean", "fudan", "exp_prior.csv")) %>% select(!X)

hierdata <- data_titre %>% select(id, age_group) %>% unique %>% mutate(age_group = as.numeric(age_group)) %>% pull(age_group) 
M <- hierdata %>% unique %>% length

# Define the biomarkers and exposure types in the model
biomarkers <- c("PreF")
exposureTypes <- c("none", "inf")
exposureFitted <- c("inf")

# Define the observational model
observationalModel <- list(
    names = c("PreF"),
    model = makeModel(
        addObservationalModel("PreF", c("sigma"), obsLogLikelihood)
        ),
    prior = bind_rows(
        addPrior("sigma", 0.0001, 4, "unif", 0.0001, 4)
        ) # observational model,
)


# Define the antibody kinetics model
abkineticsModel <- list(
    model = makeModel(
            addAbkineticsModel("none", "PreF", "none", c("wane"), noInfSerumKinetics),
            addAbkineticsModelHier("inf", "PreF", "inf", c("y1", "t1", "r"), c("y1",  "r"), hierdata, infTuenisPower2016)
        ),
    prior = bind_rows(
        addPriorHier("y1", -10, 10, "norm",  0, 1.82, "lnorm", -1, 1, M, c(0.5, 6) ), # ab kinetics
        addPrior("t1", 7, 24, "norm", 14, 3), # ab kinetics
        addPriorHier("r", -10, 10, "norm", 0, 1.82, "lnorm", -1, 1, M , c(1, 5)), # ab kinetics 
      #  add_par_df("s", 0, 1, "unif", 0, 1), # ab kinetics 
        addPrior("wane", 0, 1, "unif", 0, 1) # ab kinetics
    )
)

inf_prior_1 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K 
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) 
    logPriorExpInf
}

inf_prior_2 <- function(N, E, I, K) {
    N_adj <- N - K
    E_adj <- E - K
    logPriorExpInf <- lfactorial(E_adj) + lfactorial(N_adj - E_adj) - lfactorial(N_adj ) + dbinom(E_adj, N_adj, 0.05, log = TRUE)
    logPriorExpInf
}


seroModel_p1 <- createSeroJumpModel(
    data_sero = data_titre, 
    data_known = NULL, 
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    exposurePriorTime = exp_prior,
    exposurePriorTimeType = "empirical"
)


seroModel_p2 <- createSeroJumpModel(
    data_sero = data_titre, 
    data_known = NULL, 
    biomarkers = biomarkers,
    exposureTypes = exposureTypes,
    exposureFitted = exposureFitted,
    observationalModel = observationalModel,
    abkineticsModel = abkineticsModel,
    exposurePriorTime = exp_prior,
    exposurePriorTimeType = "empirical",
    exposurePriorPop = inf_prior_2
)


outputs_modelA <- list(seroModel_p1, seroModel_p2)


save_info_list <- list(
    list(file_name = "fudan_e3", model_name = "base_hier_2"),
    list(file_name = "fudan_e3", model_name = "inform_hier_2"),
    list(file_name = "fudan_e3_hpc", model_name = "base_hier_2"),
    list(file_name = "fudan_e3_hpc", model_name = "inform_hier_2")
)

model_run_info <- 
    list(model = outputs_modelA, save = save_info_list)
saveRDS(model_run_info, file = here::here("models", "model_base_hier_2.RDS"))


p1 <- plotSero(seroModel_p1, 1000)
ggsave(here::here("outputs", "manu", "pre", "serodata.pdf"))
p2 <- plotPriorPredictive(seroModel_p1)
ggsave(here::here("outputs", "manu", "pre", "abkin_priorpred.pdf"))
p3 <- plotPriorInfection(seroModel_p1)
ggsave(here::here("outputs", "manu", "pre", "infrate.pdf"))


