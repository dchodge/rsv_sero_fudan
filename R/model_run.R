library(furrr)
library(serojump)
library(tidyverse)
source("R/utils.R")

models <- readRDS(here::here("models", "model_base.RDS"))


plan(multisession, workers = 4)

rj_settings <- list(
        numberChainRuns = 4, 
        iterations = 500000,
        burninPosterior = 250000,
        thin = 1000,
        runParallel = TRUE,
        onDebug = FALSE
    )


save_info_base <- list(
    file_name = "fudan_e3",
    model_name = "base"
)


model_summary <- runSeroJump(models$model[[1]], rj_settings, save_info = models$save[[1]])


#model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base", "model_summary.RDS"))
 #Need to have save_info model_summary to run thes


plotMCMCDiagnosis(model_summary, save_info = models$save[[1]])
plotPostFigs(model_summary, save_info = models$save[[1]])

##############################################################################

models <- readRDS(here::here("models", "model_base_hier_2.RDS"))

#plan(multisession, workers = 4)

rj_settings <- list(
        numberChainRuns = 4, 
        iterations = 1000,
        burninPosterior = 500,
        thin = 1,
        runParallel = TRUE,
        onDebug = FALSE
    )

save_info_base_hier <- list(
    file_name = "fudan_e3",
    model_name = "base_hier_2"
)


model_summary <- runSeroJump(models$model[[1]], rj_settings, save_info = models$save[[1]])

#model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier", "model_summary.RDS"))

# Need to have save_info model_summary to run these
plotMCMCDiagnosis(model_summary, save_info = models$save[[1]])
plotPostFigs(model_summary, save_info = models$save[[1]])


models <- readRDS(here::here("models", "model_base_hier_3.RDS"))

#plan(multisession, workers = 4)

rj_settings <- list(
        numberChainRuns = 4, 
        iterations = 1000,
        burninPosterior = 500,
        thin = 1,
        runParallel = TRUE,
        onDebug = FALSE
    )

save_info_base_hier <- list(
    file_name = "fudan_e3",
    model_name = "base_hier_3"
)


model_summary <- runSeroJump(models$model[[1]], rj_settings, save_info = models$save[[1]])

#model_summary <- readRDS(here::here("outputs", "fits", "fudan_e3", "base_hier", "model_summary.RDS"))

# Need to have save_info model_summary to run these
plotMCMCDiagnosis(model_summary, save_info = models$save[[1]])
plotPostFigs(model_summary, save_info = models$save[[1]])
