run_rjmc <- function(model_info, save_info,  rj_settings, runmodel = TRUE) {

    if(runmodel) {
        model_summary <- runSeroJump(model_info, rj_settings, save_info = save_info)
    } else {
        model_summary <- readRDS(here::here("outputs", "fits", 
                save_info$file_name, save_info$model_name, paste0("model_summary.RDS")))
    }

    # Need to have save_info model_summary to run these
    plotMCMCDiagnosis(model_summary, save_info = save_info)
    plotPostFigs(model_summary, save_info = save_info)
}


theme_ft <- function() {
  theme_minimal(base_size = 14, base_family = "sans") %+replace%
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray85", size = 0.5),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "gray50"),
      axis.ticks = element_line(color = "gray50"),
      axis.text = element_text(size = 12, color = "gray30"),
      axis.title = element_text(size = 14, color = "gray30", face = "bold"),
      plot.title = element_text(size = 18, color = "gray10", face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(size = 14, color = "gray30", margin = margin(b = 15)),
      plot.caption = element_text(size = 10, color = "gray50", hjust = 0),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.spacing = unit(1, "lines")
    )
}