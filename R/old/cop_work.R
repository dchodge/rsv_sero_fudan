
#install.packages("devtools")
library(devtools)
#install(here::here("..", "serojump"))
 # install working version of this package 
library(serojump)
library(lubridate)
library(patchwork)
library(tidybayes)
library(ggdist)

library(tidyverse)

model_summary_h1n1 <- readRDS(here::here("outputs", "fits", "nih_A", "h3n2_inform", "model_summary.RDS"))
model_summary <- model_summary_h1n1



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

cop_exp_sum_plot_all$biomarker %>% unique
# A/Sydney/5/2021
# A/Sydney/5/2021e
summary_cop <- left_join(
    cop_exp_sum_plot_all %>% filter(biomarker == "A/Darwin/06/2021") %>% select(id, prop, cell = titre_val),
    cop_exp_sum_plot_all %>% filter(biomarker == "A/Darwin/09/2021e")   %>% select(id, egg = titre_val) 
)

summary_cop %>% 
    ggplot() + 
        geom_point(aes(x = cell, y = egg, size = prop), alpha = 0.5) + theme_bw()
ggsave(here::here("outputs", "manu", "cop_raw.pdf"))

library(fields)

library(caret)

fit <- Tps(summary_cop[c("cell", "egg")], summary_cop$prop)

x_pred <- seq(1, 8, length.out = 100)
y_pred <- seq(1, 9, length.out = 100)
grid <- expand.grid(x_pred, y_pred)
z_pred <- predict(fit, grid, lambda = 0.5)

interp_data <- data.frame(
  x = grid[,1],
  y = grid[,2],
  z = z_pred
)

ggplot(interp_data) +
  geom_tile(aes(x = x, y = y, fill = z)) +
  geom_contour( aes(x = x, y = y, z = z), color = "black") +
  scale_fill_viridis_c() +
  geom_point(data = summary_cop, aes(x = cell, y = egg, size = prop), alpha = 0.5, color = 'white') +
  labs(title = "Cubic Spline Interpolation (Thin-Plate Spline)") +
  labs(x = "cell", y = "egg", z = "Probabiltiy of infection")

ggsave(here::here("outputs", "manu", "cop_surface_data.pdf"))

library(plotly)

# Reshape data for plotly
z_matrix <- matrix(z_pred, nrow = 100, ncol = 100, byrow = TRUE)

p1 <- plot_ly(x = x_pred, y = y_pred, z = z_matrix, type = "surface") %>%
  layout(title = "3D Thin-Plate Spline Surface", 
        xaxis = list(title = list(text = 'cell')),
      yaxis = list(title = 'egg')
)
p1

interp_data_add <- interp_data %>% 
  mutate(
    x_group = cut(x, c(-0.1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)), 
    y_group = cut(y, c(-0.1, 0, 1, 2, 3, 5, 6, 7, 8, 9)))


p1 <- interp_data_add %>%
  ggplot() + 
    stat_lineribbon(aes(x, z, fill = y_group), alpha = 0.5, .width = 0.95) + 
    ylim(0, 1) + 
    labs(x = "cell binding titre", 
      y = "Probability of infection", 
      fill = "egg titre")

p2 <- interp_data_add %>%
  ggplot() + 
    stat_lineribbon(aes(y, z, fill = x_group), alpha = 0.5, .width = 0.95) + 
    ylim(0, 1) + 
      labs(x = "egg binding titre", 
      y = "Probability of infection", 
      fill = "cell titre")

p1 / p2 & theme_bw()
ggsave(here::here("outputs", "manu", "cop_surface_fit_split.pdf"))
##




########################################
#### CROSS VALIDATION FOR LAMBDA ####
########################################

# Example data (substitute with your own)
set.seed(123)
x <- summary_cop$spike
y <- summary_cop$iga
z <- summary_cop$prop
data <- data.frame(x = x, y = y, z = z)

# Define custom cross-validation function
cross_validate_tps <- function(data, lambda_values, k = 5) {
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))  # Randomly assign folds
  mse_values <- matrix(NA, nrow = length(lambda_values), ncol = k)
  
  for (i in seq_along(lambda_values)) {
    lambda <- lambda_values[i]
    
    for (fold in 1:k) {
      # Split into training and test sets
      train_data <- data[folds != fold, ]
      test_data <- data[folds == fold, ]
      
      # Fit Tps model on the training data
      tps_fit <- Tps(cbind(train_data$x, train_data$y), train_data$z, lambda = lambda)
      
      # Predict on the test data
      predictions <- predict(tps_fit, cbind(test_data$x, test_data$y))
      
      # Calculate Mean Squared Error (MSE)
      mse_values[i, fold] <- mean((predictions - test_data$z)^2)
    }
  }
  
  # Return the average MSE for each lambda
  rowMeans(mse_values)
}

lambda_grid <- seq(0.1, 1000, length.out = 10)

# Perform cross-validation
mse_results <- cross_validate_tps(data, lambda_grid, k = 5)

# Find the optimal lambda
optimal_lambda <- lambda_grid[which.min(mse_results)]
optimal_lambda



### PCA

data(iris)

# Remove the categorical variable (Species) since PCA works only with numeric data
iris_data <- iris[, 1:4]

# Perform PCA
pca_result <- prcomp(iris_data, center = TRUE, scale. = TRUE)

# Display summary of PCA results
summary(pca_result)
plot(pca_result, type = "l", main = "Scree Plot")
biplot(pca_result, scale = 0)


correlation_matrix <- cor(summary_cop[, 3:4], summary_cop$prop)


pca_result <- prcomp(summary_cop[2:4], center = TRUE, scale. = TRUE)
summary(pca_result)
plot(pca_result, type = "l", main = "Scree Plot")
biplot(pca_result, scale = 0)

loadings <- pca_result$rotation 
squared_loadings <- loadings^2
total_contribution <- rowSums(squared_loadings)
total_contribution

barplot(abs(loadings[, 2]), names.arg = colnames(summary_cop[, 2:4]),
        main = "Contribution to PC1", ylab = "Absolute Loadings")

# Proportion of variance explained
pca_variance <- summary(pca_result)$importance[2, ]

# Cumulative proportion of variance explained
pca_cumulative_variance <- summary(pca_result)$importance[3, ]

pca_variance
pca_cumulative_variance
