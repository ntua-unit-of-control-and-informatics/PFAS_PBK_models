# Plot comparison of observed vs predicted values for Graphene Oxide model
# Based on Liu et al. 2012 data and MCMC parameter estimation

library(deSolve)
library(nimble)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Assuming the model functions are already defined (create.params, ode.func, etc.)
# and MCMC parameter estimation has been completed

# --------------------------------
# Load data and setup parameters
# --------------------------------
Liu_1_small_tissues <- read.xlsx("Data/Liu_2012_GO_male_small_1_tissues.xlsx")
colnames(Liu_1_small_tissues)[c(2,3)] <- c("time", "mass")

# Setup parameters
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg * BW * 1e03 # ug PFOA
np_size <- 3/2 # nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list(
  'BW' = BW,
  "admin.dose" = admin.dose,
  "np_size" = np_size,
  "admin.time" = admin.time, 
  "admin.type" = admin.type,
  "sex" = sex
)

# Get initial parameters from user input
params <- create.params(user_input)

# --------------------------------
# Load MCMC results and update parameters
# --------------------------------
# If you have saved MCMC results, load them here
# Otherwise, use the posterior means from the summary_results object
# For demonstration purposes, I'll create a mock summary_results structure

# Example of how to create a mock summary_results structure
# This should be replaced with your actual MCMC results
mock_stats <- matrix(0, nrow = 13, ncol = 7)
rownames(mock_stats) <- c("coef_liver", "coef_spleen", "coef_kidney",
                          "coef_heart", "coef_lung", "coef_brain", 
                          "coef_rob", "coef_stomach", "coef_smallIn", 
                          "coef_largeIn", "CLurine", "CLfeces", "sigma")
colnames(mock_stats) <- c("Mean", "SD", "Naive SE", "Time-series SE", "2.5%", "50%", "97.5%")

# Set some example values - these should be replaced with your actual posterior means
mock_stats[, "Mean"] <- c(0.3, 0.25, 0.4, 0.35, 0.45, 0.2, 0.3, 0.25, 0.3, 0.35, 0.1, 0.15, 0.2)

summary_results <- list(statistics = mock_stats)

# Update parameters with MCMC results
updated_params <- update_params_with_mcmc(params, summary_results)

# --------------------------------
# Run model with updated parameters
# --------------------------------
time_points <- seq(0, 10, 0.05) # time points for ODE solution
updated_solution <- run_ode_model(updated_params, time_points)

# Convert solution to data frame
solution_df <- as.data.frame(updated_solution)

# --------------------------------
# Prepare experimental data
# --------------------------------
# Convert percentage to actual mass values
exp_data <- Liu_1_small_tissues
exp_data$mass_ug <- exp_data$mass * admin.dose / 100

# --------------------------------
# Generate plots comparing observed vs predicted for each tissue
# --------------------------------
# List of tissues to plot
tissues <- c("heart", "liver", "spleen", "stomach", 
             "kidneys", "lungs", "brain", 
             "small_intestine", "large_intestine")

# Create one plot for all tissues
plot_data <- data.frame()

# For each tissue, extract observed and predicted values
for (tissue in tissues) {
  # Get column name from the ODE solution
  col_name <- paste0("M", tissue)
  
  # Extract predicted values at experimental time points
  pred_at_exp_times <- solution_df[solution_df$time %in% exp_data$time, c("time", col_name)]
  names(pred_at_exp_times)[2] <- "predicted"
  
  # Filter experimental data for this tissue
  tissue_exp_data <- exp_data[exp_data$Tissue == tissue, c("time", "mass_ug")]
  names(tissue_exp_data)[2] <- "observed"
  
  # Merge observed and predicted
  merged_data <- merge(tissue_exp_data, pred_at_exp_times, by = "time")
  merged_data$Tissue <- tissue
  
  # Append to plot data
  plot_data <- rbind(plot_data, merged_data)
}

# Create plots
# 1. Individual tissue plots
individual_plots <- list()

for (tissue in tissues) {
  tissue_data <- plot_data[plot_data$Tissue == tissue, ]
  
  p <- ggplot(tissue_data, aes(x = time)) +
    geom_point(aes(y = observed), color = "blue", size = 3) +
    geom_line(aes(y = predicted), color = "red", size = 1) +
    labs(title = paste("Tissue:"),
         x = "Time (min)",
         y = "Mass (μg)") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  individual_plots[[tissue]] <- p
}

# 2. Combined tissue plot
combined_plot <- ggplot(plot_data, aes(x = time, color = Tissue)) +
  geom_point(aes(y = observed), size = 3) +
  geom_line(aes(y = predicted), linetype = "dashed", size = 1) +
  labs(title = "Observed vs Predicted Graphene Oxide Distribution",
       x = "Time (min)",
       y = "Mass (μg)") +
  facet_wrap(~ Tissue, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none")

# 3. Correlation plot
correlation_plot <- ggplot(plot_data, aes(x = observed, y = predicted, color = Tissue)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Correlation between Observed and Predicted Values",
       x = "Observed (μg)",
       y = "Predicted (μg)") +
  theme_bw()

# 4. Residual plot
plot_data$residual <- plot_data$predicted - plot_data$observed
residual_plot <- ggplot(plot_data, aes(x = observed, y = residual, color = Tissue)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals Plot",
       x = "Observed (μg)",
       y = "Residual (Predicted - Observed)") +
  theme_bw()



# Save plots
pdf("GO_model_comparison_plots.pdf", width = 12, height = 10)
print(combined_plot)
print(correlation_plot)
print(residual_plot)

# Print individual tissue plots on separate pages
for (tissue in tissues) {
  print(individual_plots[[tissue]])
}
dev.off()

# Print additional model fit statistics
cat("Model fit statistics:\n")
cat("---------------------\n")

# Calculate R-squared for each tissue
tissue_r2 <- plot_data %>%
  group_by(Tissue) %>%
  summarize(
    RMSE = sqrt(mean((predicted - observed)^2)),
    MAE = mean(abs(predicted - observed)),
    R2 = cor(predicted, observed)^2
  )

print(tissue_r2)

# Calculate overall R-squared
overall_r2 <- cor(plot_data$predicted, plot_data$observed)^2
cat(sprintf("Overall R-squared: %.4f\n", overall_r2))

# Calculate normalized RMSE (as percentage of observed mean)
cat(sprintf("Overall RMSE: %.4f\n", sqrt(mean((plot_data$predicted - plot_data$observed)^2))))
cat(sprintf("Normalized RMSE: %.2f%%\n", 
            100 * sqrt(mean((plot_data$predicted - plot_data$observed)^2)) / mean(plot_data$observed)))