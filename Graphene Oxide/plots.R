
#=====================================Prior_vs_predictive=================================
plot_prior_vs_predictive <- function(results) {
  # Extract the samples
  samples <- results$samples
  
  # Convert to a matrix for easier plotting
  sample_matrix <- as.matrix(samples)
  
  # Select parameters (excluding sigma)
  param_names <- c(
    "coef_liver", "coef_spleen", "coef_kidney", 
    "coef_heart", "coef_lung", "coef_brain", 
    "coef_rob", "coef_stomach", "coef_smallIn", 
    "coef_largeIn", "CLurine", "CLfeces"
  )
  
  # Set up the plotting area
  par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
  
  # Plot for each parameter
  for (param in param_names) {
    # Density of the parameter
    plot(density(sample_matrix[, param]), 
         main = paste("Distribution of", param),
         xlab = "Value", 
         ylab = "Density")
    
    # Add prior distribution (uniform between 0.01 and 1.0)
    curve(dunif(x, min = 0.01, max = 1.0), 
          col = "red", 
          add = TRUE, 
          lty = 2)
    
    # Add legend with reduced size and better position
    legend("topright", 
           legend = c("Posterior", "Prior"), 
           col = c("black", "red"), 
           lty = c(1, 2),
           cex = 0.7,        # Smaller text
           bty = "n",        # No box around legend
           inset = c(0.1, 0)) # Move legend inward
  }
}

# Function to generate predictive distribution with improved legend
generate_predictive_distribution <- function(results, original_params) {
  # Extract MCMC samples
  samples_matrix <- as.matrix(results$samples)
  
  # Set up time points and compartments
  time_points <- seq(0, 10, 0.05)
  compartments <- c("Mliver", "Mkidneys", "Mspleen", "Mlungs")
  
  # Prepare storage for predictions
  predictive_samples <- array(NA, 
                              dim = c(nrow(samples_matrix), 
                                      length(time_points), 
                                      length(compartments)))
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = nrow(samples_matrix), style = 3)
  
  # Generate predictions for each MCMC sample
  for (i in 1:nrow(samples_matrix)) {
    # Update parameters for this MCMC sample
    updated_params <- original_params
    
    # Update specific parameters
    param_names <- c(
      "coef_liver", "coef_spleen", "coef_kidney", 
      "coef_heart", "coef_lung", "coef_brain", 
      "coef_rob", "coef_stomach", "coef_smallIn", 
      "coef_largeIn", "CLurine", "CLfeces"
    )
    
    for (param in param_names) {
      updated_params[[param]] <- samples_matrix[i, param]
    }
    
    # Run ODE model with updated parameters
    solution <- run_ode_model(
      y = create.inits(updated_params), 
      time_points = time_points, 
      params = updated_params
    )
    
    # Store predictions
    for (j in 1:length(compartments)) {
      predictive_samples[i, , j] <- solution[, compartments[j]]
    }
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Plot predictive distributions
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  for (j in 1:length(compartments)) {
    # Prepare data for credible interval calculation
    comp_samples <- predictive_samples[, , j]
    
    # Plot median and credible intervals
    matplot(time_points, 
            t(comp_samples), 
            type = "l", 
            col = rgb(0, 0, 1, 0.1), 
            lty = 1,
            xlab = "Time (min)", 
            ylab = paste("Mass in", compartments[j], "(ug)"),
            main = paste("Predictive Distribution for", compartments[j])
    )
    
    # Add median line
    median_pred <- apply(comp_samples, 2, median)
    lines(time_points, median_pred, col = "red", lwd = 2)
    
    # Add 95% credible interval
    lower_ci <- apply(comp_samples, 2, quantile, 0.025)
    upper_ci <- apply(comp_samples, 2, quantile, 0.975)
    lines(time_points, lower_ci, col = "blue", lty = 2)
    lines(time_points, upper_ci, col = "blue", lty = 2)
    
    # Add legend with improved positioning and size
    legend(
      "topright", 
      legend = c("Median", "95% CI", "Predictions"),  # Shortened text 
      col = c("red", "blue", rgb(0, 0, 1, 0.1)), 
      lty = c(1, 2, 1),
      cex = 0.7,           # Smaller text
      bty = "n",           # No box
      inset = c(0.05, 0),  # Inset from the top right corner
      seg.len = 0.5        # Shorter line segments in legend
    )
  }
  
  # Return the predictive samples for further analysis if needed
  return(predictive_samples)
}

#==========================================Predictive_distribution plot======================================
generate_predictive_distribution <- function(results, original_params) {
  # Extract MCMC samples
  samples_matrix <- as.matrix(results$samples)
  
  # Set up time points and compartments
  time_points <- seq(0, 10, 0.05)
  compartments <- c("Mliver", "Mkidneys", "Mspleen", "Mlungs")
  
  # Prepare storage for predictions
  predictive_samples <- array(NA, 
                              dim = c(nrow(samples_matrix), 
                                      length(time_points), 
                                      length(compartments)))
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = nrow(samples_matrix), style = 3)
  
  # Generate predictions for each MCMC sample
  for (i in 1:nrow(samples_matrix)) {
    # Update parameters for this MCMC sample
    updated_params <- original_params
    
    # Update specific parameters
    param_names <- c(
      "coef_liver", "coef_spleen", "coef_kidney", 
      "coef_heart", "coef_lung", "coef_brain", 
      "coef_rob", "coef_stomach", "coef_smallIn", 
      "coef_largeIn", "CLurine", "CLfeces"
    )
    
    for (param in param_names) {
      updated_params[[param]] <- samples_matrix[i, param]
    }
    
    # Run ODE model with updated parameters
    solution <- run_ode_model(
      y = create.inits(updated_params), 
      time_points = time_points, 
      params = updated_params
    )
    
    # Store predictions
    for (j in 1:length(compartments)) {
      predictive_samples[i, , j] <- solution[, compartments[j]]
    }
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Plot predictive distributions
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  for (j in 1:length(compartments)) {
    # Prepare data for credible interval calculation
    comp_samples <- predictive_samples[, , j]
    
    # Plot median and credible intervals
    matplot(time_points, 
            t(comp_samples), 
            type = "l", 
            col = rgb(0, 0, 1, 0.1), 
            lty = 1,
            xlab = "Time (min)", 
            ylab = paste("Mass in", compartments[j], "(ug)"),
            main = paste("Predictive Distribution for", compartments[j])
    )
    
    # Add median line
    median_pred <- apply(comp_samples, 2, median)
    lines(time_points, median_pred, col = "red", lwd = 2)
    
    # Add 95% credible interval
    lower_ci <- apply(comp_samples, 2, quantile, 0.025)
    upper_ci <- apply(comp_samples, 2, quantile, 0.975)
    lines(time_points, lower_ci, col = "blue", lty = 2)
    lines(time_points, upper_ci, col = "blue", lty = 2)
    
    # Add legend with improved formatting
    legend("topright", 
           legend = c("Median", "95% CI", "Samples"), 
           col = c("red", "blue", rgb(0, 0, 1, 0.1)), 
           lty = c(1, 2, 1),
           cex = 0.7,           # Smaller text
           bty = "n",           # No box around legend
           inset = c(0.05, 0),  # Move slightly inside from the corner
           seg.len = 0.5)       # Shorter line segments
  }
  
  # Return the predictive samples for further analysis if needed
  return(predictive_samples)
}

# Prepare original parameters
original_params <- create_fixed_params(list(
  'BW' = 0.04,
  "admin.dose" = 1 * 0.04 * 1e03,
  "np_size" = 3/2,
  "admin.time" = 0, 
  "admin.type" = "iv",
  "sex" = "M"
))


# Generate and plot predictive distribution
predictive_dist <- generate_predictive_distribution(results, original_params)

#============================================Save plots==========================================================


save_parameter_plots <- function() {
  # Save prior vs posterior plots as PDF
  pdf("parameter_distributions.pdf", width = 10, height = 8)
  plot_prior_vs_predictive(results)
  dev.off()
  
  # Alternatively save as PNG with higher resolution
  png("parameter_distributions.png", width = 1200, height = 800, res = 120)
  plot_prior_vs_predictive(results)
  dev.off()
  
  cat("Parameter distribution plots saved as PDF and PNG\n")
}

# Add these lines after your generate_predictive_distribution function call
# to save the predictive distribution plots
save_predictive_plots <- function() {
  # Save predictive plots as PDF
  pdf("predictive_distributions.pdf", width = 10, height = 8)
  
  # Plot predictive distributions
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  compartments <- c("Mliver", "Mkidneys", "Mspleen", "Mlungs")
  for (j in 1:length(compartments)) {
    # Prepare data
    comp_samples <- predictive_dist[, , j]
    time_points <- seq(0, 10, 0.05)
    
    # Plot median and credible intervals
    matplot(time_points, 
            t(comp_samples), 
            type = "l", 
            col = rgb(0, 0, 1, 0.1), 
            lty = 1,
            xlab = "Time (min)", 
            ylab = paste("Mass in", compartments[j], "(ug)"),
            main = paste("Predictive Distribution for", compartments[j])
    )
    
    # Add median line
    median_pred <- apply(comp_samples, 2, median)
    lines(time_points, median_pred, col = "red", lwd = 2)
    
    # Add 95% credible interval
    lower_ci <- apply(comp_samples, 2, quantile, 0.025)
    upper_ci <- apply(comp_samples, 2, quantile, 0.975)
    lines(time_points, lower_ci, col = "blue", lty = 2)
    lines(time_points, upper_ci, col = "blue", lty = 2)
    
    # Add small legend
    legend("topright", 
           legend = c("Median", "95% CI", "Samples"), 
           col = c("red", "blue", rgb(0, 0, 1, 0.1)), 
           lty = c(1, 2, 1),
           cex = 0.7, 
           bty = "n",
           inset = c(0.05, 0),
           seg.len = 0.5)
  }
  
  dev.off()
  
  # Also save as PNG with higher resolution
  png("predictive_distributions.png", width = 1200, height = 800, res = 120)
  
  # Plot predictive distributions
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  for (j in 1:length(compartments)) {
    # Prepare data
    comp_samples <- predictive_dist[, , j]
    
    # Plot median and credible intervals
    matplot(time_points, 
            t(comp_samples), 
            type = "l", 
            col = rgb(0, 0, 1, 0.1), 
            lty = 1,
            xlab = "Time (min)", 
            ylab = paste("Mass in", compartments[j], "(ug)"),
            main = paste("Predictive Distribution for", compartments[j])
    )
    
    # Add median line
    median_pred <- apply(comp_samples, 2, median)
    lines(time_points, median_pred, col = "red", lwd = 2)
    
    # Add 95% credible interval
    lower_ci <- apply(comp_samples, 2, quantile, 0.025)
    upper_ci <- apply(comp_samples, 2, quantile, 0.975)
    lines(time_points, lower_ci, col = "blue", lty = 2)
    lines(time_points, upper_ci, col = "blue", lty = 2)
    
    # Add small legend
    legend("topright", 
           legend = c("Median", "95% CI", "Samples"), 
           col = c("red", "blue", rgb(0, 0, 1, 0.1)), 
           lty = c(1, 2, 1),
           cex = 0.7, 
           bty = "n",
           inset = c(0.05, 0),
           seg.len = 0.5)
  }
  
  dev.off()
  
  cat("Predictive distribution plots saved as PDF and PNG\n")
}

# Call these functions to save the plots
save_parameter_plots()
save_predictive_plots()
