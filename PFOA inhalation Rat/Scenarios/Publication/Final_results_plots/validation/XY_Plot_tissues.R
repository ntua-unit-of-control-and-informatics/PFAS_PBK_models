library(ggplot2)
library(dplyr)

# Set working directory
setwd("/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Publication/Final_results_plots/Validation")

# Load Results 
Cui_results <- read.csv("Cui_2008_results.csv", header = TRUE)
lupton_results <- read.csv("Lupton_2020_tissue_results.csv", header = TRUE)

# Filter by dose and add experiment label
Cui_low <- Cui_results %>% 
  filter(Dose == 5, !is.na(Observed), !is.na(Predicted), Observed > 0, Predicted > 0) %>%
  mutate(Experiment = "Cui | 5mg/kg | oral | M")

Cui_high <- Cui_results %>% 
  filter(Dose == 20, !is.na(Observed), !is.na(Predicted), Observed > 0, Predicted > 0) %>%
  mutate(Experiment = "Cui | 20mg/kg | oral | M")

lupton <- lupton_results %>% 
  filter( !is.na(Predicted), Observed > 0, Predicted > 0) %>%
  mutate(Experiment = "Lupton | 0.047mg/kg | oral | F")

# Tissue marker names (optional if used elsewhere)
Tissue_markers <- c(0:10, 15, 16)
names(Tissue_markers) <- c("Lung", "Spleen", "Liver", "Kidney", "Plasma", "Heart",
                           "Brain", "Testis", "Stomach", "Intestines", "Carcass", "Skin")

# Define color palette for experiments
Experiment <- scales::hue_pal()(3)
names(Experiment) <- unique(c(Cui_low$Experiment, Cui_high$Experiment, lupton$Experiment ))

# Manually assign shapes for up to 8 tissues
shape_values <- c(0:6, 15, 16)

# Build scatter plot
scatter_plot <- ggplot() +
  # Error range lines
  geom_abline(intercept = log10(1), slope = 1, linetype = "solid", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/2), slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "black", linewidth = 1.2) +
  
  # Data points
  geom_point(data = Cui_low, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = Cui_high, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = lupton, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  
  # Log-scale axes
  scale_x_log10() +
  scale_y_log10() +
  
  # Manual color and shape scales
  scale_color_manual(values = Experiment) +
  scale_shape_manual(values = shape_values) +
  
  # Axis labels
  labs(
    x = expression("Observed Concentration (mg PFOA / L tissue)"),
    y = expression("Predicted Concentration (mg PFOA / L tissue)"),
    color = "Experiment",
    shape = "Tissue"
  ) +
  
  # Unified theme
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(0.5, "cm")
  )

# Save figure
ggsave("validation_plot_PFOA_tissues.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)
