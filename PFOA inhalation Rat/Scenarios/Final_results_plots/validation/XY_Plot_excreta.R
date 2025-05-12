library(ggplot2)

#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")
setwd("/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Final_results_plots/Validation")

# Load Results 

# Results from Cui
Cui_2010_Excreta_results <- read.csv("Cui_2010_Excreta_results.csv", header = T)

Cui_Urine_M_5 <- Cui_2010_Excreta_results[Cui_2010_Excreta_results$Dose == 5 & 
                                            Cui_2010_Excreta_results$Tissue  == "Urine" &
                                            Cui_2010_Excreta_results$Sex == "M", ]
Cui_Urine_M_5$Experiment <- "Cui | 5mg/kg | Urine | M"


Cui_Urine_M_20 <- Cui_2010_Excreta_results[Cui_2010_Excreta_results$Dose == 20 & 
                                            Cui_2010_Excreta_results$Tissue  == "Urine" &
                                            Cui_2010_Excreta_results$Sex == "M", ]
Cui_Urine_M_20$Experiment <- "Cui | 20mg/kg | Urine | M"



Cui_Feces_M_5 <- Cui_2010_Excreta_results[Cui_2010_Excreta_results$Dose == 5 & 
                                            Cui_2010_Excreta_results$Tissue  == "Feces" &
                                            Cui_2010_Excreta_results$Sex == "M", ]
Cui_Feces_M_5$Experiment <- "Cui | 5mg/kg | Feces | M"


Cui_Feces_M_20 <- Cui_2010_Excreta_results[Cui_2010_Excreta_results$Dose == 20 & 
                                            Cui_2010_Excreta_results$Tissue  == "Feces" &
                                            Cui_2010_Excreta_results$Sex == "M", ]
Cui_Feces_M_20$Experiment <- "Cui | 20mg/kg | Feces | M"


#Results from Lupton 2020

Lupton_2020_Excreta_results <- read.csv("Lupton_2020_excreta_results.csv", header = T, colClasses = c("sex" = "character"))


Lupton_Feces_F <- Lupton_2020_Excreta_results[Lupton_2020_Excreta_results$Dose == 0.047 & 
                                              Lupton_2020_Excreta_results$Tissue  == "Feces" &
                                              Lupton_2020_Excreta_results$sex == "F", ]
Lupton_Feces_F$Experiment <- "Lupton | 0.047mg/kg | Feces | F"


Lupton_Urine_F <- Lupton_2020_Excreta_results[Lupton_2020_Excreta_results$Dose == 0.047 & 
                                              Lupton_2020_Excreta_results$Tissue  == "Urine" &
                                              Lupton_2020_Excreta_results$sex == "F", ]
Lupton_Urine_F$Experiment <- "Lupton | 0.047mg/kg | Urine | F"


Experiment <- scales::hue_pal()(6)

names(Experiment) <- unique(c( Cui_Urine_M_5$Experiment, Cui_Urine_M_20$Experiment,
                               Cui_Feces_M_5$Experiment, Cui_Feces_M_20$Experiment,
                               Lupton_Feces_F$Experiment, Lupton_Urine_F$Experiment))
# Create legend data
legend_lines <- data.frame(
  x = 1, y = 1,
  type = factor(c("1x (identity)", "±2-fold", "±3-fold", "±10-fold"),
                levels = c("1x (identity)", "±2-fold", "±3-fold", "±10-fold"))
)



scatter_plot <- ggplot() +
  # Error range lines
  geom_abline(intercept = log10(1), slope = 1, linetype = "solid", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/2), slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "black", linewidth = 1.2) +

  # Data points (original style)
  geom_point(data = Cui_Urine_M_5, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = Cui_Urine_M_20, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = Cui_Feces_M_5, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = Cui_Feces_M_20, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = Lupton_Feces_F, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  geom_point(data = Lupton_Urine_F, aes(x = Observed, y = Predicted, color = Experiment, shape = Tissue), size = 4, stroke = 1.5) +
  
  # Axes
  scale_x_log10() +
  scale_y_log10() +
  
  # Color and shape scales
  scale_color_manual(values = Experiment) +

  # Labels
  labs(
    x = expression("Observed Mass (mg PFOA)"),
    y = expression("Predicted Mass (mg PFOA )"),
    color = "Experiment",
    shape = "Tissue"
  ) +
  
  # Theme to match create.plots
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

# Save the final figure
ggsave("validation_plot_PFOA_excreta.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)
