library(ggplot2)

#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")
setwd("/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Final_results_plots/Validation")

# Load Results 

# 
Kemper_serum_results <- read.csv("Kemper_2003_serum_results.csv", header = T)
Kemper_iv_F_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "iv" &
                               Kemper_serum_results$sex == "F", ]
Kemper_iv_F_1$Experiment <- "Kemper | 1mg/kg | iv | F"

Kemper_oral_F_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "F", ]
Kemper_oral_F_1$Experiment <- "Kemper | 1mg/kg | oral | F"

Kemper_oral_F_5 <- Kemper_serum_results[Kemper_serum_results$Dose == 5 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "F", ]
Kemper_oral_F_5$Experiment <- "Kemper | 5mg/kg | oral | F"
Kemper_oral_F_25 <- Kemper_serum_results[Kemper_serum_results$Dose == 25 & Kemper_serum_results$Type == "oral" &
                                  Kemper_serum_results$sex == "F", ]
Kemper_oral_F_25$Experiment <- "Kemper | 25mg/kg | oral | F"

Kemper_iv_M_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "iv" &
                               Kemper_serum_results$sex == "M", ]
Kemper_iv_M_1$Experiment <- "Kemper | 1mg/kg | iv | M"

Kemper_oral_M_1 <- Kemper_serum_results[Kemper_serum_results$Dose == 1 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "M", ]
Kemper_oral_M_1$Experiment <- "Kemper | 1mg/kg | oral | M"

Kemper_oral_M_5 <- Kemper_serum_results[Kemper_serum_results$Dose == 5 & Kemper_serum_results$Type == "oral" &
                                 Kemper_serum_results$sex == "M", ]
Kemper_oral_M_5$Experiment <- "Kemper | 5mg/kg | oral | M"
Kemper_oral_M_25 <- Kemper_serum_results[Kemper_serum_results$Dose == 25 & Kemper_serum_results$Type == "oral" &
                                  Kemper_serum_results$sex == "M", ]
Kemper_oral_M_25$Experiment <- "Kemper | 25mg/kg | oral | M"
# 
# Kudo_high <- Kudo_results[Kudo_results$Dose == 16.560, ] 
# Kudo_high$Experiment <-  "Kudo | 16mg/kg | iv"
# 
# Kim_results <- read.csv("Kim_results.csv", header = T)
# Kim_iv <- Kim_results[Kim_results$Type == "iv",]
# Kim_iv$Experiment <-  "Kim | 1mg/kg | iv"
# 
# Kim_oral <- Kim_results[Kim_results$Type == "oral",]
# Kim_oral$Experiment <-  "Kim | 1mg/kg | oral"


Experiment <- scales::hue_pal()(8)

names(Experiment) <- unique(c( Kemper_iv_F_1$Experiment, Kemper_iv_M_1$Experiment,
                              Kemper_oral_F_1$Experiment,Kemper_oral_M_1$Experiment, Kemper_oral_F_5$Experiment, 
                               
                               Kemper_oral_M_5$Experiment,  Kemper_oral_F_25$Experiment, Kemper_oral_M_25$Experiment))



scatter_plot <- ggplot()+
  # Error range lines
  geom_abline(intercept = log10(1), slope = 1, linetype = "solid", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/2), slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "black", linewidth = 1.2) +
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "black", linewidth = 1.2) +
  
  geom_point(data = Kemper_iv_F_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_F_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_F_5, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_F_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_iv_M_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_M_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_M_5, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_oral_M_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  # Axes
  scale_x_log10() +
  scale_y_log10() +
  
  # Color and shape scales
  scale_color_manual(values = Experiment) +
  
  # Labels
  labs(
    x = expression("Observed Concentration (mg PFOA / L serum)"),
    y = expression("Predicted Concentration (mg PFOA / L serum)"),
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


#print(scatter_plot)

ggsave("validation_plot_PFOA_serum.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


