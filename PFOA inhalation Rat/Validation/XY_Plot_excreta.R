library(ggplot2)

setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")

# Load Results 

# 
Kemper_excreta_results <- read.csv("Kemper_2003_Excreta_results.csv", header = T)

Kemper_Urine_F_25 <- Kemper_excreta_results[Kemper_excreta_results$Dose == 25 & 
                                             Kemper_excreta_results$Tissue  == "Urine" &
                                          Kemper_excreta_results$sex == "F", ]
Kemper_Urine_F_25$Experiment <- "Kemper | 25mg/kg | Urine | F"


Kemper_Feces_F_25 <- Kemper_excreta_results[Kemper_excreta_results$Dose == 25 & 
                                              Kemper_excreta_results$Tissue  == "Feces" &
                                              Kemper_excreta_results$sex == "F", ]
Kemper_Feces_F_25$Experiment <- "Kemper | 25mg/kg | Feces | F"



Kemper_Urine_M_25 <- Kemper_excreta_results[Kemper_excreta_results$Dose == 25 & Kemper_excreta_results$Tissue == "Urine" &
                                             Kemper_excreta_results$sex == "M", ]
Kemper_Urine_M_25$Experiment <- "Kemper | 25mg/kg | Urine | M"



Kemper_Feces_M_25 <- Kemper_excreta_results[Kemper_excreta_results$Dose == 25 & Kemper_excreta_results$Tissue == "Feces" &
                                             Kemper_excreta_results$sex == "M", ]
Kemper_Feces_M_25$Experiment <- "Kemper | 25mg/kg | Feces | M"



Experiment <- scales::hue_pal()(4)

names(Experiment) <- unique(c( Kemper_Urine_F_25, Kemper_Feces_F_25,
                               Kemper_Urine_M_25,Kemper_Feces_M_25))



scatter_plot <- ggplot()+
  geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5, alpha = 0.7) +  # Identity line in log10 scale
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # +3-fold error line
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # -3-fold error line
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # +10-fold error line
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # -10-fold error line
  geom_point(data = Kemper_Urine_F_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_Feces_F_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_Urine_M_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kemper_Feces_M_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  scale_y_log10()+
  scale_x_log10()+
  scale_color_manual(values = Experiment)+
  theme(legend.spacing.y = unit(1, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))+
  theme_light()+
  labs(y = "Predicted PFOA (mg/L tissue)",
       x = "Observed PFOA ( mg/L tissue)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.75, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=12,  hjust = 0),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        )

print(scatter_plot)

#ggsave("validation_plot_PFOA.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


