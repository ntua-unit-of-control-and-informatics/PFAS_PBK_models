library(ggplot2)

#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")
setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_Scenario24l")

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



scatter_plot <- ggplot()+
  geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5, alpha = 0.7) +  # Identity line in log10 scale
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # +3-fold error line
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # -3-fold error line
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # +10-fold error line
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # -10-fold error line
  geom_point(data = Cui_Urine_M_5, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Cui_Urine_M_20, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Cui_Feces_M_5, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Cui_Feces_M_20, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Lupton_Feces_F, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Lupton_Urine_F, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
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

#print(scatter_plot)

ggsave("validation_plot_PFOA_excreta.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


