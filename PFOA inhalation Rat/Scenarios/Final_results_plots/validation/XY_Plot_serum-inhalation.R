library(ggplot2)

#setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results")
setwd("/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Final_results_plots/Validation")

# Load Results 

# 
Hinderliter_2006_results <- read.csv("Hinderliter_2006_results.csv", header = T)
Hinderliter_nasal_M_1 <- Hinderliter_2006_results[Hinderliter_2006_results$Dose == 1 & Hinderliter_2006_results$Type == "nasal" &
                               Hinderliter_2006_results$sex == "M", ]
Hinderliter_nasal_M_1$Experiment <- "Hinderliter | 1mg/m3 | nasal | M"

Hinderliter_nasal_F_1 <- Hinderliter_2006_results[Hinderliter_2006_results$Dose == 1 & Hinderliter_2006_results$Type == "nasal" &
                                                    Hinderliter_2006_results$sex == "F", ]
Hinderliter_nasal_F_1$Experiment <- "Hinderliter | 1mg/m3 | nasal | F"

Hinderliter_nasal_M_10 <- Hinderliter_2006_results[Hinderliter_2006_results$Dose == 10 & Hinderliter_2006_results$Type == "nasal" &
                                                    Hinderliter_2006_results$sex == "M", ]
Hinderliter_nasal_M_10$Experiment <- "Hinderliter | 10mg/m3 | nasal | M"

Hinderliter_nasal_F_10 <- Hinderliter_2006_results[Hinderliter_2006_results$Dose == 10 & Hinderliter_2006_results$Type == "nasal" &
                                                    Hinderliter_2006_results$sex == "F", ]
Hinderliter_nasal_F_10$Experiment <- "Hinderliter | 10mg/m3 | nasal | F"

Hinderliter_nasal_M_25 <- Hinderliter_2006_results[Hinderliter_2006_results$Dose == 25 & Hinderliter_2006_results$Type == "nasal" &
                                                     Hinderliter_2006_results$sex == "M", ]
Hinderliter_nasal_M_25$Experiment <- "Hinderliter | 25mg/m3 | nasal | M"

Hinderliter_nasal_F_25 <- Hinderliter_2006_results[Hinderliter_2006_results$Dose == 25 & Hinderliter_2006_results$Type == "nasal" &
                                                     Hinderliter_2006_results$sex == "F", ]
Hinderliter_nasal_F_25$Experiment <- "Hinderliter | 25mg/m3 | nasal | F"

Experiment <- scales::hue_pal()(6)

names(Experiment) <- unique(c( Hinderliter_nasal_M_1$Experiment, Hinderliter_nasal_F_1$Experiment,
                               Hinderliter_nasal_M_10$Experiment,Hinderliter_nasal_F_10$Experiment,
                               Hinderliter_nasal_M_25$Experiment, Hinderliter_nasal_F_25$Experiment))


ggplot(data = Hinderliter_2006_results[Hinderliter_2006_results$Dose == 10 & Hinderliter_2006_results$sex == "F",])+
  geom_line(aes(x =Time, y= Predicted ))+
  geom_point(aes(x =Time, y= Observed ))


scatter_plot <- ggplot()+
  geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5, alpha = 0.7) +  # Identity line in log10 scale
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # +3-fold error line
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", color = "red", linewidth = 1.5, alpha = 0.7) +  # -3-fold error line
  geom_abline(intercept = log10(10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # +10-fold error line
  geom_abline(intercept = log10(1/10), slope = 1, linetype = "dotdash", color = "blue", linewidth = 1.5, alpha = 0.7) +  # -10-fold error line
  geom_point(data = Hinderliter_nasal_M_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Hinderliter_nasal_F_1, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Hinderliter_nasal_M_10, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Hinderliter_nasal_F_10, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Hinderliter_nasal_M_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Hinderliter_nasal_F_25, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
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

ggsave("validation_plot_PFOA_nasal_serum.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


