library(ggplot2)

setwd("/Users/ptsir/Documents/GitHub/PBK_Grouping/PFOA/Validation/Validation_results")
column_names <- c("Study",	"Dose",	"Tissue",	"Type",	"Observed",	"Predicted")

# Load Results 
Cui_results <- read.csv("Cui_2008_results.csv", header = T)
Cui_low <- Cui_results[Cui_results$Dose == 5, ]
Cui_low$Experiment <- "Cui | 5mg/kg | oral"
Cui_high <- Cui_results[Cui_results$Dose == 20, ]
Cui_high$Experiment <- "Cui | 20mg/kg | oral"

Kudo_results <- read.csv("Kudo_results.csv", header = T)
Kudo_low <- Kudo_results[Kudo_results$Dose == 0.041, ] 
Kudo_low$Experiment <-  "Kudo | 0.04mg/kg | iv"

Kudo_high <- Kudo_results[Kudo_results$Dose == 16.560, ] 
Kudo_high$Experiment <-  "Kudo | 16mg/kg | iv"

Kim_results <- read.csv("Kim_results.csv", header = T)
Kim_iv <- Kim_results[Kim_results$Type == "iv",]
Kim_iv$Experiment <-  "Kim | 1mg/kg | iv"

Kim_oral <- Kim_results[Kim_results$Type == "oral",]
Kim_oral$Experiment <-  "Kim | 1mg/kg | oral"

Tissue_markers <-  c(0:10,14)
names(Tissue_markers) <- c( "Lung",   "Spleen",  "Liver",   "Kidney", "Plasma", "Heart",
                            "Brain","Testis", "Stomach", "Intestines", "Carcass")

Experiment <- scales::hue_pal()(6)
Experiment <- c("#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#000000", "#D55E00")

names(Experiment) <- unique(c(Cui_low$Experiment, Cui_high$Experiment,
                              Kudo_low$Experiment,   Kudo_high$Experiment,
                              Kim_iv$Experiment, Kim_oral$Experiment))



scatter_plot <- ggplot()+
  geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5,alpha = 0.7) +  # Identity line in log10 scale
  geom_point(data = Cui_low, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Cui_high, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kudo_low, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kudo_high, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kim_iv, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  geom_point(data = Kim_oral, aes(x=Observed, y=Predicted, color = Experiment, shape = Tissue), size=4, stroke = 1.5)+
  
  scale_y_log10()+
  scale_x_log10()+
  
  scale_color_manual(values = Experiment)+
                  
  theme(legend.spacing.y = unit(1, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))+
  
  scale_shape_manual(values = Tissue_markers)+
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

ggsave("validation_plot_PFOA.png", scatter_plot, width = 11, height = 7, units = "in", dpi = 300)


