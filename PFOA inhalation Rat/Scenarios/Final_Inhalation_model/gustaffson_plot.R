
# Set up simulations for the 8th case, i.e. Gustafsson (2022) Inhalation male blood
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
duration <- 0.375 #hours, 22.5 min
admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
depfr_head <-0
depfr_AF <-1
k = 9#partition of administration packages
admin.dose <- rep((admin.dose_per_g*BW*1000)/k, length.out = k) #ug PFOA, for 22.5 min inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "inh"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF)
params[[8]] <- create_fixed_params(user_input)


BW <- 0.5125  #kg, from Gustafsson et al., 2022
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[8]])
params <- c(fixed_params[[8]], variable_params)
params$VmLu_Oatp_ap <- 0
inits <- create.inits(params)
events <- create.events(params)

#sample_time=seq(0,0.047,0.001)
sample_time=seq(0,48,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_gus_INH_Mblood <-  solution[, c("time", "Cplasma")]
preds_gus_INH_Mtissues <-  solution[, c("time", "CalveolarLF","Cliver", "Clungtissue", "Ckidney")]


#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets
preds_gus_INH_Mblood[,2:dim(preds_gus_INH_Mblood)[2]] <- preds_gus_INH_Mblood[,2:dim(preds_gus_INH_Mblood)[2]] /1000
preds_gus_INH_Mtissues[,2:dim(preds_gus_INH_Mtissues)[2]] <- preds_gus_INH_Mtissues[,2:dim(preds_gus_INH_Mtissues)[2]] /1000


# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("PFOA concentration (" * mu* "g/g tissue)" ),
         x = "Time (hours)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("", values=cls,
                       guide = guide_legend(override.aes =
                                              list(shape = c(16,NA),
                                                   linetype = c(0,1))))+
    theme_light() + 
    theme(legend.position=c(1,1), 
          legend.justification=c(0, 1), 
          legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          axis.title=element_text(size=14),
          legend.text = element_text(size=14)
    )
  
}




# Convert Gustafsson Inhalation male blood from long to wide format using reshape
experiment_inh_8 <- reshape(gus_INH_Mblood[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_8) <- c("Time",unique(gus_INH_Mblood$Tissue))


# Convert Gustafsson Inhalation male tissues from long to wide format using reshape
experiment_inh_9 <- reshape(gus_INH_Mtissues[c("Tissue" ,"Time_hours", 
                                               "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_9) <- c("Time",unique(gus_INH_Mtissues$Tissue))

# Put the experiments in a list
experiments <- list(
                    experiment_inh_8 = experiment_inh_8,experiment_inh_9 = experiment_inh_9)



colnames(preds_gus_INH_Mblood) <- c ("Time", "Plasma")
colnames(preds_gus_INH_Mtissues) <- c ("Time", "ALF", "Liver", "Lung", "Kidney")

# Create a list containing the corresponding predictions
simulations <- list( predictions8 = preds_gus_INH_Mblood, 
                    predictions9 = preds_gus_INH_Mtissues)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment )
  })
  if(length(compartments) == 1){
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                               common.legend = TRUE, legend = "right"))
    
  }else{
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                               common.legend = TRUE, legend = "right"))
  }
  
  
  plot.margin=unit(c(0,0,0,0), "pt")
  
  print(final_plot)
}
