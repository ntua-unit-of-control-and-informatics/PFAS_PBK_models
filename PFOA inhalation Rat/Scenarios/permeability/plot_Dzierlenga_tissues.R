library(ggpubr)
library(grid)
koff_alb <- 1
f_alb_avail <- 0.04
f_lfabp_avail <- 2
# RAFOatp2_Int <- 1
# estimated_params[3] <- 1#6382.949
estimated_params[9]  <- f_alb_avail
# estimated_params[7] <- 1
# estimated_params[8] <- 1e-6#1.789484
Ka <-5e5
VmK_api_f <-0
VmK_api_m <- 0
KmK_api <- 1e20
koff_fabp <- 1
###############################################################################################
# Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[7]])
params <- c(fixed_params[[7]], variable_params)
params$koff_alb <- koff_alb
params$kon_alb <- Ka * koff_alb #1/M/s
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05
params$koff_fabp <- koff_fabp
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,864,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]

#############################################################################################
# Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
BW <- 0.2  # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[8]])
params <- c(fixed_params[[8]], variable_params)
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05
params$koff_fabp <- koff_fabp

params$koff_alb <- koff_alb
params$kon_alb <- Ka * koff_alb #1/M/s
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]


preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] <- preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000

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
    #scale_y_log10()+
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

# Convert Dzierlenga ORAL male tissues from long to wide format using reshape
experiment7 <- reshape(dzi_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(dzi_OR_Mtissues$Tissue))


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))

# Put the experiments in a list
experiments <- list(experiment7 = experiment7, 
                    experiment8 = experiment8)


colnames(preds_dzi_OR_Mtissues) <- c("Time","Liver","Kidney","Brain")
colnames(preds_dzi_OR_Ftissues) <- colnames(preds_dzi_OR_Mtissues)

# Create a list containing the corresponding predictions
simulations <- list(predictions7 = preds_dzi_OR_Mtissues, predictions8 = preds_dzi_OR_Ftissues)


all_plots <- list()

for(i in 1:length(experiments)){
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  
  compartments <- names(predictions)[2:length(predictions)]
  
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment)
  })
  
  all_plots <- c(all_plots, plots)
}

final_plot <- do.call(ggarrange, c(all_plots,
                                   ncol = 2, nrow = ceiling(length(all_plots) / 2),
                                   common.legend = TRUE, legend = "right"))

final_plot <- final_plot + theme(plot.margin = unit(c(0,0,0,0), "pt"))

print(final_plot)


