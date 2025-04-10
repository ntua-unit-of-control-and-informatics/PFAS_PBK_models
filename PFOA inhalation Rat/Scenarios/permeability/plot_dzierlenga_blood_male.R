library(ggpubr)
library(grid)
koff_alb <- 1
f_alb_avail <- 0.1
# RAFOatp2_Int <- 1
# estimated_params[3] <- 1#6382.949
estimated_params[9]  <- f_alb_avail
# estimated_params[7] <- 1
# estimated_params[8] <- 1e-6#1.789484
Ka <-2e5
VmK_api_f <-0
VmK_api_m <- 0
KmK_api <- 1e5
f_lfabp_avail <- 8

# Set up simulations for the 17th case, i.e. Dzierlenga 2021, IV male serum
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[17]])
params <- c(fixed_params[[17]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.083, 0.25, 0.5, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_IV_Mserum <-  solution[, c("time", "Cplasma")]

##############################################################################################
# Set up simulations for the 18th case, i.e. Dzierlenga 2021, ORAL male serum low
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[18]])
params <- c(fixed_params[[18]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Mserum_low <-  solution[, c("time", "Cplasma")]

########################################################################################################
# Set up simulations for the 19th case, i.e. Dzierlenga 2021, ORAL male serum medium
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[19]])
params <- c(fixed_params[[19]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Mserum_medium <-  solution[, c("time", "Cplasma")]

##############################################################################################
# Set up simulations for the 20th case, i.e. Dzierlenga 2021, ORAL male serum high
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[20]])
params <- c(fixed_params[[20]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Mserum_high <-  solution[, c("time", "Cplasma")]


preds_dzi_IV_Mserum[,2:dim(preds_dzi_IV_Mserum)[2]] <- preds_dzi_IV_Mserum[,2:dim(preds_dzi_IV_Mserum)[2]] /1000
preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] <- preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] /1000
preds_dzi_OR_Mserum_medium[,2:dim(preds_dzi_OR_Mserum_medium)[2]] <- preds_dzi_OR_Mserum_medium[,2:dim(preds_dzi_OR_Mserum_medium)[2]] /1000
preds_dzi_OR_Mserum_high[,2:dim(preds_dzi_OR_Mserum_high)[2]] <- preds_dzi_OR_Mserum_high[,2:dim(preds_dzi_OR_Mserum_high)[2]] /1000

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
    scale_y_log10()+
    
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

# Convert Dzierlenga 2021, IV male serum from long to wide format using reshape
experiment17 <- reshape(dzi_IV_Mserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment17) <- c("Time",unique(dzi_IV_Mserum$Tissue))

# Convert Dzierlenga 2021, ORAL male serum low from long to wide format using reshape
experiment18 <- reshape(dzi_OR_Mserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment18) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))

# Convert Dzierlenga 2021, ORAL male serum medium from long to wide format using reshape
experiment19 <- reshape(dzi_OR_Mserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment19) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL male serum high from long to wide format using reshape
experiment20 <- reshape(dzi_OR_Mserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment20) <- c("Time",unique(dzi_OR_Mserum_high$Tissue))




# Put the experiments in a list
experiments <- list(experiment17 = experiment17, experiment18 = experiment18, 
                    experiment19 = experiment19, experiment20 = experiment20)



colnames(preds_dzi_IV_Mserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_high) <- c ("Time", "Serum")


# Create a list containing the corresponding predictions
simulations <- list( prediction17 =preds_dzi_IV_Mserum, predictions18 =preds_dzi_OR_Mserum_low,
                     predictions19 =preds_dzi_OR_Mserum_medium, predictions20 =preds_dzi_OR_Mserum_high 
                     )


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


