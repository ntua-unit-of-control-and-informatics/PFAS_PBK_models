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

# Set up simulations for the 21th case, i.e. Dzierlenga 2021, IV female serum
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[21]])
params <- c(fixed_params[[21]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,  192, 1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_IV_Fserum <-  solution[, c("time", "Cplasma")]

#########################################################################################
# Set up simulations for the 22th case, i.e. Dzierlenga 2021, ORAL female serum low
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[22]])
params <- c(fixed_params[[22]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0, 96, 1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]

################################################################################################
# Set up simulations for the 23th case, i.e. Dzierlenga 2021, ORAL female serum medium
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[23]])
params <- c(fixed_params[[23]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0, 192, 1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Fserum_medium <-  solution[, c("time", "Cplasma")]

########################################################################################
# Set up simulations for the 24th case, i.e. Dzierlenga 2021, ORAL female serum high
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[24]])
params <- c(fixed_params[[24]], variable_params)
params$koff_alb <- koff_alb
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time <- seq(0, 96, 1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_dzi_OR_Fserum_high <-  solution[, c("time", "Cplasma")]

preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] <- preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] /1000
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] <- preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] /1000
preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] <- preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] /1000


#Convert Dzierlenga 2021, IV female serum from long to wide format using reshape
experiment21 <- reshape(dzi_IV_Fserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment21) <- c("Time",unique(dzi_IV_Fserum$Tissue))

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment22 <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment22) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL female serum medium from long to wide format using reshape
experiment23 <- reshape(dzi_OR_Fserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment23) <- c("Time",unique(dzi_OR_Fserum_medium$Tissue))

#Convert Dzierlenga 2021, ORAL female serum high from long to wide format using reshape
experiment24 <- reshape(dzi_OR_Fserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment24) <- c("Time",unique(dzi_OR_Fserum_high$Tissue))



# Put the experiments in a list
experiments <- list(experiment21 = experiment21, 
                    experiment22 = experiment22, experiment23 = experiment23, experiment24 = experiment24)

colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

# Create a list containing the corresponding predictions
simulations <- list(predictions21 =preds_dzi_IV_Fserum, predictions22 =preds_dzi_OR_Fserum_low,
                    predictions23 =preds_dzi_OR_Fserum_medium,
                    predictions24 =preds_dzi_OR_Fserum_high)


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


