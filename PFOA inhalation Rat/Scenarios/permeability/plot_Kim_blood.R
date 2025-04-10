library(ggpubr)
library(grid)
koff_alb <- 1
f_alb_avail <- 0.05
f_lfabp_avail <-0.01
# RAFOatp2_Int <- 1
# estimated_params[3] <- 1#6382.949
estimated_params[9]  <- f_lfabp_avail
estimated_params[10]  <- f_alb_avail

# estimated_params[7] <- 1
# estimated_params[8] <- 1e-6#1.789484
Ka <-5e5
VmK_api_f <-0
VmK_api_m <- 0
KmK_api <- 1e20
# Set up simulations for the 9th case, i.e. Kim (2016) ORAL male blood
BW <- 0.25  #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[9]])
params <- c(fixed_params[[9]], variable_params)
params$koff_alb <- koff_alb
params$kon_alb <- Ka * koff_alb #1/M/s
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-01, atol = 1e-01))

preds_kim_OR_Mblood <-  solution[, c("time", "Cplasma")]

######################################################################################
# Set up simulations for the 10th case, i.e. Kim (2016) IV male blood
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "M" 

variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[10]])
params <- c(fixed_params[[10]], variable_params)
params$koff_alb <- koff_alb
params$kon_alb <- Ka * koff_alb #1/M/s
params$VmK_api <- VmK_api_m
params$KmK_api <- KmK_api
#params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs

sample_time=c(0, 5/60, seq(1,288,1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-01, atol = 1e-01))

preds_kim_IV_Mblood <-  solution[, c("time", "Cplasma")]

##########################################################################################
# Set up simulations for the 25th case, i.e. Kim (2016) ORAL female blood
BW <- 0.25  #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[25]])
params <- c(fixed_params[[25]], variable_params)
params$koff_alb <- koff_alb
params$kon_alb <- Ka * koff_alb #1/M/s
params$VmK_api <- VmK_api_f
params$KmK_api <- KmK_api
#params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time= c(seq(0, 1.2, 0.2), seq(1.5,24,0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-01, atol = 1e-01))

preds_kim_OR_Fblood <-  solution[, c("time", "Cplasma")]

########################################################################################
# Set up simulations for the 26th case, i.e. Kim (2016) IV female blood
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[26]])
params <- c(fixed_params[[26]], variable_params)
params$koff_alb <- koff_alb
params$kon_alb <- Ka * koff_alb #1/M/s
params$VmK_api <- VmK_api_f
params$KmK_api <- KmK_api
#params$CFabpLT_init <- f_lfabp_avail*4.191842e-05

inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0, 25, 0.5)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-01, atol = 1e-01))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]




preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] <- preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] /1000
preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] <- preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] /1000
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000



experiment9 <- reshape(kim_OR_Mblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9) <- c("Time",unique(kim_OR_Mblood$Tissue))


# Convert Kim IV male blood from long to wide format using reshape
experiment10 <- reshape(kim_IV_Mblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment10) <- c("Time",unique(kim_IV_Mblood$Tissue))


#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment25 <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment25) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment26<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment26) <- c("Time",unique(kim_IV_Fblood$Tissue))



# Put the experiments in a list
experiments <- list(experiment9 = experiment9, 
                    experiment10 = experiment10,experiment25 = experiment25,experiment26 = experiment26)



colnames(preds_kim_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_kim_IV_Mblood) <- c ("Time", "Plasma")


colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")

# Create a list containing the corresponding predictions
simulations <- list(
  predictions9 = preds_kim_OR_Mblood,
  predictions10 = preds_kim_IV_Mblood,predictions25 = preds_kim_OR_Fblood, 
  predictions26 = preds_kim_IV_Fblood)

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
   # scale_y_log10()+
    
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


