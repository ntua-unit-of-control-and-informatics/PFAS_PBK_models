


# Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_IV_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

# Set up simulations for the 6th case, i.e. kim (2016) ORAL female tissues
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]


# Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
BW <- 0.2  # body weight (kg) not reported
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]


plot(solution$time, solution$Cliver/1000, type = "l")

#################################################################################
#--------------------------------------------------------------------------------
#                                Lupton 2020
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 11th case, i.e. Lupton (2020) ORAL female tissues
BW <- 0.184  # body weight (kg) not reported
admin.dose_per_g <- 0.047 # administered dose in mg PFOA/kg BW 
admin.dose_single <- (admin.dose_per_g*BW*1e03)/2 #ug PFOA
admin.time <- seq(0,13.5*24,12) #time when doses are administered, in hours
admin.dose <- rep(admin.dose_single, length(admin.time))

admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,324,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_Lup_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cblood", "Cskin")]

sample_time=seq(0,384,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_Lup_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cblood", "Cskin")]
preds_Lup_OR_Ffeces <-  solution[, c("time", "Mfeces")]
preds_Lup_OR_Furine <-  solution[, c("time", "Murine")]

# Convert the Lupton excreta data to cumulative masses
#Feces
exp_data <- dataset$df12 # retrieve data of Lupton (2020) ORAL female feces
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
# Convert Lupton's excreta data to cumulative mass
exp_time <- exp_data$time
#Estimate fecal dry mass 
Mfeces_wet <- (8.18/0.21) #g
Mfeces_dry <- Mfeces_wet*0.2 # Enqi et al., 2021, control rats:20.1% water content, https://doi.org/10.3389/fcimb.2020.581974

# Estimate the mass of feces by multiplying concentration by dry mass
Lupton_Ffeces <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])*Mfeces_dry
# Estimate cumulative fecal mass
obs_Lup_OR_Ffeces_cum <- cumsum(Lupton_Ffeces)

#Urine
exp_data <- dataset$df13 # retrieve data of Lupton (2020) ORAL female feces
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
# Convert Lupton's excreta data to cumulative mass
exp_time <- exp_data$time

#Estimate urine volume
Qurine_daily <- 85 * BW #  (ml/d/kg)*BW  --> mL/d, Schmidt et al., 2001, doi:10.1002/nau.1006
Lupton_urine <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])*Qurine_daily
# Estimate cumulative fecal mass
obs_Lup_OR_Furine_cum <- cumsum(Lupton_urine)

# Set up simulations for the 22nd case, i.e. Dzierlenga 2021, IV female serum 
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.083, 0.25, seq(0.5, 192, 0.5))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_IV_Fserum <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 23d case, i.e. Dzierlenga 2021, ORAL female serum low
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0, 96, 0.25)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 25th case, i.e. Dzierlenga 2021, ORAL female serum medium
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.25, seq(1, 192, 0.5))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_medium <-  solution[, c("time", "Cplasma")]



# Set up simulations for the 25th case, i.e. Dzierlenga 2021, ORAL female serum high
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 320 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0, 96, 0.25)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_high <-  solution[, c("time", "Cplasma")]



#################################################################################
#--------------------------------------------------------------------------------
#                                Kim 2016 female
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 9th case, i.e. Kim (2016) ORAL male blood
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_OR_Fblood <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 10th case, i.e. Kim (2016) IV male blood
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]




#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets
preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] <- preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] /1000
preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] <- preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000
preds_Lup_OR_Ftissues[,2:dim(preds_Lup_OR_Ftissues)[2]] <- preds_Lup_OR_Ftissues[,2:dim(preds_Lup_OR_Ftissues)[2]] /1000
preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] <- preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] /1000
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] <- preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] /1000
preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] <- preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] /1000
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000


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


# Convert Kim IV female tissues from long to wide format using reshape
experiment5 <- reshape(kim_IV_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",kim_IV_Ftissues$Tissue )

# Convert Kim ORAL female tissues from long to wide format using reshape
experiment6 <- reshape(kim_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",kim_OR_Ftissues$Tissue )


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))

# Convert Lupton ORAL female tissues from long to wide format using reshape
experiment11 <- reshape(Lup_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                          "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11) <- c("Time",unique(Lup_OR_Ftissues$Tissue))

# Convert Lupton ORAL female feces from long to wide format using reshape
experiment12 <- reshape(Lup_OR_Ffeces[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment12) <- c("Time",unique(Lup_OR_Ffeces$Tissue))
# Change original data with cumulative data
experiment12$Feces <- obs_Lup_OR_Ffeces_cum

# Convert Lupton ORAL female urine from long to wide format using reshape
experiment13 <- reshape(Lup_OR_Furine[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment13) <- c("Time",unique(Lup_OR_Furine$Tissue))
# Change original data with cumulative data
experiment13$Urine <- obs_Lup_OR_Furine_cum

#Convert Dzierlenga 2021, IV female serum from long to wide format using reshape
experiment22 <- reshape(dzi_IV_Fserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment22) <- c("Time",unique(dzi_IV_Fserum$Tissue))

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment23 <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment23) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL female serum medium from long to wide format using reshape
experiment24 <- reshape(dzi_OR_Fserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment24) <- c("Time",unique(dzi_OR_Fserum_medium$Tissue))

#Convert Dzierlenga 2021, ORAL female serum high from long to wide format using reshape
experiment25 <- reshape(dzi_OR_Fserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment25) <- c("Time",unique(dzi_OR_Fserum_high$Tissue))

#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment26 <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment26) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment27<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment27) <- c("Time",unique(kim_IV_Fblood$Tissue))



# Put the experiments in a list
experiments <- list(
                    experiment5 = experiment5, experiment6 = experiment6, experiment8 = experiment8,
                   experiment11 = experiment11, experiment12 = experiment12,
                    experiment13 = experiment13,  experiment22 = experiment22, experiment23 = experiment23, experiment24 = experiment24,
                    experiment25 = experiment25, experiment26 = experiment26, experiment27 = experiment27)


# Rename predictions so that they share the same name as the names of the experimental data dataframe


colnames(preds_kim_IV_Ftissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")
colnames(preds_kim_OR_Ftissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")

colnames(preds_dzi_OR_Ftissues) <- c("Time","Liver","Kidney","Brain")

colnames(preds_Lup_OR_Ftissues) <- c ("Time", "Liver","Kidney","Blood","Skin")
colnames(preds_Lup_OR_Ffeces) <- c ("Time", "Feces")
colnames(preds_Lup_OR_Furine) <- c ("Time", "Urine")

colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")


# Create a list containing the corresponding predictions
simulations <- list(
                    predictions5 = preds_kim_IV_Ftissues, predictions6 = preds_kim_OR_Ftissues,
                     predictions8 = preds_dzi_OR_Ftissues, predictions11 = preds_Lup_OR_Ftissues, predictions12 = preds_Lup_OR_Ffeces,
                    predictions13 = preds_Lup_OR_Furine, 
                    predictions22 =preds_dzi_IV_Fserum, predictions23 =preds_dzi_OR_Fserum_low, predictions24 =preds_dzi_OR_Fserum_medium,
                    predictions25 =preds_dzi_OR_Fserum_high, predictions26 = preds_kim_OR_Fblood, 
                    predictions27 = preds_kim_IV_Fblood)


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
  
  
  # Save the plot with dynamically adjusted dimensions
  ggsave(paste0("experiment", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}



