
# Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
BW <- 0.3  # body weight (kg) not reported
admin.dose_per_g <- 12 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "M"
estimated_params<- rep(1,10)
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=c(0,0.001,0.01,0.1, 0.5, 1,2, seq(4,864,2))
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Mtissues <-  solution[,  c("time", "Cliver","Ckidney", "Cbrain")]


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


sample_time <- c(0,0.001,0.01,0.1, 0.5, 1,2, seq(4, 96, 4))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]

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


sample_time=c(0,0.001,0.01,0.1,0.2,0.5,0.75,seq(1,24,1))
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]


# Set up simulations for the 21st case, i.e. Dzierlenga 2021, ORAL male serum low
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 6 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "M"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.01, 0.25, 1, 3, 6, seq(12, 1200, 4))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Mserum_low <-  solution[, c("time", "Cplasma")]

# Set up simulations for the 25th case, i.e. Dzierlenga 2021, ORAL female serum low
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


sample_time <- c(0,0.001,0.01,0.1,0.2,0.5,0.75,1,seq(2, 96, 4))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]

# Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL female urine LOW

sex <- "F"
BW <- 0.2 #kg
sample_time <- seq(0,192,1)
admin.type <-"oral"
admin.time <- 0

#Female, oral 1mg/kg dose
admin.dose <- 1 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
inits <- create.inits (parameters)
solution_F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))


preds_Kemp_OR_Furine_low <-  solution_F[, c("time", "Murine")]
preds_Kemp_OR_Ffeces_low <-  solution_F[, c("time", "Mfeces")]



# Set up simulations for the 17th case, i.e.Kemper 2003 (Worley) ORAL male urine LOW

sex <- "M"
BW <- 0.3 #kg
sample_time <- seq(0,673,1)
admin.type <-"oral"
admin.time <- 0

#Male, oral 1mg/kg dose
admin.dose <- 1 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution_M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Murine_low <-  solution_M[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_low <-  solution_M[, c("time", "Mfeces")]


preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] <- preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] /1000
preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] <- preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] <- preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] /1000
preds_Kemp_OR_Murine_low[,2:dim(preds_Kemp_OR_Murine_low)[2]] <- preds_Kemp_OR_Murine_low[,2:dim(preds_Kemp_OR_Murine_low)[2]] /1000
preds_Kemp_OR_Ffeces_low[,2:dim(preds_Kemp_OR_Ffeces_low)[2]] <- preds_Kemp_OR_Ffeces_low[,2:dim(preds_Kemp_OR_Ffeces_low)[2]] /1000
preds_Kemp_OR_Mfeces_low[,2:dim(preds_Kemp_OR_Mfeces_low)[2]] <- preds_Kemp_OR_Mfeces_low[,2:dim(preds_Kemp_OR_Mfeces_low)[2]] /1000


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



# Convert Dzierlenga ORAL male tissues from long to wide format using reshape
experiment7 <- reshape(dzi_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(dzi_OR_Mtissues$Tissue))
experiment7$Urine <-  rep(0,dim(experiment7)[1])
experiment7$Feces <- rep(0,dim(experiment7)[1])


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))
experiment8$Urine <- rep(0,dim(experiment8)[1])
experiment8$Feces <-  rep(0,dim(experiment8)[1])



# Convert Kemper ORAL female urine low from long to wide format using reshape
experiment14 <- reshape(Kemp_OR_Furine_low [c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment14) <- c("Time",unique(Kemp_OR_Furine_low$Tissue))
experiment14$Urine = (experiment14$Urine/100)*0.2*1



# Convert Kemper ORAL male urine low from long to wide format using reshape
experiment17 <- reshape(Kemp_OR_Murine_low [c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment17) <- c("Time",unique(Kemp_OR_Murine_low$Tissue))
experiment17$Urine = (experiment17$Urine/100)*0.3*1

# Convert Dzierlenga 2021, ORAL male serum low from long to wide format using reshape
experiment21 <- reshape(dzi_OR_Mserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment21) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))
experiment21$Urine <- rep(0,dim(experiment21)[1])
experiment21$Feces <-  rep(0,dim(experiment21)[1])

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment25 <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment25) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))
experiment25$Urine <- rep(0,dim(experiment25)[1])
experiment25$Feces <-  rep(0,dim(experiment25)[1])


# Convert Kemper ORAL female feces low from long to wide format using reshape
experiment34 <- reshape(Kemp_OR_Ffeces_low[c("Tissue" ,"Time_h", 
                                             "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment34) <- c("Time",unique(Kemp_OR_Ffeces_low$Tissue))
experiment34$Feces = (experiment34$Feces/100)*0.2*25


# Convert Kemper ORAL male feces low from long to wide format using reshape
experiment35 <- reshape(Kemp_OR_Mfeces_low[c("Tissue" ,"Time_h", 
                                             "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment35) <- c("Time",unique(Kemp_OR_Mfeces_low$Tissue))
experiment35$Feces = (experiment35$Feces/100)*0.3*25


# Put the experiments in a list
experiments <- list(experiment7 = experiment7, experiment8 = experiment8,
                    experiment14 = experiment14, experiment17 = experiment17,
                    experiment21 = experiment21, experiment25 = experiment25,
                    experiment34 = experiment34, experiment35 = experiment35)


colnames(preds_dzi_OR_Mtissues) <- c("Time","Liver","Kidney","Brain" )
colnames(preds_dzi_OR_Ftissues) <- colnames(preds_dzi_OR_Mtissues)

colnames(preds_Kemp_OR_Furine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_low) <- c ("Time", "Urine")

colnames(preds_dzi_OR_Mserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")

colnames(preds_Kemp_OR_Ffeces_low) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_low) <- c ("Time", "Feces")

# Create a list containing the corresponding predictions
simulations <- list(predictions7 = preds_dzi_OR_Mtissues, predictions8 = preds_dzi_OR_Ftissues,
                    predictions14 = preds_Kemp_OR_Furine_low, predictions17 = preds_Kemp_OR_Murine_low,
                    predictions21 = preds_dzi_OR_Mserum_low, predictions25 = preds_dzi_OR_Fserum_low,
                    predictions34 = preds_Kemp_OR_Ffeces_low, predictions35 = preds_Kemp_OR_Mfeces_low)


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



