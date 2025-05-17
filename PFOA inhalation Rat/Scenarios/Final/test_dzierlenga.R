
###############################################################################################
# Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[7]])
params <- c(fixed_params[[7]], variable_params)
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


# Convert Dzierlenga ORAL male tissues from long to wide format using reshape
experiment7 <- list()
experiment7[[1]] <- reshape(dzi_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microM")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment7[[1]]) <- c("Time",unique(dzi_OR_Mtissues$Tissue))


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment7[[2]] <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microM")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment7[[2]]) <- c("Time",unique(dzi_OR_Ftissues$Tissue))

# Put the experiments in a list
experiments <- list(experiment7 = experiment7)


colnames(preds_dzi_OR_Mtissues) <- c("Time","Liver","Kidney","Brain")
colnames(preds_dzi_OR_Ftissues) <- colnames(preds_dzi_OR_Mtissues)

# Create a list containing the corresponding predictions
simulations <- list(
                    predictions7 = list(preds_dzi_OR_Mtissues =  preds_dzi_OR_Mtissues,
                                        preds_dzi_OR_Ftissues = preds_dzi_OR_Ftissues)) 
                   


plot_names <- list(predictions7 = list(preds_dzi_OR_Mtissues =   c("Liver","Kidney","Brain"), 
                              preds_dzi_OR_Ftissues = c("Liver","Kidney","Brain")))



library(ggplot2)
library(patchwork)

create.plots <- function(predictions, observations, compartment, plot_name) {
  cls <- c("Prediction" = "#0072B2", "Observation" = "#D55E00")
  
  if(plot_name %in% c( "ELF",  "Brain",  "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                       "Heart", "Brain", "Gonads", "Stomach", "Intestine") ){
    y_name <-  expression("PFOA concentration (" * mu * "g/g tissue)")
  }else if(plot_name %in% c( "Urine Female 1 mg/kg", "Feces Female 1 mg/kg",  "Urine Female 5 mg/kg", 
                             "Feces Female 5 mg/kg", "Urine Female 25 mg/kg", "Feces Female 25 mg/kg",
                             "Urine Male 1 mg/kg","Feces Male 1 mg/kg",  "Urine Male 5 mg/kg", 
                             "Feces Male 5 mg/kg","Urine Male 25 mg/kg", "Feces Male 25 mg/kg")){
    y_name <- expression("PFOA Mass (" * mu * "g)")
  }else{
    y_name <- expression("PFOA concentration (" * mu * "g/mL)")
  }
  ggplot(data = predictions) +
    geom_line(aes_string(x = "Time", y = rlang::expr(!!compartment), color = '"Prediction"'), 
              size = 1.2, alpha = 0.9) +
    geom_point(data = observations, 
               aes_string(x = "Time", y = rlang::expr(!!compartment), color = '"Observation"'), 
               size = 3.5, shape = 21, stroke = 1, fill = "white") +
    labs(title = plot_name, 
         y = y_name,
         x = "Time (hours)",
         color = "Type") +
    scale_color_manual(values = cls) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
}

# Iterate over experiments
for (i in seq_along(experiments)) {
  all_plots <- list()
  
  for (j in seq_along(experiments[[i]])) {
    observations <- experiments[[i]][[j]]
    predictions <- simulations[[i]][[j]]
    compartments <- names(predictions)[-1]  # Skip "Time"
    
    compartment_plots <- lapply(compartments, function(compartment) {
      create.plots(
        predictions = predictions,
        observations = observations,
        compartment = compartment,
        plot_name = compartment
      )
    })
    
    all_plots <- c(all_plots, compartment_plots)  # Flatten here
  }
  
  # Combine with patchwork
  if (length(all_plots) == 1) {
    final_plot <- all_plots[[1]] + theme(legend.position = "bottom")
  } else {
    final_plot <- wrap_plots(all_plots, ncol = 3) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
  }
  
print(final_plot)
}