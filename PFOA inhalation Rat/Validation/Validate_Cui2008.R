library(deSolve)
setwd("/Users/ptsir/Documents/GitHub/PBK_Grouping/PFOA/Validation")

#===============
# Generate predictions
#===============
load("fitted_model.RData")
# Body weight 
BW <- 0.2 #kg
#5mg/kg dose
admin.time <- seq(0.01,27*24.01,24)
dose <- rep(5,28) #mg/kg
admin.dose <- dose * BW #mg
parameters <- create.params( list("BW" = BW  , sex = "M", 
                                  "admin.type" = "oral", "admin.time" = admin.time,
                                  "admin.dose" = admin.dose, "fitted_pars" = fitted_pars, 
                                  "group" = group, "N_pars" = N_pars ))
events <- create.events(parameters)
sample_time <- seq(0,28*24,0.1)
solution_5 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))


pred_comps <- c( "Cblood_art", "Cliver", "Ckidneys","Clung","Cheart",
                  "Cspleen", "Ctestis", "Cbrain" )
solution_5 <- solution_5[solution_5$time == 28*24,pred_comps]


#20mg/kg dose
dose <- rep(20,28) #mg/kg
admin.dose <- dose * BW #mg
parameters <- create.params( list("BW" = BW  , sex = "M", 
                                  "admin.type" = "oral", "admin.time" = admin.time,
                                  "admin.dose" = admin.dose, "fitted_pars" = fitted_pars, 
                                  "group" = group, "N_pars" = N_pars ))
events <- create.events(parameters)
solution_20 <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                       y = inits, parms = parameters, events = events,
                                       method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_20 <- solution_20[solution_20$time == 28*24,pred_comps]




##############################
### Load experimental data ###
##############################

# Load the mean values 
tissues <- openxlsx::read.xlsx("Raw_Data/Cui_2008/Cui_tissues.xlsx")
tissues_5 <- tissues[tissues$Dose_mg_per_kg == 5,]
tissues_20 <- tissues[tissues$Dose_mg_per_kg == 20,]

# Gather the required data for the x-y plot

results_df_5<- data.frame("Study" = "Cui_2010", "Dose" =  tissues_5$Dose_mg_per_kg,
                          "Tissue" = tissues_5$Tissue ,
                          "Type" = "oral",
                          "Observed" =tissues_5$`Concentration_ug/g`,
                          "Predicted" = unname(t(solution_5)))


results_df_20<- data.frame("Study" = "Cui_2010", "Dose" =  tissues_20$Dose_mg_per_kg,
                           "Tissue" = tissues_20$Tissue,
                           "Type" = "oral",
                           "Observed" = tissues_20$`Concentration_ug/g`,
                           "Predicted" = unname(t(solution_20)))


results_df <- rbind(results_df_5, results_df_20)

write.csv(results_df,
          "/Users/ptsir/Documents/GitHub/PBK_Grouping/PFOA/Validation/Validation_results/Cui_2008_results.csv",
          row.names =F)
