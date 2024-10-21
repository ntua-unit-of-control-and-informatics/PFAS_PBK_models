library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/Ioannis/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/full_params_PT_excreta_restricted_Ka_58")

#  absolute average fold error
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N <- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}
#===============
# Generate predictions
#===============


load("full_params_PT_excreta_restricted_Ka_58.RData")
# Body weight 


# Set up simulations for the 1st case, i.e. Cui (2010) ORAL male urine low
BW <-  0.200 #kg, before quarantine
# Each day the male rat body weight increases by 5.9 g
for (i in 1:7){
  BW <- BW + 5.9/1000
}
BW_init <- BW #Body weight at the beginning of the experiment
#Initialize a vector of daily body weight
BW <- c(BW_init, rep(NA, 27))
# Estimate the BW each day
for (i in 2:length(BW)){
  BW[i] <- BW[i-1] + 5.9/1000
} 
admin.dose_per_BW <- 5 # administered dose in mg PFOA/kg BW 
admin.time <- seq(0,27*24,24) #time when doses are administered, in hours
admin.dose <- admin.dose_per_BW*BW*1e03#ug PFOA
admin.type <- "oral"
sex <- "M" 


user_input <- list('BW'=BW_init,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,672,2)
solution_M_oral_5_urine <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-07, atol = 1e-07))


# Set up simulations for the 2nd case, i.e. Cui (2010) ORAL male urine high
BW <-  0.210 #kg, before quarantine
# Each day the male rat body weight increases by 5.9 g
for (i in 1:7){
  BW <- BW + 5.9/1000
}
BW_init <- BW #Body weight at the beginning of the experiment
#Initialize a vector of daily body weight
BW <- c(BW_init, rep(NA, 27))
# Estimate the BW each day
for (i in 2:length(BW)){
  BW[i] <- BW[i-1] + 5.9/1000
} 
admin.dose_per_BW <- 20 # administered dose in mg PFOA/kg BW 
admin.time <- seq(0,27*24,24) #time when doses are administered, in hours
admin.dose <- admin.dose_per_BW*BW*1e03#ug PFOA
admin.type <- "oral"
sex <- "M" 


user_input <- list('BW'=BW_init,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,672,2)
solution_M_oral_20_urine <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))



# Set up simulations for the 3d case, i.e. Cui (2010) ORAL male feces low
BW <-  0.200 #kg, before quarantine
# Each day the male rat body weight increases by 5.9 g
for (i in 1:7){
  BW <- BW + 5.9/1000
}
BW_init <- BW #Body weight at the beginning of the experiment
#Initialize a vector of daily body weight
BW <- c(BW_init, rep(NA, 27))
# Estimate the BW each day
for (i in 2:length(BW)){
  BW[i] <- BW[i-1] + 5.9/1000
} 
admin.dose_per_BW <- 5 # administered dose in mg PFOA/kg BW 
admin.time <- seq(0,27*24,24) #time when doses are administered, in hours
admin.dose <- admin.dose_per_BW*BW*1e03#ug PFOA
admin.type <- "oral"
sex <- "M" 


user_input <- list('BW'=BW_init,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,672,2)
solution_M_oral_5_feces <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))



# Set up simulations for the 4th case, i.e. Cui (2010) ORAL male feces high
BW <-  0.210 #kg, before quarantine
# Each day the male rat body weight increases by 5.9 g
for (i in 1:7){
  BW <- BW + 5.9/1000
}
BW_init <- BW #Body weight at the beginning of the experiment
#Initialize a vector of daily body weight
BW <- c(BW_init, rep(NA, 27))
# Estimate the BW each day
for (i in 2:length(BW)){
  BW[i] <- BW[i-1] + 5.9/1000
} 
admin.dose_per_BW <- 20 # administered dose in mg PFOA/kg BW 
admin.time <- seq(0,27*24,24) #time when doses are administered, in hours
admin.dose <- admin.dose_per_BW*BW*1e03#ug PFOA
admin.type <- "oral"
sex <- "M" 


user_input <- list('BW'=BW_init,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,672,2)
solution_M_oral_20_feces <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))




##############################
### Load experimental data ###
##############################

score <- rep(NA, 4)
df_excreta <- openxlsx::read.xlsx("Raw_Data/Cui_2008/Cui_excreta.xlsx")
obs_M_oral_5_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 5 & 
                                     df_excreta$Tissue == "Urine" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_oral_5_urine$Time_day*24)
rounded_soltime <- round(solution_M_oral_5_urine$time)
preds_M_oral_5_urine <- solution_M_oral_5_urine[rounded_soltime %in% rounded_time, "Murine"]/1000
score[1] <- AAFE(preds_M_oral_5_urine, (obs_M_oral_5_urine$Mass_mg))
results_urine_M_oral_5<- data.frame("Study" = "Cui_2010_excreta", "Dose" =  obs_M_oral_5_urine$Dose_mg_per_kg,
                          "Tissue" = obs_M_oral_5_urine$Tissue ,
                          "Type" = obs_M_oral_5_urine$Type, "sex" = obs_M_oral_5_urine$Sex,
                          "Observed" = (obs_M_oral_5_urine$Mass_mg),
                          "Predicted" = preds_M_oral_5_urine,  "Time" = obs_M_oral_5_urine$Time_day*24)


obs_M_oral_20_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 20 & 
                                   df_excreta$Tissue == "Urine" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_oral_20_urine$Time_day*24)
rounded_soltime <- round(solution_M_oral_20_urine$time)
preds_M_oral_20_urine <- solution_M_oral_20_urine[rounded_soltime %in% rounded_time, "Murine"]/1000
score[2] <- AAFE(preds_M_oral_20_urine, (obs_M_oral_20_urine$Mass_mg))
results_urine_M_oral_20<- data.frame("Study" = "Cui_2010_excreta", "Dose" =  obs_M_oral_20_urine$Dose_mg_per_kg,
                                    "Tissue" = obs_M_oral_20_urine$Tissue ,
                                    "Type" = obs_M_oral_20_urine$Type, "sex" = obs_M_oral_20_urine$Sex,
                                    "Observed" = (obs_M_oral_20_urine$Mass_mg),
                                    "Predicted" = preds_M_oral_20_urine,  "Time" = obs_M_oral_20_urine$Time_day*24)



obs_M_oral_5_feces <- df_excreta[df_excreta$Dose_mg_per_kg == 5 & 
                                   df_excreta$Tissue == "Feces" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_oral_5_feces$Time_day*24)
rounded_soltime <- round(solution_M_oral_5_feces$time)
preds_M_oral_5_feces <- solution_M_oral_5_feces[rounded_soltime %in% rounded_time, "Mfeces"]/1000
score[3] <- AAFE(preds_M_oral_5_feces, (obs_M_oral_5_feces$Mass_mg))
results_feces_M_oral_5<- data.frame("Study" = "Cui_2010_excreta", "Dose" =  obs_M_oral_5_feces$Dose_mg_per_kg,
                                    "Tissue" = obs_M_oral_5_feces$Tissue ,
                                    "Type" = obs_M_oral_5_feces$Type, "sex" = obs_M_oral_5_feces$Sex,
                                    "Observed" = (obs_M_oral_5_feces$Mass_mg),
                                    "Predicted" = preds_M_oral_5_feces,  "Time" = obs_M_oral_5_feces$Time_day*24)


obs_M_oral_20_feces <- df_excreta[df_excreta$Dose_mg_per_kg == 20 & 
                                   df_excreta$Tissue == "Feces" & df_excreta$Sex == "M",]
rounded_time <- round(obs_M_oral_20_feces$Time_day*24)
rounded_soltime <- round(solution_M_oral_20_feces$time)
preds_M_oral_20_feces <- solution_M_oral_20_feces[rounded_soltime %in% rounded_time, "Mfeces"]/1000
score[4] <- AAFE(preds_M_oral_20_feces, (obs_M_oral_20_feces$Mass_mg))
results_feces_M_oral_20<- data.frame("Study" = "Cui_2010_excreta", "Dose" =  obs_M_oral_20_feces$Dose_mg_per_kg,
                                    "Tissue" = obs_M_oral_20_feces$Tissue ,
                                    "Type" = obs_M_oral_20_feces$Type, "sex" = obs_M_oral_20_feces$Sex,
                                    "Observed" = (obs_M_oral_20_feces$Mass_mg),
                                    "Predicted" = preds_M_oral_20_feces,  "Time" = obs_M_oral_20_feces$Time_day*24)

results_df <- rbind(results_urine_M_oral_5, results_urine_M_oral_20, results_feces_M_oral_5, results_feces_M_oral_20)


AAFE_Cui <- mean(score)
print(paste0("The AAFE on the excreta data of Cui et al. (2010) is ", AAFE_Cui))


write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Kemper_2003_excreta_Loccisano_results.csv",
          "C:/Users/Ioannis/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/full_params_PT_excreta_restricted_Ka_58/Cui_2010_Excreta_results.csv",
          row.names =F)
