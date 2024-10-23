library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/Ioannis/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/full_params_PT_excreta_restricted_Ka_58_no_papp")

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


load("full_params_PT_excreta_restricted_Ka_58_no_papp.RData")

# Set up simulations for the 1st case, i.e. Lupton (2020) ORAL female feces
BW <- 0.184  # body weight (kg) 
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

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,384,1)

# ode(): The solver of the ODEs
solution_F_oral_feces <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))



# Set up simulations for the 1st case, i.e. Lupton (2020) ORAL female urine
BW <- 0.184  # body weight (kg) 
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

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,384,1)

# ode(): The solver of the ODEs
solution_F_oral_urine <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                                 y = inits, parms = params,
                                                 events = events,
                                                 method="lsodes",rtol = 1e-05, atol = 1e-05))




##############################
### Load experimental data ###
##############################

score <- rep(NA, 2)
df_excreta <- openxlsx::read.xlsx("Raw_Data/Lupton_2020/Lupton_2020_excreta.xlsx")
obs_F_oral_feces <- df_excreta[df_excreta$Dose_mg_per_kg == 0.047 & 
                                     df_excreta$Tissue == "Feces" & df_excreta$Sex == "F",]

#Estimate fecal dry mass 
Mfeces_wet <- (8.18/0.21)*BW #g
Mfeces_dry <- Mfeces_wet*0.8 # Enqi et al., 2021, control rats:20.1% water content, https://doi.org/10.3389/fcimb.2020.581974


Lupton_Ffeces <- c(df_excreta[df_excreta$Tissue == "Feces", "Concentration_microg_per_g_organ"])*Mfeces_dry
# Estimate cumulative fecal mass
obs_Lup_OR_Ffeces_cum <- cumsum(Lupton_Ffeces)

rounded_time <- round(obs_F_oral_feces$Time_hours)
rounded_soltime <- round(solution_F_oral_feces$time)
preds_F_oral_feces <- solution_F_oral_feces[rounded_soltime %in% rounded_time, "Mfeces"]


score[1] <- AAFE(preds_F_oral_feces, obs_Lup_OR_Ffeces_cum)

results_feces_F_oral<- data.frame("Study" = "Lupton_2020", "Dose" =  obs_F_oral_feces$Dose_mg_per_kg,
                          "Tissue" = obs_F_oral_feces$Tissue ,
                          "Type" = obs_F_oral_feces$Type, sex = "F",
                          "Observed" = obs_Lup_OR_Ffeces_cum,
                          "Predicted" = preds_F_oral_feces,  "Time" = obs_F_oral_feces$Time_hours)




obs_F_oral_urine <- df_excreta[df_excreta$Dose_mg_per_kg == 0.047 & 
                                 df_excreta$Tissue == "Urine" & df_excreta$Sex == "F",]

Qurine_daily <- 85 * BW#  (ml/d/kg)*BW  --> mL/d, Schmidt et al., 2001, doi:10.1002/nau.1006
Lupton_Furine <- c(df_excreta[df_excreta$Tissue == "Urine", "Concentration_microg_per_g_organ"])*Qurine_daily
# Estimate cumulative fecal mass
obs_Lup_OR_Furine_cum <- cumsum(Lupton_Furine)


rounded_time <- round(obs_F_oral_urine$Time_hours)
rounded_soltime <- round(solution_F_oral_urine$time)
preds_F_oral_urine <- solution_F_oral_urine[rounded_soltime %in% rounded_time, "Murine"]


score[2] <- AAFE(preds_F_oral_urine, obs_Lup_OR_Furine_cum)

results_urine_F_oral<- data.frame("Study" = "Lupton_2020", "Dose" =  obs_F_oral_urine$Dose_mg_per_kg,
                                  "Tissue" = obs_F_oral_urine$Tissue,
                                  "Type" = obs_F_oral_urine$Type, sex = "F",
                                  "Observed" = obs_Lup_OR_Furine_cum,
                                  "Predicted" = preds_F_oral_urine,  "Time" = obs_F_oral_urine$Time_hours)




results_df <- rbind(results_feces_F_oral, results_urine_F_oral)


AAFE_Lupton <- mean(score)
print(paste0("The AAFE on the excreta data of Lupton et al. (2020) is ", AAFE_Lupton))


write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Kemper_2003_excreta_Loccisano_results.csv",
          "C:/Users/Ioannis/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/full_params_PT_excreta_restricted_Ka_58_no_papp/Lupton_2020_excreta_results.csv",
          row.names =F)
