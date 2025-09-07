library(deSolve)
setwd("/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Publication/Final_results_plots/Validation")

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
load("PBK_validation.RData")
# Body weight 
BW <- 0.21  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "M"
BFn = 166*60# 1/h
TVn = 1.71*1e-3# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 1 # administered dose in mg/m^3
depfr_head <- 0.3057
depfr_AF <- (0.1195+0.0243)
k = 1*duration
dose <- admin.dose_mg_per_m3*duration*BFn*TVn/k#ug PFOA, hour dose of the 6h inhalation
admin.time <- c(seq(0,5/24,1/24), seq(1,1+5/24,1/24), seq(2,2+5/24,1/24), seq(3,3+5/24,1/24), seq(4,4+5/24,1/24),
              seq(7,7+5/24,1/24), seq(8,8+5/24,1/24),  seq(9,9+5/24,1/24), seq(10,10+5/24,1/24), seq(11,11+5/24,1/24),
              seq(14,14+5/24,1/24), seq(15,15+5/24,1/24),  seq(16,16+5/24,1/24), seq(17,17+5/24,1/24), seq(18,18+5/24,1/24))*24
admin.dose <- rep(dose, length(admin.time))
admin.type <- "nasal"

user_input <-  list('BW'=BW,"admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                    "depfr_head" = depfr_head, "depfr_AF" = depfr_AF,
                                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0,22*24,2)
solution_1M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-2, atol = 1e-2))




#10mg/m^3 dose male
admin.dose_mg_per_m3 <- 10 # administered dose in mg/m^3
dose <- admin.dose_mg_per_m3*duration*BFn*TVn/k#ug PFOA, hour dose of the 6h inhalation
admin.dose <- rep(dose, length(admin.time))
user_input <-  list('BW'=BW,"admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                    "depfr_head" = depfr_head, "depfr_AF" = depfr_AF,
                                   "sex" = sex)
params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)



solution_10M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                       y = inits, parms = params, events = events,
                                       method="lsodes",rtol = 1e-2, atol = 1e-2))



#25mg/m^3 dose male
admin.dose_mg_per_m3 <- 25 # administered dose in mg/m^3
dose <- admin.dose_mg_per_m3*duration*BFn*TVn/k#ug PFOA, hour dose of the 6h inhalation
admin.dose <- rep(dose, length(admin.time))
user_input <-  list('BW'=BW,"admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "depfr_head" = depfr_head, "depfr_AF" = depfr_AF,
                    "sex" = sex)
params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

solution_25M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-2, atol = 1e-2))



#1mg/m^3 dose female
BW <- 0.15  #kg, not reported in the study - 180-240 g average BW of female CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "F"
BFn = 166*60# 1/h
TVn = 1.05*1e-3# L
admin.dose_mg_per_m3 <- 1 # administered dose in mg/m^3
dose <- admin.dose_mg_per_m3*duration*BFn*TVn/k#ug PFOA, hour dose of the 6h inhalation
admin.dose <- rep(dose, length(admin.time))
user_input <-  list('BW'=BW,"admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "depfr_head" = depfr_head, "depfr_AF" = depfr_AF,
                    "sex" = sex)
params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


solution_1F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-2, atol = 1e-2))




#10mg/m^3 dose male
admin.dose_mg_per_m3 <- 10 # administered dose in mg/m^3
dose <- admin.dose_mg_per_m3*duration*BFn*TVn/k#ug PFOA, hour dose of the 6h inhalation
admin.dose <- rep(dose, length(admin.time))
user_input <-  list('BW'=BW,"admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "depfr_head" = depfr_head, "depfr_AF" = depfr_AF,
                    "sex" = sex)
params <- create.params(user_input)

inits <- create.inits(params)
events <- create.events(params)
solution_10F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-2, atol = 1e-2))



#25mg/m^3 dose male
admin.dose_mg_per_m3 <- 25 # administered dose in mg/m^3
dose <- admin.dose_mg_per_m3*duration*BFn*TVn/k#ug PFOA, hour dose of the 6h inhalation
admin.dose <- rep(dose, length(admin.time))
user_input <-  list('BW'=BW,"admin.dose"= admin.dose,
                    "admin.time" = admin.time, 
                    "admin.type" = admin.type,
                    "depfr_head" = depfr_head, "depfr_AF" = depfr_AF,
                    "sex" = sex)
params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

solution_25F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-2, atol = 1e-2))



##############################
### Load experimental data ###
##############################

# Custom rounding function
custom_round <- function(time) {
  if (time < 1) {
    return(round(time, 1))
  } else {
    return(round(time, 0))
  }
}

# Load the mean values 

score <- rep(NA,6)

# Gather the required data for the x-y plot
df <- openxlsx::read.xlsx("Raw_Data/Hinderliter_2006/All_data_Hinderliter_2006_repeated.xlsx")
obs_plasma_1M <- df[df$Dose_mg_per_m3 == 1 & df$Sex == "M",]
rounded_time <- sapply(obs_plasma_1M$Time_h, custom_round)
rounded_soltime <- sapply(solution_1M$time, custom_round)
preds_plasma_1M <- unique(solution_1M[rounded_soltime %in% rounded_time, "Cplasma"]/1000)
score[1] <- AAFE(preds_plasma_1M,obs_plasma_1M$Concentration_microg_per_g_organ)

results_df_1M<- data.frame("Study" = "Hinderliter_2006", "Dose" =  obs_plasma_1M$Dose_mg_per_m3,
                          "Tissue" = obs_plasma_1M$Tissue ,
                          "Type" = obs_plasma_1M$Type, sex = obs_plasma_1M$Sex,
                          "Observed" =obs_plasma_1M$Concentration_microg_per_g_organ,
                          "Predicted" = preds_plasma_1M, "Time" = obs_plasma_1M$Time_h, time = rounded_time)



df <- openxlsx::read.xlsx("Raw_Data/Hinderliter_2006/All_data_Hinderliter_2006_repeated.xlsx")
obs_plasma_10M <- df[df$Dose_mg_per_m3 == 10 & df$Sex == "M",]
rounded_time <- sapply(obs_plasma_10M$Time_h, custom_round)
rounded_soltime <- sapply(solution_10M$time, custom_round)
preds_plasma_10M <- unique(solution_10M[rounded_soltime %in% rounded_time, "Cplasma"]/1000)
score[2] <- AAFE(preds_plasma_10M,obs_plasma_10M$Concentration_microg_per_g_organ)

results_df_10M<- data.frame("Study" = "Hinderliter_2006", "Dose" =  obs_plasma_10M$Dose_mg_per_m3,
                           "Tissue" = obs_plasma_10M$Tissue ,
                           "Type" = obs_plasma_10M$Type, sex = obs_plasma_10M$Sex,
                           "Observed" =obs_plasma_10M$Concentration_microg_per_g_organ,
                           "Predicted" = preds_plasma_10M, "Time" = obs_plasma_10M$Time_h, time = rounded_time )




df <- openxlsx::read.xlsx("Raw_Data/Hinderliter_2006/All_data_Hinderliter_2006_repeated.xlsx")
obs_plasma_25M <- df[df$Dose_mg_per_m3 == 25 & df$Sex == "M",]
rounded_time <- sapply(obs_plasma_25M$Time_h, custom_round)
rounded_soltime <- sapply(solution_25M$time, custom_round)
preds_plasma_25M <- unique(solution_25M[rounded_soltime %in% rounded_time, "Cplasma"]/1000)
score[3] <- AAFE(preds_plasma_25M,obs_plasma_25M$Concentration_microg_per_g_organ)

results_df_25M<- data.frame("Study" = "Hinderliter_2006", "Dose" =  obs_plasma_25M$Dose_mg_per_m3,
                            "Tissue" = obs_plasma_25M$Tissue ,
                            "Type" = obs_plasma_25M$Type, sex = obs_plasma_25M$Sex,
                            "Observed" =obs_plasma_25M$Concentration_microg_per_g_organ,
                            "Predicted" = preds_plasma_25M, "Time" = obs_plasma_25M$Time_h , time = rounded_time)





df <- openxlsx::read.xlsx("Raw_Data/Hinderliter_2006/All_data_Hinderliter_2006_repeated.xlsx")
obs_plasma_1F <- df[df$Dose_mg_per_m3 == 1 & df$Sex == "F",]
rounded_time <- sapply(obs_plasma_1F$Time_h, custom_round)
rounded_soltime <- sapply(solution_1F$time, custom_round)
preds_plasma_1F <- unique(solution_1F[rounded_soltime %in% rounded_time, "Cplasma"]/1000)
score[4] <- AAFE(preds_plasma_1F,obs_plasma_1F$Concentration_microg_per_g_organ)

results_df_1F<- data.frame("Study" = "Hinderliter_2006", "Dose" =  obs_plasma_1F$Dose_mg_per_m3,
                           "Tissue" = obs_plasma_1F$Tissue ,
                           "Type" = obs_plasma_1F$Type, sex = obs_plasma_1F$Sex,
                           "Observed" =obs_plasma_1F$Concentration_microg_per_g_organ,
                           "Predicted" = preds_plasma_1F, "Time" = obs_plasma_1F$Time_h, time = rounded_time )


df <- openxlsx::read.xlsx("Raw_Data/Hinderliter_2006/All_data_Hinderliter_2006_repeated.xlsx")
obs_plasma_10F <- df[df$Dose_mg_per_m3 == 10 & df$Sex == "F",]
rounded_time <- sapply(obs_plasma_10F$Time_h, custom_round)
rounded_soltime <- sapply(solution_10F$time, custom_round)
preds_plasma_10F <- unique(solution_10F[rounded_soltime %in% rounded_time, "Cplasma"]/1000)
score[5] <- AAFE(preds_plasma_10F,obs_plasma_10F$Concentration_microg_per_g_organ)

results_df_10F<- data.frame("Study" = "Hinderliter_2006", "Dose" =  obs_plasma_10F$Dose_mg_per_m3,
                            "Tissue" = obs_plasma_10F$Tissue ,
                            "Type" = obs_plasma_10F$Type, sex = obs_plasma_10F$Sex,
                            "Observed" =obs_plasma_10F$Concentration_microg_per_g_organ,
                            "Predicted" = preds_plasma_10F, "Time" = obs_plasma_10F$Time_h, time = rounded_time )



df <- openxlsx::read.xlsx("Raw_Data/Hinderliter_2006/All_data_Hinderliter_2006_repeated.xlsx")
obs_plasma_25F <- df[df$Dose_mg_per_m3 == 25 & df$Sex == "F",]
rounded_time <- sapply(obs_plasma_25F$Time_h, custom_round)
rounded_soltime <- sapply(solution_25F$time, custom_round)
preds_plasma_25F <- unique(solution_25F[rounded_soltime %in% rounded_time, "Cplasma"]/1000)
score[6] <- AAFE(preds_plasma_25F,obs_plasma_25F$Concentration_microg_per_g_organ)

results_df_25F<- data.frame("Study" = "Hinderliter_2006", "Dose" =  obs_plasma_25F$Dose_mg_per_m3,
                            "Tissue" = obs_plasma_25F$Tissue ,
                            "Type" = obs_plasma_25F$Type, sex = obs_plasma_25F$Sex,
                            "Observed" =obs_plasma_25F$Concentration_microg_per_g_organ,
                            "Predicted" = preds_plasma_25F, "Time" = obs_plasma_25F$Time_h, time = rounded_time )


results_df <- rbind(results_df_1M, results_df_10M, results_df_25M, results_df_1F, results_df_10F, results_df_25F)

AAFE_Hinderliter <- mean(score)
print(paste0("The AAFE on the Plasma data of Hinderliter et al. (2006) was ", AAFE_Hinderliter))

write.csv(results_df,
          '/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Publication/Final_results_plots/Validation/Hinderliter_2006_results.csv',
          row.names =F)


