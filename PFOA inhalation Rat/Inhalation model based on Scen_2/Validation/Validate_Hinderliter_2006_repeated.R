library(deSolve)
#setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Inhalation model based on Scen_2/Validation")

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
load("Inhalation model_Scen_2.RData")
# Body weight 
BW <- 0.25  #kg, not reported in the study
sex <- "M"
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 1 # administered dose in mg/m^3
depfr_head <- 0.3057
depfr_AF <- (0.1195+0.0243)
dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn))#ug PFOA, for 6h inhalation
admin.time <- sort(c(seq(0,4), seq(7,11), seq(14,18))*24) #time when doses are administered, in hours
admin.dose <- rep(dose, length(admin.time))
admin.type <- "nasal"

parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
sample_time <- seq(0,22*24,2)
solution_1M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

# Carefully select the order so that it matches the experimental data presented in the xlsx file
pred_comps <- c( "Cplasma" )
solution_1M <- solution_1M[solution_1M$time %in% sample_time, pred_comps] / 1000  #[ug/L]/1000-->[ug/g]


#10mg/m^3 dose male
admin.dose_mg_per_m3 <- 10 # administered dose in mg/m^3
dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn))#ug PFOA, for 6h inhalation
admin.dose <- rep(dose, length(admin.time))
parameters <- create.params( list('BW'=BW,
                                  "admin.dose"= admin.dose,
                                  "admin.time" = admin.time, 
                                  "admin.type" = admin.type,
                                  "estimated_params" = estimated_params,
                                  "sex" = sex))
events <- create.events(parameters)
solution_10M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                       y = inits, parms = parameters, events = events,
                                       method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_10M <- solution_10M[solution_10M$time %in% sample_time, pred_comps] / 1000 #[ug/L]/1000-->[ug/g]

#25mg/m^3 dose male
admin.dose_mg_per_m3 <- 25 # administered dose in mg/m^3
dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn))#ug PFOA, for 6h inhalation
admin.dose <- rep(dose, length(admin.time))
parameters <- create.params( list('BW'=BW,
                                  "admin.dose"= admin.dose,
                                  "admin.time" = admin.time, 
                                  "admin.type" = admin.type,
                                  "estimated_params" = estimated_params,
                                  "sex" = sex))
events <- create.events(parameters)
solution_25M <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_25M <- solution_25M[solution_25M$time %in% sample_time, pred_comps] / 1000 #[ug/L]/1000-->[ug/g]

#1mg/m^3 dose female
sex <- "F"
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# 
admin.dose_mg_per_m3 <- 1 # administered dose in mg/m^3
dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn))#ug PFOA, for 6h inhalation
admin.dose <- rep(dose, length(admin.time))
parameters <- create.params( list('BW'=BW,
                                  "admin.dose"= admin.dose,
                                  "admin.time" = admin.time, 
                                  "admin.type" = admin.type,
                                  "estimated_params" = estimated_params,
                                  "sex" = sex))
events <- create.events(parameters)
solution_1F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_1F <- solution_1F[solution_1F$time %in% sample_time, pred_comps] / 1000 #[ug/L]/1000-->[ug/g]


#10mg/m^3 dose male
admin.dose_mg_per_m3 <- 10 # administered dose in mg/m^3
dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn))#ug PFOA, for 6h inhalation
admin.dose <- rep(dose, length(admin.time))
parameters <- create.params( list('BW'=BW,
                                  "admin.dose"= admin.dose,
                                  "admin.time" = admin.time, 
                                  "admin.type" = admin.type,
                                  "estimated_params" = estimated_params,
                                  "sex" = sex))
events <- create.events(parameters)
solution_10F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_10F <- solution_10F[solution_10F$time %in% sample_time, pred_comps] / 1000 #[ug/L]/1000-->[ug/g]

#25mg/m^3 dose male
admin.dose_mg_per_m3 <- 25 # administered dose in mg/m^3
dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn))#ug PFOA, for 6h inhalation
admin.dose <- rep(dose, length(admin.time))
parameters <- create.params( list('BW'=BW,
                                  "admin.dose"= admin.dose,
                                  "admin.time" = admin.time, 
                                  "admin.type" = admin.type,
                                  "estimated_params" = estimated_params,
                                  "sex" = sex))
events <- create.events(parameters)
solution_25F <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = parameters, events = events,
                                      method="lsodes",rtol = 1e-7, atol = 1e-7))

solution_25F <- solution_25F[solution_25F$time %in% sample_time, pred_comps] / 1000 #[ug/L]/1000-->[ug/g]

##############################
### Load experimental data ###
##############################

# Load the mean values 
df <- openxlsx::read.xlsx("Raw_Data/All_data_Hinderliter_2006_repeated.xlsx")
plasma_1M <- df[df$Dose_mg_per_m3 == 1 & df$Sex == "M",]
plasma_10M <- df[df$Dose_mg_per_m3 == 10 & df$Sex == "M",]
plasma_25M <- df[df$Dose_mg_per_m3 == 25 & df$Sex == "M",]
plasma_1F <- df[df$Dose_mg_per_m3 == 1 & df$Sex == "F",]
plasma_10F <- df[df$Dose_mg_per_m3 == 10 & df$Sex == "F",]
plasma_25F <- df[df$Dose_mg_per_m3 == 25 & df$Sex == "F",]

score <- rep(NA,6)
# Gather the required data for the x-y plot
score[1] <- AAFE(solution_1M,plasma_1M$`Concentration_microg_per_g_organ`)

results_df_1M<- data.frame("Study" = "Hinderliter_2006", "Dose" =  plasma_1M$Dose_mg_per_m3,
                          "Tissue" = plasma_1M$Tissue ,
                          "Type" = plasma_1M$Type, sex = plasma_1M$Sex,
                          "Observed" =plasma_1M$`Concentration_microg_per_g_organ`,
                          "Predicted" = unname(t(solution_1M)), "Time" = plasma_1M$Time_h )


score[2] <- AAFE(solution_10M,plasma_1M$`Concentration_microg_per_g_organ`)

results_df_10M<- data.frame("Study" = "Hinderliter_2006", "Dose" =  plasma_10M$Dose_mg_per_m3,
                            "Tissue" = plasma_10M$Tissue ,
                            "Type" = plasma_10M$Type, sex = plasma_10M$Sex,
                            "Observed" =plasma_10M$`Concentration_microg_per_g_organ`,
                            "Predicted" = unname(t(solution_10M)), "Time" = plasma_10M$Time_h )

score[3] <- AAFE(solution_25M,plasma_25M$`Concentration_microg_per_g_organ`)

results_df_25M<- data.frame("Study" = "Hinderliter_2006", "Dose" =  plasma_25M$Dose_mg_per_m3,
                            "Tissue" = plasma_25M$Tissue ,
                            "Type" = plasma_25M$Type, sex = plasma_25M$Sex,
                            "Observed" =plasma_25M$`Concentration_microg_per_g_organ`,
                            "Predicted" = unname(t(solution_25M)), "Time" = plasma_25M$Time_h )


score[4] <- AAFE(solution_1F,plasma_1F$`Concentration_microg_per_g_organ`)

results_df_1F<- data.frame("Study" = "Hinderliter_2006", "Dose" =  plasma_1F$Dose_mg_per_m3,
                            "Tissue" = plasma_1F$Tissue ,
                            "Type" = plasma_1F$Type, sex = plasma_1F$Sex,
                            "Observed" =plasma_1F$`Concentration_microg_per_g_organ`,
                            "Predicted" = unname(t(solution_1F)), "Time" = plasma_1F$Time_h )


score[5] <- AAFE(solution_10F,plasma_10F$`Concentration_microg_per_g_organ`)

results_df_10F<- data.frame("Study" = "Hinderliter_2006", "Dose" =  plasma_10F$Dose_mg_per_m3,
                            "Tissue" = plasma_10F$Tissue ,
                            "Type" = plasma_10F$Type, sex = plasma_10F$Sex,
                            "Observed" =plasma_10F$`Concentration_microg_per_g_organ`,
                            "Predicted" = unname(t(solution_10F)), "Time" = plasma_10F$Time_h )


score[6] <- AAFE(solution_25F,plasma_25F$`Concentration_microg_per_g_organ`)

results_df_25F<- data.frame("Study" = "Hinderliter_2006", "Dose" =  plasma_25F$Dose_mg_per_m3,
                            "Tissue" = plasma_25F$Tissue ,
                            "Type" = plasma_25F$Type, sex = plasma_25F$Sex,
                            "Observed" =plasma_25F$`Concentration_microg_per_g_organ`,
                            "Predicted" = unname(t(solution_25F)), "Time" = plasma_25F$Time_h )

results_df <- rbind(results_df_1M, results_df_10M, results_df_25M, results_df_1F, results_df_10F, results_df_25F)

AAFE_Hinderliter <- mean(score)
print(paste0("The AAFE on the Plasma data of Hinderliter et al. (2006) was ", AAFE_Hinderliter))

write.csv(results_df,
          #"C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/Validation_results/Cui_2008_results.csv",
          "C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Inhalation model based on Scen_2/Validation_results/Hinderliter_2006_results.csv",
          row.names =F)
