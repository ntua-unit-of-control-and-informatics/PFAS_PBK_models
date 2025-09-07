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

##########################
#-------------------------
# Lupton ORAL female tissues
#-------------------------
##########################
# Set up simulations for the 11th case, i.e. Lupton (2020) ORAL female tissues
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
                   
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,384,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

##############################
### Load experimental data ###
##############################

# Load the mean values 
df <- openxlsx::read.xlsx("Raw_Data/Lupton_2020/PFOA_female_tissues_Lupton_2020.xlsx")
colnames(df)[c(2,3)] <- c("time", "concentration")
preds_Lup_OR_Ftissues <- solution[round(solution$time,1) %in% unique(round(df$time,1)), c("Cliver","Ckidney", 
                                                                              "Cblood", "Cskin")]/1000



obs_Lup_OR_Ftissues <- list(df[df$Tissue == "Liver", "concentration"],
                            df[df$Tissue == "Kidney", "concentration"],
                            df[df$Tissue == "Blood", "concentration"],
                            df[df$Tissue == "Skin", "concentration"])

# Gather the required data for the x-y plot
results_df<- data.frame("Study" = "Lupton_2020", "Dose" =  df$Dose_mg_per_kg,
                          "Tissue" = df$Tissue ,
                          "Type" = "oral",
                          "Observed" =df$concentration,
                          "Predicted" = unname(t(preds_Lup_OR_Ftissues)), "Time" = df$time )




AAFE_Lupton <-  AAFE(preds_Lup_OR_Ftissues,obs_Lup_OR_Ftissues)
print(paste0("The AAFE on the tissue data of Lupton et al. (2020) was ", AAFE_Lupton))

write.csv(results_df,
          "/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Publication/Final_results_plots/Validation/Lupton_2020_tissue_results.csv",
          row.names =F)
