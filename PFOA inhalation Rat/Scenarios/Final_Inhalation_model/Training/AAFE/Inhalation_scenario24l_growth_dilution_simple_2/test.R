##########################
#-------------------------
# Hinderliter Inhalation female single high
#-------------------------
##########################
# Set up simulations for the 15th case, i.e. Hinderliter Inhalation female single high
BW <- 0.13  #kg, not reported in the study - 180-240 g average BW of female CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "F"
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 27 # administered dose in mg/m^3
depfr_head <- 0.3372
depfr_AF <- (0.1327+0.0177)
k = 1*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user.input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )


params <- create_params(user.input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time= sort(c(0.5,union(events$data$time, seq(0,24,1))))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df15=========================================================

exp_data <- dataset$df15 # retrieve data of Hinderliter Inhalation male single high
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
column_names <- c("Cplasma")

preds_hind_INH_Fblood_high <- list()
# loop over compartments with available data
for (i in 1:length(unique(exp_data$Tissue))) {
  compartment <- unique(exp_data$Tissue)[i]
  #Retrieve time points at which measurements are available for compartment i
  exp_time <- exp_data[exp_data$Tissue == compartment, 2]
  
  preds_hind_INH_Fblood_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
}

obs_hind_INH_Fblood_high <- list (exp_data[exp_data$Tissue == "Plasma", "concentration"])



print(AAFE(predictions = preds_hind_INH_Fblood_high, observations = obs_hind_INH_Fblood_high)) 

df <- data.frame(predictions = solution$Cplasma/1000, solution_time = solution$time)
df2 <- data.frame(observations = obs_hind_INH_Fblood_high[[1]], observation_time = exp_time)

library(ggplot2)
ggplot(data = df)+
  geom_line(aes(x =solution_time, y= predictions ))+
  geom_point(data = df2, aes(x =exp_time, y= observations ))
