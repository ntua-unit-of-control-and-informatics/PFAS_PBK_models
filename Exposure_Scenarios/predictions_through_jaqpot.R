# Set the model ID of the uploaded Jaqpot model
modelID <- "pYkcPk31g8OtDCF4ZJLd"
# Log in the Jaqpot server
#jaqpotr::login.cred()
jaqpotr::login.api()
#Get model features
jaqpotr::get.model.feats(modelID)
# Create a prediction instance

# PFOA ORAL exposure scenario 

daily_intake <- (18.92+0.06+0.056)*70 # ng of PFAS per day

df = data.frame("substance"= "PFOA",  
                "f_unabs"= 0.47, #range (0,1)
                "BW"= 70,  # kg
                "admin.dose"= rep(daily_intake, 40*365 ), # ng
                "admin.time"=seq(0, 40*365*24-1, 24), # hours
                "sim.start"=0, # hours
                "sim.end"=40*365*24-1, # hours
                "sim.step"=24)# hours

# Acquire model predictions
start.time <- Sys.time()
predictions <- jaqpotr::jaqpot.predict( df = df, modelID = modelID)
end.time <- Sys.time()
total_duration <- end.time - start.time
print(total_duration)

predictions$predictions
# Plot an output
plot(predictions$predictions$time,predictions$predictions$CLu)
