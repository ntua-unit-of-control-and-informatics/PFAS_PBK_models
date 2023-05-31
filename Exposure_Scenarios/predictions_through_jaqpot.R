# Set the model ID of the uploaded Jaqpot model
modelID <- "JwEkaXHV2sDZVpEUgQdE"
# Log in the Jaqpot server
jaqpotr::login.cred()
#Get model features
jaqpotr::get.model.feats(modelID)
# Create a prediction instance
df = data.frame("admin.type" = "iv",
                "admin.dose" = 10, 
                "admin.time" = 0,
                "BW"=0.25, "BW.times" = 0,
                "F_unabs" = 0, "sex" = "M",  check.names = FALSE)
# Acquire model predictions
predictions <- jaqpotr::jaqpot.predict( df = df, modelID = modelID)
# Plot an output
plot(predictions$predictions$time,predictions$predictions$Blood_total)