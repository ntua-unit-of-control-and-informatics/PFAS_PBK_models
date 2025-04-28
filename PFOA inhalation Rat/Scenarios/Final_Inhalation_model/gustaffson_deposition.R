

setwd( "C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Final_Inhalation_model/Training/AAFE/inhalation_permeability")
inh_df <- read.csv("gustaffson.csv")
mean(inh_df$concentration)
C_PFOA <- 0.333 #ug PFOA/mg dust
TV <- 3.89*1E-3 #L
BF <- 166 #breaths per minute
df <- 0.67
dose <- rep(NA, dim(inh_df)[1])
inh_df$time[1] <- 0.05
inh_df$time[2] <- 0.1
inh_df$time[3] <- 0.15

for (i in 1:(dim(inh_df)[1]-1)){
  duration <- (inh_df$time[i+1] - inh_df$time[i]) #min
  aerosol_con <-  1000*(inh_df$concentration[i]+inh_df$concentration[i+1])/2
  if (aerosol_con < 0){
    aerosol_con <- 0
  }
  dose[i] <- duration* aerosol_con *0.34*df#*C_PFOA
}
print(sum(dose, na.rm = TRUE))

0.51