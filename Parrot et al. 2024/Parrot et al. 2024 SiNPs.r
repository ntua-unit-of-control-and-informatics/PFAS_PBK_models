#THis code replicates the model in "A Minimal PBPK Model Describes 
#the Differential Disposition of Silica Nanoparticles In Vivo "

library(tidyverse)
library(deSolve)  


#=========================
# 1. Parameters of the model
#=========================
create.params <- function(BW) {
  #Organ volumes for 70 kg human
  Vp70= 3 #Plasma Volume in L;  Davis and Morris et al.1993
  Vmps70= 2.125 #MPS (mononuclear phagocytic system) consists of liver and spleen in L; Davis and Morris et al.1993
  Vlu70= 1.170 #Lung in L; Davis and Morris et al.1993
  Vk70= 0.28 #Kidneys in L;Davis and Morris et al.1993
  Vo70= 38.4539 #Other Body Organs in L; Davis and Morris et al.1993
  
  #Organ volumes scaling with body weight
  Vp= Vp70*(BW/70)^0.75
  Vmps=Vmps70*(BW/70)^0.75
  Vlu=Vlu70*(BW/70)^0.75
  Vk=Vk70*(BW/70)^0.75
  Vo=Vo70*(BW/70)^0.75
  
  #Blood flows through organs for 70 kg human
  Qlu70= 330 #Blood flow through lungs in L/h; Ivanov et al, 2013
  Qmps70= 87.0#80.8 #Blood flow through MPS in L/h; Davis and Morris et al.1993
  Qk70= 74.4 #Blood flow through Kidneys in L/h; Davis and Morris et al.1993
  Qo70= 167.4 #Blood flow through other organs in L/h; Tian et al. 2012
  
  #Blood flows through organs scaling with body weight
  Qlu=Qlu70*(BW/70)^0.75
  Qmps=Qmps70*(BW/70)^0.75
  Qk=Qk70*(BW/70)^0.75
  Qo=Qo70*(BW/70)^0.75
  
  #Chemical Specific parameters
  #Fraction unbound of NPs in compartments; Parrot et al. 2024
  flu=0.889
  fmps=0.155
  fk=0.078
  fo=0.964
  
  ke_mps=0.452 #Excrection rate constant in MPS in 1/h; Parrot et al. 2024
  ke_k= 2.67e-4 #Excrection rate constant in Kidneys in 1/h; Parrot et al. 2024
  
  return(list(
    "Vp"= Vp, "Vmps"=Vmps, "Vlu"=Vlu, "Vk"=Vk, "Vo"= Vo,
    "Qlu"=Qlu, "Qmps"=Qmps, "Qk"=Qk, "Qo"=Qo, 
    "flu"= flu, "fmps"= fmps, "fk"= fk, "fo" = fo,
    "ke_mps"= ke_mps, "ke_k"=ke_k))
}

#===============================================
# 2. Function to create initial values for ODEs
#===============================================
create.inits<-function(parameters){
  with(as.list(parameters),{
    return(c(
      "Cp"= 0, "Cmps"= 0, "Clu"=0, "Ck"=0, "Co"=0, 
      "Mu"=0, "Mf"=0
    ))
  })
}

#===================
# 3. Events function
#===================
create.events <- function(admin_type,admin_dose_iv,admin_time_iv ) {
  if (admin_type == "iv") {
    ldose <- length(admin_dose_iv)
    ltimes <- length(admin_time_iv)
    if (ltimes != ldose) {
      stop("The times of administration should be equal in number to the doses")
    }
    events <- list(data = data.frame(
      var = "Cp", time = admin_time_iv,
      value = admin_dose_iv, method = "rep"))}
  else{
    print("Not available administration type.")
  }
  return(events)
}

#===============
# 5. ODEs System
#===============
ode.func<- function(time, inits, params) {
  with(as.list(c(inits, params)), {
    
    #SiNP concentration kinetics in plasma
    dCp=((Qlu*flu*Clu+Qmps*fmps*Cmps+Qk*fk*Ck+Qo*fo*Co)-(Qlu+Qmps+Qk+Qo)*Cp)/Vp
    
    #SiNP concentration kinetics in MPS
    dCmps=Qmps/Vmps*(Cp-fmps*Cmps)-(1-fmps)*ke_mps*Cmps
    
    #SiNP concentration kinetics in lungs
    dClu=Qlu/Vlu*(Cp-flu*Clu)
    
    #SiNP concentration kinetics in kidneys
    dCk=Qk/Vk*(Cp-fk*Ck)-(1-fk)*ke_k*Ck
    
    #SiNP concentration kinetics in others
    dCo=Qo/Vo*(Cp-fo*Co)
    
    #SiNP mass kinetics in urine
    dMu=(1-fk)*ke_k*Ck*Vk
    
    #SiNP mass kinetics in feces
    dMf=(1-fmps)*ke_mps*Cmps*Vmps
    
    list(c(dCp, dCmps, dClu, dCk, dCo, dMu, dMf))
  })
}

# ================================================================================
# 6. SIMULATION SETUP
# ================================================================================

# --- Dosing and Administration Settings ---
BW= 70 # Body weight (kg)
admin_type = "iv"  # administration type: "iv"

# Consider molar dose that ranges from 3.4 to 6.7 nmol
admin_dose_iv = 6.7  # administered dose through IV (nmol)
admin_time_iv = 4

# Define simulation time points
sample_time <- seq(4, 72, 0.01)  # 0 to 168 hours (1 week), with 0.1 hour steps

params <- create.params(BW)
inits <- create.inits(params)
events <- create.events(admin_type,admin_dose_iv/params$Vp,admin_time_iv)

solution <- data.frame(deSolve::ode(
  times = sample_time, func = ode.func, y = inits, parms = params,
  events = events,
  method = "lsodes", rtol = 1e-05, atol = 1e-05
))


time<-solution$time
plasma_pred<-(solution$Cp/1000)/admin_dose_iv*100
mps_pred<-(solution$Cmps/1000)/admin_dose_iv*100
kidney_pred<-(solution$Ck/1000)/admin_dose_iv*100
lung_pred<-(solution$Clu/1000)/admin_dose_iv*100
other_pred<-(solution$Co/1000)/admin_dose_iv*100
urine_pred<-(solution$Mu)/admin_dose_iv*100
preds<-data.frame(time,plasma_pred,mps_pred,kidney_pred,lung_pred,other_pred, urine_pred)

# Validation Data

speciment_data<-read.csv("C_dot_biological_specimens_blood_urine.csv")
organ_data<-read.csv("C_dot_biological_specimens_organs.csv")

# Create the plot
#Plasma plot
plot_plasma <- ggplot(preds) +
  geom_line(aes(x = time, y = plasma_pred)) +
  geom_point(
    data = speciment_data,
    aes(
      x = time,  
      y = Blood.Median
    #   ymin = Blood.min,
    #   ymax = Blood.max
    )
  ) +
  labs(
    title = "Plasma Concentration of SiNPs \n as a fraction of injected dose ",
    x = "Time (hours)",
    y = "%ID/g"
  ) +
  theme_minimal()
print(plot_plasma)

#Urine plot
plot_urine <- ggplot(preds) +
  geom_line(aes(x = time, y = urine_pred+0.04)) +
  geom_point(
    data = speciment_data,
    aes(
      x = time,
      y = Urine.Median
    #   ymin = Urine.min,
    #   ymax = Urine.max
    )
  ) +
  labs(
    title = "Urine Concentration of SiNPs \n as a fraction of injected dose ",
    x = "Time (hours)",
    y = "%ID/g"
  ) +
  theme_minimal()

print(plot_urine)


#MPS plot
plot_mps <- ggplot(preds) +
  geom_line(aes(x = time, y =  (mps_pred)))+
  geom_point(
    data = organ_data,
    aes(
      x = time,
      y = MPS.Median
    #   ymin = MPS.min,
    #   ymax = MPS.max
    )
  ) +
  labs(
    title = "MPS Concentration of SiNPs \n as a fraction of injected dose ",
    x = "Time (hours)",
    y = "%ID/g"
  ) +
  theme_minimal()
print(plot_mps)


#Lung plot
plot_lung <- ggplot(preds) +
  geom_line(aes(x = time, y =  (lung_pred)))+
  geom_point(
    data = organ_data,
    aes(
      x = time,
      y = Lung.Median
    #   ymin = Lung.min,
    #   ymax = Lung.max
    )
  ) +
  labs(
    title = "lung Concentration of SiNPs \n as a fraction of injected dose ",
    x = "Time (hours)",
    y = "%ID/g"
  ) +
  theme_minimal()
print(plot_lung)



#Kidney plot
plot_kidney <- ggplot(preds) +
  geom_line(aes(x = time, y =  (kidney_pred))) +
  geom_point(
    data = organ_data,
    aes(
      x = time,
      y = Kidney.Median
    #   ymin = Kidney.min,
    #   ymax = Kidney.max
    )
  ) +
  labs(
    title = "kidney Concentration of SiNPs \n as a fraction of injected dose ",
    x = "Time (hours)",
    y = "%ID/g"
  ) +
  theme_minimal()
print(plot_kidney)