#THis code replicates the model in "A Minimal PBPK Model Describes 
#the Differential Disposition of Silica Nanoparticles In Vivo "

library(tidyverse)
library(deSolve)  


#=========================
# 1. Parameters of the model
#=========================
create.params <- function(user_ipnuts){
  with(as.list(user_ipnuts),{
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
    "ke_mps"= ke_mps, "ke_k"=ke_k, "admin_dose_iv"= admin_dose_iv, "admin_time_iv"= admin_time_iv
  ))
})
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
create.events <- function(parameters){
  with(as.list(parameters),{
    if(length(admin_time_iv) != length(admin_dose_iv)) {
      stop("The times of administration should be equal in number to the doses")
    }else{  
      events <- list(data = data.frame(
      var = "Cp", time = admin_time_iv,
      value = admin_dose_iv, method = "rep"))
      }
  
    return(events)
  })
}
#===================
# 4. Custom function
#===================
custom.func <- function(){
 }
#===============
# 5. ODEs System
#===============
ode.func<- function(time, inits, params,custom.func) {
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

user_inputs <- list(
BW= 70 ,
admin_dose_iv = 6.7,
admin_time_iv = 4
)


params <- create.params(user_inputs)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0, 48, by = 1.0)
# (admin_type,admin_dose_iv/params$Vp,admin_time_iv)
solution <- data.frame(deSolve::ode(
  times = sample_time, func = ode.func, y = inits, parms = params,
  events = events,
  method = "lsodes", rtol = 1e-05, atol = 1e-05
))

predicted.feats<-c("Cp","Cmps","Clu","Ck","Co","Mu","Mf")

jaqpotr::deploy.pbpk(user.input = user_inputs,out.vars = predicted.feats,
create.params = create.params,  create.inits = create.inits,
create.events = create.events, custom.func = custom.func,
envFile = "C:\\Users\\fotis\\Desktop\\KEY.env")
