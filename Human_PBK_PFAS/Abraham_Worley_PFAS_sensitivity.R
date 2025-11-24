# The code replicates the Loccisano et al.(2012) PBK model
# and simulates IV or Oral PFOA distribution for female or male rats of constant or varying
# body weight. The Bernstein et al. (2021) code was used for replicating the
# model and the parameter set was expanded to also describe female rats
# (the original code only described male rats). 
# We did not include simultaneous IV and oral because there are no such experiments

library(deSolve)
library(ggplot2)
library(tidyverse)
library(sensitivity)


#=========================
#1. Parameters of the model
#=========================



create.params  <- function(user_input){
  with( as.list(user_input),{
    # User input: BW(kg)
    
    #Cardiac Output and Bloodflow (as fraction of cardiac output)
    QCC = 12.5*24 #cardiac output in L/day/kg^0.75; Brown 1997		 
    QLC = 0.25 #fraction blood flow to liver; Brown 1997	
    QKC = 0.175 #fraction blood flow to kidney; Brown 1997.
    Htc = 0.467 #reported in Abraham. 0.44 hematocrit for the rat; Davies 1993
    
    #Tissue Volumes 
    VplasC = 0.0428 #fraction vol. of plasma (L/kg BW); Davies 1993
    VLC = 0.026 #fraction vol. of liver (L/kg BW); Brown 1997
    VKC = 0.004 #fraction vol. of kidney (L/kg BW); Brown 1997
    VfilC = 4e-4	#fraction vol. of filtrate (L/kg BW)
    VPTCC = 1.35e-4 #vol. of proximal tubule cells (L/g kidney) (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
    
    
    #Kidney Transport Parameters
    Vmax_baso_invitro = 439.2 #Vmax of basolateral transporter (pmol/mg protein/min); averaged in vitro value of OAT1 and OAT3 from Nakagawa, 2007
    Km_baso = 20100 #Km of basolateral transporter (ug/L) Average of OAT1 and OAT3 from Nakagawa et. al, 2007
    Vmax_apical_invitro = 37400 #Vmax of apical transporter (pmol/mg protein/min); invitro value for OAT4 from Yang et al, 2010
    Km_apical = 77500#Km of apical transporter (ug/L), in vitro value for OAT4 and URAT1 from Yang et al, 2010.
    protein = 2.0e-6	#amount of protein in proximal tubule cells (mg protein/proximal tubule cell)
    GFRC = 24.19*24	#glomerular filtration rate (L/day/kg kidney) (male); Corley, 2005
    
    
    #rate constants
    kdif = 0.001*24	#diffusion rate from proximal tubule cells (L/day)
    GEC = 3.5*24#gastric emptying time (1/(day*BW^-0.25)); from Yang, 2013
    keffluxc =0.1*24 #rate of clearance of PFOA from proximal tubule cells into blood (1/(day*BW^-0.25))
    kvoid = 0.06974*24  #daily urine volume rate (L/day); Van Haarst, 2004 
    
    #Scaled Parameters
    #Cardiac output and blood flows
    QC = QCC*(BW^0.75)*(1-Htc)	#cardiac output in L/day; adjusted for plasma
    QK = (QKC*QC)	#plasma flow to kidney (L/day)
    QL = (QLC*QC)	#plasma flow to liver (L/day)
    QR = QC - QK - QL 	#plasma flow to rest of body (L/day)
    QBal = QC - (QK + QL + QR) #Balance check of blood flows; should equal zero
    
    #Tissue Volumes
    
    VPlas = VplasC*BW 	#volume of plasma (L) 
    VK = VKC*BW 	#volume of kidney (L)
    MK = VK*1.0*1000	#mass of the kidney (g)
    VKb = VK*0.16	#volume of blood in the kidney (L); fraction blood volume of kidney (0.16) from Brown, 1997
    Vfil = VfilC*BW	#volume of filtrate (L)
    VL = VLC*BW	#volume of liver (L)
    ML = VL*1.05*1000	#mass of the liver (g)
    
    #Kidney Parameters
    
    PTC = VKC*1000*6e7	#number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney, assuming density of 1 kg/L)
    VPTC = VK*1000*VPTCC	#volume of proximal tubule cells (L)	
    MPTC = VPTC*1000 #mass of the proximal tubule cells (g) (assuming density 1 kg/L)	
    VR = (0.93*BW) - VPlas - VPTC - Vfil - VL	#volume of remaining tissue (L); 
    VBal = (0.93*BW) - (VR + VL + VPTC + Vfil + VPlas)	#Balance check of tissue volumes; should equal zero 
    
    Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1e6)*24#Vmax of basolateral transporters (ug/day/kg BW)
    Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1e6)*24 #Vmax of basolateral transporters (ug/day/kg BW)
    Vmax_baso = Vmax_basoC*BW^0.75	#(ug/day)
    Vmax_apical = Vmax_apicalC*BW^0.75	#(ug/day)
    kbile = kbilec*BW^(-0.25)	#biliary elimination; liver to feces storage (/day)
    kurine = kurinec*BW^(-0.25)	#urinary elimination, from filtrate (/day)
    kefflux = keffluxc*BW^(-0.25)	#efflux clearance rate, from PTC to blood (/day)
    GFR = 163.65#glomerular filtration rate,(L/day).Calculated 190 L/Day from GFR= GFRC*VK scaled to mass of kidney(in kg)
    
    #GI Tract Parameters
    kabs = kabsc*BW^(-0.25)	#rate of absorption of chemical from small intestine to liver (/day)
    kunabs = kunabsc*BW^(-0.25)	#rate of unabsorbed dose to appear in feces (/day)
    GE = GEC*BW^(-0.25)	#gastric emptying time (/day)
    k0 = k0C*BW^(-0.25) 	#rate of uptake from the stomach into the liver (/day)
    
    water_consumption <- 1.36# L/day
    
    return(list( "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR, 
                 "VPlas" = VPlas,
                 "VKb" = VKb, "Vfil" = Vfil, "VL" = VL, "VR" = VR, "ML" = ML,
                 "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
                 'kdif' = kdif, "Km_baso" = Km_baso, "Km_apical" = Km_apical,
                 "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
                 "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs, "GE" = GE, "k0" = k0,
                 "kvoid" = kvoid, "admin_type" = admin_type,
                 "admin_time_iv" = admin_time_iv, "admin_dose_iv" = admin_dose_iv,
                 "admin_time_bolus" = admin_time_bolus, "admin_dose_bolus"= admin_dose_bolus,
                 "Cwater" = Cwater, "Cwater_time" = Cwater_time,
                 "ingestion" = ingestion,"ingestion_time" = ingestion_time,
                 "water_consumption" = water_consumption 
                 
    ))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
    "AR" = 0; "Adif" = 0; "A_baso" = 0; "AKb" = 0;
    "ACl" = 0; "Aefflux" = 0;
    "A_apical" = 0; "APTC" = 0; "Afil" = 0;
    "Aurine" = 0; "AST" = 0;
    "AabsST" = 0; "ASI" = 0; "AabsSI" = 0; "Afeces" = 0;
    "AL" = 0; "Abile" = 0; "Aplas_free" = 0;
    "ingestion" = 0; "Cwater" = 0
    
    return(c("AR" = AR, "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
             "ACl" = ACl, "Aefflux" = Aefflux,
             "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
             "Aurine" = Aurine, "AST" =AST,
             "AabsST" = AabsST, "ASI" = ASI, "AabsSI" = AabsSI, "Afeces" = Afeces, 
             "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free,
             "ingestion" = ingestion, "Cwater" = Cwater))
  })
}

#===================
#3. Events function
#===================

create.events <- function(parameters){
  with(as.list(parameters), {
    if (admin_type == "iv"){
      # Calculate number of administrated doses and corresponding administration time for IV
      ldose <- length(admin_dose_iv)
      ltimes <- length(admin_time_iv)
      # If not equal, then stop 
      if (ltimes != ldose){
        stop("The times of administration should be equal in number to the doses")
      }else{
        
        events <- list(data = rbind(data.frame(var = c("Aplas_free"),  time = admin_time_iv, 
                                               value = admin_dose_iv, method = c("add")) ))
      }
    }else if (admin_type == "oral"){
      # Calculate number of administrated doses and corresponding administration time
      lcwater <- length(Cwater)
      lcwatertimes <- length(Cwater_time)
      lingest <- length(ingestion)
      lingesttimes <- length(ingestion_time)
      # If not equal, then stop 
      if (lcwater != lcwatertimes){
        stop("The times of water concentration change should be equal in vector of Cwater")
      }else if (lingest != lingesttimes){
        stop("The times of ingestion rate change should be equal in vector of Cwater")
      }else{
        events <- list(data = rbind(data.frame(var = c("Cwater"),  time = Cwater_time, 
                                               value = Cwater, method = c("rep")),
                                    
                                    data.frame(var = c("ingestion"),  time = ingestion_time, 
                                               value = ingestion, method = c("rep"))))
      }
      
    } else if (admin_type == "bolus"){
      # Calculate number of administrated doses and corresponding administration time for oral bolus
      ldose <- length(admin_dose_bolus)
      ltimes <- length(admin_time_bolus)
      # If not equal, then stop 
      if (ltimes != ldose){
        stop("The times of administration should be equal in number to the doses")
      }else{
        
        events <- list(data = rbind(data.frame(var = c("AST"),  time = admin_time_bolus, 
                                               value = admin_dose_bolus, method = c("add"))
        ))
      }
    }
    
    
    return(events)
  })
}

#==================
#4. Custom function 
#==================

#Fitting functions
mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}

AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}


#Cumulative mass function aggragating experimental data
#Used for urine and feces samples
cumulative_exp_data<-function(df,time_col,concentration_col, multiply_col){
  
  working_df <- df[, c(time_col, concentration_col, multiply_col)]
  complete_cases <- complete.cases(working_df)
  clean_df <- working_df[complete_cases, ]
  
  if (nrow(clean_df) == 0) {
    stop("No complete cases found after removing NA values")
  }
  
  # Multiply and calculate cumulative sum
  multiplied_concentration <- clean_df[[concentration_col]] * clean_df[[multiply_col]]
  
  # Create new dataframe with time and cumulative concentration only
  result_df <- data.frame(
    time = clean_df[[time_col]],
    cumulative_mass = cumsum(multiplied_concentration)
  )
  
  return(result_df)
}

trapz <- function(x, y) {
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
}

#==============
#5. ODEs System
#==============

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    
    CR = AR/VR #concentration in rest of body (ug/L)
    CVR = CR/PR	#concentration in venous blood leaving the rest of the body (ug/L)
    CKb = AKb/VKb	#concentration in kidney blodd (ug/L) 
    CVK = CKb #/PK	#concentration in venous blood leaving kidney (ug/L)
    CPTC = APTC/VPTC	#concentration in PTC (ug/L)
    Cfil = Afil/Vfil	#concentration in filtrate (ug/L)
    CL = AL/VL	#concentration in the liver (ug/L)
    CLiver = AL/ML #	concentration in the liver (ug/g)
    CVL = CL/PL	#concentration in the venous blood leaving the liver (ug/L)
    CA_free = Aplas_free/VPlas		#concentration in plasma (ug)
    CA = CA_free/Free	#concentration of total PFOA in plasma (ug/L)
    Curine = Aurine/kvoid
    
    # Rest of Body (Tis)
    dAR = QR*(CA-CVR)*Free	#rate of change in rest of body (ug/day)
    
    #Kidney 
    #Kidney Blood (Kb)
    dAdif = kdif*(CKb - CPTC)	#rate of diffusion from into the PTC (ug/day)
    dA_baso = (Vmax_baso*CKb)/(Km_baso + CKb)	
    dAKb = QK*(CA-CVK)*Free - CA*GFR*Free - dAdif - dA_baso #rate of change in kidney blood (ug/day).
    dACl = CA*GFR*Free	#rate of clearance via glormerular filtration (ug/day)
    
    #Proximal Tubule Cells (PTC)
    dAefflux = kefflux*APTC
    dA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil)
    dAPTC =  dAdif + dA_apical + dA_baso - dAefflux #rate of change in PTC(ug/day)
    
    #Filtrate (Fil)
    dAfil = CA*GFR*Free - dA_apical - Afil*kurine	#rate of change in filtrate (ug/day)
    
    #Urinary elimination 
    dAurine = kurine*Afil	#rate of change in urine (ug/day)
    
    #GI Tract (Absorption site of oral dose)
    #Stomach
    dAST=  ingestion + Cwater*water_consumption - k0*AST - GE*AST	#rate of change in the stomach (ug/day)
    dAabsST = k0*AST	#rate of absorption in the stomach (ug/day)
    
    #Small Intestine
    dASI = GE*AST - kabs*ASI - kunabs*ASI	#rate of change in the small intestine (ug/day)
    dAabsSI = kabs*ASI	#rate of absorption inthe small intestine (ug/day)
    
    total_oral_uptake = AabsSI + AabsST	#total oral uptake in the GI tract (ug) 
    
    #Feces compartment
    dAfeces = kbile*AL + kunabs*ASI #rate of change in the feces compartment (ug/day)
    
    #Liver
    dAL = QL*(CA-CVL)*Free - kbile*AL + kabs*ASI + k0*AST #rate of change in the liver (ug/day)
    dAbile = kbile*AL
    amount_per_gram_liver = CLiver	#amount of PFOA in liver per gram liver (ug/g)
    
    #Plasma compartment
    dAplas_free = (QR*CVR*Free) + (QK*CVK*Free) + (QL*CVL*Free) - 
      (QC*CA*Free) + dAefflux  #rate of change in the plasma (ug/day) 
    
    dCwater = 0
    dingestion = 0
    
    #Mass Balance Check
    Atissue = Aplas_free + AR + AKb + Afil + APTC + AL + AST + ASI 	#sum of mass in all compartments (ug)
    Aloss = Aurine + Afeces #sum of mass lost through urinary and fecal excretion (ug)
    Atotal = Atissue + Aloss 	#total mass; should equal total dose
    
    list(c("dAR" = dAR, "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
           "dACl" = dACl, "dAefflux" = dAefflux,
           "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
           "dAurine" = dAurine, "dAST" =dAST,
           "dAabsST" = dAabsST, "dASI" = dASI, "dAabsSI" = dAabsSI, "dAfeces" = dAfeces, 
           "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
           "dCwater" = dCwater, "dingestion" = dingestion), 
         "total_oral_uptake" = total_oral_uptake, "amount_per_gram_liver" = amount_per_gram_liver,
         "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal, "CR" =CR, "CVR" = CVR, "CKb" = CKb, 
         "CVK" = CVK, "CPTC" = CPTC,
         "Cfil" = Cfil, "CL" = CL, "CVL" = CVL, "CA_free" = CA_free, "CA" = CA)
    
  })
}

#=============
#6. User input 
#=============
# PFOA Abraham Trial Runs

#--------------------#
# Changing parameter #
#--------------------#

BW <- 82

PL = 0.434698291544763 #liver:blood
PK = 0.413283707125888 #kidney:blood
PR = 0.534206103089267 #rest of body:blood
#Chemical Specific Parameters
MW = 414.07	#PFOA molecular mass (g/mol)
Free = 0.001 #value 0.001 Smeltz 2023 alternative value 0.00245 from Ryu 2024

admin_type <-  "bolus" # administration type values: iv/oral/oral bolus
admin_dose_bolus <- c(3.96) # administered dose through bolus in ug
admin_time_bolus <- c(0) # time when bolus doses are administered, in days
admin_time_iv=0 # administered dose through IV
admin_dose_iv=0 # time when IV doses are administered, in days

Cwater <- 0.00 #ug/L
Cwater_time <- 0
ingestion <- 0 #ug/day
ingestion_time <- c(0)

#Parameter for estimation
RAFbaso = 1	#relative activity factor, basolateral transporters (male) (fit to data)	
RAFapi = 0.0007	#relative activity factor, apical transporters (male) (fit to data)
kabsc = 2.12*24	#rate of absorption of chemical from small intestine to liver (1/(day*BW^-0.25))(fit to data)
kunabsc = 7.06e-5*24	#rate of unabsorbed dose to appear in feces (1/(day*BW^-0.25))(fit to data)
k0C = 1.0*24	#rate of uptake from the stomach into the liver (1/(day*BW^-0.25)) (fit to data)
kbilec = 0.0001*24 #biliary elimination rate ((male); liver to feces storage (1/(day*BW^-0.25)) (fit to data)
kurinec = 0.063*24 #rate of urine elimination from urine storage (male) (1/(day*BW^-0.25))(fit to data)


user_input <- list( "admin_type" = admin_type,
                    "admin_dose_bolus" = admin_dose_bolus, 
                    "admin_time_bolus" = admin_time_bolus,
                    "BW"=BW, "Cwater" = Cwater, 
                    "Cwater_time" = Cwater_time, "ingestion" = ingestion,
                    "ingestion_time" = ingestion_time,"PL"=PL, "PK"=PK,
                    "PR"=PR, "MW"=MW, "RAFbaso"=RAFbaso, "RAFapi"=RAFapi,
                    "kabsc"=kabsc, "kunabsc"=kunabsc, "k0C"=k0C, "kbilec"=kbilec,
                    "kurinec"=kurinec)

obj.func<-function(x){
  
 
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  sample_time=seq(0,60,0.001)
  
  solution <-as.data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                               events = events, method="bdf",rtol = 1e-05, atol = 1e-05))
  
  result <- c(
    Cmax = max(solution[, "CA"]),
    AUC = trapz(sample_time, solution[, "CA"]),
    C_24h = solution[which.min(abs(sample_time - 12)), "CA"],
    C_48h = solution[which.min(abs(sample_time - 6)), "CA"],
    Total_Excreted = max(solution[, "Aurine"]) + max(solution[, "Afeces"])
  )
  
  return(result)
  
}

parameters <- c(
  "Free",
  "kbilec",
  "kurinec",
  "Km_baso",
  "Km_apical",
  "kabsc",
  "kunabsc",
  "k0c",
  "RAFbaso",
  "RAFapi")

# Define parameter ranges (min, max)
# You'll need to set appropriate ranges based on literature or uncertainty
param_ranges <- list(
  Free = c(0.0005, 0.005),      
  kbilec = c(0.00005, 0.0005),   
  kurinec = c(0.05, 0.15),                
  Km_baso=c(20100/2,20100*2),
  Km_apical=c(37400/2,37400*2),
  kabsc=c(2.12*24/2,2.12*24*2),
  kunabsc=c(7.06e-5*24/2,7.06e-5*24*2),
  k0c=c(1.0*24/2,1.0*24*2),
  RAFbaso=c(0.0001,1),
  RAFapi =c(0.0001,1)             
)

#Morris screening

# Morris screening - efficient for identifying important parameters

set.seed(123)
morris_result <- morrisMultOut(
  model = obj.func,
  factors = parameters,
  r = 4,           # Number of trajectories
  design = list(type = "oat", levels = 10, grid.jump = 5),
  binf = sapply(param_ranges, function(x) x[1]),
  bsup = sapply(param_ranges, function(x) x[2]),scale=FALSE
)

# Plot Morris results
plot(morris_result)
print(morris_result)
