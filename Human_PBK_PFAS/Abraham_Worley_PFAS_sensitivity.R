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
    Htc = 0.467 #hematocrit for the rat; Davies 1993
    
    #Tissue Volumes 
    VplasC = 0.0428 #fraction vol. of plasma (L/kg BW); Davies 1993
    VLC = 0.026 #fraction vol. of liver (L/kg BW); Brown 1997
    VKC = 0.004 #fraction vol. of kidney (L/kg BW); Brown 1997
    VfilC = 4e-4	#fraction vol. of filtrate (L/kg BW)
    VPTCC = 1.35e-4 #vol. of proximal tubule cells (L/g kidney) (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
    
    #Chemical Specific Parameters
    MW = 414.07	#PFOA molecular mass (g/mol)
    
    #Free = 0.001#value 0.001 Smeltz 2023 alternative value 0.00245 from Ryu 2024
    
    #Kidney Transport Parameters
    Vmax_baso_invitro = 439.2 #Vmax of basolateral transporter (pmol/mg protein/min); averaged in vitro value of OAT1 and OAT3 from Nakagawa, 2007
    Km_baso = 20100 #Km of basolateral transporter (ug/L) Average of OAT1 and OAT3 from Nakagawa et. al, 2007
    Vmax_apical_invitro = 37400 #Vmax of apical transporter (pmol/mg protein/min); invitro value for OAT4 from Yang et al, 2010
    Km_apical = 77500#Km of apical transporter (ug/L), in vitro value for OAT4 and URAT1 from Yang et al, 2010.
    #RAFbaso = 1	#relative activity factor, basolateral transporters (male) (fit to data)	
    #RAFapi = 0.0007	#relative activity factor, apical transporters (male) (fit to data)	
    protein = 2.0e-6	#amount of protein in proximal tubule cells (mg protein/proximal tubule cell)
    GFRC = 24.19*24	#glomerular filtration rate (L/day/kg kidney) (male); Corley, 2005
    
    #Partition Coefficients (from Allendorf 2021)
    #PL = 0.434698291544763 #liver:blood
    #PK = 0.413283707125888 #kidney:blood
    #PR = 0.534206103089267 #rest of body:blood
    
    #rate constants
    #kdif = 0.001*24	#diffusion rate from proximal tubule cells (L/day)
    #kabsc = 2.12*24	#rate of absorption of chemical from small intestine to liver (1/(day*BW^-0.25))(fit to data)
    #kunabsc = 7.06e-5*24	#rate of unabsorbed dose to appear in feces (1/(day*BW^-0.25))(fit to data)
    #GEC = 3.5*24#gastric emptying time (1/(day*BW^-0.25)); from Yang, 2013
    #k0C = 1.0*24	#rate of uptake from the stomach into the liver (1/(day*BW^-0.25)) (fit to data)
    
    #keffluxc =0.1*24 #rate of clearance of PFOA from proximal tubule cells into blood (1/(day*BW^-0.25))
    #kbilec = 0.0001*24 #biliary elimination rate ((male); liver to feces storage (1/(day*BW^-0.25)) (fit to data)
    #kurinec = 0.063*24 #rate of urine elimination from urine storage (male) (1/(day*BW^-0.25))(fit to data)
    #kvoid = 0.06974*24  #daily urine volume rate (L/day); Van Haarst, 2004                                                   
    
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
                 "PL" = PL, "PK" = PK, "PR" = PR, "kvoid" = kvoid, "admin_type" = admin_type,
                 "admin_time_iv" = admin_time_iv, "admin_dose_iv" = admin_dose_iv,
                 "Cwater" = Cwater, "Cwater_time" = Cwater_time,
                 "ingestion" = ingestion,"ingestion_time" = ingestion_time,
                 "water_consumption" = water_consumption,
                 "RAFbaso"=RAFbaso,"RAFapi"=RAFapi ,"kabsc"=kabsc ,
                 "kunabsc"=kunabsc,"GEC"=GEC ,"k0C"=k0C ,"keffluxc"=keffluxc,
                 "kbilec"=kbilec ,"kurinec"=kurinec, "kvoid"=kvoid)) 
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

AUC <- function(x, y){
  individual_auc <- c()
  for (i in 1:(length(x)-1)){
    individual_auc[i] <- (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return(sum(individual_auc))
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
    CVL = CL/PL#concentration in the venous blood leaving the liver (ug/L)
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
BW <- 82

admin_type <-  "bolus" # administration type values: iv/oral/oral bolus
admin_dose_bolus <- c(3.96) # administered dose through bolus in ug
admin_time_bolus <- c(0) # time when bolus doses are administered, in days
admin_time_iv=0 # administered dose through IV
admin_dose_iv=0 # time when IV doses are administered, in days

Cwater <- 0.00 #ug/L
Cwater_time <- 0
ingestion <- 0 #ug/day
ingestion_time <- c(0)

# The following parameters have been commented out of
# create.params in ordered to be changed
Free  = 0.001
RAFbaso	= 1
RAFapi = 0.0007	
PL = 0.434698291544763
PK = 0.413283707125888
PR = 0.534206103089267 
kdif = 0.001*24
kabsc = 2.12*24
kunabsc = 7.06e-5*24
GEC = 3.5*24
k0C = 1.0*24
keffluxc =0.1*24
kbilec = 0.0001*24
kurinec =0.063*24 
kvoid = 0.06974*24

 user_input <- list( "admin_type" = admin_type,
                    "admin_dose_bolus" = admin_dose_bolus, 
                    "admin_time_bolus" = admin_time_bolus,
                    "BW"=BW, "Cwater" = Cwater, 
                    "Cwater_time" = Cwater_time, "ingestion" = ingestion,
                    "ingestion_time" = ingestion_time, "Free" =Free,
                    "RAFbaso" = RAFbaso, "RAFapi" = RAFapi, "PL"=PL,
                    "PK"=PK, "PR" = PR, "kdif"=kdif, "kabsc" = kabsc,
                    "kunabsc" = kunabsc, "GEC" = GEC, "k0C" = k0C,
                    "keffluxc" =keffluxc, "kbilec" = kbilec, "kurinec" =kurinec,
                    "kvoid" = kvoid)

#================
#7. Wrap function
#================
obj.func<-function(x,dp){
 

 if (is.na(x)){
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  sample_time=seq(0,6,0.01)
  
  solution <-as.data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                               events = events, method="bdf",rtol = 1e-05, atol = 1e-05))
  
  result <- c(
    Cmax = max(solution[, "CA"]),
    AUC = AUC(sample_time, solution[, "CA"]),
    C_6h = solution[which.min(abs(sample_time - 0.25)), "CA"],
    C_6d = solution[which.min(abs(sample_time - 6)), "CA"],
    Feces_Excreted = sum(solution[, "Afeces"]), 
    Urine_Excreted = sum(solution[, "Aurine"]))  
  
  return(result)
  } else{
      loop_user_input<-list()
      loop_user_input[[x]]<-user_input[[x]]*(1+dp)
      loopprams <- create.params(loop_user_input) # modify user_input to apply the pertribution 
      inits <- create.inits(loopprams)
      events <- create.events(loopprams)
      sample_time=seq(0,6,0.01)
      solution <-as.data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = loopprams,
                               events = events, method="bdf",rtol = 1e-05, atol = 1e-05))
      result <- c(
      Cmax = max(solution[, "CA"]),
      AUC = AUC(sample_time, solution[, "CA"]),
      C_6h = solution[which.min(abs(sample_time - 0.25)), "CA"],
      C_6d = solution[which.min(abs(sample_time - 6)), "CA"],
      Feces_Excreted = sum(solution[, "Afeces"]), 
      Urine_Excreted = sum(solution[, "Aurine"]))  
      
  }
}

#=================================
#8. Sensitivity Analysis Parameter
#=================================
thetas<-list(
"Free" ,
"RAFbaso",	
"RAFapi" ,
"PL" , 
"PK" , 
"PR" ,
"kdif" ,
"kabsc" ,
"kunabsc",
"GEC" ,
"k0C" ,
"keffluxc",
"kbilec" ,
"kurinec",
"kvoid")


#=======================
#9. Sensitivity analysis
#=======================

# dp is the percentile change of a parameter 
dp=0.5

results<- as.data.frame(
  matrix(NA, 
         nrow = length(c("Cmax", "AUC", "C_6h", "C_6d", "Feces_Excreted", 
                         "Urine_Excreted" )), 
         ncol = length(thetas),
         dimnames = list(c("Cmax", "AUC", "C_6h", "C_6d", "Feces_Excreted", 
                         "Urine_Excreted" ), thetas)))

for(i in 1:length(thetas)){
 results[,i]<-obj.func(thetas[[i]],dp)
 print(thetas[i])
}

#The Sensitivity Index table 
SI<- as.data.frame(
  matrix(NA, 
         nrow = length(c("Cmax", "AUC", "C_6h", "C_6d", "Feces_Excreted", 
                         "Urine_Excreted" )), 
         ncol = length(thetas),
         dimnames = list(c("Cmax", "AUC", "C_6h", "C_6d", "Feces_Excreted", 
                         "Urine_Excreted" ), thetas)))

results0<-obj.func(NA,dp)

for(i in 1:length(thetas)){
  SI[,i]<-(abs((results[,i]-results0)/results0)/dp)
}
SIt<-t(SI)


plot_list <- list()

# Loop through each column (each output variable)
for(j in 1:ncol(SIt)) {
  
  # Get column name for title
  col_name <- colnames(SIt)[j]
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Parameter = rownames(SIt),
    Sensitivity = SIt[, j]
  )

  plot_data$Parameter <- factor(plot_data$Parameter, 
                                levels = plot_data$Parameter[order(plot_data$Sensitivity, decreasing = FALSE)])
  
  # Create the plot
  plot_list[[j]] <- ggplot(plot_data, aes(x = Parameter, y = Sensitivity)) +
    geom_bar(stat = "identity", width = 0.6, fill = "#2E86AB") +
    labs(
      title = paste("Sensitivity Analysis:", col_name),
      x = "Parameters",
      y = "Sensitivity Index"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    # Add value labels
    geom_text(aes(label = format(Sensitivity, scientific = TRUE, digits = 2)), 
              hjust = -0.1, size = 3)
}

# Name the list elements
names(plot_list) <- colnames(SIt)

grid.arrange(grobs = plot_list, ncol = 2)

# Option 2: Display plots one by one (useful for detailed inspection)
for(j in 1:length(plot_list)) {
  print(plot_list[[j]])
}

# =========================
# Save Plots to Files
# =========================

# Create directory for plots
if(!dir.exists("sensitivity_plots")) {
  dir.create("sensitivity_plots")
}

# Save each plot as PNG
for(j in 1:length(plot_list)) {
  col_name <- colnames(SIt)[j]
  
  png_file <- paste0("sensitivity_plots/", col_name, "_sensitivity.png")
  
  ggsave(
    filename = png_file,
    plot = plot_list[[j]],
    width = 10,
    height = 8,
    dpi = 300
  )
  print(paste("Saved:", png_file))
}
