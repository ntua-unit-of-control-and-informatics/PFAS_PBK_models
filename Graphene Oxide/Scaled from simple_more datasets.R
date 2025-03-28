# Graphene oxide model with population modeling using nimble
# Real implementation using actual Liu et al. 2012 dataset

library(deSolve)
library(nimble)
library(openxlsx)
library(coda)

setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide")
# Define molecular weight
MW <- 124.91 # g/mol

# Keep your original functions unchanged
create.params <- function(user.input){
  with(as.list(user.input),{
    
    #Volumes of organs as percent of BW
    #blood, kidneys, liver, stomach, Small intestine, large intestine, lungs, 
    #spleen, heart, brain, RoB  
    
    #Blood
    PVB <- 1.7e-3/0.02 #Davies et al. 1993, 1.7 for BW = 0.02 kg
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 1e-3 /0.02 #Davies et al. 1993, 1.0 for BW = 0.02 kg
    Vplasma <- PVplasma * BW #plasma volume kg=L
    VBven <- BW*0.5/20 	#volume of venous plasma (L); from doi:https://doi.org/10.1007/bf02353860
    VBart <- BW*0.22/20	#volume of arterial plasma (L); from doi:https://doi.org/10.1007/bf02353860
    
    #Kidney
    PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
    VKi <- PVKi * BW #kidney volume kg=L
    
    #Liver
    PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
    VLi <- PVLi * BW #liver volume kg=L
    V_macro_Li = 27.5/100*VLi # https://doi.org/10.3892/etm.2019.7450
    
    #Stomach
    PVSt <- 6e-3 #Brown et al. 1997, Table 4
    VSt <- PVSt * BW #Stomach volume kg=L
    
    #Small intestine
    PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
    VSIn <- PVSIn * BW #Small intestine volume kg=L
    
    #Large intestine
    PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
    VLIn <- PVLIn * BW #Large intestine volume kg=L
    
    #Lung
    PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
    VLn <- PVLn * BW #Lung volume kg=L
    V_macro_Ln <- 7.5/100 * VLn #random
    
    #Spleen
    PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
    VSpl <- PVSpl * BW #Spleen volume kg=L
    V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
    
    #Heart
    PVH <- 5e-3 #Brown et al. 1997, Table 4
    VH <- PVH * BW #Heart volume kg=L
    
    #Brain
    PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
    VBr <- PVBr * BW #Brain volume kg=L
    
    #RoB
    PVRe <- 1 - PVB - PVKi - PVLi - PVSt - PVSIn - PVLIn - PVLn - PVSpl - PVH - PVBr
    VRe <- PVRe * BW #volume of the rest of the body kg=L
    
    #Capillary surface area for each tissue (Ai) as percentage of body weight (m^2/kg), values from pkSim "Endothelial Surface area"
    
    PAKi <- 33.92e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AKi <- PAKi * BW #the surface area of kidney (m^2)
    
    PALi <- 142.0e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALi <- PALi * BW #liver surface area (m^2)
    
    PASt <- 3.34e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASt <- PASt * VSt * 1e3 #stomach surface area (m^2)
    
    PASIn <- 9.62e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASIn <- PASIn * VLIn * 1e3 #Small intestine surface area (m^2)
    
    PALIn <- 5.88e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALIn <- PALIn * VLIn * 1e3 #Small intestine surface area (m^2)
    
    PLn <- 59.47e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALn <- PLn* BW #lung surface area (m^2)
    
    PSpl <- 26.79e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASpl <- PSpl* BW #spleen surface area (m^2), same as muscle #assumption
    
    PH <-  23.65e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AH <- PH* VH #heart surface area (m^2)
    
    PBr <- 5.98e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ABr <- PBr* VBr *1e03 #brain surface area (m^2)
    
    ARe <- (AKi+ALi+ASt+ASIn+ALIn+ALn+ASpl+AH+ABr)/9 #assumption ???
    
    
    np_size_small <- 148 #nm,  #3/2
    np_size_large <- 556 #nm, #300/2
    
    
    #(QBi, in L/min) to different tissues (i=L, K, G, A, M, R)
    #Hall et al., 2012, https://doi.org/10.1002/jps.22811
    #PQi, fraction of cardiac output
    
    BW_ref <- 0.03 #kg
    Cardiac_output <- 11.4/1000 * BW/BW_ref #mL/min --> L/min
    
    QKi <- 1.3/1000 * BW/BW_ref # mL/min --> L/min
    PQKi <- QKi/Cardiac_output
    QLi <- 1.8/1000 * BW/BW_ref # mL/min --> L/min
    PQLi <- QLi/Cardiac_output
    QLn <- 1
    QSpl <- 0.09/1000 * BW/BW_ref # mL/min --> L/min
    PQSpl <- QSpl/Cardiac_output
    QH <- 0.28/1000 * BW/BW_ref # mL/min --> L/min
    PQH <- QH/Cardiac_output
    QBr <- 0.26/1000 * BW/BW_ref # mL/min --> L/min
    PQBr <- QBr/Cardiac_output
    PQSt <-((0.53+1.45)/2)/100 #Table 1 (mean of fore and glandular), https://doi.org/10.1002/jat.2550030607
    QSt <- PQSt*Cardiac_output 
    PSIn <- 13.5/100 #Table 1, https://doi.org/10.1002/jat.2550030607
    QSIn <- PSIn*Cardiac_output 
    PLIn <- 3.67/100 #Table 1, https://doi.org/10.1002/jat.2550030607
    QLIn <- PLIn*Cardiac_output 
    PQRe <- 1 - PQKi - PQLi - PQSpl - PQH - PQBr - PQSt - PSIn - PLIn
    QRe <- PQRe*Cardiac_output
    
    
    Qtotal <- QKi+QLi+QRe+QSpl+QH+QBr+QSIn+QSt+QLIn
    QLitot <- QLi+QSpl+QSIn+QSt+QLIn
    
    coef_liver <- 0.5     # Initial values - these will be updated by MCMC
    coef_spleen <- 0.5    
    coef_kidney <- 0.5    
    coef_heart <- 0.5     
    coef_lung <- 0.5       
    coef_brain <- 0.5      
    coef_rob <- 0.5        
    coef_stomach <- 0.5   
    coef_smallIn <- 0.5   
    coef_largeIn <- 0.5    
    CLurine <- 0.05        # Initial elimination rate constants
    CLfeces <- 0.05
    
    return(list('VB'=VB,'Vplasma'=Vplasma,'VBven'=VBven,'VBart'=VBart,
                'VKi'=VKi,'VLi'=VLi, 'V_macro_Li'=V_macro_Li,'VSt'=VSt,'VSIn'=VSIn,
                'VLIn'=VLIn,'VLn'=VLn, 'V_macro_Ln'=V_macro_Ln,
                'VSpl'=VSpl,'VH'=VH,'VBr'=VBr,'VRe'=VRe,
                'V_macro_Spl'=V_macro_Spl,
                'AKi'=AKi,'ALi'=ALi,'ASt'=ASt,'ASIn'=ASIn,'ALIn'=ALIn,
                'ALn'=ALn,'ASpl'=ASpl,'AH'=AH,'ABr'=ABr,'ARe'=ARe,
                
                'QKi'=QKi, 'QLi'=QLi, 'QRe'=QRe, 'QLn'=QLn, 'QSpl'=QSpl, 'QH'=QH,
                'QBr'=QBr, 'QSt'=QSt,'QSIn'=QSIn, 'QLIn'=QLIn,
                
                "coef_liver"=coef_liver, "coef_spleen"=coef_spleen,"coef_kidney"=coef_kidney,
                "coef_heart"=coef_heart, "coef_lung"=coef_lung, "coef_brain"=coef_brain,
                "coef_rob"=coef_rob,"coef_stomach"=coef_stomach, "coef_smallIn"=coef_smallIn,
                "coef_largeIn"=coef_largeIn,"CLurine"=CLurine,
                "CLfeces"=CLfeces,
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW, "np_size"=np_size,
                "np_size_small"=np_size_small, "np_size_large"=np_size_large,
                "Qtotal"=Qtotal, "QLitot"=QLitot,"sex"=sex
                
                
                
    ))
  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
    # Blood concentration
    CBven <- MBven/VBven
    CBart <- MBart/VBart
    
    # Kidney 
    CKi = MKi/VKi # tissue concentration
    
    #Liver
    CLi = MLi/VLi # tissue concentration
    
    #Stomach
    CSt = MSt/VSt # tissue concentration
    
    #Small Intestine
    CSIn = MSIn/VSIn # tissue concentration
    
    #Large Intestine
    CLIn = MLIn/VLIn # tissue concentration
    
    #Lungs
    CLn = MLn/VLn # tissue concentration
    
    #Spleen
    CSpl = MSpl/VSpl # tissue concentration
    
    #Heart
    CH = MH/VH # tissue concentration
    
    #Brain
    CBr = MBr/VBr # tissue concentration
    
    #Rest-of-the-body
    CRe = MRe/VRe # tissue concentration
    
    #========================================================================================================
    #Arterial Blood
    dMBart =  QLn*(CLn/coef_lung) -QKi*CBart -QLi*CBart -QSt*CBart -QSIn*CBart -QLIn*CBart -
      QSpl*CBart -QH*CBart -QBr*CBart -QRe*CBart
    
    #Venous Blood
    dMBven = -QLn*CBven +QKi*(CKi/coef_kidney) +QLitot*(CLi/coef_liver) +
      QH*(CH/coef_heart) +QBr*(CBr/coef_brain) +QRe*(CRe/coef_rob)
    
    #Kidney
    dMKi = QKi*(CBart - (CKi/coef_kidney)) - CLurine*MKi
    
    #Liver
    dMLi =  QLi*CBart - QLitot*(CLi/coef_liver) + QSt*(CSt/coef_stomach) + QSpl*(CSpl/coef_spleen) + 
      QSIn*(CSIn/coef_smallIn) + QLIn*(CLIn/coef_largeIn) 
    
    #Stomach
    dMSt = QSt*(CBart - (CSt/coef_stomach))
    
    #Small Intestine
    dMSIn = QSIn*(CBart - (CSIn/coef_smallIn)) 
    
    #Large Intestine
    dMLIn = QLIn*(CBart - (CLIn/coef_largeIn)) - CLfeces*CLIn
    
    #Lung 
    dMLn = QLn*(CBven - (CLn/coef_lung))
    
    #Spleen
    dMSpl = QSpl*(CBart - (CSpl/coef_spleen)) 
    
    #Heart
    dMH = QH*(CBart - (CH/coef_heart))
    
    #Brain
    dMBr = QBr*(CBart - (CBr/coef_brain))
    
    #Rest of body
    dMRe = QRe*(CBart - (CRe/coef_rob))
    
    # Urine
    dMurine = CLurine*MKi
    # Feces
    dMfeces = CLfeces*CLIn
    
    dVurine = CLurine
    dVfeces = CLfeces
    
    #Concentration calculation in each compartment 
    Cven <- CBven
    Cart <- CBart
    Cblood <- (MBven + MBart)/ (VBven + VBart)
    Ckidneys <- MKi/VKi
    Cliver <- MLi /VLi
    Cstomach <-  MSt/VSt
    Csmall_intestine <-  MSIn/VSIn
    Clarge_intestine <- MLIn/VLIn
    Clungs <-  MLn /VLn
    Crest <-  MRe/VRe
    Cspleen <- MSpl /VSpl
    Cheart <- MH/VH
    Cbrain <-  MBr/VBr
    
    Mven <- MBven
    Mart <- MBart
    Mblood <- MBven + MBart
    Mkidneys <- MKi
    Mliver <- MLi
    Mstomach <- MSt
    Msmall_intestine <- MSIn
    Mlarge_intestine <- MLIn
    Mlungs <- MLn 
    Mrest <- MRe
    Mfeces <- Mfeces
    Murine <- Murine
    Mspleen <-  MSpl
    Mheart <-  MH
    Mbrain <-  MBr
    
    list(c( 
      "dMBart"=dMBart, "dMBven"=dMBven,"dMKi"=dMKi,
      "dMLi"=dMLi, "dMSt"=dMSt, "dMSIn"=dMSIn, "dMLIn"=dMLIn, 
      "dMLn"=dMLn, "dMSpl"=dMSpl, "dMH"=dMH, "dMBr"=dMBr,
      "dMRe"=dMRe, "dMurine"=dMurine, 
      "dMfeces"=dMfeces,  
      "dVurine"=dVurine, "dVfeces"=dVfeces), 
      
      'Cblood'=Cblood, 'Ckidneys'=Ckidneys, 'Cliver'=Cliver, 'Cstomach'=Cstomach,
      'Csmall_intestine'=Csmall_intestine, 'Clarge_intestine'=Clarge_intestine,
      'Clungs'=Clungs, 'Crest'=Crest,
      'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain,
      
      'Mblood'=Mblood, 'Mkidneys'=Mkidneys, 'Mliver'=Mliver, 'Mstomach'=Mstomach,
      'Msmall_intestine'=Msmall_intestine, 'Mlarge_intestine'=Mlarge_intestine,
      'Mlungs'=Mlungs, 'Mrest'=Mrest, 
      'Mspleen'=Mspleen, 'Mheart'=Mheart, 'Mbrain'=Mbrain,
      
      'CBven'=CBven, 'CBart'=CBart,'CKi'=CKi, 
      'CLi'=CLi, 'CSt'=CSt, 'CSIn'=CSIn, 
      'CLIn'=CLIn, 'CLn'=CLn, 'CSpl'=CSpl, 
      'CH'=CH, 'CBr'=CBr,'CRe'=CRe)
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MBart<-0; MBven<-0;MKi<-0; MLi<-0; MSt<-0;
    MSIn<-0; MLIn<-0; MLn<-0; MSpl<-0; MH<-0; MBr<-0;
    MRe<-0; Murine<-0; Mfeces<-0; Vurine <-0; Vfeces <-0
    
    return(c(
      "MBart"=MBart, "MBven"=MBven,"MKi"=MKi,
      "MLi"=MLi, "MSt"=MSt, "MSIn"=MSIn,
      "MLIn"=MLIn,"MLn"=MLn,"MSpl"=MSpl, 
      "MH"=MH, "MBr"=MBr,"MRe"=MRe,
      "Murine"=Murine, "Mfeces"=Mfeces, 
      "Vurine"=Vurine, "Vfeces"=Vfeces
    ))
  })
}

create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "iv"){
        events <- list(data = rbind(data.frame(var = c("MBven"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("MSt"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }
    }
    return(events)
  })
}

run_ode_model <- function(params, time_points) {
  
  
  inits <- create.inits(params)
  events <- create.events(params)
  
  # Run ODE model 
 
    ode_result <- ode(
      y = inits,
      times = time_points,
      func = ode.func,
      parms = params,
      events = events,
      method = "lsodes"
    )
    
    return(ode_result)
    
}

update_params_with_mcmc <- function(params, mcmc_samples) {
  # Update partition coefficients with posterior means from MCMC
  params$coef_liver <- mcmc_samples$statistics["coef_liver", "Mean"]
  params$coef_spleen <- mcmc_samples$statistics["coef_spleen", "Mean"]
  params$coef_kidney <- mcmc_samples$statistics["coef_kidney", "Mean"]
  params$coef_heart <- mcmc_samples$statistics["coef_heart", "Mean"]
  params$coef_lung <- mcmc_samples$statistics["coef_lung", "Mean"]
  params$coef_brain <- mcmc_samples$statistics["coef_brain", "Mean"]
  params$coef_rob <- mcmc_samples$statistics["coef_rob", "Mean"]
  params$coef_stomach <- mcmc_samples$statistics["coef_stomach", "Mean"]
  params$coef_smallIn <- mcmc_samples$statistics["coef_smallIn", "Mean"]
  params$coef_largeIn <- mcmc_samples$statistics["coef_largeIn", "Mean"]
  params$CLurine <- mcmc_samples$statistics["CLurine", "Mean"]
  params$CLfeces <- mcmc_samples$statistics["CLfeces", "Mean"]
  
  return(params)
}

pbkCode <- nimbleCode({
  
  coef_liver ~ dunif(0.01, 1.0)     
  coef_spleen ~ dunif(0.01, 1.0)    
  coef_kidney ~ dunif(0.01, 1.0)    
  coef_heart ~ dunif(0.01, 1.0)     
  coef_lung ~ dunif(0.01, 1.0)       
  coef_brain ~ dunif(0.01, 1.0)      
  coef_rob ~ dunif(0.01, 1.0)        
  coef_stomach ~ dunif(0.01, 1.0)   
  coef_smallIn ~ dunif(0.01, 1.0)   
  coef_largeIn ~ dunif(0.01, 1.0)    
  CLurine ~ dunif(0.01, 1.0)
  CLfeces ~ dunif(0.01, 1.0)  

  # Measurement error model - simplified
  sigma ~ dunif(0.01, 1.0)  
  
  
  # Likelihood - observations with log-normal error
  for(i in 1:N) {
    # Log-normal observation model
    y[i] ~ dlnorm(log(pred[i]), sd = sigma)
  }
})
  
  Liu_1_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues.xlsx")
  dataset <- list("df1" = Liu_1_small_tissues)
  ##########################
  #-------------------------
  # Liu et al., 2012
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 3/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  
  tissue_cols <- c("Mheart", "Mliver", "Mspleen", "Mstomach", 
                   "Mkidneys", "Mlungs", "Mbrain", 
                   "Msmall_intestine", "Mlarge_intestine")
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df1=========================================================
  
  
  N <- 9
  preds <- solution[solution$time %in% unique(exp_data$time), c("time", tissue_cols)]
  # Print column names to debug
  print("Columns in preds dataframe:")
  print(colnames(preds))
  
  # Flatten the predictions to match observed data in long format
  library(tidyr)
  preds_long <- gather(preds, key = "Tissue", value = "predicted", -time)
  
  # Print to verify structure
  print("Structure of preds_long:")
  print(head(preds_long))
  
  # Remove the M prefix to match tissue names
  preds_long$tissue <- gsub("^M", "", preds_long$Tissue)
  
  # Print to verify the tissue column
  print("Unique values in preds_long$tissue:")
  print(unique(preds_long$Tissue))
  
  # Prepare data for MCMC
  # Convert from wide to long format for your observed data
  exp_data_long <- exp_data
  exp_data_long$mass <- exp_data_long$mass * admin.dose / 100 # Convert to mass units
  
  # Print to verify structure
  print("Structure of exp_data_long:")
  print(head(exp_data_long))
  print("Unique values in exp_data_long$tissue:")
  print(unique(exp_data_long$Tissue))
  
  # Merge predictions with experimental data
  comparison <- merge(exp_data_long, preds_long, by = c( "Tissue","time"), all.x = TRUE, all.y = FALSE)
  
  # Print to verify the merge worked
  print("Number of rows in comparison:")
  print(nrow(comparison))
  print("Head of comparison dataframe:")
  print(head(comparison))
  
  
  # # Create data and constants for the model
  # pbkData <- list(
  #   y = y_obs,  # Your observed concentration data
  #   t = t_obs,  # Your observed time points
  #   N = length(t_obs)
  # )
  
  # Prepare Nimble data
  nimble_data <- list(
    y = comparison$mass,
    pred = comparison$predicted,
    N = nrow(comparison)
  )
  
  # Set up initial values for each chain with different starting points
  inits <- list()
  for(i in 1:4) {
    inits[[i]] <- list(
      
      coef_liver = runif(1, 0.05, 0.5),     
      coef_spleen = runif(1, 0.05, 0.5),  
      coef_kidney = runif(1, 0.05, 0.5),    
      coef_heart = runif(1, 0.05, 0.5),    
      coef_lung = runif(1, 0.05, 0.5),       
      coef_brain = runif(1, 0.05, 0.5),       
      coef_rob = runif(1, 0.05, 0.5),         
      coef_stomach = runif(1, 0.05, 0.5),    
      coef_smallIn = runif(1, 0.05, 0.5),    
      coef_largeIn = runif(1, 0.05, 0.5),     
      CLurine = runif(1, 0.05, 0.5), 
      CLfeces = runif(1, 0.05, 0.5), 
      
      sigma = runif(1, 0.1, 0.3)
      
    )
  }
  
  
  # Create and compile the model
  pbkModel <- nimbleModel(pbkCode, data = nimble_data, inits = inits[[1]])
  cPbkModel <- compileNimble(pbkModel)
  
  # Configure MCMC
  mcmcConf <- configureMCMC(pbkModel, monitors = c("coef_liver", "coef_spleen","coef_kidney",
                                    "coef_heart", "coef_lung", "coef_brain", "coef_rob", 
                                    "coef_stomach", "coef_smallIn", "coef_largeIn",
                                    "CLurine", "CLfeces","sigma"))
  # Build the MCMC from the configuration
  mcmc <- buildMCMC(mcmcConf)
  
  # Compile the MCMC
  cmcmc <- compileNimble(mcmc)
  
  # Run MCMC with multiple chains
  niter <- 10000
  nburnin <- 5000
  samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, 
                     nchains = 4, inits = inits, 
                     samplesAsCodaMCMC = TRUE)
  

  # Summary of posterior distributions
  summary_results <- summary(samples)
  print(summary_results)
  
  # Update parameters with MCMC results
  updated_params <- update_params_with_mcmc(params, summary_results)
  new_inits <- create.inits(updated_params)
  new_events <- create.events(updated_params)
  
  # Run model with updated parameters
  updated_solution <- data.frame(deSolve::ode(times = sample_time, 
                                              func = ode.func,
                                              y = new_inits, 
                                              parms = updated_params,
                                              events = new_events,
                                              method="lsodes", 
                                              rtol = 1e-05, 
                                              atol = 1e-05))
  
  # Extract updated predictions for comparison
  updated_preds <- updated_solution[updated_solution$time %in% unique(exp_data$time), 
                                    c("time", tissue_cols)]
  
  # Compare original and updated predictions
  plot(exp_data$time, exp_data$mass * admin.dose / 100, 
       xlab = "Time", ylab = "Mass", pch = 19,
       main = "Model fit before and after MCMC calibration")
  lines(preds$time, preds$Mliver, col = "blue", lty = 2)
  lines(updated_preds$time, updated_preds$Mliver, col = "red")
  legend("topright", legend = c("Observed", "Initial model", "Calibrated model"),
         col = c("black", "blue", "red"), pch = c(19, NA, NA), lty = c(NA, 2, 1))