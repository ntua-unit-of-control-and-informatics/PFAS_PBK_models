#Graphene oxide model

library(deSolve)
setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide")


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
    # Qurine = (5.4+3.7+2.7)/3/1000/60/0.02 #mL/h/kg--> L/min/kg, unknown BW, https://doi.org/10.1007/BF02035147
    # Qfeces = 6.68/20.6 #mg/g --> g/kg BW, https://doi.org/10.3390/toxins10050204
     
    
    return(list('VB'=VB,'Vplasma'=Vplasma,'VBven'=VBven,'VBart'=VBart,
                'VKi'=VKi,'VLi'=VLi, 'V_macro_Li'=V_macro_Li,'VSt'=VSt,'VSIn'=VSIn,
                'VLIn'=VLIn,'VLn'=VLn, 'V_macro_Ln'=V_macro_Ln,
                'VSpl'=VSpl,'VH'=VH,'VBr'=VBr,'VRe'=VRe,
                'V_macro_Spl'=V_macro_Spl,
                'AKi'=AKi,'ALi'=ALi,'ASt'=ASt,'ASIn'=ASIn,'ALIn'=ALIn,
                'ALn'=ALn,'ASpl'=ASpl,'AH'=AH,'ABr'=ABr,'ARe'=ARe,
                
                'QKi'=QKi, 'QLi'=QLi, 'QRe'=QRe, 'QLn'=QLn, 'QSpl'=QSpl, 'QH'=QH,
                'QBr'=QBr, 'QSt'=QSt,'QSIn'=QSIn, 'QLIn'=QLIn,
                
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW, "np_size"=np_size,
                "np_size_small"=np_size_small, "np_size_large"=np_size_large,
                "Qtotal"=Qtotal, "QLitot"=QLitot,
                #"Qurine"=Qurine,"Qfeces"=Qfeces, 
                "sex"=sex, "estimated_params"=estimated_params
          
                
    
  ))

  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
#coef_i is the ratio of permeability coefficient (xi) between capillary blood and organ
#and partition coefficient of nanoparticles between tissue and blood (Pi), coef_i=xi/Pi
# based on the work of Li et al., 2013, https://doi.org/10.3109/17435390.2013.863406
    
    if(np_size == "np_size_small"){
      
      coef_liver <- estimated_params[1]
      coef_spleen <- estimated_params[2]
      coef_kidney <- estimated_params[3]
      CLurine <- estimated_params[4]
      
    }else{
      
      coef_liver <- estimated_params[5]
      coef_spleen <- estimated_params[6]
      coef_kidney <- estimated_params[7]
      CLurine <- estimated_params[8]
      
    }
    
    coef_heart <- estimated_params[9]
    coef_lung <- estimated_params[10]
    coef_rob <- estimated_params[11]
    coef_stomach <- estimated_params[12]
    coef_smallIn <- estimated_params[13]
    coef_largeIn <- estimated_params[14]
    coef_brain <- estimated_params[15]
    
    CLfeces <- estimated_params[16]
    
   

    # Blood concentration
    CBven <- MBven/VBven
    CBart <- MBart/VBart
    #CB <- MB/VB
    
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
    
    #Blood
    #dMB = QLn*(CLn/coef_lung) -QKi*CB -QLi*CB -QSt*CB-QSIn*CB -QLIn*CB -
          # QSpl*CB -QH*CB-QBr*CB -QRe*CB -QLn*CB +QKi*(CKi/coef_kidney) +QLi*(CLi/coef_liver) +QSt*(CSt/coef_stomach) +
          # QSIn*(CSIn/coef_smallIn) +QLIn*(CLIn/coef_largeIn) +QSpl*(CSpl/coef_spleen) +
          # QH*(CH/coef_heart) +QBr*(CBr/coef_brain) +QRe*(CRe/coef_rob)
      
    #Kidney
    #tissue subcompartment
    dMKi = QKi*(CBart - (CKi/coef_kidney)) - CLurine*MKi #Qurine*MKi
    
    #Liver
    #tissue subcompartment
    dMLi =  QLi*CBart - QLitot*(CLi/coef_liver) + QSt*(CSt/coef_stomach) + QSpl*(CSpl/coef_spleen) + 
            QSIn*(CSIn/coef_smallIn) + QLIn*(CLIn/coef_largeIn) 
    
    #Stomach
    #tissue subcompartment 
    dMSt = QSt*(CBart - (CSt/coef_stomach))
    
    #Small Intestine
    #tissue subcompartment 
    dMSIn = QSIn*(CBart - (CSIn/coef_smallIn)) 
    
    #Large Intestine
    #tissue subcompartment 
    dMLIn = QLIn*(CBart - (CLIn/coef_largeIn)) - CLfeces*CLIn #Qfeces*CLIn 
    
    #Lung 
    #tissue subcompartment
    dMLn = QLn*(CBven - (CLn/coef_lung))
    
    
    #Spleen
    #tissue subcompartment 
    dMSpl = QSpl*(CBart - (CSpl/coef_spleen)) 
    
    
    #Heart
    #tissue subcompartment 
    dMH = QH*(CBart - (CH/coef_heart))
    
    #Brain
    #Tissue subcompartment 
    dMBr = QBr*(CBart - (CBr/coef_brain))
    
    #Rest of body
    #Tissues ubcompartment 
    dMRe = QRe*(CBart - (CRe/coef_rob))
    
    # Urine
    dMurine = CLurine*MKi #Qurine*MKi
    # Feces
    dMfeces = CLfeces*CLIn #Qfeces*CLIn
    
    
    dVurine = CLurine #Qurine
    dVfeces = CLfeces #Qfeces
    
   
    
    #Concentration calculation in each compartment 
    Cven <- CBven
    Cart <- CBart
    #CB <- CB
    Cblood <- (MBven + MBart)/ (VBven + VBart)
    #Cblood <- MB/VB
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

obj.func <- function(x, dataset){
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)

  ##########################
  #-------------------------
  # Liu et al., 2012
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  preds_Liu_1_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                        "Mliver","Mspleen", "Mstomach",
                                                                                        "Mkidneys", "Mlungs", "Mbrain",
                                                                                        "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_1_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                        exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100)
  
  score[1] <- AAFE(predictions = preds_Liu_1_small_tissues, observations = obs_Liu_1_small_tissues)
  
  
  # Set up simulations for the 2nd case, Liu (2012) 1 mg/kg small_tissues_diff_time_points
  
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mheart", "Mliver","Mspleen", "Mstomach","Mkidneys", "Mlungs", "Mbrain",
                    "Msmall_intestine", "Mlarge_intestine")
  preds_Liu_1_small_diftp_tissues <- list()
  

  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_small_diftp_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  obs_Liu_1_small_diftp_tissues <- list( exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                               exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  
  score[2] <- AAFE(predictions = preds_Liu_1_small_diftp_tissues, observations = obs_Liu_1_small_diftp_tissues)
  
  
  
  # Set up simulations for the 3rd case, i.e. Liu (2012) 1 mg/kg large_tissues_diff_time_points
  
  
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_large <- 556/2
  np_size <- np_size_large #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of Liu et al. 2012 tissues large p.s., 1 mg/kg diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mheart", "Mliver","Mspleen", "Mstomach","Mkidneys", "Mlungs", "Mbrain",
                    "Msmall_intestine", "Mlarge_intestine")
  preds_Liu_1_large_diftp_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_large_diftp_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  obs_Liu_1_large_diftp_tissues <- list( exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  score[3] <- AAFE(predictions = preds_Liu_1_large_diftp_tissues, observations = obs_Liu_1_large_diftp_tissues)
  
  
  # Set up simulations for the 4th case, i.e. Liu (2012) 2 mg/kg small_tissues
 
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 2 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
   user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of Liu et al. 2012 tissues Small p.s., 2 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_Liu_2_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                  "Mliver","Mspleen", "Mstomach",
                                                                                                  "Mkidneys", "Mlungs", "Mbrain",
                                                                                                  "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_2_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Liver", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Spleen", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Stomach", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Kidneys", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Lungs", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Brain", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Small_intestine", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Large_intestine", "concentration"]*admin.dose/100)
  
  score[4] <- AAFE(predictions = preds_Liu_2_small_tissues, observations = obs_Liu_2_small_tissues)
  
  
  # # Set up simulations for the 5th case, Liu (2012) 10 mg/kg small_tissues
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of Liu et al. 2012 tissues Small p.s., 10 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_Liu_10_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                  "Mliver","Mspleen", "Mstomach",
                                                                                                  "Mkidneys", "Mlungs", "Mbrain",
                                                                                                  "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_10_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Liver", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Spleen", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Stomach", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Kidneys", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Lungs", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Brain", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Small_intestine", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Large_intestine", "concentration"]*admin.dose/100)
  
  score[5] <- AAFE(predictions = preds_Liu_10_small_tissues, observations = obs_Liu_10_small_tissues)
    
 
  # Set up simulations for the 6th case, i.e. Liu (2012) 1 mg/kg small_blood diff_time_points
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  
  preds_Liu_1_small_diftp_blood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_small_diftp_blood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  
  obs_Liu_1_small_diftp_blood <- list(exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100)
  
  score[6] <- AAFE(predictions = preds_Liu_1_small_diftp_blood, observations = obs_Liu_1_small_diftp_blood)
  
  
  # Set up simulations for the 7th case, i.e. Liu (2012) 1 mg/kg large_blood_diff_time_points
  
 
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_large <- 556/2
  np_size <- np_size_large #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df7=========================================================
  
  exp_data <- dataset$df7 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  
  preds_Liu_1_large_diftp_blood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_large_diftp_blood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  
  obs_Liu_1_large_diftp_blood <- list(exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100)
  
  score[7] <- AAFE(predictions = preds_Liu_1_large_diftp_blood, observations = obs_Liu_1_large_diftp_blood)
  
  # Estimate final score
  if (sum(is.na(score))>0){
    final_score <- 100
    
  }else{
    final_score <- mean(score)
    
  }
  return(final_score)
  
}

################################################################################


setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide")

MW <- 124.91 #g/mol
source("Goodness-of-fit-metrics.R")
# Read data
Liu_1_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues.xlsx")
Liu_1_small_diftp_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues_dif_times.xlsx")
Liu_1_large_diftp_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_large_1_tissues_dif_times.xlsx")
Liu_2_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_2_tissues.xlsx")
Liu_10_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_10_tissues.xlsx")
Liu_1_small_diftp_blood <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_blood_dif_times.xlsx")
Liu_1_large_diftp_blood <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_large_1_blood_dif_times.xlsx")


setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide/Training/AAFE/GO_model")


dataset <- list("df1" = Liu_1_small_tissues,"df2" = Liu_1_small_diftp_tissues,
                "df3" = Liu_1_large_diftp_tissues,"df4" = Liu_2_small_tissues,
                "df5" = Liu_10_small_tissues,"df6" = Liu_1_small_diftp_blood,
                "df7" = Liu_1_large_diftp_blood)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA"
              "xtol_rel" = 1e-03,
              "ftol_rel" = 0.0,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0, 
              "maxeval" = 1500, 
              "print_level" = 1)

# Create initial conditions (zero initialisation)
#Parameter names:

N_pars <- 16 # Number of parameters to be fitted
fit <-  c(rep(log(1), 16))

lb = c(rep(log(1e-10),16))
ub = c(rep(log(1e10),16))

# lb = c(rep(log(1e-3),13), rep(log(1e-6),5))
# ub = c(rep(log(1e3),13), rep(log(1e6),5))

# Run the optimization algorithm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             lb	= lb,
                             ub = ub,
                             opts = opts,
                             dataset = dataset)


estimated_params <- exp(optimizer$solution)
save.image("GO_model.RData")



# Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df1=========================================================

exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_tissues <- solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                          "Mkidneys", "Mlungs", "Mbrain",
                                          "Msmall_intestine", "Mlarge_intestine")]



# Set up simulations for the 2nd case, Liu (2012) 1 mg/kg small_tissues_diff_time_points
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
#======================================df2=========================================================

exp_data <- dataset$df2 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_diftp_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                                 "Mkidneys", "Mlungs", "Mbrain",
                                                 "Msmall_intestine", "Mlarge_intestine")]



# Set up simulations for the 3rd case, i.e. Liu (2012) 1 mg/kg large_tissues_diff_time_points

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_large <- 556/2
np_size <- np_size_large #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df3=========================================================

exp_data <- dataset$df3 # retrieve data of Liu et al. 2012 tissues large p.s., 1 mg/kg diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_large_diftp_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                                 "Mkidneys", "Mlungs", "Mbrain",
                                                 "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 4th case, i.e. Liu (2012) 2 mg/kg small_tissues

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 2 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))


#======================================df4=========================================================

exp_data <- dataset$df4 # retrieve data of Liu et al. 2012 tissues Small p.s., 2 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
preds_Liu_2_small_tissues <- solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                          "Mkidneys", "Mlungs", "Mbrain",
                                          "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 5th case, Liu (2012) 10 mg/kg small_tissues
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

#======================================df5=========================================================

exp_data <- dataset$df5 # retrieve data of Liu et al. 2012 tissues Small p.s., 10 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
preds_Liu_10_small_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                            "Mkidneys", "Mlungs", "Mbrain",
                                            "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 6th case, i.e. Liu (2012) 1 mg/kg small_blood diff_time_points

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #h

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df6=========================================================

exp_data <- dataset$df6 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_diftp_blood <- solution[, c("time", "Mblood")]


# Set up simulations for the 7th case, i.e. Liu (2012) 1 mg/kg large_blood_diff_time_points

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_large <- 556/2
np_size <- np_size_large #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df7=========================================================

exp_data <- dataset$df7 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_large_diftp_blood <- solution[, c("time", "Mblood")]



# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("GO mass (" * mu* "g)" ),
         x = "Time (min)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("", values=cls,
                       guide = guide_legend(override.aes =
                                              list(shape = c(16,NA),
                                                   linetype = c(0,1))))+
    theme_light() + 
    theme(legend.position=c(1,1), 
          legend.justification=c(0, 1), 
          legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          axis.title=element_text(size=14),
          legend.text = element_text(size=14)
    )
  
}


# Convert Liu 2012, male small_1_tissues from long to wide format using reshape
experiment1 <- reshape(Liu_1_small_tissues[c("Tissue" ,"Time_min", 
                                        "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",Liu_1_small_tissues$Tissue )

# Convert Liu 2012, male small_1_tissues_dif_times from long to wide format using reshape
experiment2 <- reshape(Liu_1_small_diftp_tissues[c("Tissue" ,"Time_min", 
                                        "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",unique(Liu_1_small_diftp_tissues$Tissue))

# Convert Liu 2012, male large_1_tissues_dif_times from long to wide format using reshape
experiment3 <- reshape(Liu_1_large_diftp_tissues[c("Tissue" ,"Time_min", 
                                        "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment3) <- c("Time",unique(Liu_1_large_diftp_tissues$Tissue))

# Convert Liu 2012, male small_2_tissues from long to wide format using reshape
experiment4 <- reshape(Liu_2_small_tissues[c("Tissue" ,"Time_min", 
                                                   "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment4) <- c("Time",Liu_2_small_tissues$Tissue)

# Convert Liu 2012, male small_10_tissues from long to wide format using reshape
experiment5 <- reshape(Liu_10_small_tissues[c("Tissue" ,"Time_min", 
                                             "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",Liu_10_small_tissues$Tissue)

# Convert Liu 2012, male small_1_blood_dif_times from long to wide format using reshape
experiment6 <- reshape(Liu_1_small_diftp_blood[c("Tissue" ,"Time_min", 
                                              "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",unique(Liu_1_small_diftp_blood$Tissue))

# Convert Liu 2012, male large_1_blood_dif_times from long to wide format using reshape
experiment7 <- reshape(Liu_1_large_diftp_blood[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(Liu_1_large_diftp_blood$Tissue))



# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment3,
                    experiment4 = experiment4, experiment5 = experiment5, experiment6 = experiment6,
                    experiment7 = experiment7)


# Rename predictions so that they share the same name as the names of the experimental data dataframe
colnames(preds_Liu_1_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                          "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine")
colnames(preds_Liu_1_small_diftp_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                                "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_1_large_diftp_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                          "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_2_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                    "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_10_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                    "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_1_small_diftp_blood) <- c("Time", "Blood")
colnames(preds_Liu_1_large_diftp_blood) <- c("Time", "Blood") 


# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_Liu_1_small_tissues, predictions2 = preds_Liu_1_small_diftp_tissues,
                    predictions3 = preds_Liu_1_large_diftp_tissues, predictions4 = preds_Liu_2_small_tissues,
                    predictions5 = preds_Liu_10_small_tissues, predictions6 = preds_Liu_1_small_diftp_blood,
                    predictions7 = preds_Liu_1_large_diftp_blood)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment )
  })
  if(length(compartments) == 1){
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                               common.legend = TRUE, legend = "right"))
    
  }else{
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                               common.legend = TRUE, legend = "right"))
  }
  
  
  plot.margin=unit(c(0,0,0,0), "pt")
  
  
  # Save the plot with dynamically adjusted dimensions
  ggsave(paste0("experiment", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}


