setwd('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK')
main_func <- function(substance){
  create.params <- function(user.input){
    # Transform input temperature into Kelvin scale
    Texp <- 273 + Texp # K
    
    # Cardiac output reference value at T = 6 C (Barron et al. 1987, Table II)
    F_card_ref_6 <- 1.188 # ml/h/g 
    # Cardiac output reference value at T = 12 C (Barron et al. 1987, Table II)
    F_card_ref_12 <- 2.322 # ml/h/g 
    # Cardiac output reference value at T = 18 C (Barron et al. 1987, Table II)
    F_card_ref_18 <- 3.75 # ml/h/g 
    
    if(Texp <= 273+6){
      F_card_ref <- F_card_ref_6
    }else if (Texp >= 273+18){
      F_card_ref <- F_card_ref_18
    }else if(Texp >= 273+6 & Texp <= 273+12){
      F_card_ref <- approx(x=c(273+6, 273+12), y=c(F_card_ref_6, F_card_ref_12), xout = Texp)$y
    }else if(Texp >= 273+12 & Texp <= 273+18){
      F_card_ref <- approx(x=c(273+12, 273+18), y=c(F_card_ref_12, F_card_ref_18), xout = Texp)$y
    }
    
    # Body weight reference value at T = 6 C (Barron et al. 1987, Table II)
    BW_ref_6 <- 270.1 # g 
    # Body weight reference value at T = 12 C (Barron et al. 1987, Table II)
    BW_ref_12 <- 296.4 # g 
    # Body weight reference value at T = 18 C (Barron et al. 1987, Table II)
    BW_ref_18 <- 414.5 # g
    
    if(Texp <= 273+6){
      BW_ref <- BW_ref_6
    }else if (Texp >= 273+18){
      BW_ref <- BW_ref_18
    }else if(Texp >= 273+6 & Texp <= 273+12){
      BW_ref <- approx(x=c(273+6, 273+12), y=c(BW_ref_6, BW_ref_12), xout = Texp)$y
    }else if(Texp >= 273+12 & Texp <= 273+18){
      BW_ref <- approx(x=c(273+12, 273+18), y=c(BW_ref_12, BW_ref_18), xout = Texp)$y
    }
    
    # Arrhenius Temperature function 
    TA <- 6930 # Arrhenius Temperature K - Grech et al.2018
    Tref <- 273 + c(6,12,18) # Reference Temperature K - Grech et al.2018
    Tr <- Tref[which.min(abs(Tref - Texp))]
    KT <- exp(TA/Tr - TA/Texp)
    
    
    
    # Load the xlsx file with the physiological params pf rainbow trout
    phys.params <- openxlsx::read.xlsx('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 1)
    
    # Keep only the physiological parameters from the paper of Vidal et al. 2019
    # fw are the fractions of tissue_weight/total_weight (unitless)
    fw <- phys.params[phys.params$Source=='Vidal et al. 2019', c('Liver', 'Blood', 'Skin',
                                                                 'Muscle', 'Gills', 'Kidney', 'Viscera')]
    fw_Liver <- fw$Liver
    fw_Blood <- fw$Blood
    fw_Skin <- fw$Skin
    fw_Muscle <- fw$Muscle
    fw_Gills <- fw$Gills
    fw_Kidney <- fw$Kidney
    fw_Viscera <- fw$Viscera
    fw_lumen <- 0.012
    
    # Load the xlsx file with the physiological params pf rainbow trout
    phys.params <- openxlsx::read.xlsx('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 2)
    
    # Keep only the physiological parameters from the paper of Vidal et al. 2019
    # fb are the fractions of blood flow of each tissue (unitless)
    fb <- phys.params[phys.params$Source=='Vidal et al. 2019', c('Liver', 'Skin', 'Muscle',
                                                                 'Gills', 'Kidney', 'Viscera')]
    
    fb_Liver <- fb$Liver
    fb_Skin <- fb$Skin
    fb_Muscle <- fb$Muscle
    fb_Gills <- fb$Gills
    fb_Kidney <- fb$Kidney
    fb_Viscera <- fb$Viscera
    
    # Cl_bile <- 2.16e-06 # L/h 
    # Cl_feces <- 2.62e-07 # L/h 
    # Cl_urine <- 7e-05 # L/h 
    Ku <- 0.13 # 1/h
    # P_liver <- 2.05
    # P_muscle <- 0.15
    # P_kidney <- 0.58
    # P_viscera <- 0.87
    # P_brain <- 0.63
    # P_skin <- 0.25
    # P_gills <- 0.40
    # P_carcass <- 0.00129             
    Free = 3.2e-02
    
    a_skin <- 0.9 # 90% of venous blood of skinwas assumed to flow directly to kidney (Nischols et al.1996)
    a_muscle <- 0.6 # 90% of venous blood of skinwas assumed to flow directly to kidney (Nischols et al.1996)
    
    P_viscera = 0.87
    P_brain = 0.63
    
    plasma <- 0.7
    
    return(list('F_card_ref'=F_card_ref, 'BW_ref'=BW_ref, 'KT'=KT, 
                'admin.dose_dietary'=admin.dose_dietary, 
                'admin.time_dietary'=admin.time_dietary,
                
                'fw_Liver'=fw_Liver, 'fw_Blood'=fw_Blood, 'fw_Skin'=fw_Skin, 
                'fw_Muscle'=fw_Muscle, 'fw_Gills'=fw_Gills, 'fw_Kidney'=fw_Kidney,
                'fw_Viscera'=fw_Viscera, 'fw_lumen'=fw_lumen, 
                
                'fb_Liver'=fb_Liver, 'fb_Skin'=fb_Skin, 'fb_Muscle'=fb_Muscle,
                'fb_Gills'=fb_Gills, 'fb_Kidney'=fb_Kidney, 'fb_Viscera'=fb_Viscera,
                
                #'Cl_bile'=Cl_bile, 'Cl_feces'=Cl_feces, 'Cl_urine'=Cl_urine,
                'Ku'=Ku, 'Free'=Free, 'a_skin'=a_skin, 'a_muscle'=a_muscle,
                'P_viscera'=P_viscera, 'P_brain'=P_brain,
                'plasma'=plasma))
  }
  
  create.inits <- function(parameters){
    with(as.list(parameters),{
      M_art<-0; M_venous<-0;
      M_gills<-0; M_lumen=0; M_bile<-0; M_viscera<-0; M_liver<-0; M_kidney<-0;
      M_muscle<-0; M_skin<-0; M_carcass<-0; M_urine<-0; M_feces<-0; M_input<-0 
      
      return(c('M_art'=M_art, 'M_venous'=M_venous, 'M_gills'=M_gills, 'M_lumen'=M_lumen,
               'M_bile'=M_bile, 'M_viscera'=M_viscera, 'M_liver'=M_liver, 'M_kidney'=M_kidney, 
               'M_muscle'=M_muscle, 'M_skin'=M_skin, 'M_carcass'=M_carcass,
               'M_urine'=M_urine, 'M_feces'=M_feces, 'M_input'=M_input))
    })
  }
  
  create.events <- function(parameters){
    with(as.list(parameters),{
      
      # Calculate number of administrated doses and corresponding administration time
      ldose_dietary <- length(admin.dose_dietary)
      ltimes_dietary <- length(admin.time_dietary)
      
      # If not equal, then stop 
      if (ltimes_dietary != ldose_dietary){
        stop("The times of administration should be equal in number to the doses")
      }else{
        events <- data.frame(var = c(rep(c('M_lumen', 'M_input'), ltimes_dietary)), 
                             time = sort(rep(admin.time_dietary,2)),
                             value = rep(admin.dose_dietary,each=2),
                             method = 'add')
      }
      return(list(data=events))
    })
  }
  
  
  
  ode.func <- function(time, inits, params){
    with(as.list(c(inits, params)),{
      
      
      # Body weight - g
      #BW <- 1000*a_bio*L^b_bio
      BW <- fish_weight(time)
      
      # Total cardiac output ml/h
      Q_total <- F_card_ref*KT*(BW/BW_ref)^(-0.1)*BW*plasma  
      
      # Calculate the mass of each tissue - g
      w_blood <- fw_Blood*BW*plasma     # Blood mass - g
      w_liver <- fw_Liver*BW     # Liver mass - g
      w_skin <- fw_Skin*BW       # Skin weight - g
      w_muscle <- fw_Muscle*BW   # Muscle weight - g
      w_gills <- fw_Gills*BW     # Gills weight - g
      w_kidney <- fw_Kidney*BW   # Kidney weight - g
      w_viscera <- fw_Viscera*BW # Viscera weight - g
      w_lumen <- fw_lumen*BW
      w_bile <- w_lumen
      w_art <- 1/3*w_blood
      w_venous <- 2/3*w_blood
      w_carcass <- BW - (w_blood + w_liver + w_skin + w_muscle +
                           w_gills + w_kidney + w_viscera + w_lumen + w_bile)
      
      # Calculate the regional blood flows - ml/h
      Q_liver <- fb_Liver*BW     # Liver blood flow - ml/h
      Q_skin <- fb_Skin*BW       # Skin blood flow - ml/h
      Q_muscle <- fb_Muscle*BW   # Muscle blood flow - ml/h
      Q_gills <- Q_total #fb_Gills*BW     # Gills blood flow - ml/h
      Q_kidney <- fb_Kidney*BW   # Kidney blood flow - ml/h
      Q_viscera <- fb_Viscera*BW # Viscera blood flow - ml/h
      Q_carcass <- Q_total -(Q_liver + Q_skin + Q_muscle + #Q_gills + 
                               Q_kidney + Q_viscera)
      
      # Tissue concentrations ug PFAS/g tissue
      C_gills <- M_gills/w_gills
      C_bile <- M_bile/w_bile
      C_viscera <- M_viscera/w_viscera
      C_liver <- M_liver/w_liver
      C_kidney <- M_kidney/w_kidney
      C_muscle <- M_muscle/w_muscle 
      C_skin <- M_skin/w_skin
      C_carcass <- M_carcass/w_carcass
      C_lumen <- M_lumen/w_lumen
      C_art <- M_art/w_art
      C_venous <- M_venous/w_venous
      C_blood <- (M_art + M_venous)/w_blood
      
      # Arterial Blood
      dM_art <- Free*Q_gills*C_gills/P_gills - 
        (Q_viscera + Q_liver + Q_kidney +
           Q_muscle + Q_skin + Q_carcass)*Free*C_art
      
      dM_venous <- - Free*Q_total*C_venous + ((Q_liver + Q_viscera)*C_liver/P_liver + 
                                                Q_kidney*C_kidney/P_kidney + (1-a_muscle)*Q_muscle*C_muscle/P_muscle +
                                                (1-a_skin)*Q_skin*C_skin/P_skin + Q_carcass*C_carcass/P_carcass)*Free
      
      # Gills
      dM_gills <- Q_gills*Free*(C_venous - C_gills/P_gills)
      
      dM_input=0
      
      # Viscera lumen - Available PFAS for absorption and elimination
      dM_lumen = - Ku*M_lumen - Cl_feces*C_lumen 
      
      # Bile - Billiary eliminated PFAS (Considering no reabsoprtion)
      dM_bile = Free*Cl_bile*C_liver - Cl_feces*C_bile
      
      # Viscera tissue
      dM_viscera <- Q_viscera*Free*(C_art - C_viscera/P_viscera) + Ku*M_lumen 
      
      # Liver
      dM_Liver <- Q_liver*Free*C_art + Q_viscera*Free*C_viscera/P_viscera - 
        (Q_liver + Q_viscera)*Free*C_liver/P_liver - Free*Cl_bile*C_liver
      
      # Kidney
      dM_kidney <- Q_kidney*Free*(C_art - C_kidney/P_kidney) + 
        a_muscle*Q_muscle*Free*C_muscle/P_muscle + 
        a_skin*Q_skin*Free*C_skin/P_skin - Cl_urine*C_kidney
      
      # Muscle
      dM_muscle <- Q_muscle*Free*(C_art - C_muscle/P_muscle)
      
      # Skin
      dM_skin <- Q_skin*Free*(C_art - C_skin/P_skin)
      
      # Carcass 
      dM_carcass <- Q_carcass*Free*(C_art - C_carcass/P_carcass)
      
      # Urine
      dM_urine <- Cl_urine*C_kidney
      
      # Feces
      dM_feces <- Cl_feces*C_lumen + Cl_feces*C_bile
      
      Mass_balance <- M_input - (M_art + M_venous + M_gills + M_lumen + M_bile + 
                                   M_viscera + M_liver + M_kidney + M_muscle + 
                                   M_skin + M_carcass + M_urine + M_feces)
      
      return(list(c('dM_art'=dM_art, 'dM_venous'=dM_venous, 
                    'dM_gills'=dM_gills, 'dM_lumen'=dM_lumen, 'dM_bile'=dM_bile,
                    'dM_viscera'=dM_viscera, 'dM_Liver'=dM_Liver, 
                    'dM_kidney'=dM_kidney, 'dM_muscle'=dM_muscle,
                    'dM_skin'=dM_skin, 'dM_carcass'=dM_carcass,
                    'dM_urine'=dM_urine, 'dM_feces'=dM_feces, 'dM_input'=dM_input),
                  'C_Gills'=C_gills, 'C_Bile'=C_bile, 'C_Viscera'=C_viscera,
                  'C_Liver'=C_liver, 'C_Kidney'=C_kidney, 'C_Muscle'=C_muscle,
                  'C_Skin'=C_skin, 'C_Carcass'=C_carcass, 'C_Lumen'=C_lumen,
                  'C_Blood'=C_blood,
                  'Mass_balance'=Mass_balance, 'BW'=BW))
    })
  }
  
  # Experimental data from Falk et al.2015
  #---------------------------------------
  # The concentrations in the data are given in ug PFAS/kg tissue units. 
  # The time is given in days and will be transformed in hours, to be compatible
  # with the model
  
  # Directory of folder with saved data files
  data_dir <- 'C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015'
  
  # Load PFOS data
  #---------------
  PFOS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFOS.xlsx'))
  PFOS_data$Time <- PFOS_data$Time*24
  
  # Load PFOA data
  #---------------
  PFOA_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFOA.xlsx'))
  PFOA_data$Time <- PFOA_data$Time*24
  
  # Load PFBS data
  #---------------
  PFBS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFBS.xlsx'))
  PFBS_data$Time <- PFBS_data$Time*24
  
  # Load PFHxS data
  #----------------
  PFHxS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFHxS.xlsx'))
  PFHxS_data$Time <- PFHxS_data$Time*24
  
  # Load PFNA data
  #----------------
  PFNA_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFNA.xlsx'))
  PFNA_data$Time <- PFNA_data$Time*24
  
  # Put all data in a list
  data_list <- list('PFOS'=PFOS_data, 'PFOA'=PFOA_data,
                    'PFBS'=PFBS_data, 'PFHxS'=PFHxS_data,
                    'PFNA'=PFNA_data)
  
  # Absolute average fold error (AAFE)
  #--------------------------------
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
  
  rmse <- function(observed, predicted){
    y_obs <- unlist(observed)
    y_pred <- unlist(predicted)
    return(sqrt(mean((y_obs-y_pred)^2)) )
  }
  
  fitness.metric <- function(observed, predicted, comp.names =NULL){
    # Check if the user provided the correct input format
    if (!is.list(observed) || !is.list(predicted)){
      stop(" The observations and predictions must be lists")
    }
    # Check if the user provided equal length lists
    if (length(observed) != length(predicted)){
      stop(" The observations and predictions must have the same compartments")
    }
    Ncomp <- length(observed) # Number of compartments
    I <- rep(NA, Ncomp) # Compartment discrepancy index
    N_obs <- rep(NA, Ncomp) #Number of observations per compartment
    #loop over the compartments
    for (i in 1:Ncomp){
      Et <- 0 #relative error with observations
      St <- 0  #relative error with simulations
      N <- length(observed[[i]]) # number of observations for compartment i
      # Check if observations and predictions have equal length
      if(N != length(predicted[[i]])){
        stop(paste0("Compartment ",i," had different length in the observations and predictions"))
      }
      N_obs[i] <- N # populate the N_obs vector
      for (j in 1:N){
        # sum of relative squared errors (error = observed - predicted)
        Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
        St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
      }
      
      # root mean of the square of observed values
      RMEt <- sqrt(Et/N)
      # root mean of the square of simulated values
      RMSt <- sqrt( St/N)
      
      I[i] <- (RMEt + RMSt)/2   
    }
    # Total number of observations
    Ntot <- sum(N_obs)
    # Initialise the consolidated discrepancy index
    Ic <-0
    for (i in 1:Ncomp){
      # Give weight to compartments with more observations (more information)
      Ic <- Ic +  I[i]* N_obs[i]/Ntot
    }
    # Name the list of compartment discrepancy indices
    if ( !is.null(comp.names)){
      names(I) <- comp.names
    }else if (!is.null(names(observed))){
      names(I) <- names(observed)
    } else if (!is.null(names(predicted)) && is.null(comp.names) ){
      names(I) <- names(predicted)
    }
    return(Ic)
    #return(list(Total_index = Ic, Compartment_index= I))
  }
  
  # Objective function
  #-------------------
  
  obj.func <- function(x, user.input, substance){
    names(x) <- c('P_liver', 'P_muscle', 'P_kidney', 
                  'P_skin', 'P_gills', 'P_carcass',
                  'Cl_bile', 'Cl_feces', 'Cl_urine')
    
    # lets try just for pfos
    exp_data <- data_list[[substance]]
    
    params <- create.params(user.input)
    # inits <- create.inits(params)
    # events <- create.events(params)
    
    solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = c(x,params),
                                        events = events,
                                        method="lsodes",rtol = 1e-05, atol = 1e-05))
    
    # Keep the predictions only for the time points at which there are available data
    predictions_df <- solution[solution$time %in% exp_data$Time,c('time', 'C_Liver', 'C_Blood', 
                                                                  'C_Skin', 'C_Muscle', 'C_Gills',
                                                                  'C_Kidney', 'C_Carcass')]
    # Transform predictions and exp_data in list data type
    observations <- list()
    predictions <- list()
    for (i in 2:dim(exp_data)[2]) {
      observations[[i-1]] <- exp_data[,i]
      names(observations)[i-1] <- colnames(exp_data)[i] 
      predictions[[i-1]] <- predictions_df[,i]*1000 # transform from ug/g to ug/kg
      names(predictions)[i-1] <- colnames(predictions_df)[i]
    }
    
    score <- fitness.metric(observations,predictions)
    
    return(score)
  }
  
  
  ################################################################
  #x0 <- c(1.17613431, 0.08398678, 0.07102088, 0.33793084, 0.39039634, 0.23137421, 0.78637215, 0.63154606, 0.70663080)
  x0 <- runif(9)
  names(x0) <- c('P_liver', 'P_muscle', 'P_kidney', 
                 'P_skin', 'P_gills', 'P_carcass',
                 'Cl_bile', 'Cl_feces', 'Cl_urine')
  # Total body weight of fish
  Texp <- 15 # C
  
  # Fish feeding
  # The fish were fed once per day. The daily food intake was 2.6% of the mean body weight.
  # The amount of food was  consumed immediately and completely. The uptake phase
  # lasted 28 days. Nominal concentration of PFAS in spiked food is 500 ug/kg.
  
  # We estimate the mean body weight of the fish from day 0 till day 28 in order to 
  # estimate the added mount of food everyday.
  
  fish_weight <- function(time){
    x <- c(0,28,56)*24 
    y <- c(314, 655, 808)
    
    if(time <= x[1]){
      w = y[1]
    }else if(time >= x[3]){
      w = y[3]
    }else if(time >= x[1] & time < x[2]){
      w = approx(x=x[1:2],y=y[1:2], xout = time)$y
    }else if(time >= x[2] & time < x[3]){
      w = approx(x=x[2:3],y=y[2:3], xout = time)$y
    }
    return(w)
  }
  
  # Time points of added food
  admin.time_dietary <- seq(0,27*24,24)
  # Calculate fish weight over time (g)
  fish_weights <- unlist(lapply(admin.time_dietary, fish_weight))
  # Multiply fish_weights * g daily_food_intake/g of BW * Concentration (ug/g of food)
  admin.dose_dietary <- fish_weights*2.6/100*500/1000
  
  user.input <- list('Texp'=Texp,
                     'admin.dose_dietary'=admin.dose_dietary,
                     'admin.time_dietary'=admin.time_dietary)
  
  params <- create.params(user.input)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- seq(0,56*24,1)
  
  N_iter <- 1000
  opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",  #"NLOPT_LN_SBPLX" ,
                "xtol_rel" = 0.0,
                "ftol_rel" = 0.0,
                "ftol_abs" = 0.0,
                "xtol_abs" = 0.0 ,
                "maxeval" = N_iter,
                "print_level" = 1 )
  
  optimization <- nloptr::nloptr( x0 = x0,
                                  eval_f = obj.func,
                                  lb	= rep(1e-08, length(x0)),
                                  ub = rep(20, length(x0)),
                                  opts = opts,
                                  user.input=user.input,
                                  substance = substance)
  
  
  
  
  # Plot Concentration - Time profiles
  #------------------------------------
  
  library(ggplot2)
  
  x_opt <- optimization$solution
  names(x_opt) <- names(x0)
  
  user.input <- list('Texp'=Texp,
                     'admin.dose_dietary'=admin.dose_dietary,
                     'admin.time_dietary'=admin.time_dietary)
  
  params <- c(create.params(user.input), x_opt)
  inits <- create.inits(params)
  events <- create.events(params)
  
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = c(x_opt, params),
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # Keep the predictions only for the time points at which there are available data
  predictions_df <- solution[,c('time' , 'C_Liver', 'C_Blood',
                                'C_Skin', 'C_Muscle', 'C_Gills',
                                'C_Kidney', 'C_Carcass')]
  predictions_df[,-1] <- predictions_df[,-1]*1000
  exp_data <- data_list[[substance]]
  
  compartments <- colnames(exp_data)[2:8]
  color_codes <- scales::hue_pal()(length(compartments))
  names(color_codes) <-  colnames(exp_data)[2:8]
  
  plot <- ggplot()+
    geom_line(data = predictions_df, aes(x = time, y = C_Liver, color='Liver'), size=1.3)+
    geom_line(data = predictions_df, aes(x = time, y = C_Blood, color='Blood'), size=1.3)+
    geom_line(data = predictions_df, aes(x = time, y = C_Skin, color='Skin'), size=1.3)+
    geom_line(data = predictions_df, aes(x = time, y = C_Muscle, color='Muscle'), size=1.3)+
    geom_line(data = predictions_df, aes(x = time, y = C_Gills, color='Gills'), size=1.3)+
    geom_line(data = predictions_df, aes(x = time, y = C_Kidney, color='Kidney'), size=1.3)+
    geom_line(data = predictions_df, aes(x = time, y = C_Carcass, color='Carcass'), size=1.3)+
    
    geom_point(data = exp_data, aes(x = Time, y = Liver, color = 'Liver'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Blood, color = 'Blood'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Skin, color = 'Skin'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Muscle, color = 'Muscle'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Gills, color = 'Gills'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Kidney, color = 'Kidney'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Carcass, color = 'Carcass'), size=5)+
    #scale_y_log10(limits = c(1, 600))+
    #ylim(c(1, 600))+
    geom_vline(xintercept = 28*24, size=1.0)+ # end of uptake
    
    
    
    labs(title = paste0("Tissues Predicted vs Observed Concentrations of ", substance),
         y = 'Concentration (ug/kg)' , x = "Time (hours)")+
    theme(plot.title = element_text(hjust = 0.5,size=30), 
          axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.y=element_text(size=22),
          axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.x=element_text(size=22),
          legend.title=element_text(hjust = 0.5,size=25), 
          legend.text=element_text(size=22),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0)) + 
    
    scale_color_manual("Tissues", values=color_codes)+
    theme(legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          axis.text = element_text(size = 14))
  
  #plot
  
  return(list('substance'=substance,
              'optimizations'=optimization,
              'plot'=plot))
}

#===============================================================================
substances <- c("PFOS",  "PFOA",  "PFBS",  "PFHxS", "PFNA" )
library(parallel)
cores <- detectCores()-2

start.time <- Sys.time()
cl = makeCluster(cores)
#clusterExport(cl=cl)
output <- parLapply(cl, substances, main_func)
stopCluster(cl)
total.duration <- Sys.time() - start.time
print(total.duration)

params_names <- c('P_liver', 'P_muscle', 'P_kidney', 
                  'P_skin', 'P_gills', 'P_carcass',
                  'Cl_bile', 'Cl_feces', 'Cl_urine')
optimized_params <- data.frame(matrix(NA, nrow = 5, ncol = length(params_names)))
rownames(optimized_params) <- substances
colnames(optimized_params) <- params_names
for (i in 1:5) {
  optimized_params[i,] <- output[[i]]$optimizations$solution
}
