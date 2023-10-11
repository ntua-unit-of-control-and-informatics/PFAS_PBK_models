create.params <- function(user.input){
  with(as.list(user.input),{
    # Transform input temperature into Kelvin scale
    Texp <- 273 + Texp # K
    
    Tref <- 273 + c(6,12,18) # Reference Temperature K - Grech et al.2018
    keep_ref_value <- which.min(abs(Tref - Texp))
    
    # Cardiac output reference value at T = 6 C (Barron et al. 1987, Table II)
    F_card_ref_6 <- 1.188 # ml/h/g 
    # Cardiac output reference value at T = 12 C (Barron et al. 1987, Table II)
    F_card_ref_12 <- 2.322 # ml/h/g 
    # Cardiac output reference value at T = 18 C (Barron et al. 1987, Table II)
    F_card_ref_18 <- 3.75 # ml/h/g 
    F_card_ref_values <- c(F_card_ref_6,F_card_ref_12,F_card_ref_18)
    F_card_ref <- F_card_ref_values[keep_ref_value]
    
    # Body weight reference value at T = 6 C (Barron et al. 1987, Table II)
    BW_ref_6 <- 270.1 # g 
    # Body weight reference value at T = 12 C (Barron et al. 1987, Table II)
    BW_ref_12 <- 296.4 # g 
    # Body weight reference value at T = 18 C (Barron et al. 1987, Table II)
    BW_ref_18 <- 414.5 # g
    BW_ref_values <- c(BW_ref_6,BW_ref_12,BW_ref_18)
    BW_ref <- BW_ref_values[keep_ref_value]
    
    # Arrhenius Temperature function 
    TA <- 6930 # Arrhenius Temperature K - Grech et al.2018
    Tr <- Tref[which.min(abs(Tref - Texp))]
    KT <- exp(TA/Tr - TA/Texp)
    
    # Load the xlsx file with the physiological params pf rainbow trout
    phys.params <- openxlsx::read.xlsx('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 1)
    
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
    phys.params <- openxlsx::read.xlsx('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 2)
    
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
    
    #Ku <- 0.13 # 1/h
    
    # Reabsorption coefficients from bile to intestine
    # estimated by Cao et al., 2022
    # K_urine = Cl_urine/f_reab_urine estimated by Ng et al., 2013 (unitless)
    if(substance=='PFOA'){
      a <- 0.138 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.30 # Cao et al., 2022
      K_urine <- 2.08 
      Cl_urine <- 0.029*3600 # 1/h (Sun et al., 2022)
    }else if(substance=='PFNA'){
      a <- 0.522 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.34 # Cao et al., 2022
      K_urine <- 1.35
      Cl_urine <- 0.050*3600 # 1/h (Sun et al., 2022)
    }else if(substance=='PFBS'){
      a <- 0.0598 # Goeritz et al.2013
      f_reab_hep <- 0.23 # Cao et al., 2022
      K_urine <- 5.88
      Cl_urine <- 0.023*3600 # 1/h (Sun et al., 2022) # Assumed equal to PFHxS
    }else if(substance=='PFHxS'){
      a <- 0.558 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.30 # Cao et al., 2022
      K_urine <- 5.88
      Cl_urine <- 0.023*3600 # 1/h (Sun et al., 2022)
    }else if(substance=='PFOS'){
      a <- 0.721 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.42 # Cao et al., 2022
      K_urine <- 1.35
      Cl_urine <- 0.050*3600 # 1/h (Sun et al., 2022)
    }
    
    # Bile flow coefficient
    Q_bile_coef <- 7.5e-05 # ml/g BW/h Grosell et al., 2000
    Q_urine_coef <- 2.755e-03 # ml/h/g of BW Urinary flow rate (Curtis et al., 1981)
    V_urine_coef <- 2.2e-03 # ml/g of BW Urine volume inside urinary bladder (Curtis et al., 1981)
    
    a_skin <- 0.9 # 90% of venous blood of skin was assumed to flow directly to kidney (Nichols et al.1996)
    a_muscle <- 0.6 # 60% of venous blood of muscle was assumed to flow directly to kidney (Nichols et al.1996)
    
    plasma <- 0.7
    
    return(list('F_card_ref'=F_card_ref, 'BW_ref'=BW_ref, 'KT'=KT, 
                
                'fw_Liver'=fw_Liver, 'fw_Blood'=fw_Blood, 'fw_Skin'=fw_Skin, 
                'fw_Muscle'=fw_Muscle, 'fw_Gills'=fw_Gills, 'fw_Kidney'=fw_Kidney,
                'fw_Viscera'=fw_Viscera, 'fw_lumen'=fw_lumen, 
                
                'fb_Liver'=fb_Liver, 'fb_Skin'=fb_Skin, 'fb_Muscle'=fb_Muscle,
                'fb_Gills'=fb_Gills, 'fb_Kidney'=fb_Kidney, 'fb_Viscera'=fb_Viscera,
                
                'a_skin'=a_skin, 'a_muscle'=a_muscle,
                'Q_bile_coef'=Q_bile_coef,
                'Q_urine_coef'=Q_urine_coef, 'V_urine_coef'=V_urine_coef,
                'K_urine'=K_urine, 'Cl_urine'=Cl_urine,
                'f_reab_hep'=f_reab_hep, 'plasma'=plasma,"a"=a))
  })
}
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


user.input <- list('substance'='PFOS',
                   'Texp'=15)

time_specific_params <- function(time, params){
  with(as.list(params),{
    BW = fish_weight(time)
    # Calculate the mass of each tissue - g
    w_blood <- fw_Blood*BW*plasma     # Blood mass - g
    w_liver <- fw_Liver*BW     # Liver mass - g
    w_skin <- fw_Skin*BW       # Skin weight - g
    w_muscle <- fw_Muscle*BW   # Muscle weight - g
    w_gills <- fw_Gills*BW     # Gills weight - g
    w_kidney <- fw_Kidney*BW   # Kidney weight - g
    w_viscera <- fw_Viscera*BW # Viscera weight - g
    w_lumen <- fw_lumen*BW
    w_art <- 1/3*w_blood
    w_venous <- 2/3*w_blood
    w_carcass <- BW - (w_blood/plasma + w_liver + w_skin + w_muscle +
                         w_gills + w_kidney + w_viscera + w_lumen)
    
    return(list('w_blood'=w_blood, 'w_liver'=w_liver, 'w_skin'=w_skin,
                'w_muscle'=w_muscle, 'w_gills'=w_gills, 'w_kidney'=w_kidney,
                'w_viscera'=w_viscera, 'w_lumen'=w_lumen, 'w_art'=w_art,
                'w_venous'=w_venous, 'w_carcass'=w_carcass))
  })
}
params <- create.params(user.input)
time_168 = time_specific_params(168, params)
time_336 = time_specific_params(336, params)
time_672 = time_specific_params(672, params)
time_744 = time_specific_params(744, params)
time_840 = time_specific_params(840, params)
time_1008 = time_specific_params(1008, params)
time_1344 = time_specific_params(1344, params)


estimate_scores <- function(id){
  # Estimate metrics for the best PINN 
  model_id <- paste0('w_data_', id)
  # Load the predictions of the PINN on experimental time points
  pinn_preds = read.csv(file = paste0('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/Hyperparameters_tuning/Loss_weights_tuning/', model_id, '/model_w', id, '_preds.csv' ))
  # Keep only the time points and columns of interest
  pinn_preds <- pinn_preds[1:7,c('Time', 'M_liver', 'M_blood', 'M_skin', 'M_muscle', 'M_gills', 
                                 'M_kidney', 'M_carcass')]
  
  # Tranform the predicted variables from mass into concentrations (ug -> ug PFOS/kg tissue)
  pinn_preds_conc <- pinn_preds
  pinn_preds_conc[1, 2:8] = pinn_preds[1, 2:8]/c(time_168$w_liver, time_168$w_blood, time_168$w_skin, time_168$w_muscle, time_168$w_gills, time_168$w_kidney, time_168$w_carcass)*1000
  pinn_preds_conc[2, 2:8] = pinn_preds[2, 2:8]/c(time_336$w_liver, time_336$w_blood, time_336$w_skin, time_336$w_muscle, time_336$w_gills, time_336$w_kidney, time_336$w_carcass)*1000
  pinn_preds_conc[3, 2:8] = pinn_preds[3, 2:8]/c(time_672$w_liver, time_672$w_blood, time_672$w_skin, time_672$w_muscle, time_672$w_gills, time_672$w_kidney, time_672$w_carcass)*1000
  pinn_preds_conc[4, 2:8] = pinn_preds[4, 2:8]/c(time_744$w_liver, time_744$w_blood, time_744$w_skin, time_744$w_muscle, time_744$w_gills, time_744$w_kidney, time_744$w_carcass)*1000
  pinn_preds_conc[5, 2:8] = pinn_preds[5, 2:8]/c(time_840$w_liver, time_840$w_blood, time_840$w_skin, time_840$w_muscle, time_840$w_gills, time_840$w_kidney, time_840$w_carcass)*1000
  pinn_preds_conc[6, 2:8] = pinn_preds[6, 2:8]/c(time_1008$w_liver, time_1008$w_blood, time_1008$w_skin, time_1008$w_muscle, time_1008$w_gills, time_1008$w_kidney, time_1008$w_carcass)*1000
  pinn_preds_conc[7, 2:8] = pinn_preds[7, 2:8]/c(time_1344$w_liver, time_1344$w_blood, time_1344$w_skin, time_1344$w_muscle, time_1344$w_gills, time_1344$w_kidney, time_1344$w_carcass)*1000
  
  RSQUARE = function(y_actual,y_predict){
    cor(y_actual,y_predict)^2
  }
  
  #Load experimental data (Concentrations - ug PFOS/kg tissue)
  experimental_df <- openxlsx::read.xlsx('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015/PFOS.xlsx')
  
  #Create lists with predicted and observed values
  observations <- list()
  predictions <- list()
  
  for (z in 2:dim(pinn_preds_conc)[2]) {
    observations[[z-1]] <- experimental_df[,z]
    names(observations)[z-1] <- colnames(experimental_df)[z] 
    predictions[[z-1]] <- pinn_preds_conc[,z]
    names(predictions)[z-1] <- colnames(pinn_preds_conc)[z]
  }
  
  
  pbkof_score <- PBKtools::PBKOF(observations, predictions)
  
  aafe_score <- PBKtools::AAFE(predictions, observations)

  rmse_score <- PBKtools::rmse(unlist(observations), unlist(predictions))
  
  R2 <- RSQUARE(unlist(observations), unlist(predictions))
  
  return(c('PBKOF'=pbkof_score, 'AAFE'=aafe_score, 'RMSE'= rmse_score, 'R2'=R2))
}

scores_df <- data.frame(matrix(NA, ncol = 5))
colnames(scores_df) <- c('Model', 'PBKOF', 'AAFE', 'RMSE', 'R2')
counter <- 1
for (i in c(0.01, 0.1,1,2,5,10,100)) {
  scores_df[counter,1] <- paste0('model_',i)
  scores_df[counter,2:5] <-  estimate_scores(i)
  counter = counter+1
}

csv_df <- scores_df
csv_df[,-1] <- round(csv_df[,-1], 4)

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/Hyperparameters_tuning/Loss_weights_tuning')
write.csv(csv_df, 'w_results.csv', row.names = F)
