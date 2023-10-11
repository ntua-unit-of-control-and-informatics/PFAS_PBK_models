# merge and create the dataset for the PINN implementation

added_pfas <- read.csv('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/NN_events.csv',)[,c(2,3)]
measurements_times <- c(168, 336, 672, 744, 840, 1008, 1344)

cumulative_pfas <- data.frame(matrix(data=c(0, 0,
                                            168, sum(added_pfas[added_pfas$time<168,'value']),
                                            336, sum(added_pfas[added_pfas$time<336,'value']),
                                            672, sum(added_pfas[,'value']),
                                            744, sum(added_pfas[,'value']),
                                            840, sum(added_pfas[,'value']),
                                            1008, sum(added_pfas[,'value']),
                                            1344, sum(added_pfas[,'value'])), nrow=8, byrow = T))
colnames(cumulative_pfas) <- c('Time', 'Cumulative_added_pfas')

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015')

Feeding_period <- c(0, 168, 336, 672, 0, 0, 0, 0)
Depuration_period <- c(0, 0, 0, 0, 72, 168, 336, 672)

PFBS <- openxlsx::read.xlsx('PFBS.xlsx')
PFBS <- rbind(rep(0,8), PFBS)
PFBS$Time <- PFBS$Time*24
PFBS <- cbind(rep('PFBS',8), PFBS[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFBS[,-c(1)] )
colnames(PFBS)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFHxS <- openxlsx::read.xlsx('PFHxS.xlsx')
PFHxS <- rbind(rep(0,8), PFHxS)
PFHxS$Time <- PFHxS$Time*24
PFHxS <- cbind(rep('PFHxS',8), PFHxS[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFHxS[,-c(1)] )
colnames(PFHxS)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFOS <- openxlsx::read.xlsx('PFOS.xlsx')
PFOS <- rbind(rep(0,8), PFOS)
PFOS$Time <- PFOS$Time*24
PFOS <- cbind(rep('PFOS',8), PFOS[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFOS[,-c(1)] )
colnames(PFOS)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFOA <- openxlsx::read.xlsx('PFOA.xlsx')
PFOA <- rbind(rep(0,8), PFOA)
PFOA$Time <- PFOA$Time*24
PFOA <- cbind(rep('PFOA',8), PFOA[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFOA[,-c(1)] )
colnames(PFOA)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

PFNA <- openxlsx::read.xlsx('PFNA.xlsx')
PFNA <- rbind(rep(0,8), PFNA)
PFNA$Time <- PFNA$Time*24
PFNA <- cbind(rep('PFNA',8), PFNA[,1], cumulative_pfas[,2], Feeding_period, Depuration_period,PFNA[,-c(1)] )
colnames(PFNA)[c(1:3)] <- c('Substance', 'Time', 'Cumulative_added_PFAS') 

df <- rbind(PFBS,PFHxS,PFOS,PFOA,PFNA)

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN')
#write.csv(df, 'dataset.csv', row.names = F)


dummy_df <- PFOS[-1,]
#write.csv(dummy_df, 'PFOS_dataset.csv')


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
mass_dummy_df <- dummy_df
mass_dummy_df[1, 6:12] = dummy_df[1, 6:12]*c(time_168$w_liver, time_168$w_blood, time_168$w_skin, time_168$w_muscle, time_168$w_gills, time_168$w_kidney, time_168$w_carcass)/1000
mass_dummy_df[2, 6:12] = dummy_df[2, 6:12]*c(time_336$w_liver, time_336$w_blood, time_336$w_skin, time_336$w_muscle, time_336$w_gills, time_336$w_kidney, time_336$w_carcass)/1000
mass_dummy_df[3, 6:12] = dummy_df[3, 6:12]*c(time_672$w_liver, time_672$w_blood, time_672$w_skin, time_672$w_muscle, time_672$w_gills, time_672$w_kidney, time_672$w_carcass)/1000
mass_dummy_df[4, 6:12] = dummy_df[4, 6:12]*c(time_744$w_liver, time_744$w_blood, time_744$w_skin, time_744$w_muscle, time_744$w_gills, time_744$w_kidney, time_744$w_carcass)/1000
mass_dummy_df[5, 6:12] = dummy_df[5, 6:12]*c(time_840$w_liver, time_840$w_blood, time_840$w_skin, time_840$w_muscle, time_840$w_gills, time_840$w_kidney, time_840$w_carcass)/1000
mass_dummy_df[6, 6:12] = dummy_df[6, 6:12]*c(time_1008$w_liver, time_1008$w_blood, time_1008$w_skin, time_1008$w_muscle, time_1008$w_gills, time_1008$w_kidney, time_1008$w_carcass)/1000
mass_dummy_df[7, 6:12] = dummy_df[7, 6:12]*c(time_1344$w_liver, time_1344$w_blood, time_1344$w_skin, time_1344$w_muscle, time_1344$w_gills, time_1344$w_kidney, time_1344$w_carcass)/1000

write.csv(mass_dummy_df, file = '/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/PFOS_Mass_dataset.csv')
