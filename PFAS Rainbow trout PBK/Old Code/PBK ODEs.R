library(deSolve)

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
  
  # Total cardiac output ml/h
  Q_total <- F_card_ref*KT*(BW/BW_ref)^(-0.1)*BW  
  
  # Load the xlsx file with the physiological params pf rainbow trout
  phys.params <- openxlsx::read.xlsx('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 1)
  
  # Keep only the physiological parameters from the paper of Vidal et al. 2019
  # fw are the fractions of tissue_weight/total_weight (unitless)
  fw <- phys.params[phys.params$Source=='Vidal et al. 2019', c('Liver', 'Blood', 'Skin',
                                                               'Muscle', 'Gills', 'Kidney', 'Viscera')]
  fw_lumen <- 0.012
  # Calculate the mass of each tissue - g
  w_blood <- fw$Blood*BW     # Blood mass - g
  w_liver <- fw$Liver*BW     # Liver mass - g
  w_skin <- fw$Skin*BW       # Skin weight - g
  w_muscle <- fw$Muscle*BW   # Muscle weight - g
  w_gills <- fw$Gills*BW     # Gills weight - g
  w_kidney <- fw$Kidney*BW   # Kidney weight - g
  w_viscera <- fw$Viscera*BW # Viscera weight - g
  w_lumen <- fw_lumen*BW
  w_bile <- w_lumen
  w_art <- 1/3*w_blood
  w_venous <- 2/3*w_blood
  w_carcass <- BW - (w_blood + w_liver + w_skin + w_muscle +
                 w_gills + w_kidney + w_viscera + w_lumen + w_bile)
  
  # Load the xlsx file with the physiological params pf rainbow trout
  phys.params <- openxlsx::read.xlsx('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 2)
  
  # Keep only the physiological parameters from the paper of Vidal et al. 2019
  # fb are the fractions of blood flow of each tissue (unitless)
  fb <- phys.params[phys.params$Source=='Vidal et al. 2019', c('Liver', 'Skin', 'Muscle',
                                                               'Gills', 'Kidney', 'Viscera')]
  
  # Calculate the regional blood flows - ml/h
  Q_liver <- fb$Liver*BW     # Liver blood flow - ml/h
  Q_skin <- fb$Skin*BW       # Skin blood flow - ml/h
  Q_muscle <- fb$Muscle*BW   # Muscle blood flow - ml/h
  Q_gills <- fb$Gills*BW     # Gills blood flow - ml/h
  Q_kidney <- fb$Kidney*BW   # Kidney blood flow - ml/h
  Q_viscera <- fb$Viscera*BW # Viscera blood flow - ml/h
  Q_carcass <- Q_total -(Q_liver + Q_skin + Q_muscle + Q_gills + 
                       Q_kidney + Q_viscera)
  
  Cl_bile <- 2.16e-06 # L/h 
  Cl_feces <- 2.62e-07 # L/h 
  Cl_urine <- 7e-05 # L/h 
  Ku <- 0.13 # 1/h
  P_liver <- 2.05
  P_muscle <- 0.15
  P_kidney <- 0.58
  P_viscera <- 0.87
  P_brain <- 0.63
  P_skin <- 0.25
  P_gills <- 0.40
  P_carcass <- 0.00129             
  Free = 3.2e-02
  
  
  return(list('admin.dose_dietary'=admin.dose_dietary, 
              'admin.time_dietary'=admin.time_dietary,
              'BW'=BW, 'Q_total'=Q_total,
              'w_art'=w_art, 'w_venous'=w_venous,  
              'w_blood'=w_blood, 'w_liver'=w_liver, 'w_skin'=w_skin,
              'w_muscle'=w_muscle, 'w_gills'=w_gills, 'w_kidney'=w_kidney,
              'w_viscera'=w_viscera, 'w_lumen'=w_lumen, 'w_bile'=w_bile, 'w_carcass'=w_carcass,
              'Q_liver'=Q_liver, 'Q_skin'=Q_skin, 'Q_muscle'=Q_muscle,
              'Q_gills'=Q_gills, 'Q_kidney'=Q_kidney, 'Q_viscera'=Q_viscera,
              'Q_carcass'=Q_carcass,
              'Cl_bile'=Cl_bile, 'Cl_feces'=Cl_feces, 'Cl_urine'=Cl_urine,
              'Ku'=Ku, 'P_liver'=P_liver, 'P_muscle'=P_muscle, 'P_kidney'=P_kidney,
              'P_viscera'=P_viscera, 'P_brain'=P_brain, 'P_skin'=P_skin, 
              'P_gills'=P_gills, 'P_carcass'=P_carcass, 'Free'=Free))
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
    
    # Arterial Blood
    dM_art <- Free*Q_total*C_venous - 
      (Q_gills + Q_viscera + Q_liver + Q_kidney +
         Q_muscle + Q_skin + Q_carcass)*Free*C_art
    
    dM_venous <- - Free*Q_total*C_venous + (Q_gills*C_gills/P_gills +
      (Q_liver + Q_viscera)*C_liver/P_liver + 
      Q_kidney*C_kidney/P_kidney + Q_muscle*C_muscle/P_muscle +
      Q_skin*C_skin/P_skin + Q_carcass*C_carcass/P_carcass)*Free
    
    # Gills
    dM_gills <- Q_gills*Free*(C_art - C_gills/P_gills)
    
    dM_input=0
    
    # Viscera lumen 1 - Available PFAS for absorption and elimination
    dM_lumen = - Ku*M_lumen - Cl_feces*C_lumen 
    
    # Bile - Billiary eliminated PFAS (Considering no reabsoprtion)
    dM_bile = Free*Cl_bile*C_liver - Cl_feces*C_bile
    
    # Viscera tissue
    dM_viscera <- Q_viscera*Free*(C_art - C_viscera/P_viscera) + Ku*M_lumen 
    
    # Liver
    dM_Liver <- Q_liver*Free*C_art + Q_viscera*Free*C_viscera/P_viscera - 
      (Q_liver + Q_viscera)*Free*C_liver/P_liver - Free*Cl_bile*C_liver
    
    # Kidney
    dM_kidney <- Q_kidney*Free*(C_art - C_kidney/P_kidney) - Cl_urine*C_kidney
    
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
                'C_gills'=C_gills, 'C_bile'=C_bile, 'C_viscera'=C_viscera,
                'C_liver'=C_liver, 'C_kidney'=C_kidney, 'C_muscle'=C_muscle,
                'C_skin'=C_skin, 'C_carcass'=C_carcass, 'C_lumen'=C_lumen,
                'Mass_balance'=Mass_balance))
  })
}

# Total body weight of fish
BW <- 314 # g 
Texp <- 15 # C
admin.dose_dietary <- c(1, 5, 0)
admin.time_dietary <- c(0, 5, 7)

user.input <- list('BW'=BW,
                   'Texp'=Texp,
                   'admin.dose_dietary'=admin.dose_dietary,
                   'admin.time_dietary'=admin.time_dietary)

params <- create.params(user.input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,10,1)

solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-05, atol = 1e-05))
solution$Mass_balance
