library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    T_C <- 18.6 # Temperature in Celsium 
    T_K <- T_C + 273 # Temperature in Kelvin
    # Biometric parameters 
    Lm <- 59.68 # Maximum length - cm 
    f <- 0.90 # Ad-libitum fraction - unitless
    k <- 8.73e-03 # Growth rate - cm/h
    a_bio <- 1.08e-05 # Constant in the allometric growth equation - unitelss
    b_bio <- 3.03 # Constant in the allometric growth equation - unitelss
    
    TA <- 6930 # Arrhenius temperature - K
    T_ref <- 279.2 #  Temperature at which F_cardref was recorded - K
    # Correction by Arrhenius function for F_card
    KT <- exp((TA/T_K) - (TA/T_ref))
    
    F_card_ref <- 1.19 # Cardiac output reference - L/h/kg
    BW_f_card_ref <- 0.27 # Body weight reference at which F_cardref was - kg
    plasma <- 0.7 # Fraction of plasma in trout
    
    Free <- 3.2e-02 # Free PFOS in plasma
    
    # frac_i: Relative arterial blood fraction to tissue i - unitless
    frac_liver <- 0.04
    frac_muscle <- 0.655
    frac_kidney <- 0.071
    frac_viscera <- 0.069
    frac_brain <- 0.055
    frac_skin <- 0.073
    frac_gill <- 0.0002
    
    # Sc_i: Fraction of body weight - unitless
    sc_blood <- 0.045
    sc_liver <- 0.012
    sc_muscle <- 0.566
    sc_kidney <- 0.016
    sc_viscera <- 0.051
    sc_brain <- 0.001
    sc_skin <- 0.064
    sc_gill <- 0.020
    sc_lumen <- 0.001
    sc_carcass <- 1 - (sc_blood + sc_liver + sc_muscle + sc_kidney + sc_viscera + sc_brain + 
                         sc_skin + sc_gill + sc_lumen)
    
    # Partition coefficients
    P_liver <- 2.09
    P_muscle <- 0.15
    P_kidney <- 0.60
    P_viscera <- 0.75
    P_brain <- 0.64
    P_skin <- 0.40
    P_gill <- 0.25
    P_carcass <- 0.00129
    
    Cl_urine <- 0.18e-05 # Urinary excretion L/h
    Cl_bile <- 2.16e-03 # Biliary excretion L/h
    Cl_feces <- 2.62e-04 # Fecal excretion L/h
    Ku <- 0.13 # Absorption rate constant 1/h
    
    a <- 0.6 #Fraction of muscle going to kidney
    b <- 0.9 #Fraction of skin going to kidney
    
    return(list('Lm'=Lm, 'f'=f, 'k'=k, 'a_bio'=a_bio, 'b_bio'=b_bio, 'a'=a, 'b'=b,
                'KT'=KT, 'F_card_ref'=F_card_ref, 'BW_f_card_ref'=BW_f_card_ref,
                'plasma'=plasma, 'Free'=Free, 'frac_liver'=frac_liver,
                'frac_muscle'=frac_muscle, 'frac_kidney'=frac_kidney,
                'frac_viscera'=frac_viscera, 'frac_brain'=frac_brain,
                'frac_skin'=frac_skin, 'frac_gill'=frac_gill,
                'sc_blood'=sc_blood, 'sc_liver'=sc_liver,
                'sc_muscle'=sc_muscle, 'sc_kidney'=sc_kidney,
                'sc_viscera'=sc_viscera, 'sc_brain'=sc_brain,
                'sc_skin'=sc_skin, 'sc_gill'=sc_gill,
                'sc_lumen'=sc_lumen, 'sc_carcass'=sc_carcass,
                'P_liver'=P_liver, 'P_muscle'=P_muscle, 'P_kidney'=P_kidney,
                'P_viscera'=P_viscera, 'P_brain'=P_brain, 'P_skin'=P_skin, 
                'P_gill'=P_gill, 'P_carcass'=P_carcass, 
                'Cl_urine'=Cl_urine, 'Cl_bile'=Cl_bile,
                'Cl_feces'=Cl_feces, 'Ku'=Ku, 
                'L_0'=L_0, 'BW_0'=BW_0,
                'food_frac'=food_frac, 'C_food'=C_food))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    Q_art<-0; Q_ven<-0; Q_kidney<-0;
    Q_excret_urine<-0; Q_liver<-0;
    Q_excret_bile<-0; Q_viscera<-0;
    Q_lumen_1<-0; Q_lumen_2<-0;
    Q_excret_feces<-0; Q_brain<-0;
    Q_gill<-0; Q_muscle<-0; Q_skin<-0;
    Q_carcass<-0; L<-L_0; Q_food <- BW_0*food_frac*C_food*1000/24
    
    return(c('Q_art'=Q_art, 'Q_ven'=Q_ven, 'Q_kidney'=Q_kidney,
             'Q_excret_urine'=Q_excret_urine, 'Q_liver'=Q_liver,
             'Q_excret_bile'=Q_excret_bile, 'Q_viscera'=Q_viscera,
             'Q_lumen_1'=Q_lumen_1, 'Q_lumen_2'=Q_lumen_2,
             'Q_excret_feces'=Q_excret_feces, 'Q_brain'=Q_brain,
             'Q_gill'=Q_gill, 'Q_muscle'=Q_muscle, 'Q_skin'=Q_skin,
             'Q_carcass'=Q_carcass,'L'=L, 'Q_food'=Q_food))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    
  })
}

custom.func <- function(parameters){
  with(as.list(parameters),{
    
  })
}

ode.func <- function(time,inits, params){
  with(as.list(c(inits, params)),{
    
    # Rainbow trout length in cm
    dL <- k*f*(1-L/Lm) 
    # Calculate the total body weight in kg
    BW = a_bio*L^b_bio
    
    # F_card: Cardiac output - L/h
    F_card <- F_card_ref * KT *(BW/BW_f_card_ref)^(-0.1)*BW*plasma
    
    # Calculate the blood flows of the compartment - L/h
    F_liver <- F_card*frac_liver
    F_muscle <- F_card*frac_muscle
    F_kidney <- F_card*frac_kidney
    F_viscera <- F_card*frac_viscera
    F_brain <- F_card*frac_brain
    F_skin <- F_card*frac_skin
    F_gill <- F_card*frac_gill
    F_carcass <- F_card - (F_liver+F_muscle+F_kidney+F_viscera+F_brain+F_skin+F_gill)

    # Calculate the weight of the organs - kg
    w_blood <- sc_blood*BW
    w_art <- w_blood*1/3
    w_ven <- w_blood*2/3
    w_liver <- sc_liver*BW
    w_muscle <- sc_muscle*BW
    w_kidney <- sc_kidney*BW
    w_viscera <- sc_viscera*BW
    w_brain <- sc_brain*BW
    w_skin <- sc_skin*BW
    w_gill <- sc_gill*BW
    w_lumen <- sc_lumen*BW
    w_carcass <- sc_carcass*BW
    
    # Calculate concentrations - ng/kg
    C_art <- Q_art/w_art
    C_venous <- Q_ven/w_ven
    C_kidney <- Q_kidney/w_kidney
    C_liver <- Q_liver/w_liver
    C_viscera <- Q_viscera/w_viscera
    C_brain <- Q_brain/w_brain
    C_gill <- Q_gill/w_gill
    C_muscle <- Q_muscle/w_muscle
    C_skin <- Q_skin/w_skin
    C_carcass <- Q_carcass/w_carcass
    C_lumen_1 <- Q_lumen_1/w_lumen
    C_lumen_2 <- Q_lumen_2/w_lumen
    
    
    # Arterial blood
    dQ_art <- Free*F_card*C_venous - F_liver*C_art - F_muscle*C_art - F_brain*C_art -
              F_viscera*C_art - F_kidney*C_art - F_skin*C_art - F_gill*C_art - F_carcass*C_art
    
    # Venous blood
    dQ_ven <-  - Free*F_card*C_venous + #dQ_admin_water - dQ_excret_water
      (F_liver + F_viscera)*C_liver/P_liver + (1-a)*F_muscle*C_muscle/P_muscle +
      F_brain*C_brain/P_brain + (a*F_muscle+b*F_skin+F_kidney)*C_kidney/P_kidney + F_gill*C_gill/P_gill + 
      (1-b)*F_skin*C_skin/P_skin + F_carcass*C_carcass/P_carcass
    
    # Kidney 
    dQ_kidney <- F_kidney*Free*C_art + a*F_muscle*C_muscle/P_muscle + 
      b*F_skin*C_skin/P_skin - (a*F_muscle + b*F_skin + F_kidney)*C_kidney/P_kidney - Cl_urine*Free*C_kidney/P_kidney 
    
    # Amount excreted via urine
    dQ_excret_urine <- Cl_urine*Free*C_kidney/P_kidney
                  
    # Liver
    dQ_liver <- F_liver*Free*C_art + F_viscera*C_viscera/P_viscera - 
      (F_liver + F_viscera)*C_liver/P_liver - Cl_bile*Free*C_liver
    
    # Bile elimination 
    dQ_excret_bile <- Cl_bile*Free*C_liver
    
    # Viscera
    dQ_viscera <- F_viscera*(C_art-C_viscera/P_viscera) + Ku*Q_lumen_1
    
    dQ_food <- BW*food_frac*C_food*1000/24
    
    # Gut lumen 1
    if(time < 42*24){
      dQ_lumen_1 <- Q_food - Ku*Q_lumen_1 -Cl_feces*C_lumen_1
    }else{
      dQ_lumen_1 <- - Ku*Q_lumen_1 -Cl_feces*C_lumen_1
    }
    
    # Gut lumen 2
    dQ_lumen_2 <- Cl_bile*Free*C_liver - Cl_feces*C_lumen_2
    
    # Amount excreted via feces
    dQ_excret_feces <- Cl_feces*(C_lumen_1 + C_lumen_2)
    
    # Brain
    dQ_brain <- F_brain*Free*(C_art-C_brain/P_brain)
    
    # Gill
    dQ_gill <- F_gill*Free*(C_art-C_gill/P_gill)
    
    # Muscle
    dQ_muscle <- F_muscle*Free*(C_art-C_muscle/P_muscle)
    
    # Skin
    dQ_skin <- F_skin*Free*(C_art-C_skin/P_skin)
    
    # Carcass
    dQ_carcass <- F_carcass*Free*(C_art-C_carcass/P_carcass)
    
    return(list(c('dQ_art'=dQ_art, 'dQ_ven'=dQ_ven, 'dQ_kidney'=dQ_kidney,
                  'dQ_excret_urine'=dQ_excret_urine, 'dQ_liver'=dQ_liver,
                  'dQ_excret_bile'=dQ_excret_bile, 'dQ_viscera'=dQ_viscera,
                  'dQ_lumen_1'=dQ_lumen_1, 'dQ_lumen_2'=dQ_lumen_2,
                  'dQ_excret_feces'=dQ_excret_feces, 'dQ_brain'=dQ_brain,
                  'dQ_gill'=dQ_gill, 'dQ_muscle'=dQ_muscle, 'dQ_skin'=dQ_skin,
                  'dQ_carcass'=dQ_carcass, 'dL'=dL, 'dQ_food'=dQ_food), 
                'C_art'=C_art, 'C_venous'=C_venous, 'C_kidney'=C_kidney,
                'C_liver'=C_liver, 'C_brain'=C_brain, 'C_gill'=C_gill,
                'C_muscle'=C_muscle, 'C_skin'=C_skin, 'C_carcass'=C_carcass))
  })
}


BW_0 <- 0.32 # initial body weight in kg  
L_0 <- mean(c(280,328,328,328,280))/10 # average initial length in cm
C_food <- 462.5 # ng/g food
food_frac <- 1/100 # of total BW once per day
sample_time <- seq(0, (42+37)*24, 1) 


user.input <- list('BW_0'=BW_0,
                   'L_0'=L_0,
                   'food_frac'=food_frac,
                   'C_food'=C_food)

params <- create.params(user.input)
inits <- create.inits(params)
#events <- create.events(params)

solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           #events = events,
                           method="lsodes",rtol = 1e-05, atol = 1e-05))

plot(solution$time, solution$C_brain)
