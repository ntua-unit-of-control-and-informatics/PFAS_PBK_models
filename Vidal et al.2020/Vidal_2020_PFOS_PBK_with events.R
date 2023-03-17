library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    #============================================
    # Units 
    #============================================
    # Volume           : L 
    # Time             : hours
    # Flow             : L/h
    # Concentration    : ng/kg
    # Masses           : kg
    # Temperature      : C
    # Ventilation rate : L/h
    # VO2		   : mg(O2)/h
    # Cox	           : mg(O2)/L
    # Length           : cm
    # Clearance 	   : L/h
    #============================================
    
    # Plasma fraction 
    #--------------------------------------------
    plasma = 0.7 # Plasmatic fraction = 1 - Haematocrit
    
    # Free fraction in plasma
    #--------------------------------------------
    Free = 0.0315
    
    # Reabsorbed fraction - entero-hepatic cycle 
    #--------------------------------------------
    f_reab = 0.0
    
    # Fraction entrant par les ouies
    #--------------------------------------------
    C_permeab = 1.0
    
    # Biometric data
    #--------------------------------------------
    Lm = 73.91;			# Optimised value RB revision (cm)
    kappa =  0.0096 # Optimised value RB revision (cm/day)
    a_bio = 9.57767e-06 
    b_bio = 3.05088 
    
    # volume scaling factor : fraction of BW
    # sc_organ have been calculated for the 7C experiment fish : mean mass of fish and organ at the start of experiment
    #--------------------------------------------
    sc_blood    = 0.0449 
    sc_liver    = 0.01155	
    sc_muscle   = 0.566	
    sc_skin     = 0.0638
    sc_gill     = 0.0196
    sc_kidney   = 0.01643
    sc_viscera  = 0.0514
    sc_brain    = 0.001
    sc_lumen    = 0.012 #nichols et al., 2004
    
    # Fraction of arterial blood flow
    # frac_organ = %blood_flow_rate * Qc / organ_weight
    #--------------------------------------------
    frac_liver    = 0.0035		#Barron et al., 1987, 18C
    frac_muscle_ref   = 0.655	#Barron et al., 1987, 18C
    frac_skin      = 0.0728		#Nichols et al., 1996
    frac_gill      = 0.0021		#Barron et al., 1987, 18C
    frac_kidney    = 0.0712		#Barron et al., 1987, 18C
    frac_viscera   = 0.069		#Barron et al., 1987, 18C
    frac_brain     = 0.055		#Pery et al., 2014
    frac_carcass = 1 - (frac_liver+ frac_muscle_ref + frac_skin  + frac_gill + frac_kidney+ frac_viscera  + frac_brain)
    
    #--------------------------------------------
    a = 0.6	# Muscle fraction of flow to the kidney (Nichols et al., 1990)
    b = 0.9	# Skin fraction of flow to the kidney  (Nichols et al., 1996)
    
    # Partition coefficient  
    #--------------------------------------------
    PC_liver_ref = 2.05	
    PC_muscle_ref = 0.15	
    PC_skin      = 0.289531 
    PC_gill      = 0.355483 
    PC_kidney    = 0.58
    PC_viscera   = 0.87
    PC_brain     = 0.63
    PC_carcass   = 0.00129
    
    PC_blood_water_ref = 4239.0
    
    # Rate constants (absorption and elimination) 
    #--------------------------------------------
    Ku_ref 	 =  0.07	
    Cl_urine_ref =  0.0000179 #Consoer et al., 2016 - BW_urine_ref = 1.037 ;# Kg 
    Cl_feces_ref =  0.00019
    Cl_bile_ref  =  0.00083
    
    # Physico-chemical parameters 
    #--------------------------------------------
    T_ref = 291.75 
    TA_VO2 = 6930
    TA_Qc = 6930
    TA_perf = 6930
    TA_PC = 5664
    TA_ku = 5423
    TA_clearance = 8267
    
    
    VO2_ref = 135.8	# mgO2/h/kg Elliott 1969 - 3.26 mgO2/d/g --> 3.26/24*1000 = 135.8 mgO2/h/kg
    BW_VO2_ref = 1.0 
    T_VO2_ref = 283.15
    
    Qc_ref     = 1.188   	# Barron et al., 1987 - 12C : Qc = 19.8 mL/min/kg --> 38.7*60/1000 L/h/kg
    BW_Qc_ref  = 0.2701  	# kg 
    T_Qc_ref   = 279.15  	# Kelvin
    
    return(list('L0'=L0, 'Texp'=Texp, 'Cox'=Cox, 'Concentration_water'=Concentration_water,
                'admin.dose_dietary'=admin.dose_dietary, 'admin.time_dietary'=admin.time_dietary,
                'plasma'=plasma, 'Free'=Free, 'f_reab'=f_reab, 'C_permeab'=C_permeab,
                'Lm'=Lm, 'kappa'=kappa, 'a_bio'=a_bio, 'b_bio'=b_bio, 'f'=f, 
                'sc_blood'=sc_blood, 'sc_liver'=sc_liver, 'sc_muscle'=sc_muscle,
                'sc_skin'=sc_skin, 'sc_gill'=sc_gill, 'sc_kidney'=sc_kidney,
                'sc_viscera'=sc_viscera, 'sc_brain'=sc_brain, 'sc_lumen'=sc_lumen,
                'frac_liver'=frac_liver, 'frac_muscle_ref'=frac_muscle_ref,
                'frac_skin'=frac_skin, 'frac_gill'=frac_gill, 'frac_kidney'=frac_kidney,
                'frac_viscera'=frac_viscera, 'frac_brain'=frac_brain, 'frac_carcass'=frac_carcass,
                'a'=a, 'b'=b, 'PC_liver_ref'=PC_liver_ref, 'PC_muscle_ref'=PC_muscle_ref,
                'PC_skin'=PC_skin, 'PC_gill'=PC_gill, 'PC_kidney'=PC_kidney, 
                'PC_viscera'=PC_viscera, 'PC_brain'=PC_brain, 'PC_carcass'=PC_carcass,
                'PC_blood_water_ref'=PC_blood_water_ref, 'Ku_ref'=Ku_ref, 
                'Cl_urine_ref'=Cl_urine_ref, 'Cl_feces_ref'=Cl_feces_ref,
                'Cl_bile_ref'=Cl_bile_ref, 'T_ref'=T_ref, 'TA_VO2'=TA_VO2,
                'TA_Qc'=TA_Qc, 'TA_perf'=TA_perf, 'TA_PC'=TA_PC, 'TA_ku'=TA_ku,
                'TA_clearance'=TA_clearance, 'VO2_ref'=VO2_ref, 'BW_VO2_ref'=BW_VO2_ref,
                'T_VO2_ref'=T_VO2_ref, 'Qc_ref'=Qc_ref, 'BW_Qc_ref'=BW_Qc_ref, 
                'T_Qc_ref'=T_Qc_ref))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    Qadmin_water=0; Qadmin_food=0; 
    L = L0; Q_lumen_1 = 0; Q_lumen_2 = 0; Q_art= 0; Q_ven= 0; Q_liver= 0;
    Q_muscle= 0; Q_brain= 0; Q_viscera= 0; Q_kidney= 0; Q_skin= 0;
    Q_gill= 0; Q_carcass= 0; Qexcret_gill= 0; Qexcret_bile= 0; Qexcret_urine= 0;
    Qexcret_feces<-0
    
    return(c('Qadmin_food'=Qadmin_food, 'Qadmin_water'=Qadmin_water,
             'L'=L, 'Q_lumen_1'=Q_lumen_1, 'Q_lumen_2'=Q_lumen_2,
             'Q_art'=Q_art, 'Q_ven'=Q_ven, 'Q_liver'=Q_liver,
             'Q_muscle'=Q_muscle, 'Q_brain'=Q_brain, 'Q_viscera'=Q_viscera,
             'Q_kidney'=Q_kidney, 'Q_skin'=Q_skin, 'Q_gill'=Q_gill,
             'Q_carcass'=Q_carcass, 'Qexcret_gill'=Qexcret_gill,
             'Qexcret_bile'=Qexcret_bile, 'Qexcret_urine'=Qexcret_urine,
             'Qexcret_feces'=Qexcret_feces))
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
      events <- data.frame(var = c(rep(c('Q_lumen_1', 'Qadmin_food'), ltimes_dietary)), 
                           time = sort(rep(admin.time_dietary,2)),
                           value = sort(rep(admin.dose_dietary,2)),
                           method = 'add')
    }
    
    #events <- events[order(events$time),] 
    return(list(data=events))
  })
}

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Arrhenius physiological processes
    KT_VO2= exp((TA_VO2/T_VO2_ref) - (TA_VO2/Texp))
    KT_Qc = exp((TA_Qc/T_Qc_ref)-(TA_Qc/Texp))
    KT_PC = exp((TA_PC/T_ref)-(TA_PC/Texp))
    KT_ku = exp((TA_ku/T_ref)-(TA_ku/Texp))
    KT_clearance = exp((TA_clearance/T_ref)-(TA_clearance/Texp))
    
    # Body weight - kg
    BW <- a_bio*L^b_bio
    
    KT = exp((TA_perf/T_ref)-(TA_perf/Texp)) ; 	
    frac_muscle   = frac_muscle_ref * KT ;		# KT : effect of temperature
    Delta_muscle = frac_muscle_ref * (1-KT) ;
    
    a_liver = frac_liver /(1 -frac_muscle_ref)		# 1 = sum of frac_i
    a_skin = frac_skin /(1 -frac_muscle_ref)
    a_gill = frac_gill /(1 -frac_muscle_ref)
    a_kidney = frac_kidney /(1 -frac_muscle_ref)
    a_viscera = frac_viscera/(1 -frac_muscle_ref)
    a_brain = frac_brain /(1 -frac_muscle_ref)
    a_carcass = frac_carcass /(1 -frac_muscle_ref)
    
    #Volumes (L or Kg) of the organs changing with the time
    # Organ density considered equal to 1 (Kg/L)
    #--------------------------------------------
    V_art = sc_blood * BW  * 1/3.0 * plasma  # Kg
    V_ven = sc_blood * BW  * 2/3.0 * plasma 
    
    V_liver = sc_liver * BW 
    V_muscle = sc_muscle * BW 
    V_brain = sc_brain * BW 
    V_viscera = sc_viscera* BW 
    V_skin = sc_skin * BW 
    V_gill = sc_gill * BW 
    V_kidney = sc_kidney * BW 
    V_lumen = sc_lumen * BW 
    V_carcass = (1 - sc_blood - sc_liver - sc_muscle - sc_brain - sc_viscera - sc_kidney - sc_skin - sc_gill - sc_lumen) * BW  
    
    # Calculation of cardiac output and flow  -- correction for plasma
    #--------------------------------------------    
    Qc = Qc_ref*KT_Qc*(BW/BW_Qc_ref)^(-0.1)*BW*plasma # L/h
    
    # Partition Coefficients
    PC_liver = PC_liver_ref*KT_PC
    PC_muscle = PC_muscle_ref*(1/KT_PC)
    PC_blood_water = PC_blood_water_ref*KT_PC
    
    # Flows to tissues
    F_liver    = (frac_liver + a_liver*Delta_muscle) * Qc * Free
    F_muscle   = frac_muscle * Qc * Free
    F_brain    = (frac_brain + a_brain * Delta_muscle) * Qc * Free
    F_viscera  = (frac_viscera + a_viscera * Delta_muscle) * Qc * Free
    F_kidney   = (frac_kidney + a_kidney * Delta_muscle) * Qc * Free 
    F_skin     = (frac_skin + a_skin * Delta_muscle) * Qc * Free 
    F_gill     = (frac_gill + a_gill * Delta_muscle) * Qc * Free 
    #F_carcass  = (1 - frac_liver - frac_muscle - frac_brain - frac_viscera - frac_kidney - frac_skin - frac_gill) * Qc * Free 
    F_carcass = (frac_carcass + a_carcass * Delta_muscle) * Qc * Free
    
    #Ventilation rate 
    #--------------------------------------------
    VO2_arr  =  VO2_ref * KT_VO2 * (BW/BW_Qc_ref)^(-0.1) #  mgO2/h/kg   GRECH
    VO2 =  VO2_arr * BW #  mg O2/h
    
    Qw =  VO2/(Cox-0.2*Cox) # Effective respiratory volume (L/h)
    
    # Chemical flux at fish gills (general expression)
    #--------------------------------------------
    #tmp = Qc* PC_blood_water ;
    #Kx = (tmp < Qw ? tmp : Qw);						# Kx equals the flow term that is limiting (L/d)
    Kx = Qw
    
    # Concentrations in tissues (ng/kg) 
    #--------------------------------------------
    C_liver   = Q_liver  / V_liver
    C_muscle  = Q_muscle / V_muscle	
    C_brain   = Q_brain  / V_brain	
    C_viscera = Q_viscera/ V_viscera 
    C_kidney  = Q_kidney / V_kidney	
    C_skin    = Q_skin   / V_skin	
    C_gill    = Q_gill   / V_gill	
    C_carcass = Q_carcass/ V_carcass
    C_lumen_1= Q_lumen_1 / V_lumen 
    C_lumen_2= Q_lumen_2 / V_lumen 
    C_lumen_viscera = (Q_lumen_1 + Q_lumen_2) / V_lumen 
    
    # because equal to 1 (PC_ven and PC_art) 
    C_art     = Q_art    / V_art	
    C_ven     = Q_ven    / V_ven	
    C_blood_total <- (Q_art+Q_ven)/(V_art+V_ven)
    
    # Kinetics
    #-------------------------------------------	
    Ku    = Ku_ref * KT_ku	      		
    Cl_urine = Cl_urine_ref * KT_clearance 		
    Cl_feces = Cl_feces_ref * KT_clearance 			
    Cl_bile  = Cl_bile_ref  * KT_clearance
    
    # Differentials
    dL = kappa *f* (1 - (L/Lm))
    
    # Absorption
    #--------------------------------------------
    dQadmin_food = 0	# ng/h
    
    dQadmin_water = C_permeab * Kx * Concentration_water
    
    # Excretion
    #--------------------------------------------
    dQexcret_gill =  Kx * ( Free * C_ven/PC_blood_water )		
    dQexcret_bile =  Free * Cl_bile * C_liver      			
    dQexcret_urine =  Cl_urine * (Free * C_kidney / PC_kidney)   	
    dQexcret_feces =  Cl_feces * ( C_lumen_2 + C_lumen_1) 
    
    # Distribution
    #--------------------------------------------
    dQ_carcass = F_carcass * (C_art - (C_carcass/PC_carcass))
    
    dQ_gill = F_gill    * (C_art - C_gill/PC_gill)		
    
    dQ_muscle  = F_muscle  * (C_art - C_muscle/PC_muscle)		  # F_muscle = a*F_muscle + (1-a)*F_muscle
    
    dQ_skin = F_skin  * (C_art - C_skin/PC_skin)		  # F_skin = b*F_skin + (1-b)*F_skin
    
    dQ_brain = F_brain * (C_art - C_brain/PC_brain)		
    
    dQ_kidney = F_kidney  * C_art + a * F_muscle * C_muscle/PC_muscle + 
      b * F_skin * C_skin/PC_skin - (a * F_muscle + b * F_skin + F_kidney) * C_kidney/PC_kidney -
      Cl_urine * (Free * C_kidney / PC_kidney)
    
    dQ_lumen_1 = ( - Ku * Q_lumen_1 ) - ( Cl_feces * C_lumen_1) + ( f_reab * (Free * Cl_bile * C_liver))
    
    dQ_lumen_2 = (1 - f_reab ) * (Free * Cl_bile * C_liver) - ( Cl_feces * C_lumen_2)
    
    dQ_viscera =  F_viscera*(C_art - C_viscera/PC_viscera) + Ku * Q_lumen_1
    
    dQ_liver = F_liver * C_art + F_viscera *(C_viscera/PC_viscera) - 
      (F_liver + F_viscera) * (C_liver/PC_liver) - (Free * Cl_bile * C_liver)
    
    dQ_art = Qc * C_ven * Free - F_liver * C_art - F_muscle * C_art - F_brain * C_art -
      F_viscera * C_art - F_kidney  * C_art - F_skin * C_art - F_gill * C_art -
      F_carcass * C_art
    
    dQ_ven = (C_permeab * Kx * Concentration_water) - (Kx * ( Free * C_ven/PC_blood_water)) -
      Qc * C_ven * Free + (F_liver + F_viscera)* C_liver/PC_liver +
      (1-a) * F_muscle  * C_muscle/PC_muscle +  F_brain * C_brain/PC_brain +
      (a* F_muscle + b * F_skin + F_kidney ) * C_kidney/PC_kidney +
      F_gill * C_gill/PC_gill +  (1-b) * F_skin * C_skin/PC_skin +  
      F_carcass * C_carcass/PC_carcass
    
    Total_excreted <- Qexcret_urine + Qexcret_gill + Qexcret_feces
    Total_admin <- Qadmin_water + Qadmin_food
    Total_accumulated <- Q_art + Q_ven + Q_liver + Q_muscle + Q_brain  + 
      Q_viscera + Q_kidney + Q_skin + Q_gill + Q_carcass + Q_lumen_1 + Q_lumen_2
    Mass_equilibrium = (Total_admin - Total_excreted - Total_accumulated)
    
    return(list(c('dQadmin_food'=dQadmin_food, 'dQadmin_water'=dQadmin_water, 
                  'dL'=dL, 'dQ_lumen_1'=dQ_lumen_1, 'dQ_lumen_2'=dQ_lumen_2,
                  'dQ_art'=dQ_art, 'dQ_ven'=dQ_ven, 'dQ_liver'=dQ_liver,
                  'dQ_muscle'=dQ_muscle, 'dQ_brain'=dQ_brain, 'dQ_viscera'=dQ_viscera,
                  'dQ_kidney'=dQ_kidney, 'dQ_skin'=dQ_skin, 'dQ_gill'=dQ_gill,
                  'dQ_carcass'=dQ_carcass, 'dQexcret_gill'=dQexcret_gill,
                  'dQexcret_bile'=dQexcret_bile, 'dQexcret_urine'=dQexcret_urine, 
                  'dQexcret_feces'=dQexcret_feces),
                'C_liver'=C_liver, 'C_muscle'=C_muscle, 'C_brain'=C_brain,
                'C_viscera'=C_viscera, 'C_kidney'=C_kidney, 'C_skin'=C_skin, 'C_gill'=C_gill,
                'C_carcass'=C_carcass, 'C_lumen_1'=C_lumen_1, 'C_lumen_2'=C_lumen_2, 
                'C_lumen_viscera'=C_lumen_viscera, 'C_art'=C_art, 'C_ven'=C_ven,
                'C_blood_total'=C_blood_total,
                'Total_excreted'=Total_excreted, 'Total_admin'=Total_admin,
                'Total_accumulated'=Total_accumulated, 'Mass_equilibrium'=Mass_equilibrium))
  })
}


################################################################################

# Reproduce results for T_exp <- 19 # C

growth_inits <- c("L"=30.88)
gworth_params <- c(Lm = 73.91,			# Optimised value RB revision (cm)
                   kappa =  0.0096, # Optimised value RB revision (cm/day)
                   a_bio = 9.57767e-06, 
                   b_bio = 3.05088,
                   f = 0.616003)

growth_function <- function(time, inits, params){
  with(as.list(c(inits,params)),{
    dL = kappa *f* (1 - (L/Lm))
    # Body weight - kg
    BW <- a_bio*L^b_bio
    return(list(c("dL"=dL),"BW"=BW))
  })
}
growth_times <- seq(0,41*24,24)
growth_solution <- data.frame(ode(times = growth_times,  func = growth_function, y = growth_inits, parms = gworth_params,
                                  #events = events,
                                  method="lsodes",rtol = 1e-05, atol = 1e-05))

admin.dose_dietary <- 0.01*growth_solution$BW*1000*500 # ng PFOS
admin.time_dietary <- growth_solution$time

user.input_19 <- list(Texp = 19 + 273,      # C
                      Cox = 8.05,   # mg(O2)/L
                      L0 = 30.88,   # cm
                      f = 0.616003, # unitless
                      admin.dose_dietary=admin.dose_dietary, # ng PFOS/g dw of food
                      admin.time_dietary=admin.time_dietary,
                      Concentration_water =	0.4	# PFOS concentration in water (ng/L)
)
params <- create.params(user.input_19)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0, (42+35)*24, 1)

solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-05, atol = 1e-05))


plot(solution$time, solution$C_liver/1000)
plot(solution$time, solution$C_muscle/1000)
plot(solution$time, solution$C_blood_total/1000)
plot(solution$time, solution$C_brain/1000)
plot(solution$time, solution$C_kidney/1000)
plot(solution$time, solution$C_viscera/1000)
