library(deSolve)

create.params <- funtion(user.input){
  with(as.list(user.input),{
    
    # Ng and Hungerbuhler 2013 Table S1
    # Tissue Sub-compartments volumes as fractions of total body weight (unitless) 
    f_LT <- 0.0114
    f_KT <- 0.0084
    f_MT <- 0.465
    f_AT <- 0.0089
    # Blood volume coefficient (ml/kg)
    f_blood <- (33+40)/2
    
    # Tissue Sub-compartments volumes - ml
    V_LT <- f_LT*BW # Liver tissue
    V_KT <- f_KT*BW # Kidney tissue
    V_MT <- f_MT*BW # Muscle tissue
    V_AT <- f_AT*BW # Adipose tissue
    V_B <- f_blood*BW/1000 # Blood volume
    
    # Ng and Hungerbuhler 2013 Table S2
    # Interstitial fluid Sub-compartments volumes as fractions of tissue weight (unitless) 
    f_LF <- 0.283
    f_KF <- 0.672
    f_MF <- 0.054
    f_AF <- 0.174
    
    # Interstitial fluid Sub-compartments volumes - ml
    V_LF <- f_LF*V_LT # Liver interstitial fluid
    V_KF <- f_KF*V_KT # Kidney interstitial fluid
    V_MF <- f_MF*V_MT # Muscle interstitial fluid
    V_AF <- f_AF*V_AT # Adipose interstitial fluid
    
    # Ng and Hungerbuhler 2013 Table S3
    # Total capillary surface area - cm^2 
    # Calculated for 
    A_G <- 71 # Total capillary surface area of gill
    A_L <- 16 # Total capillary surface area of liver
    A_K <- 22 # Total capillary surface area of kidney
    A_M <- 34 # Total capillary surface area of muscle
    A_A <- 23 # Total capillary surface area of adipose
    
    # Ng and Hungerbuhler 2013 Table S4
    # Fluid flow rates for 8 g rainbow trout - ml/day
    Q_W <- 9500  
    Q_G <- 511 # Gills fluid flow rate
    Q_L <- 4.6 # Liver fluid flow rate
    Q_K <- 21 # Kidney fluid flow rate
    Q_M <- 62 # Muscle fluid flow rate
    Q_A <- 0.5 # Adipose fluid flow rate
    Q_Ur <- 0.53 # Urine flow rates
    
    # Ng and Hungerbuhler 2013 Table S5
    # Albumin and liver fatty acid binding protein concentrations - mmol/L
    C_Alb_B <- 0.2 # Albumin binding protein concentration in blood  
    C_Alb_LF <- 0.1 # Albumin binding protein concentration in liver fluid
    C_Alb_KF <- 0.1 # Albumin binding protein concentration in kidney fluid
    C_Alb_MF <- 0.06 # Albumin binding protein concentration in muscle fluid
    C_Alb_AF <- 0.03 # Albumin binding protein concentration in adipose fluid
    C_Fapb_LT <- 0.05 # Fatty acid binding protein concentration in liver tissue
  })
}

create.inits <- funtion(parameters){
  with(as.list(parameters),{
    
  })
}

create.events <- funtion(parameters){
  with(as.list(parameters),{
    
  })
}

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    # Blood - Free
    dM_B_free <- b_W_B*C_W_free - b_B_W*M_B_free - 
      (b_B_AF*M_B_free + b_B_MF*M_B_free + b_B_LF*M_B_free + b_B_KF*M_B_free) +
      (b_AF_B*M_AF_free + b_MF_B*M_MF_free + b_LF_B*M_LF_free + b_KF_B*M_KF_free) -
      b_B_on*M_B_free + b_B_off*M_B_bound
    
    # Blood - Bound
    dM_B_bound <- b_B_on*M_B_free - b_B_off*M_B_bound
    
    # Adipose fluid - Free
    dM_AF_free <- b_B_AF*M_B_free + b_AT_AF*M_AT_free - (b_AF_B + b_AF_AT)*M_AF_free - 
      b_AF_on*M_AF_free + b_AF_off*M_AF_bound
    
    # Adipose fluid - Bound
    dM_AF_bound <- b_AF_on*M_AF_free - b_AF_off*M_AF_bound
    
    # Adipose tissue - Free
    dM_AT_free <- b_AF_AT*M_AF_free - b_AT_AF*M_AT_free
    
    # Muscle Fluid - Free
    dM_MF_free <- b_B_MF*M_B_free + b_MT_MF*M_MT_free - (b_MF_B + b_MF_MT)*M_MF_free - 
      b_MF_on*M_MF_free + b_MF_off*M_MF_bound
    
    # Muscle Fluid - Bound
    dM_MF_bound <- b_MF_on*M_MF_free - - b_MF_off*M_MF_bound
    
    # Muscle Tissue - Free
    dM_MT_free <- b_MF_MT*M_MF_free - b_MT_MF*M_MT_free
    
    # Liver Fluid - Free
    dM_LF_free <- b_B_LF*M_B_free + b_LT_LF*M_LT_free - (b_LF_B + b_LF_LT)*M_LF_free - 
      b_LF_on*M_LF_free + b_LF_off*M_LF_bound
    
    # Liver Fluid - Bound
    dM_LF_bound <- b_LF_on*M_LF_free - b_LF_off*M_LF_bound
    
    # Liver Tissue - Free
    dM_LT_free <- b_LF_LT*M_LF_free - b_LT_LF*M_LT_free - b_LT_on*M_LT_free + b_LT_off*M_LT_bound
    
    # Liver Tissue - Bound
    dM_LT_bound <- b_LT_on*M_LT_free - b_LT_off*M_LT_bound
    
    # Kidney Fluid - Free
    dM_KF_free <- b_B_KF*M_B_free + b_KT_KF*M_KT_free - (b_KF_B + b_KF_KT)*M_KF_free - 
      b_KF_on*M_KF_free +b_KF_off*M_KF_bound
    
    # Kidney Fluid - bound
    d_M_KF_bound <- b_KF_on*M_KF_free - b_KF_off*M_KF_bound
    
    # Kidney Tissue - Free
    d_M_KT_free <- b_KF_KT*M_KF_free - (b_KT_KF + b_KT_Ur)*M_KT_free + b_Ur_KT*M_Ur_free -
      b_clear*M_KT_free + b_reab*M_Ur_free
    
    # Urine 
    dM_Ur_free <- b_KT_Ur*M_KT_free - b_Ur_KT*M_Ur_free + b_clear*M_KT_free - b_reab_M_Ur_free - 
      (Q_Ur/V_Ur)*M_Ur_free
    
    return(list('dM_B_free'=dM_B_free, 'dM_B_bound'=dM_B_bound, 'dM_AF_free'=dM_AF_free,
                'dM_AF_bound'=dM_AF_bound, 'dM_AT_free'=dM_AT_free, 'dM_MF_free'=dM_MF_free,
                'dM_MF_bound'=dM_MF_bound, 'dM_MT_free'=dM_MT_free, 'dM_LF_free'=dM_LF_free,
                'dM_LF_bound'=dM_LF_bound, 'dM_LT_free'=dM_LT_free, 'dM_LT_bound'=dM_LT_bound,
                'dM_KF_free'=dM_KF_free, 'd_M_KF_bound'=d_M_KF_bound, 'd_M_KT_free'=d_M_KT_free,
                'dM_Ur_free'=dM_Ur_free))
    
  })
}

################################################################################
BW = 8 # Total body weight - g
