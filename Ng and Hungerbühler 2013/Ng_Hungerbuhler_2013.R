library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    
    # Ng and Hungerbuhler 2013 Table S1
    # Tissue Sub-compartments volumes as fractions of total body weight (unitless) 
    f_LT <- 0.0114
    f_KT <- 0.0084
    f_MT <- 0.465
    f_AT <- 0.0089
    # Blood volume coefficient (ml/kg)
    f_blood <- (33+40)/2
    
    # Tissue Sub-compartments volumes - L transformed to m^3
    V_LT <- f_LT*BW*1e-03 # Liver tissue
    V_KT <- f_KT*BW*1e-03 # Kidney tissue
    V_MT <- f_MT*BW*1e-03 # Muscle tissue
    V_AT <- f_AT*BW*1e-03 # Adipose tissue
    # ml transformed to m^3
    V_B <- f_blood*BW*1e-06 # Blood volume
    
    # Ng and Hungerbuhler 2013 Table S2
    # Interstitial fluid Sub-compartments volumes as fractions of tissue weight (unitless) 
    f_LF <- 0.283
    f_KF <- 0.672
    f_MF <- 0.054
    f_AF <- 0.174
    
    # Interstitial fluid Sub-compartments volumes - m^3
    V_LF <- f_LF*V_LT # Liver interstitial fluid
    V_KF <- f_KF*V_KT # Kidney interstitial fluid
    V_MF <- f_MF*V_MT # Muscle interstitial fluid
    V_AF <- f_AF*V_AT # Adipose interstitial fluid
    
    # Urine Volume - m^3
    V_Ur <- 1.8e-08 # m^3 given in Table s9
    
    # Ng and Hungerbuhler 2013 Table S3
    # Total capillary surface area - cm^2 transformed to m^2
    # Calculated for 
    A_B_W <- 71*1e-04 # Total capillary surface area of gill
    A_L <- 16*1e-04 # Total capillary surface area of liver
    A_K <- 22*1e-04 # Total capillary surface area of kidney
    A_M <- 34*1e-04 # Total capillary surface area of muscle
    A_A <- 23*1e-04 # Total capillary surface area of adipose
    
    # Ng and Hungerbuhler 2013 Table S4
    # Fluid flow rates for 8 g rainbow trout - ml/day transformed to m^3/s
    Q_B_W <- 9500*(1e-06/86400)  
    Q_W_B <- 511*(1e-06/86400) # Gills fluid flow rate
    Q_L <- 4.6*(1e-06/86400) # Liver fluid flow rate
    Q_K <- 21*(1e-06/86400) # Kidney fluid flow rate
    Q_M <- 62*(1e-06/86400) # Muscle fluid flow rate
    Q_A <- 0.5*(1e-06/86400) # Adipose fluid flow rate
    Q_Ur <- 0.53*(1e-06/86400) # Urine flow rates
    
    # Ng and Hungerbuhler 2013 Table S5
    # Albumin and liver fatty acid binding protein concentrations - mmol/L transformed to mol/m^3
    C_Alb_B <- 0.2 # Albumin binding protein concentration in blood  
    C_Alb_LF <- 0.1 # Albumin binding protein concentration in liver fluid
    C_Alb_KF <- 0.1 # Albumin binding protein concentration in kidney fluid
    C_Alb_MF <- 0.06 # Albumin binding protein concentration in muscle fluid
    C_Alb_AF <- 0.03 # Albumin binding protein concentration in adipose fluid
    C_Fapb_LT <- 0.05 # Fatty acid binding protein concentration in liver tissue
    
    if(substance == 'PFOA'){
      
      # Ng and Hungerbuhler 2013 Table S6
      # Parameters associated with PFAA uptake and loss via the gills.
      P_eff <- 1.7e-08 # m/s
      CR_ss_cw <- 1.62 # unitless
      k_W_B <- 1.238*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      k_B_W <- 0.59*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      
      # Ng and Hungerbuhler 2013 Table S7
      # Flux (J) for protein-facilitated uptake and 
      # calculated rate constants (b) for clearance and reabsorption used in our model.
      # J_Oat1 <- 0.46    # nmol/mg/protein/min
      # J_Oat3 <-  0.60   # nmol/mg/protein/min
      # J_Oatp1a1 <- 0.50 # nmol/mg/protein/min
      b_clear <- 0.029  # 1/s
      b_reab <- 0.14    # 1/s
      
      # Association Constant KA for Albumin - Table 2
      KA_ALB <- 3.7e+03 # 1/M
      k_off <- 0.01 # 1/s
      k_ALB_on <- KA_ALB*k_off
      
      # Association Constant KA for FAPB - Table 2
      KA_FAPB <- 5.6e+04 # 1/M
      k_FAPB_on <- KA_FAPB*k_off
      
    }else if(substance == 'PFDA'){
      
      # Ng and Hungerbuhler 2013 Table S6
      # Parameters associated with PFAA uptake and loss via the gills.
      P_eff <- 9.5e-08 # m/s
      CR_ss_cw <- 15.95 # unitless
      k_W_B <- 6.72*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      k_B_W <- 0.30*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      
      # Ng and Hungerbuhler 2013 Table S7
      # Flux (J) for protein-facilitated uptake and 
      # calculated rate constants (b) for clearance and reabsorption used in our model.
      # J_Oat1 <- 0.76    # nmol/mg/protein/min
      # J_Oat3 <-  1.04   # nmol/mg/protein/min
      # J_Oatp1a1 <- 1.56 # nmol/mg/protein/min
      b_clear <- 0.049  # 1/s
      b_reab <- 0.042   # 1/s
      
      # Association Constant KA for Albumin - Table 2
      KA_ALB <- 4.7e+04 # 1/M
      k_off <- 0.01 # 1/s
      k_ALB_on <- KA_ALB*k_off
      
      # Association Constant KA for FAPB - Table 2
      KA_FAPB <- 6.4+05 # 1/M
      k_FAPB_on <- KA_FAPB*k_off
      
    }else if(substance == 'PFUnA'){
      
      # Ng and Hungerbuhler 2013 Table S6
      # Parameters associated with PFAA uptake and loss via the gills.
      P_eff <- 1.2e-07 # m/s
      CR_ss_cw <- 52.9 # unitless
      k_W_B <- 8.68*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      k_B_W <- 0.12*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      
      # Ng and Hungerbuhler 2013 Table S7
      # Flux (J) for protein-facilitated uptake and 
      # calculated rate constants (b) for clearance and reabsorption used in our model.
      # J_Oat1 <- 0.89    # nmol/mg/protein/min
      # J_Oat3 <-  1.42   # nmol/mg/protein/min
      # J_Oatp1a1 <- 2.19 # nmol/mg/protein/min
      b_clear <- 0.062  # 1/s
      b_reab <- 0.059    # 1/s
      
      # Association Constant KA for Albumin - Table 2
      KA_ALB <- 4.3e+04 # 1/M
      k_off <- 0.01 # 1/s
      k_ALB_on <- KA_ALB*k_off
      
      # Association Constant KA for FAPB - Table 2
      KA_FAPB <- 2.2+06 # 1/M
      k_FAPB_on <- KA_FAPB*k_off
      
    }else if(substance == 'PFDoA'){
      
      # Ng and Hungerbuhler 2013 Table S6
      # Parameters associated with PFAA uptake and loss via the gills.
      P_eff <- 1.5e-07 # m/s
      CR_ss_cw <- 175  # unitless
      k_W_B <- 10.9*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      k_B_W <- 0.04*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      
      # Ng and Hungerbuhler 2013 Table S7
      # Flux (J) for protein-facilitated uptake and 
      # calculated rate constants (b) for clearance and reabsorption used in our model.
      # J_Oat1 <- 1.01    # nmol/mg/protein/min
      # J_Oat3 <-  1.66   # nmol/mg/protein/min
      # J_Oatp1a1 <- 2.70 # nmol/mg/protein/min
      b_clear <- 0.072  # 1/s
      b_reab <- 0.073   # 1/s
      
      # Association Constant KA for Albumin - Table 2
      KA_ALB <- 1.2e+06 # 1/M
      k_off <- 0.01 # 1/s
      k_ALB_on <- KA_ALB*k_off
      
      # Association Constant KA for FAPB - Table 2
      KA_FAPB <- 7.4+06 # 1/M
      k_FAPB_on <- KA_FAPB*k_off
      
    }else if(substance == 'PFHxS'){
      
      # Ng and Hungerbuhler 2013 Table S6
      # Parameters associated with PFAA uptake and loss via the gills.
      P_eff <- 5.9e-09 # m/s
      CR_ss_cw <- 0.39 # unitless
      k_W_B <- 0.42*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      k_B_W <- 0.84*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      
      # Ng and Hungerbuhler 2013 Table S7
      # Flux (J) for protein-facilitated uptake and 
      # calculated rate constants (b) for clearance and reabsorption used in our model.
      # J_Oat1 <- 0.41    # nmol/mg/protein/min
      # J_Oat3 <-  0.44   # nmol/mg/protein/min
      # J_Oatp1a1 <- 0.13 # nmol/mg/protein/min
      b_clear <- 0.023  # 1/s
      b_reab <- 0.0035   # 1/s
      
      # Association Constant KA for Albumin - Table 2
      KA_ALB <- 1.2e+04 # 1/M
      k_off <- 0.01 # 1/s
      k_ALB_on <- KA_ALB*k_off
      
      # Association Constant KA for FAPB - Table 2
      KA_FAPB <- 1.7+04 # 1/M
      k_FAPB_on <- KA_FAPB*k_off
      
    }else if(substance == 'PFOS'){
      
      # Ng and Hungerbuhler 2013 Table S6
      # Parameters associated with PFAA uptake and loss via the gills.
      P_eff <- 6.1e-08 # m/s
      CR_ss_cw <- 4.30 # unitless
      k_W_B <- 4.37*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      k_B_W <- 0.79*(1.157*1e-08) # L/kg/day transformed to m^3/kg/s
      
      # Ng and Hungerbuhler 2013 Table S7
      # Flux (J) for protein-facilitated uptake and 
      # calculated rate constants (b) for clearance and reabsorption used in our model.
      # J_Oat1 <- 0.66    # nmol/mg/protein/min
      # J_Oat3 <-  1.20   # nmol/mg/protein/min
      # J_Oatp1a1 <- 0.13 # nmol/mg/protein/min
      b_clear <- 0.050  # 1/s
      b_reab <- 0.0035    # 1/s
      
      # Association Constant KA for Albumin - Table 2
      KA_ALB <- 8.9e+03 # 1/M
      k_off <- 0.01 # 1/s
      k_ALB_on <- KA_ALB*k_off
      
      # Association Constant KA for FAPB - Table 2
      KA_FAPB <- 1.9+05 # 1/M
      k_FAPB_on <- KA_FAPB*k_off
    }
    
    k_ALB_off <- k_off
    k_FAPB_off <- k_off
    
    b_KT_Ur <- 0
    b_Ur_KT <- 0
    
    return(list('V_LT'=V_LT, 'V_KT'=V_KT, 'V_MT'=V_MT, 'V_AT'=V_AT, 'V_B'=V_B,
                'V_LF'=V_LF, 'V_KF'=V_KF, 'V_MF'=V_MF, 'V_AF'=V_AF, 'V_Ur'=V_Ur,
                'A_B_W'=A_B_W, 'A_L'=A_L, 'A_K'=A_K, 'A_M'=A_M, 'A_A'=A_A,
                'Q_B_W'=Q_B_W, 'Q_W_B'=Q_W_B, 'Q_L'=Q_L, 'Q_K'=Q_K, 'Q_M'=Q_M, 'Q_A'=Q_A, 'Q_Ur'=Q_Ur,
                'C_Alb_B'=C_Alb_B, 'C_Alb_LF'=C_Alb_LF, 'C_Alb_KF'=C_Alb_KF, 
                'C_Alb_MF'=C_Alb_MF, 'C_Alb_AF'=C_Alb_AF, 'C_Fapb_LT'=C_Fapb_LT,
                'P_eff'=P_eff, 'CR_ss_cw'=CR_ss_cw, 'k_W_B'=k_W_B, 'k_B_W'=k_B_W,
                #'J_Oat1'=J_Oat1, 'J_Oat3'=J_Oat3, 'J_Oatp1a1'=J_Oatp1a1,
                'b_clear'=b_clear, 'b_reab'=b_reab,
                'KA_ALB'=KA_ALB, 'k_off'=k_off, 'k_ALB_on'=k_ALB_on,
                'k_ALB_off'=k_ALB_off, 'k_FAPB_off'=k_FAPB_off,
                'KA_FAPB'=KA_FAPB, 'k_FAPB_on'=k_FAPB_on,
                'b_KT_Ur'=b_KT_Ur, 'b_Ur_KT'=b_Ur_KT,
                "admin.dose" = admin.dose, 
                "admin.time" = admin.time
                ))
    
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    C_W_free<-0;
    M_B_free<-0; M_B_bound<-0;
    M_AF_free<-0; M_AF_bound<-0; M_AT_free<-0;
    M_MF_free<-0; M_MF_bound<-0; M_MT_free<-0; 
    M_LF_free<-0; M_LF_bound<-0; M_LT_free<-0; M_LT_bound<-0;
    M_KF_free<-0; M_KF_bound<-0; M_KT_free<-0;
    M_Ur_free<-0; 
    C_B_ALB<-0; C_AF_ALB<-0; C_MF_ALB<-0; C_LF_ALB<-0;
    C_LT_FAPB<-0; C_KF_ALB<-0
    
    return(c('C_W_free'=C_W_free,
            'M_B_free'=M_B_free, 'M_B_bound'=M_B_bound,
            'M_AF_free'=M_AF_free, 'M_AF_bound'=M_AF_bound, 'M_AT_free'=M_AT_free,
            'M_MF_free'=M_MF_free, 'M_MF_bound'=M_MF_bound, 'M_MT_free'=M_MT_free, 
            'M_LF_free'=M_LF_free, 'M_LF_bound'=M_LF_bound, 'M_LT_free'=M_LT_free, 'M_LT_bound'=M_LT_bound,
            'M_KF_free'=M_KF_free, 'M_KF_bound'=M_KF_bound, 'M_KT_free'=M_KT_free,
            'M_Ur_free'=M_Ur_free, 
            'C_B_ALB'=C_B_ALB, 'C_AF_ALB'=C_AF_ALB, 'C_MF_ALB'=C_MF_ALB, 'C_LF_ALB'=C_LF_ALB,
            'C_LT_FAPB'=C_LT_FAPB, 'C_KF_ALB'=C_KF_ALB))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    ltimes <- length(admin.time)
    lexposure <- length(admin.dose)
    
    events <- data.frame(var = rep('C_W_free', ltimes),
                         time = admin.time,
                         value = admin.dose, method = rep('rep',ltimes))
    return(list(data=events))
  })
}

custom.func <- function(){
  return()
}

ode.func <- function(time, inits, params,custom.func){
  with(as.list(c(inits, params)),{
    
    # Uptake and loss rate constant from gills
    b_W_B <- k_W_B     # m^3/s
    b_B_W <- k_B_W/V_B # 1/s
    
    # Mass Transfer Coefficients - MTC 
    # For diffusion from the interstitial fluid back to blood, the same overall resistance is
    # assumed (that is, kiF-B = kB-iF)
    k_B_AF = (1/Q_A + 1/(P_eff*A_A))^(-1) # adipose MTC
    k_B_MF = (1/Q_M + 1/(P_eff*A_M))^(-1) # muscle MTC
    k_B_LF = (1/Q_L + 1/(P_eff*A_L))^(-1) # liver MTC
    k_B_KF = (1/Q_K + 1/(P_eff*A_K))^(-1) # kidney <TC

    # Rate constants 
    # Adipose
    b_B_AF = k_B_AF/V_B
    b_AF_B = k_B_AF/V_AF
    # Muscle
    b_B_MF = k_B_MF/V_B
    b_MF_B = k_B_MF/V_MF
    # Liver
    b_B_LF = k_B_LF/V_B
    b_LF_B = k_B_LF/V_LF
    # Kidney
    b_B_KF = k_B_KF/V_B
    b_KF_B = k_B_KF/V_KF
    
    # MTC from Interstitial fluid to Tissue sub-compartments
    k_AF_AT <- P_eff*A_A # Adipose
    k_MF_MT <- P_eff*A_M # Muscle
    k_LF_LT <- P_eff*A_L # Liver
    k_KF_KT <- P_eff*A_K # Kidney
    
    # First-order rate constants for diffusion 
    # Adipose
    b_AF_AT <- k_AF_AT/V_AF # from fluid to issue
    b_AT_AF <- k_AF_AT/V_AT # from tissue to fluid
    # Muscle
    b_MF_MT <- k_MF_MT/V_MF # from fluid to issue
    b_MT_MF <- k_MF_MT/V_MT # from tissue to fluid
    # Liver
    b_LF_LT <- k_LF_LT/V_LF # from fluid to issue
    b_LT_LF <- k_LF_LT/V_LT # from tissue to fluid
    # Kidney
    b_KF_KT <- k_KF_KT/V_KF # from fluid to issue
    b_KT_KF <- k_KF_KT/V_KT # from tissue to fluid
    
    # Urinary removal rate constant  - 1/??
    b_Ur <- Q_Ur/V_Ur

    # Concentration
    C_B_bound <- M_B_bound/V_B
    C_B_free <- M_B_free/V_B
    C_AF_bound <- M_AF_bound/V_AF
    C_AF_free <- M_AF_free/V_AF
    C_MF_bound <- M_MF_bound/V_MF
    C_MF_free <- M_MF_free/V_MF
    C_LF_bound <- M_LF_bound/V_LF
    C_LF_free <- M_LF_free/V_LF
    C_LT_bound <- M_LT_bound/V_LT
    C_LT_free <- M_LT_free/V_LT
    C_KF_bound <- M_KF_bound/V_KF
    C_KF_free <- M_KF_free/V_KF
    
    
    # Blood - Unoccupied Albumin binding cites concentration 
    dC_B_ALB <- k_ALB_off*C_B_bound - k_ALB_on*C_B_ALB*C_B_free
    
    # Binding rate constant in Blood - 1/s
    b_B_on <- k_ALB_on*C_B_ALB
    b_B_off <- k_off
    
    # Water concentration
    dC_W_free <- 0 
    
    # Blood - Free
    dM_B_free <- b_W_B*C_W_free - b_B_W*M_B_free - 
      (b_B_AF*M_B_free + b_B_MF*M_B_free + b_B_LF*M_B_free + b_B_KF*M_B_free) +
      (b_AF_B*M_AF_free + b_MF_B*M_MF_free + b_LF_B*M_LF_free + b_KF_B*M_KF_free) -
      b_B_on*M_B_free + b_B_off*M_B_bound
    
    # Blood - Bound
    dM_B_bound <- b_B_on*M_B_free - b_B_off*M_B_bound
    
    # Adipose - Unoccupied Albumin binding cites concentration in fluid
    dC_AF_ALB <- k_ALB_off*C_AF_bound - k_ALB_on*C_AF_ALB*C_AF_free
    
    # Binding rate constant in Adipose fluid - 1/s
    b_AF_on <- k_ALB_on*C_AF_ALB
    b_AF_off <- k_off
    
    # Adipose fluid - Free
    dM_AF_free <- b_B_AF*M_B_free + b_AT_AF*M_AT_free - (b_AF_B + b_AF_AT)*M_AF_free - 
      b_AF_on*M_AF_free + b_AF_off*M_AF_bound
    
    # Adipose fluid - Bound
    dM_AF_bound <- b_AF_on*M_AF_free - b_AF_off*M_AF_bound
    
    # Adipose tissue - Free
    dM_AT_free <- b_AF_AT*M_AF_free - b_AT_AF*M_AT_free
    
    # Muscle - Unoccupied Albumin binding cites concentration in fluid
    dC_MF_ALB <- k_ALB_off*C_MF_bound - k_ALB_on*C_MF_ALB*C_MF_free
    
    # Binding rate constant in Muscle fluid - 1/s
    b_MF_on <- k_ALB_on*C_MF_ALB
    b_MF_off <- k_off
    
    # Muscle Fluid - Free
    dM_MF_free <- b_B_MF*M_B_free + b_MT_MF*M_MT_free - (b_MF_B + b_MF_MT)*M_MF_free - 
      b_MF_on*M_MF_free + b_MF_off*M_MF_bound
    
    # Muscle Fluid - Bound
    dM_MF_bound <- b_MF_on*M_MF_free - - b_MF_off*M_MF_bound
    
    # Muscle Tissue - Free
    dM_MT_free <- b_MF_MT*M_MF_free - b_MT_MF*M_MT_free
    
    # Liver - Unoccupied Albumin binding cites concentration in fluid 
    dC_LF_ALB <- k_ALB_off*C_LF_bound - k_ALB_on*C_LF_ALB*C_LF_free
    # Liver - Unoccupied FAPB binding cites concentration in tissue 
    dC_LT_FAPB <- k_FAPB_off*C_LT_bound - k_FAPB_on*C_LT_FAPB*C_LT_free
    
    # Binding rate constant in Liver fluid - 1/s
    b_LF_on <- k_ALB_on*C_LF_ALB
    b_LF_off <- k_off
    
    # Binding rate constant in Liver Tissue - 1/s
    b_LT_on <- k_FAPB_on*C_LT_FAPB
    b_LT_off <- k_off
    
    # Liver Fluid - Free
    dM_LF_free <- b_B_LF*M_B_free + b_LT_LF*M_LT_free - (b_LF_B + b_LF_LT)*M_LF_free - 
      b_LF_on*M_LF_free + b_LF_off*M_LF_bound
    
    # Liver Fluid - Bound
    dM_LF_bound <- b_LF_on*M_LF_free - b_LF_off*M_LF_bound
    
    # Liver Tissue - Free
    dM_LT_free <- b_LF_LT*M_LF_free - b_LT_LF*M_LT_free - b_LT_on*M_LT_free + b_LT_off*M_LT_bound
    
    # Liver Tissue - Bound
    dM_LT_bound <- b_LT_on*M_LT_free - b_LT_off*M_LT_bound
    
    # Kidney - Unoccupied Albumin binding cites concentration in fluid 
    dC_KF_ALB <- k_ALB_off*C_KF_bound - k_ALB_on*C_KF_ALB*C_KF_free
    
    # Binding rate constant in Kidney fluid - 1/s
    b_KF_on <- k_ALB_on*C_KF_ALB
    b_KF_off <- k_off
    
    # Kidney Fluid - Free
    dM_KF_free <- b_B_KF*M_B_free + b_KT_KF*M_KT_free - (b_KF_B + b_KF_KT)*M_KF_free - 
      b_KF_on*M_KF_free +b_KF_off*M_KF_bound
    
    # Kidney Fluid - bound
    dM_KF_bound <- b_KF_on*M_KF_free - b_KF_off*M_KF_bound
    
    # Kidney Tissue - Free
    dM_KT_free <- b_KF_KT*M_KF_free - (b_KT_KF + b_KT_Ur)*M_KT_free + b_Ur_KT*M_Ur_free -
      b_clear*M_KT_free + b_reab*M_Ur_free

    # Urine 
    dM_Ur_free <- b_KT_Ur*M_KT_free - b_Ur_KT*M_Ur_free + b_clear*M_KT_free - b_reab*M_Ur_free - 
      (Q_Ur/V_Ur)*M_Ur_free
    
    return(list(c('dC_W_free'=dC_W_free,
                  'dM_B_free'=dM_B_free, 'dM_B_bound'=dM_B_bound,
                  'dM_AF_free'=dM_AF_free, 'dM_AF_bound'=dM_AF_bound, 'dM_AT_free'=dM_AT_free,
                  'dM_MF_free'=dM_MF_free, 'dM_MF_bound'=dM_MF_bound, 'dM_MT_free'=dM_MT_free, 
                  'dM_LF_free'=dM_LF_free, 'dM_LF_bound'=dM_LF_bound, 'dM_LT_free'=dM_LT_free, 'dM_LT_bound'=dM_LT_bound,
                  'dM_KF_free'=dM_KF_free, 'dM_KF_bound'=dM_KF_bound, 'dM_KT_free'=dM_KT_free,
                  'dM_Ur_free'=dM_Ur_free, 
                  'dC_B_ALB'=dC_B_ALB, 'dC_AF_ALB'=dC_AF_ALB, 'dC_MF_ALB'=dC_MF_ALB, 'dC_LF_ALB'=dC_LF_ALB,
                  'dC_LT_FAPB'=dC_LT_FAPB, 'dC_KF_ALB'=dC_KF_ALB)))
      
  })
}

################################################################################
BW <- 8/1000 # Total body weight - kg
substance <- 'PFOA'
admin.dose <- c(1000) # administered dose in ug/L or ug
admin.time <- c(0) # time when doses are administered, in hours
user_input <- list('BW'=BW,
                   'substance'=substance,
                   "admin.dose"=admin.dose,
                   "admin.time"= admin.time)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

sample_time <- seq(0,40,1)
solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           events = events, 
                           method="lsodes",rtol = 1e-05, atol = 1e-05)) 

#====================
# Upload on Jaqpot 
#===================
# Subset of features to be displayed on the user interface
predicted.feats <- c('C_W_free',
                     'M_B_free', 'M_B_bound',
                     'M_AF_free', 'M_AF_bound', 'M_AT_free',
                     'M_MF_free', 'M_MF_bound', 'M_MT_free', 
                     'M_LF_free', 'M_LF_bound', 'M_LT_free', 'M_LT_bound',
                     'M_KF_free', 'M_KF_bound', 'M_KT_free',
                     'M_Ur_free', 
                     'C_B_ALB', 'C_AF_ALB', 'C_MF_ALB', 'C_LF_ALB',
                     'C_LT_FAPB', 'C_KF_ALB') 

# Log in Jaqpot server
jaqpotr::login.api()
# Deploy the model on the Jaqpot server to create a web service
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func, 
                     method = "bdf",url = "https://api.jaqpot.org/jaqpot/")
