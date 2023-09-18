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
    #Free = 3.2e-02
    
    # Reabsorption coefficients from bile to intestine
    # estimated by Cao et al., 2022
    # K_urine = Cl_urine/f_reab_urine estimated by Ng et al., 2013 (unitless)
    if(substance=='PFOA'){
      a <- 0.138 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.30 # Cao et al., 2022
      K_urine <- 2.08 
      Cl_urine <- 0.029*3600 # 1/h (Sun et al., 2022)
      Free <- 0.385
    }else if(substance=='PFNA'){
      a <- 0.522 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.34 # Cao et al., 2022
      K_urine <- 1.35
      Cl_urine <- 0.050*3600 # 1/h (Sun et al., 2022)
      Free <- 0.622
    }else if(substance=='PFBS'){
      a <- 0.0598 # Goeritz et al.2013
      f_reab_hep <- 0.23 # Cao et al., 2022
      K_urine <- 5.88
      Cl_urine <- 0.023*3600 # 1/h (Sun et al., 2022) # Assumed equal to PFHxS
      Free <- 0.1 # assumed
    }else if(substance=='PFHxS'){
      a <- 0.558 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.30 # Cao et al., 2022
      K_urine <- 5.88
      Cl_urine <- 0.023*3600 # 1/h (Sun et al., 2022)
      Free <- 0.217
    }else if(substance=='PFOS'){
      a <- 0.721 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.42 # Cao et al., 2022
      K_urine <- 1.35
      Cl_urine <- 0.050*3600 # 1/h (Sun et al., 2022)
      Free <- 0.819
    }
    
    # Bile flow coefficient
    Q_bile_coef <- 7.5e-05 # ml/g BW/h Grosell et al., 2000
    Q_urine_coef <- 2.755e-03 # ml/h/g of BW Urinary flow rate (Curtis et al., 1981)
    V_urine_coef <- 2.2e-03 # ml/g of BW Urine volume inside urinary bladder (Curtis et al., 1981)
    
    a_skin <- 0.9 # 90% of venous blood of skin was assumed to flow directly to kidney (Nichols et al.1996)
    a_muscle <- 0.6 # 60% of venous blood of muscle was assumed to flow directly to kidney (Nichols et al.1996)
    
    plasma <- 0.7
    
    return(list('F_card_ref'=F_card_ref, 'BW_ref'=BW_ref, 'KT'=KT, 
                'admin.dose_dietary'=admin.dose_dietary, 
                'admin.time_dietary'=admin.time_dietary,
                
                'fw_Liver'=fw_Liver, 'fw_Blood'=fw_Blood, 'fw_Skin'=fw_Skin, 
                'fw_Muscle'=fw_Muscle, 'fw_Gills'=fw_Gills, 'fw_Kidney'=fw_Kidney,
                'fw_Viscera'=fw_Viscera, 'fw_lumen'=fw_lumen, 
                
                'fb_Liver'=fb_Liver, 'fb_Skin'=fb_Skin, 'fb_Muscle'=fb_Muscle,
                'fb_Gills'=fb_Gills, 'fb_Kidney'=fb_Kidney, 'fb_Viscera'=fb_Viscera,
                
                'a_skin'=a_skin, 'a_muscle'=a_muscle,
                'Q_bile_coef'=Q_bile_coef,
                'Q_urine_coef'=Q_urine_coef, 'V_urine_coef'=V_urine_coef,
                'K_urine'=K_urine, 'Cl_urine'=Cl_urine,
                'f_reab_hep'=f_reab_hep, 'plasma'=plasma, "Free"=1, "a"=a))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    M_art<-0; M_venous<-0;
    M_gills<-0; M_lumen=0; M_lumen_2=0; M_viscera<-0; M_liver<-0; M_kidney<-0;
    M_muscle<-0; M_skin<-0; M_carcass<-0; M_storage<-0; M_urine<-0; M_feces<-0; M_input<-0 
    
    return(c('M_art'=M_art, 'M_venous'=M_venous, 'M_gills'=M_gills, 'M_lumen'=M_lumen,
             'M_lumen_2'=M_lumen_2, 'M_viscera'=M_viscera, 'M_liver'=M_liver, 'M_kidney'=M_kidney, 
             'M_muscle'=M_muscle, 'M_skin'=M_skin, 'M_carcass'=M_carcass,
             'M_storage'=M_storage,
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
    
    # Total cardiac output ml/h considered as plasma flow
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
    w_art <- 1/3*w_blood
    w_venous <- 2/3*w_blood
    w_carcass <- BW - (w_blood/plasma + w_liver + w_skin + w_muscle +
                         w_gills + w_kidney + w_viscera + w_lumen)
    
    # Calculate the regional blood flows - ml/h
    Q_liver <- fb_Liver*Q_total     # Liver blood flow - ml/h
    Q_skin <- fb_Skin*Q_total      # Skin blood flow - ml/h
    Q_muscle <- fb_Muscle*Q_total   # Muscle blood flow - ml/h
    Q_gills <- Q_total #fb_Gills*BW     # Gills blood flow - ml/h
    Q_kidney <- fb_Kidney*Q_total   # Kidney blood flow - ml/h
    Q_viscera <- fb_Viscera*Q_total # Viscera blood flow - ml/h
    Q_carcass <- Q_total - (Q_liver + Q_skin + Q_muscle + 
                              Q_kidney + Q_viscera)
    
    # Calculate the absolute bile flow rate - ml/h
    Q_bile <- Q_bile_coef*BW
    # Calculate Urinary flow rate - ml/h
    Q_urine <- Q_urine_coef*BW
    
    # Calculate urine volume  - ml 
    v_urine <- V_urine_coef*BW
    
    # Calculate f_reab_urine based on Cl_urine and K_urine - 1/h
    f_reab_urine <- Cl_urine*CLU_coef/K_urine
    
    # Tissue concentrations ug PFAS/g tissue
    C_gills <- M_gills/w_gills
    C_viscera <- M_viscera/w_viscera
    C_liver <- M_liver/w_liver
    C_kidney <- M_kidney/w_kidney
    C_muscle <- M_muscle/w_muscle 
    C_skin <- M_skin/w_skin
    C_carcass <- M_carcass/w_carcass
    C_lumen <- (M_lumen+M_lumen_2)/w_lumen
    C_art <- M_art/w_art
    C_venous <- M_venous/w_venous
    C_blood <- (M_art + M_venous)/w_blood
    C_storage <- M_storage/v_urine
    
    # Arterial Blood
    dM_art <- Free*Q_gills*C_gills/P_gills - 
      (Q_viscera + Q_liver + Q_kidney +
         Q_muscle + Q_skin + Q_carcass)*Free*C_art
    
    dM_venous <- - Free*Q_total*C_venous + 
      ((Q_liver + Q_viscera)*C_liver/P_liver +
         (Q_kidney + a_muscle*Q_muscle + a_skin*Q_skin)*C_kidney/P_kidney +
         (1-a_muscle)*Q_muscle*C_muscle/P_muscle +
         (1-a_skin)*Q_skin*C_skin/P_skin + Q_carcass*C_carcass/P_carcass)*Free
    
    # Gills
    dM_gills <- Q_gills*Free*(C_venous - C_gills/P_gills)
    
    dM_input=0
    
    # Viscera lumen - Available PFAS for absorption and elimination
    dM_lumen = - Ku*a*M_lumen - Cl_feces*(1-a)*M_lumen 
    
    # Viscera lumen_2- Unavailable PFAS for absorption. Can be only eliminated.
    dM_lumen_2 = (1-f_reab_hep)*Q_bile*C_liver - Cl_feces*M_lumen_2 
    
    # Viscera tissue
    dM_viscera <- Q_viscera*Free*(C_art - C_viscera/P_viscera) + Ku*a*M_lumen +
      f_reab_hep*Q_bile*C_liver
    
    # Liver
    dM_Liver <- Q_liver*Free*C_art + Q_viscera*Free*C_viscera/P_viscera - 
      (Q_liver + Q_viscera)*Free*C_liver/P_liver - Q_bile*C_liver
    
    # Kidney
    dM_kidney <- Q_kidney*Free*C_art - 
      (Q_kidney + a_muscle*Q_muscle + a_skin*Q_skin)*Free*C_kidney/P_kidney + 
      a_muscle*Q_muscle*Free*C_muscle/P_muscle + 
      a_skin*Q_skin*Free*C_skin/P_skin - Cl_urine*CLU_coef*M_kidney + f_reab_urine*M_storage
    
    # Muscle
    dM_muscle <- Q_muscle*Free*(C_art - C_muscle/P_muscle)
    
    # Skin
    dM_skin <- Q_skin*Free*(C_art - C_skin/P_skin)
    
    # Carcass 
    dM_carcass <- Q_carcass*Free*(C_art - C_carcass/P_carcass)
    
    # Urine storage
    dM_storage <- Cl_urine*CLU_coef*M_kidney - f_reab_urine*M_storage - Q_urine*C_storage
    
    # Urine
    dM_urine <- Q_urine*C_storage
    
    # Feces
    dM_feces <- Cl_feces*((1-a)*M_lumen + M_lumen_2)
    
    Mass_balance <- M_input - (M_art + M_venous + M_gills + M_lumen + M_lumen_2 + 
                                 M_viscera + M_liver + M_kidney + M_muscle + 
                                 M_skin + M_carcass + M_storage + M_urine + M_feces)
    
    return(list(c('dM_art'=dM_art, 'dM_venous'=dM_venous, 
                  'dM_gills'=dM_gills, 'dM_lumen'=dM_lumen, 'dM_lumen_2'=dM_lumen_2,
                  'dM_viscera'=dM_viscera, 'dM_Liver'=dM_Liver, 
                  'dM_kidney'=dM_kidney, 'dM_muscle'=dM_muscle,
                  'dM_skin'=dM_skin, 'dM_carcass'=dM_carcass, 'dM_storage'=dM_storage,
                  'dM_urine'=dM_urine, 'dM_feces'=dM_feces, 'dM_input'=dM_input),
                'C_Gills'=C_gills, 'C_Viscera'=C_viscera,
                'C_Liver'=C_liver, 'C_Kidney'=C_kidney, 'C_Muscle'=C_muscle,
                'C_Skin'=C_skin, 'C_Carcass'=C_carcass, 'C_Lumen'=C_lumen,
                'C_Blood'=C_blood*plasma,
                'Mass_balance'=Mass_balance, 'BW'=BW))
  })
}


obj_f <- function(x, constant_theta, constant_theta_names, params_names, constant_params=NULL,
                  data_df, errors_df){
  
  if(!is.null(constant_theta)){
    if(length(constant_theta_names) != length(constant_theta)){
      stop("The constant_theta_names vector must be of equal length with the constant_theta vector")
    }
    for (j in 1:length(constant_theta)){
      assign(constant_theta_names[j], 10^constant_theta[[j]])
    }  
  }
  # Assign the values of the x vector to the corresponding parameters
  if(length(x) != length(params_names)){
    stop("The params_names must be of equal length with the x vector")
  }
  for (k in 1:length(x)) {
    assign(params_names[k], 10^x[k])
  }
  
  if(!is.null(constant_params)){
    for (k in 1:length(constant_params)) {
      assign(names(constant_params)[k], constant_params[[k]])
    }
  }
  
  params <- create.params(user.input)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- seq(0,56*24,2)
  # Time of measurement of selected PFAS 
  exp_time <- data_df[,1]
  
  pbk_params <- c('P_liver'=P_liver[[1]], 'P_muscle'=P_muscle[[1]], 'P_kidney'=P_kidney[[1]], 
                  'P_skin'=P_skin[[1]], 'P_gills'=P_gills[[1]], 'P_carcass'=P_carcass[[1]],
                  'P_viscera'=P_viscera[[1]], 'Cl_feces'=Cl_feces[[1]], 'CLU_coef'=CLU_coef[[1]],
                  'Ku'=Ku[[1]])
  
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = c(pbk_params,params),
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  
  
  output_variables <- c('time','C_Liver', 'C_Blood', 
                        'C_Skin', 'C_Muscle', 'C_Gills',
                        'C_Kidney', 'C_Carcass')
  
  if(sum(solution$time %in% exp_time) == length(exp_time)){
    results <- solution[which(solution$time %in% exp_time), output_variables]
    results[,-1] <- results[,-1]
  }else{
    stop(print("Length of predictions is not equal to the length of data"))
  }
  
  # Transform predictions and data_df in list data type
  observations <- list()
  predictions <- list()
  weights_list <- list()
  scores <- c()
  for (z in 2:dim(data_df)[2]) {
    observations[[z-1]] <- data_df[,z]
    names(observations)[z-1] <- colnames(data_df)[z] 
    predictions[[z-1]] <- results[,z]*1000 # transform from ug/g to ug/kg
    names(predictions)[z-1] <- colnames(results)[z]
    weights_list[[z-1]] <- errors_df[,z]
    names(weights_list)[z-1] <- colnames(errors_df)[z]
  }
  #score <- PBKtools::PBKOF(observations,predictions)
  score <- PBKtools::WSSR(observations,predictions,weights_list)
  
  return(score)
}

profile_likelihood_plots <- function(analysis_results, thetas, global_optimum, alpha = 0.95,
                                     df = 1){
  output <- analysis_results$Likelihood_profiles
  plot_list <- list()
  for (i in 1:length(output)) {
    data_to_plot <- output[[i]]$plik
    current_param <- names(data_to_plot)[1]
    names(data_to_plot)[1] <- "Parameter"
    optimal_value <- data.frame(log10(thetas[i]), global_optimum)
    names(optimal_value) <- c("Parameter", "Likelihood")
    data_to_plot[,1] <- data_to_plot[,1]
    
    plot <- ggplot()+
      geom_hline(yintercept=global_optimum + qchisq(alpha,df), linetype="dashed", color = "red", size=1)+
      #geom_hline(yintercept=global_optimization$objective + qchisq(0.95,1), linetype="dashed", color = "green", size=1)+
      geom_hline(yintercept=global_optimum , linetype="dashed", color = "blue", size=1)+
      geom_line(data = data_to_plot,  aes(x=Parameter, y=Likelihood), color = 'black', size=2)+
      #geom_smooth(data = data_to_plot,  aes(x=Parameter, y=Likelihood), method = "loess", span = 0.5, se =0, color = 'black', size=2)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=11)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, colour="pink", size=10)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=5)+
      
      #scale_y_log10()+
      #ylim(c(6,NA))+
      
      labs(#title = paste0( current_param), #"Profile Likelihood of ",
        y = expression(paste(chi^2, "(", theta, ")")) , x = current_param)+
      theme(plot.title = element_text(hjust = 0.5,size=30), 
            axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.y=element_text(size=22),
            axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.x=element_text(size=22),
            legend.title=element_text(hjust = 0.5,size=25), 
            legend.text=element_text(size=22),
            panel.border = element_rect(colour = "black", fill=NA, size=1.0))
    
    #print(plot)
    plot_list[[i]] <- plot
  }
  # Arrange and print the plots in a grid
  grid.arrange(grobs = plot_list, nrow = 3)
}

################################################################################
################################################################################
################################################################################
#                           END OF FUNCTIONS                                   #
################################################################################
################################################################################
################################################################################


################################
#   Identifiability Analysis   #
################################

# Experimental data from Falk et al.2015
#---------------------------------------
# The concentrations in the data are given in ug PFAS/kg tissue units. 
# The time is given in days and will be transformed in hours, to be compatible
# with the model

# Directory of folder with saved data files
data_dir <- '/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015'

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

# Generate the errors
# Consider a Coefficient of Variation of the data points (CV = sd/mean)
errors_list <- list()
CV <- 100/100
for (k in 1:5) {
  data_df <- data_list[[k]]
  errors_list[[k]] <- data.frame(matrix(NA, nrow = nrow(data_df), ncol = ncol(data_df)))
  for (i in 1:nrow(data_df)) {
    for (j in 2:ncol(data_df)) {
      set.seed(100)
      errors_list[[k]][i,j] <- data_df[i,j]*CV #abs(rnorm(1, data_df[i,j]*CV, 1))
    }
  }
  errors_list[[k]][,1] <- data_df[,1]
  colnames(errors_list[[k]]) <- colnames(data_list[[1]]) 
}
names(errors_list) <- names(data_list)

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Identifiability_Analysis')
# Those are the PFOS parameters
x_opt <- c(0.773479754, 0.002495272, 0.680588549, 1.4247478, 0.09900998, 0.12376605,
           0.2668910, 0.2168347, 0.10476366, 3.802160470)
names(x_opt) <- c('Ku', 'CLU_coef', 'Cl_feces', 'P_liver', 'P_muscle', 'P_kidney', 
                  'P_skin', 'P_gills', 'P_carcass', 'P_viscera')

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
user.input <- list('substance'='PFOS',
                   'Texp'=15,
                   'admin.dose_dietary'=admin.dose_dietary,
                   'admin.time_dietary'=admin.time_dietary)


thetas <- log10(x_opt[-c(1:3)])
thetas_names <- c('P_liver', 'P_muscle', 'P_kidney', 
                  'P_skin', 'P_gills', 'P_carcass', 'P_viscera')

names(thetas) <- thetas_names

# lower bounds of parameters
lb <- log10(rep(1e-07,length(thetas)))
# upper bounds of parameters
ub <- log10(rep(1e05, length(thetas)))

names(thetas) <- thetas_names

constant_params <- c('Ku'=x_opt[[1]], 'CLU_coef'=x_opt[[2]], 'Cl_feces'=x_opt[[3]])
substance <- "PFOS" 
Texp <- 15
exported_to_cluster = list("fish_weight"=fish_weight,
                           "create.inits"=create.inits,
                           "create.params"=create.params,
                           "create.events"=create.events,
                           "ode.func"=ode.func,
                           "WSSR"=WSSR,
                           "substance"=substance, 
                           "Texp"=Texp,
                           "admin.dose_dietary"=admin.dose_dietary,
                           "admin.time_dietary"=admin.time_dietary,
                           "user.input"=user.input)


global_optimum = 5.38859
test <- Identifiability_analysis(obj_f = obj_f,
                                 thetas=thetas,
                                 thetas_names=thetas_names,
                                 data_df=data_list[[1]] ,
                                 errors_df=errors_list[[1]],
                                 lb=lb ,
                                 ub=ub,
                                 N_samples = 50,
                                 alpha = 0.95 ,
                                 df = 1,
                                 q = 0.20,
                                 global_optimum = global_optimum ,
                                 min_step_coef = 0.01 ,
                                 max_step_coef = 0.25,
                                 N_cores = 6,
                                 constant_params = constant_params,
                                 exported_to_cluster = exported_to_cluster,
                                 break_at_bounds = TRUE,
                                 # nlopt settings for the main optimization problem
                                 opts = list("algorithm" = "NLOPT_LN_SBPLX",
                                             "xtol_rel" = 1e-03, 
                                             "ftol_rel" = 1e-03,
                                             "ftol_abs" = 0.0,
                                             "xtol_abs" = 0.0 ,
                                             "maxeval" = 500,
                                             "print_level" = 0),
                                 # nlopt settings for the estimation of theta_step
                                 opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                        "xtol_rel" = 1e-03,
                                                        "ftol_rel" = 1e-03,
                                                        "ftol_abs" = 0.0,
                                                        "xtol_abs" = 0.0 ,
                                                        "maxeval" = 50,
                                                        "print_level" = 0),
                                 create_txt = TRUE)

profile_likelihood_plots(analysis_results = test, 10^thetas, global_optimum, alpha = 0.95,
                         df = 1)