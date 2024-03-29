library(ggplot2)
#===============================================================#
# create.params: Function that estimates the model's parameters #
#===============================================================#

create.params <- function(BW, T_water, sex, pH_w, C_water){
  # BW: total body weight of fish (g)
  # T_water: Temperature of water (K)
  # sex: Sex of fish
  # pH_w: pH of water
  
  F_card <- exp(25.5)*exp(-5790/T_water)*BW^0.75 # Cardiac output (ml/day) - Wang et al. (2022)
  VO2 <- exp(19.1)*exp(-5040/T_water)*BW^0.75 # Oxygen Consumption Rate (mg O2/day) - Wang et al. (2022)
  Cox <- exp(-10.6)*exp(1743.3/T_water) # Dissolved oxygen concentration (mg O2 /ml) - Wang et al. (2022)
  
  # Estimate tissues volumes
  V_ad <- 0.020*BW^1.24 # Adipose volume (ml): Wang et al. (2022) Table S2
  V_bl <- 0.030*BW^1.02 # Blood volume (ml): Wang et al. (2022) Table S2
  V_br <- 0.017*BW^0.61 # Brain volume (ml): Wang et al. (2022) Table S2
  V_ki <- 0.014*BW^0.87 # Kidney volume (ml): Wang et al. (2022) Table S2
  V_li <- 0.018*BW^0.95 # Liver volume (ml): Wang et al. (2022) Table S2
  V_git <- 0.068*BW^0.96 # GIT volume (ml): Wang et al. (2022) Table S2
  if(sex== "Male"){
    V_go <- 0.012*BW^1.06 # Male Gonads volume (ml): Wang et al. (2022) Table S2  
  }else if (sex == "Female"){
    V_go <- 0.13*BW^1.04 # Female Gonads volume (ml): Wang et al. (2022) Table S2  
  }
  V_rp <- 0.0025*BW^1.04 # Richly Perfused tissues volume (ml): Wang et al. (2022) Table S2  
  V_skin <- 0.1*BW^0.67 # Skin volume (ml): Wang et al. (2022) Table S2  
  V_skeleton <- 0.033*BW^1.03 # Skeleton volume (ml): Wang et al. (2022) Table S2 
  
  # Data for gill mass estimation --> Zhang et al. 2019
  f_gill <- mean(c(0.0078,	0.0086,	0.01,	0.0114,	0.0129,	0.0102)) /mean(c(0.5469,	0.5378,	0.666,	0.6587,	0.5857,	0.5521))
  V_gill <- f_gill*BW
  
  V_muscle <- BW - (V_ad + V_bl + V_br + V_ki + V_li + V_git + V_go + V_rp +
                      V_skin + V_skeleton + V_gill)
  
  # As poorly perfused compartment is considered the sum of skeleton and muscle
  V_pp <- V_skeleton + V_muscle
  
  # Estimate tissue blood flows
  # Weighting factors - Wang et al. (2022) Table 1
  wf_ad <- 0.626
  wf_br <- 3.600
  wf_git <- 4.846
  wf_go <- 3.367
  wf_ki <- 7.005
  wf_li <- 2.230
  wf_pp <- 0.733  
  wf_rp <- 3.600
  wf_skin <- 0.570
  
  volumes.vector <- c(V_ad, V_br, V_ki, V_li, V_git, V_go, V_rp, V_skin, V_pp)
  wf.vector <- c(wf_ad, wf_br, wf_ki, wf_li, wf_git, wf_go, wf_rp, wf_skin, wf_pp)
  
  F_ad <- ((wf_ad * V_ad/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Adipose blood flow (ml/day), Wang et al. (2022)
  F_br <- ((wf_br * V_br/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Brain blood flow (ml/day), Wang et al. (2022)
  F_git <- ((wf_git * V_git/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # GIT blood flow (ml/day), Wang et al. (2022)
  F_go <- ((wf_go * V_go/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Gonads blood flow (ml/day), Wang et al. (2022)
  F_ki <- ((wf_ki * V_ki/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Kidney blood flow (ml/day), Wang et al. (2022)
  F_li <- ((wf_li * V_li/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Liver blood flow (ml/day), Wang et al. (2022)
  F_pp <- ((wf_pp * V_pp/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Poorly Perfused blood flow (ml/day), Wang et al. (2022)
  F_rp <- ((wf_rp * V_rp/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Richly Perfused lood flow (ml/day), Wang et al. (2022)
  F_skin <- ((wf_skin * V_skin/BW)/(sum(wf.vector * volumes.vector/BW)))*F_card # Skin blood flow (ml/day), Wang et al. (2022)
  #F_card - sum(c(F_ad,F_br,F_git,F_go,F_ki,F_li,F_pp,F_rp,F_skin))
  
  pKa = 0.5 # PFOA dissociation constant - Golosovskaia et al.(2024) Table 1
  logKow = 4.35 # PFOA Kow - Golosovskaia et al.(2024) Table 1
  Kow = 10^logKow
  delta_logKow = 3.1  # Golosovskaia et al.(2024) 
  logKow_ion = logKow - delta_logKow # Golosovskaia et al.(2024)
  Kow_ion = 10^logKow_ion 
  
  # Fraction of the neutral molecule at gill surface - Wang et al.(2022)
  fn_gill <- 1/(1 + 10^(1*(pH_w - pKa)))
  
  # Fraction of the neutral molecule in fish body - Wang et al.(2022)
  fn_fish <- 1/(1 + 10^(1*(7.4 - pKa)))
  
  # Octanol-Water distribution coefficient at gill surface pH - Wang et al.(2022)
  Dow_gill <- fn_gill * Kow + (1-fn_gill)*Kow_ion
  
  # Octanol-Water distribution coefficient at fish body pH 7.4 - Wang et al.(2022)
  Dow <- fn_fish * Kow +(1-fn_fish)*Kow_ion
  
  # 1st column: Fraction of neutral lipid content (fnl)
  # 2nd column: Fraction of polar lipid content (fpl)
  # 3nd column: Fraction of non-lipid organic matter (NLOM) content (fnlom)
  # 4th column: Fraction of water (fw)
  # 5th column: Weighting factor (wf)
  tissue_contents <- data.frame(matrix(data = c(
    88.4, 	0.9, 	5.9, 	  4.8,
    #0.8, 	  0.7, 	85.1, 	13.4, # blood contents, not needed
    4.7, 	  4.3, 	68.2, 	22.8,	
    5.2, 	  2.3, 	67.0, 	25.5,	
    4.2, 	  2.4, 	67.2, 	26.2,	
    5.9, 	  2.3, 	70.3, 	21.5,	
    3.7, 	  2.5, 	70.4, 	23.4,	
    2.4, 	  0.9, 	76.9, 	19.9,	
    6.6, 	  0.8, 	58.7, 	33.8,	
    3.7, 	  1.2, 	69.2, 	25.9
  )/100, byrow = T, ncol = 4))
  colnames(tissue_contents) <- c("fnl", "fpl", "fnlom", "fw")
  rownames(tissue_contents) <- c("ad", "br", "git", "go", "ki", "li",
                                 "pp", "rp", "skin")
  
  
  # Estimate Blood-Water Partition coefficient 
  Pbw <- 0.008*0.3*Dow + 0.007*2.0*(Dow^0.94) + 0.134*2.9*(Dow^0.63) + 0.851
  
  # Function to estimate the Partition coefficient of each tissue - Wang et al. (2022)
  P_ib <- function(fnl, fpl, fnlom, fw, Pbw){
    P_iw <- fnl*0.3*Dow + fpl*2.0*(Dow^0.94) + fnlom*2.9*(Dow^0.63)+fw
    P_ib <- P_iw/Pbw
    return(P_ib)
  }
  
  PCs_vector <- c()
  # Estimate tissue:blood partition coefficient of each compartment
  for (i in 1:dim(tissue_contents)[1]) { #loop over all compartments
    current_input <- tissue_contents[i,1:4]
    PCs_vector[i] <- P_ib(current_input[1],current_input[2], current_input[3],current_input[4], Pbw)
    names(PCs_vector)[i] <- paste0("P_",rownames(tissue_contents)[i]) 
    assign(paste0("P_",rownames(tissue_contents)[i]), PCs_vector[[i]])
  }
  
  # Gill ventilation coefficient (L/kg^0.75/day) - Wang et al. (2022)
  y_water <- (VO2/(0.71*Cox))*(1/((1000^0.25) * (BW^0.75)))
  
  # Blood perfusion coefficient (L/kg^0.75/day) - Wang et al. (2022)
  y_blood <- F_card*Pbw*(1/(1000^0.25 * (BW^0.75)))
  
  # Exchange coefficient between blood and water (ml/day) - Wang et al. (2022)
  kx <- (((BW/1000)^0.75) / (2.8 * 10^-3 + (68/Dow_gill) + (1/y_water) + (1/y_blood)))*1000
  
  k_out <- 1/(68*(Dow_gill-1)+1) * kx # Golosovskaia et al. (2024) - Hendriks et al. (2001)
  
  # Plasma:water distribution ratio - Wang et al. (2022)
  D_plasma_w <- 10^(0.75*logKow_ion + 0.58)
  
  # Unbound fraction in blood - Wang et al. (2022)
  UF = 1/(1 + D_plasma_w)
  # Fraction of PPT blood flow to venous blood - Wang et al. (2022)
  a_fpp = 0.4 
  
  # Fraction of skin blood flow to venous blood - Wang et al. (2022)
  a_fs = 0.1
  
  # Fraction of arterial and venous blood volume to total blood - Wang et al. (2022)
  f_art <- 1/3
  f_ven <- 2/3
  
  ku=0
  ke_feces=0
  cl_feces = 0.001 # 1/day random initial value
  CalbB = 0 #0.008 # M
  Ka = 2.2e05  # 1/M
  
  return(list(#Volume of tissues
    "V_ad"=V_ad, "V_bl"=V_bl, "V_br"=V_br, "V_git"=V_git, "V_go"=V_go,
    "V_ki"=V_ki, "V_li"=V_li, "V_pp"=V_pp, "V_rp"=V_rp, 
    "V_skin"=V_skin, "V_skeleton"=V_skeleton, "V_muscle"=V_muscle,
    "V_gill"=V_gill,
    # Blood flows of tissues
    "F_card"=F_card,
    "F_ad"=F_ad, "F_br"=F_br, "F_git"=F_git, "F_go"=F_go,
    "F_ki"=F_ki, "F_li"=F_li, "F_pp"=F_pp, "F_rp"=F_rp, 
    "F_skin"=F_skin,
    # Partition Coefficients
    "P_ad"=P_ad, "P_br"=P_br, "P_git"=P_git, "P_go"=P_go, "P_ki"=P_ki,
    "P_li"=P_li, "P_pp"=P_pp, "P_rp"=P_rp, "P_skin"=P_skin, "Pbw"=Pbw,
    
    "y_water"=y_water, "y_blood"=y_blood, "kx"=kx, "D_plasma_w"=D_plasma_w,
    "UF"=UF, "a_fpp"=a_fpp, "a_fs"=a_fs, "k_out"=k_out,
    
    "f_art"=f_art, "f_ven"=f_ven,
    
    "C_water"=C_water,
    "ku"=ku, "ke_feces"=ke_feces, "cl_feces"=cl_feces, "BW"=BW,
    "CalbB"=CalbB, "Ka"=Ka
  ))  
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    M_gill=0; M_excr=0;
    M_ven=0; M_art=0; M_lumen=0; M_git=0; M_li=0; M_ki=0; M_br=0; M_go=0;
    M_rp=0; M_pp=0; M_ad=0; M_skin=0; C_water=C_water
    
    return(c("M_gill"=M_gill, "M_excr"=M_excr, "M_ven"=M_ven, "M_art"=M_art, "M_lumen"=M_lumen, 
             "M_git"=M_git, "M_li"=M_li, "M_ki"=M_ki, "M_br"=M_br, "M_go"=M_go, 
             "M_rp"=M_rp, "M_pp"=M_pp, "M_ad"=M_ad, "M_skin"=M_skin, 
             "C_water"=C_water))
  })
}

#=====================================#
# ode.func: The ODEs of the PBK model #
#=====================================#

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    # Concentrations (ug/ml)
    C_art <- M_art/(V_bl*f_art)
    C_ven <- M_ven/(V_bl*f_ven)
    C_git <- (M_lumen+M_git)/V_git
    C_li <- M_li/V_li
    C_ki <- M_ki/V_ki
    C_br <- M_br/V_br
    C_go <- M_go/V_go
    C_rp <- M_rp/V_rp
    C_pp <- M_pp/V_pp
    C_ad <- M_ad/V_ad
    C_skin <- M_skin/V_skin
    C_gill <- M_gill/V_gill
    
    # Concentrations of unbound pfoa in the blood of each compartment
    CF_art <- C_art * 1.0 / (1.0 + CalbB * Ka)
    CF_ven <- C_ven * 1.0 / (1.0 + CalbB * Ka)
    CF_git <- C_git * 1.0 / (1.0 + CalbB * Ka)
    CF_li <- C_li * 1.0 / (1.0 + CalbB * Ka)
    CF_ki <- C_ki * 1.0 / (1.0 + CalbB * Ka)
    CF_br <- C_br * 1.0 / (1.0 + CalbB * Ka)
    CF_go <- C_go * 1.0 / (1.0 + CalbB * Ka)
    CF_rp <- C_rp * 1.0 / (1.0 + CalbB * Ka)
    CF_pp <- C_pp * 1.0 / (1.0 + CalbB * Ka)
    CF_ad <- C_ad * 1.0 / (1.0 + CalbB * Ka)
    CF_skin <- C_skin * 1.0 / (1.0 + CalbB * Ka)
    CF_gill <- C_gill * 1.0 / (1.0 + CalbB * Ka)
    
    #dM_i/dt: ug/d
    #dM_gill <- kx*C_water - kx/Pbw*CF_gill + F_card*(CF_ven - CF_gill)
    dM_gill <- kx*C_water - k_out*CF_gill + F_card*(CF_ven - CF_gill)
    dM_excr <- k_out*CF_gill + cl_feces*M_lumen
    # Water concentration: ug/ml
    dC_water = 0
    
    # Venous blood
    dM_ven <- (F_li + F_rp + F_git + F_go)*(CF_li/P_li) + F_br*(CF_br/P_br) +
      F_ad*(CF_ad/P_ad) + (F_ki + (1-a_fpp)*F_pp + (1-a_fs)*F_skin)*(CF_ki/P_ki) +
      a_fpp*F_pp*(CF_pp/P_pp) + a_fs*F_skin*(CF_skin/P_skin) - 
      F_card*CF_ven 
    
    # Arterial blood
    dM_art <- F_card*(CF_gill - CF_art) 
    
    # GIT lumen
    dM_lumen <- -(ku + ke_feces) * M_git
    
    # GIT
    dM_git <- ku*M_lumen + F_git*(CF_art - CF_git/P_git) - cl_feces*M_git
    
    # Liver
    dM_li <- F_li*(CF_art - CF_li/P_li) + F_rp*(CF_rp/P_rp - CF_li/P_li) +
      F_git*(CF_git/P_git - CF_li/P_li) + F_go*(CF_go/P_go - CF_li/P_li)
    
    # Kidney
    dM_ki <- F_ki*(CF_art - CF_ki/P_ki) + (1-a_fpp)*F_pp*(CF_pp/P_pp - CF_ki/P_ki) +
      (1 - a_fs)*F_skin*(CF_skin/P_skin - CF_ki/P_ki) 
    
    # Brain
    dM_br <- F_br*(CF_art - CF_br/P_br) # ug/day
    
    # Gonads
    dM_go <- F_go*(CF_art - CF_go/P_go) # ug/day
    
    # Richly perfused tissues
    dM_rp <- F_rp*(CF_art - CF_rp/P_rp) # ug/day
    
    # Poorly perfused tissues
    dM_pp <- F_pp*(CF_art - CF_pp/P_pp) # ug/day
    
    # Adipose
    dM_ad <- F_ad*(CF_art - CF_ad/P_ad) # ug/day
    
    # Skin
    dM_skin <- F_skin*(CF_art - CF_skin/P_skin) # ug/day
    
    return(list(c("dM_gill"=dM_gill, "dM_excr"=dM_excr, "dM_ven"=dM_ven, "dM_art"=dM_art, 
                  "dM_lumen"=dM_lumen, "dM_git"=dM_git,
                  "dM_li"=dM_li, "dM_ki"=dM_ki, "dM_br"=dM_br, "dM_go"=dM_go, "dM_rp"=dM_rp,
                  "dM_pp"=dM_pp, "dM_ad"=dM_ad, "dM_skin"=dM_skin, "dC_water"=dC_water),
                "C_gill"=C_gill,
                "C_art"=C_art, "C_ven"=C_ven, "C_git"=C_git, "C_li"=C_li,
                "C_ki"=C_ki, "C_br"=C_br, "C_go"=C_go, "C_rp"=C_rp, 
                "C_pp"=C_pp, "C_ad"=C_ad, "C_skin"=C_skin,
                "CF_gill"=CF_gill,
                "CF_art"=CF_art, "CF_ven"=CF_ven, "CF_git"=CF_git, "CF_li"=CF_li,
                "CF_ki"=CF_ki, "CF_br"=CF_br, "CF_go"=CF_go, "CF_rp"=CF_rp, 
                "CF_pp"=CF_pp, "CF_ad"=CF_ad, "CF_skin"=CF_skin
    ))
    
  })
}

#====================#
# Objective Function #
#====================#

obj.func <- function(x, compartments, compartments_last_point, params, rest, 
                     grouped_parameters, fitted_params_names){
  
  mass_names_list<- list("Liver" = "M_li", "Brain" = "M_br", "Gill" = "M_gill",
                         "Git" = "M_git", "Kidney" = "M_ki", "RPT" = "M_rp",
                         "PPT"= "M_pp")
  
  conc_names_list<- list("Liver" = "C_li", "Brain" = "C_br", "Gill" = "C_gill",
                         "Git" = "C_git", "Kidney" = "C_ki", "RPT" = "C_rp",
                         "PPT"= "C_pp")
  
  tweights_names_list <- list("Liver" = "V_li", "Brain" = "V_br", "Gill" = "V_gill",
                              "Git" = "V_git", "Kidney" = "V_ki", "RPT" = "V_rp",
                              "PPT"= "V_pp")
  
  if(length(x) != length(fitted_params_names)){
    stop("Length of x must be equal to length of fitted_params_names")
  }
  
  if(any(compartments %in% compartments_last_point)){
    stop("Objects compartments  and compartments_last_point are not allowed to have common values")
  }
  
  # Distinguish and assign values to correction factors 
  # according to the number of selected groups of parameters
  cf_vector <- c()
  for (i in 1:length(fitted_params_names)) {
    if( !(fitted_params_names[i] %in% names(params)) ){
      cf_vector[i] <- x[i]
    }else{
      params[[fitted_params_names[i]]] <- x[i]
    }
  }
  
  # Apply the correction factors to the corresponding parameters
  for (i in 1:length(grouped_parameters)) {
    for (j in 1:length(grouped_parameters[[i]])) {
      params[[ grouped_parameters[[i]][j] ]] <- params[[ grouped_parameters[[i]][j] ]] * cf_vector[i]
    }
  }
  
  
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  
  keep_results <- solution[solution$time %in% data_df$Time, c("time", "C_br", "C_gill", "C_git",
                                                              "C_ki", "C_li", "C_rp" , "C_pp")]
  colnames(keep_results) <- c("Time", "Brain", "Gill", "Git", "Kidney", "Liver", "RPT", "PPT")
  
  results <- keep_results[compartments]
  exp_data <- data_df[compartments]
  # Transform predictions and data_df in list data type
  observations <- list()
  predictions <- list()
  for (z in 1:dim(exp_data)[2]) {
    observations[[z]] <- exp_data[,z]
    names(observations)[z] <- colnames(exp_data)[z]
    predictions[[z]] <- results[,z]
    names(predictions)[z] <- colnames(results)[z]
  }
  
  if(!is.null(compartments_last_point)){
    time_point <- max(data_df[,1])
    last_point_obs <- data_df[data_df$Time == time_point, which(colnames(data_df) %in% compartments_last_point)]
    last_point_preds <- solution[solution$time == time_point, unname(unlist(conc_names_list[compartments_last_point]))]
    names(last_point_preds) <- names(last_point_obs)
    
    for (i in 1:length(last_point_obs)) {
      current_name <- names(last_point_obs)[i]
      
      observations[[current_name]] <- last_point_obs[[i]]
      predictions[[current_name]] <- last_point_preds[[i]]
      
    }
  }
  
  if(rest){
    # Define the compartments not explicitly fitted on data
    rest_compartments_names <- colnames(keep_results)[- which(colnames(keep_results) %in% c("Time", compartments))]
    
    # Define columns to keep from solution
    solution_cols <- unname(unlist(mass_names_list[rest_compartments_names]))
    
    # PFOA mass predictions of compartments of interest
    mass_pred <- solution[solution$time %in% data_df$Time, colnames(solution) %in% solution_cols ]
    # Rowsums() of mass_pred (The total mass of PFOA in the rest compartments)
    total_rest_mass_pred <- rowSums(mass_pred)
    # Estimate the equivalent concentration
    total_rest_conc_pred <- total_rest_mass_pred / sum(unlist(params[unname(unlist(tweights_names_list[rest_compartments_names]))]))
    
    # Estimate the total concentration of the rest compartments 
    # (excluding the compartments that are explicitly fitted).
    # Keep only the observation (concentrations) for the rest compartments
    rest_data_conc <- data_df[which(colnames(data_df) %in% rest_compartments_names)]
    tissue_weights <- unlist(params[ unname(unlist(tweights_names_list[rest_compartments_names])) ])
    # for (i in 1:length(tissue_weights)) {
    #   names(tissue_weights)[i] <- tweights_names_list
    # }
    #Initialize the mass dataframe
    rest_data_mass <- rest_data_conc
    # Transform concentrations to mass by multiplying each column
    # with the corresponding tissu weight
    for (i in 1:dim(rest_data_conc)[2]) {
      rest_data_mass[,i] <- rest_data_conc[,i] * unlist(params[unname(unlist(tweights_names_list[colnames(rest_data_conc)[i]]))])
    }
    total_rest_mass_obs <- rowSums(rest_data_mass)
    # Estimate total concentration of the rest compartments
    total_rest_conc_obs <- total_rest_mass_obs / sum(unlist(params[unname(unlist(tweights_names_list[rest_compartments_names]))]))
    
    # Add the total_rest_conc_pred and total_rest_conc_obs to 
    # predictions and observations lists
    
    observations[["Rest"]] = total_rest_conc_obs
    predictions[["Rest"]] = total_rest_conc_pred
    
  }
  
  compartments_scores <- c()
  for (i in 1:length(observations)) {
    compartments_scores[i] <- PBKtools::AAFE(observations[i],predictions[i])
  }
  
  return(mean(compartments_scores))
}

prepare_plot<-function(Test){
  # Part 1: plot the training data
  # Training input
  L = 3 # cm
  BW <- 0.0096 * L^3.55 # g Grech et al. (2019)
  T_water = 23 + 273 #K 
  sex = "Male"
  pH_w = 8.0
  C_water = 25 # ug/ml
  params <- create.params(BW, T_water, sex, pH_w, C_water)
  inits <- create.inits(params)
  sample_time <- seq(0,30,0.5)
  
  # Distinguish and assign values to correction factors 
  # according to the number of selected groups of parameters
  cf_vector <- c()
  for (i in 1:length(Test$fitted_params_names)) {
    if( !(Test$fitted_params_names[i] %in% names(params)) ){
      cf_vector[i] <- Test$optimization$solution[i]
    }else{
      params[[Test$fitted_params_names[i]]] <- Test$optimization$solution[i]
    }
  }
  
  # Apply the correction factors to the corresponding parameters
  for (i in 1:length(Test$grouped_parameters)) {
    for (j in 1:length(Test$grouped_parameters[[i]])) {
      params[[ Test$grouped_parameters[[i]][j] ]] <- params[[ Test$grouped_parameters[[i]][j] ]] * cf_vector[i]
    }
  }
  
  final_solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                            method="lsodes",rtol = 1e-05, atol = 1e-05))
  predictions <- final_solution[,c("time", "C_br", "C_gill", "C_git", "C_ki", "C_li", "C_rp", "C_pp")]
  colnames(predictions) <- colnames(data_df)
  
  #color_codes <- scales::hue_pal()(2)
  #names(color_codes) <-  c("Predicted", "Observed")
  color_codes <- c("Predictions" = "#000000", "Observations" = "#E69F00")
  
  p1 <- create.plot("Brain", predictions, data_df, color_codes)
  p2 <- create.plot("Gill", predictions, data_df, color_codes)
  p3 <- create.plot("Git", predictions, data_df, color_codes)
  p4 <- create.plot("Kidney", predictions, data_df, color_codes)
  p5 <- create.plot("Liver", predictions, data_df, color_codes)
  p6 <- create.plot("RPT", predictions, data_df, color_codes)
  p7 <- create.plot("PPT", predictions, data_df, color_codes)
  
  grid_plot <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7,
                                 ncol=4, nrow=2, 
                                 common.legend = TRUE, legend="bottom")
  #print(grid_plot)
  
  ##############################
  
  # Part 2: plot the validation data (Ulhaq et al. (2015))
  # Validation input
  L = 3 # cm
  BW <- 0.0096 * L^3.55 # g Grech et al. (2019)
  T_water = 25 + 273 #K 
  pH_w = 8.0
  # The water concentration is reported in ug/L and is transformed to 
  # ug/ml
  C_water = c(0.3, 1, 3, 10, 30)/1000 # ug/ml
  sample_time <- seq(0,40,1) # days
  
  # Initiate a dataframe to keep the predictions for each
  # simulated concentration
  
  predictions_male <- data.frame(matrix(NA, ncol=5))
  colnames(predictions_male) <- c("Concentration", "Brain", "Git", "Liver", "Gonads")
  predictions_female <- predictions_male
  T_water<- 25 + 273 #K
  
  
  for (i in 1:length(C_water)) {
    params_male <- create.params(BW, T_water, "Male", pH_w, C_water[i])
    params_female <- create.params(BW, T_water, "Female", pH_w, C_water[i])
    
    # Distinguish and assign values to correction factors 
    # according to the number of selected groups of parameters
    cf_vector <- c()
    for (k in 1:length(Test$fitted_params_names)) {
      if( !(Test$fitted_params_names[k] %in% names(params_male)) ){
        cf_vector[k] <- Test$optimization$solution[k]
      }else{
        params_male[[Test$fitted_params_names[k]]] <- Test$optimization$solution[k]
        params_female[[Test$fitted_params_names[k]]] <- Test$optimization$solution[k]
      }
    }
    
    # Apply the correction factors to the corresponding parameters
    for (z in 1:length(Test$grouped_parameters)) {
      for (j in 1:length(Test$grouped_parameters[[z]])) {
        params_male[[ Test$grouped_parameters[[z]][j] ]] <- params_male[[ Test$grouped_parameters[[z]][j] ]] * cf_vector[z]
        params_female[[ Test$grouped_parameters[[z]][j] ]] <- params_female[[ Test$grouped_parameters[[z]][j] ]] * cf_vector[z]
        
      }
    }

    inits <- create.inits(params_male)
    validation_solution_male <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params_male,
                                                        method="lsodes",rtol = 1e-05, atol = 1e-05))
    inits <- create.inits(params_female)
    validation_solution_female <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params_female,
                                                          method="lsodes",rtol = 1e-05, atol = 1e-05))
    
    predictions_male[i,] <- c(C_water[i]*1000, tail(validation_solution_male[c("C_br", "C_git", "C_li", "C_go")], 1) )  
    predictions_female[i,] <- c(C_water[i]*1000, tail(validation_solution_female[c("C_br", "C_git", "C_li", "C_go")], 1) )  
    
  }
  
  # Prepare a dataframe with the male observations and predictions
  results_male <- data.frame(matrix(NA, ncol = 4))
  colnames(results_male) <- c("Concentration", "Tissue", "Observed", "Predicted")
  
  #Loop over all tissues - Male data
  df <- validation_data$male_data
  for (j in 2:dim(df)[2]) {
    added_data <- cbind.data.frame(df$Concentration, colnames(df)[j] ,df[,j], predictions_male[, colnames(df)[j]])
    colnames(added_data) <- c("Concentration", "Tissue", "Observed", "Predicted")
    
    results_male <- rbind(results_male, added_data)
    
  }
  # Remove the first row of NA values
  results_male<-results_male[-1,]
  results_male$Concentration <- as.character(results_male$Concentration)
  
  
  # Prepare a dataframe with the female observations and predictions
  
  results_female <- data.frame(matrix(NA, ncol = 4))
  colnames(results_female) <- c("Concentration", "Tissue", "Observed", "Predicted")
  
  for (j in 2:dim(validation_data$female_data)[2]) {
    df <- validation_data$female_data
    
    added_data <- cbind.data.frame(df$Concentration, colnames(df)[j] ,df[,j], predictions_female[, colnames(df)[j]])
    colnames(added_data) <- c("Concentration", "Tissue", "Observed", "Predicted")
    
    results_female <- rbind(results_female, added_data)
    
  }
  # Remove the first row of NA values
  results_female<-results_female[-1,]
  results_female$Concentration <- as.character(results_female$Concentration)
  
  Tissue_markers <-  c(0,1,2,5)
  names(Tissue_markers) <- c( "Brain", "Git", "Liver", "Gonads")
  
  Concentration <- c("#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#CD3333") #scales::hue_pal()(5)
  names(Concentration) <- as.character(sort(unique(results_female$Concentration)))
  
  male_validation <- validation_plot(results_male, "Male Zebrafish", Concentration, Tissue_markers)
  female_validation <- validation_plot(results_female, "Female Zebrafish", Concentration, Tissue_markers)
  
  grid_plot_validation <- ggpubr::ggarrange(female_validation, male_validation, 
                                            ncol=2, nrow=1, 
                                            common.legend = TRUE, legend="bottom")
  
  return(list("Training_plot" = grid_plot, 
              "Validation_plot"= grid_plot_validation))
  
}



create.plot <- function(compartment, predictions, data_df, color_codes){
  plot <- ggplot()+
    geom_line(data = predictions, aes(x = Time, y = !!as.name(compartment), color='Predictions'), size=1.3)+
    geom_point(data=data_df, aes(x=Time, y= !!as.name(compartment), color='Observations'), size=4)+
    labs(title = compartment, 
         y = "PFOA (ug/ml)",
         x = "Time (d)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("Data", values=color_codes)+
    theme_light()+
    theme(legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          axis.text = element_text(size = 14))
}


validation_plot <- function(xy_data,set_title, Concentration, Tissue_markers){
  ggplot()+
    geom_abline(intercept = log10(1), slope = 1, linetype = "dashed", color = "black", linewidth = 1.5,alpha = 0.7)+  # Identity line in log10 scale
    geom_point(data = xy_data, aes(x=Observed, y=Predicted, color = Concentration, shape = Tissue), size=4, stroke = 1.5)+
    
    
    scale_y_log10(limits=c(1e-04,1e01))+
    scale_x_log10(limits=c(1e-04,1e01))+
    
    scale_color_manual(values = Concentration, name = "Concentrations (ug/L):")+
    
    theme(legend.spacing.y = unit(1, 'cm')) +
    guides(fill = guide_legend(byrow = TRUE))+
    
    scale_shape_manual(values = Tissue_markers, name= "Tissues:")+
    theme_light()+
    labs(title = set_title ,
         y = expression("Predicted PFOA (" * mu* "g/g tissue)"),
         x = expression("Observed PFOA (" * mu* "g/g tissue)"))+
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(1.0, 'cm'),  
          legend.title = element_text(size=14),
          legend.text = element_text(size=12,  hjust = 0),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
    )
}




#=============================#
#.      END OF FUNCTIONS      # 
#=============================#







setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/Zebrafish_PFAS_Experiment_Design')

#============#
# User Input #
#============#

# Training input
L = 3 # cm
BW <- 0.0096 * L^3.55 # g Grech et al. (2019)
T_water = 23 + 273 #K 
sex = "Male"
pH_w = 8.0
C_water = 25 # ug/ml
params <- create.params(BW, T_water, sex, pH_w, C_water)
inits <- create.inits(params)
sample_time <- seq(0,30,0.5)

# solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
#                                     method="lsodes",rtol = 1e-05, atol = 1e-05))





#===========#
# Load Data #
#===========#
brain_data <- read.csv("Data/Bian_2022_data/csv_files/brain.csv")
#gallbladder_data <- read.csv("Data/Bian_2022_data/csv_files/gallbladder.csv")
gill_data <- read.csv("Data/Bian_2022_data/csv_files/gill.csv")
heart_data <- read.csv("Data/Bian_2022_data/csv_files/heart.csv")
intestine_data <- read.csv("Data/Bian_2022_data/csv_files/intestines.csv")
kidney_data <- read.csv("Data/Bian_2022_data/csv_files/kidney.csv")
liver_data <- read.csv("Data/Bian_2022_data/csv_files/liver.csv")
muscle_data <- read.csv("Data/Bian_2022_data/csv_files/muscle.csv")
spinal_cord_data <- read.csv("Data/Bian_2022_data/csv_files/spinal_cord.csv")
#swim_bladder_data <- read.csv("Data/Bian_2022_data/csv_files/swim_bladder.csv")

# Estimate the total concentration of muscle and spinal_cord as
# they are modeled through the poorly perfused tissues
PPT_data <- data.frame(cbind(muscle_data[,1], (muscle_data[,2]*params$V_muscle + spinal_cord_data[,2]*params$V_skeleton)/params$V_pp ))

Time <- c(1,2,4,8,12,16,20,24,30)
data_df <- data.frame(cbind(Time, brain_data[,2], gill_data[,2], heart_data[,2],
                            intestine_data[,2], kidney_data[,2], liver_data[,2], PPT_data[,2])) 
colnames(data_df) <- c("Time", "Brain", "Gill", "RPT", "Git",
                       "Kidney", "Liver", "PPT")

# Transform the intensity units to concentration based on figure 2(d) of Bian et al. (2022)
intensity_to_concentrtion <-data.frame(concentration = c(3.991297134, 26.03656992,
                                                         50.23768957, 100.6318519,
                                                         499.4784241),
                                       intensity = c(1402.641217, 12473.87966, 16992.55354,
                                                     43119.9451, 153203.5173))

# Concentration is estimated in ug PFOA/ml tissue
for (i in 2:dim(data_df)[2]) {
  data_df[,i] <- Hmisc::approxExtrap(intensity_to_concentrtion$intensity, intensity_to_concentrtion$concentration, data_df[,i])$y
}
data_df <- data_df[c("Time", "Brain", "Gill", "Git", "Kidney", "Liver", "RPT", "PPT")]


### Load and process the validation data from Ulaq et al. (2015) (Used for validation)
brain_female <- read.csv("Data/Ulhaq_2015_data/csv_files/brain_female.csv")
intestine_female <- read.csv("Data/Ulhaq_2015_data/csv_files/intestine_female.csv")
liver_female <- read.csv("Data/Ulhaq_2015_data/csv_files/liver_female.csv")
ovaries_female <- read.csv("Data/Ulhaq_2015_data/csv_files/ovaries.csv")

brain_male <- read.csv("Data/Ulhaq_2015_data/csv_files/brain_male.csv")
intestine_male <- read.csv("Data/Ulhaq_2015_data/csv_files/intestine_male.csv")
liver_male <- read.csv("Data/Ulhaq_2015_data/csv_files/liver_male.csv")

# The observed concentrations are reported in ng/g
# So they are transformed to ug/g by multiplying with 1/1000

male_data <- data.frame(cbind(brain_male[,1], brain_male[,2]/1000, 
                              intestine_male[,2]/1000, liver_male[,2]/1000))
colnames(male_data) <- c("Concentration", "Brain", "Git", "Liver")

female_data <- data.frame(cbind(brain_female[,1], brain_female[,2]/1000, 
                                intestine_female[,2]/1000, liver_female[,2]/1000, ovaries_female[,2]/1000))
colnames(female_data) <- c("Concentration", "Brain", "Git", "Liver", "Gonads")


validation_data <- list("male_data" = male_data,
                        "female_data" = female_data)

#====================#
# Optimization Tests #
#====================#

N_iter <- 5000
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",  #"NLOPT_LN_SBPLX" ,
              "xtol_rel" = 1e-07,
              "ftol_rel" = 1e-07,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = N_iter,
              "print_level" = 1 )


Test_1 <- list("compartments" = c("Liver"),
               "compartments_last_point" = NULL,#c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
               "rest" = T,
             "grouped_parameters" = list(c("P_li"),
                                         c("P_rp", "P_ki", "P_br", "P_git", "P_pp")),
               "fitted_params_names" = c("cf_1", "cf_2", "cl_feces",  "CalbB"),
               "x0" = c(1, 1, 1, 0.0001),
               "lb"	= c(1e-10, 1e-10, 1e-10, 1e-10),
               "ub" = c(1e05, 1e05, 1e05, 1e05)
)

Test_2 <- list("compartments" = c("Liver"),
               "compartments_last_point" = NULL,#c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
               "rest" = T,
               "grouped_parameters" = list(c("P_li"),
                                           c("P_rp", "P_ki", "P_br", "P_git", "P_pp")),
               "fitted_params_names" = c("cf_1", "cf_2", "cl_feces"),
               "x0" = c(1, 1, 1),
               "lb"	= c(1e-10, 1e-10, 1e-10),
               "ub" = c(1e05, 1e05, 1e05)
)

Test_3 <- list("compartments" = c("Liver"),
               "compartments_last_point" = NULL,#c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
               "rest" = T,
               "grouped_parameters" = list(c("P_li"),
                                           c("P_rp", "P_ki", "P_br", "P_git", "P_pp")),
               "fitted_params_names" = c("cf_1", "cf_2", "cl_feces"),
               "x0" = c(1, 1, 1),
               "lb"	= c(1e-10, 1e-10, 1e-10),
               "ub" = c(1e05, 1e05, 1e05)
)


# Test_1 <- list("compartments" = c("Liver"),
#                "compartments_last_point" = c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
#                "rest" = F,
#                "grouped_parameters" = list(c("P_li", "P_br", "P_git", "P_ki", "P_rp", "P_pp")),
#                "fitted_params_names" = c("cf_1", "cl_feces",  "CalbB"),
#                "x0" = c(1, 1, 0.001),
#                "lb"	= c( 1e-10, 1e-10, 1e-05),
#                "ub" = c(1e03, 1e03, 10)
# )
# 
# Test_2 <- list("compartments" = c("Liver"),
#                "compartments_last_point" = c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
#                "rest" = T,
#                "grouped_parameters" = list(c("P_li", "P_br", "P_git", "P_ki", "P_rp", "P_pp")),
#                "fitted_params_names" = c("cf_1", "cl_feces",  "CalbB"),
#                "x0" = c(1, 1, 0.001),
#                "lb"	= c(1e-10, 1e-10, 1e-05),
#                "ub" = c(1e03, 1e03, 10)
# )
# 
# 
# Test_3 <- list("compartments" = c("Liver"),
#                "compartments_last_point" = c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
#                "rest" = F,
#                "grouped_parameters" = list(c("P_li", "P_rp"),
#                                            c("P_ki", "P_br", "P_git", "P_pp")),
#                "fitted_params_names" = c("cf_1", "cf_2", "cl_feces",  "CalbB"),
#                "x0" = c(1, 1, 1, 0.001),
#                "lb"	= c(1e-10, 1e-10, 1e-10, 1e-05),
#                "ub" = c(1e03, 1e03, 1e03, 10)
# )
# 
# Test_4 <- list("compartments" = c("Liver"),
#                "compartments_last_point" = c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
#                "rest" = T,
#                "grouped_parameters" = list(c("P_li", "P_rp"),
#                                            c("P_ki", "P_br", "P_git", "P_pp")),
#                "fitted_params_names" = c("cf_1", "cf_2", "cl_feces",  "CalbB"),
#                "x0" = c(1, 1, 1, 0.001),
#                "lb"	= c(1e-10, 1e-10, 1e-10, 1e-05),
#                "ub" = c(1e03, 1e03, 1e03, 10)
# )
# 
# Test_5 <- list("compartments" = c("Liver"),
#                "compartments_last_point" = c("Brain", "Gill",  "Git", "Kidney", "RPT", "PPT"),
#                "rest" = F,
#                "grouped_parameters" = list(c("P_li", "P_rp"),
#                                            c( "P_git", "P_pp"),
#                                            c("P_br", "P_ki")),
#                "fitted_params_names" = c("cf_1", "cf_2", "cf_3", "cl_feces",  "CalbB"),
#                "x0" = c(1, 1, 1, 1, 0.001),
#                "lb"	= c(1e-10, 1e-10, 1e-10, 1e-10, 1e-05),
#                "ub" = c(1e03, 1e03, 1e03, 1e03, 10)
# )



# Tests_list <-list(Test_1, Test_2, Test_3, Test_4, Test_5)

Tests_list <- list(Test_1, Test_2)

for (i in 1:length(Tests_list)) {
  Test <- Tests_list[[i]]
  params <- create.params(BW, T_water, sex, pH_w, C_water)
  optimization <- nloptr::nloptr( x0 = Test$x0,
                                  eval_f = obj.func,
                                  lb	= Test$lb,
                                  ub = Test$ub,
                                  opts = opts,
                                  compartments = Test$compartments,
                                  compartments_last_point=Test$compartments_last_point,
                                  params = params,
                                  rest = Test$rest,
                                  grouped_parameters = Test$grouped_parameters,
                                  fitted_params_names = Test$fitted_params_names)
  
  Test[["optimization"]] <- optimization
  Test[["Plots"]] <- prepare_plot(Test)
  Tests_list[[i]] <- Test
  
  ggsave(paste0("Results/Correction_factor_results/Train/", "Test_", i, "_Train.png" ), Test$Plots$Training_plot, width = 12, height = 8, units = "in")
  ggsave(paste0("Results/Correction_factor_results/Validation/", "Test_", i, "_Validation.png" ), Test$Plots$Validation_plot, width = 12, height = 8, units = "in")
}


#save.image("Results/Tests.RData")
#ggsave("Results/Test_1.png", grid_plot, width = 12, height = 8, units = "in")
