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
  Dow_gill <- fn_gill * Kow *(1-fn_gill)*Kow_ion
  
  # Octanol-Water distribution coefficient at fish body pH 7.4 - Wang et al.(2022)
  Dow <- fn_fish * Kow *(1-fn_fish)*Kow_ion
  
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
  
  # Plasma:water distribution ratio - Wang et al. (2022)
  D_plasma_w <- 10^(0.75*Kow_ion + 0.58)
  
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
    "UF"=UF, "a_fpp"=a_fpp, "a_fs"=a_fs,
    
    "f_art"=f_art, "f_ven"=f_ven,
    
    "C_water"=C_water,
    "ku"=ku, "ke_feces"=ke_feces
  ))  
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    M_gill=0;
    M_ven=100; M_art=0; M_lumen=0; M_git=0; M_li=0; M_ki=0; M_br=0; M_go=0;
    M_rp=0; M_pp=0; M_ad=0; M_skin=0; M_excreted = 0; C_water=C_water
    
    return(c("M_gill"=M_gill, "M_ven"=M_ven, "M_art"=M_art, "M_lumen"=M_lumen, 
             "M_git"=M_git, "M_li"=M_li, "M_ki"=M_ki, "M_br"=M_br, "M_go"=M_go, 
             "M_rp"=M_rp, "M_pp"=M_pp, "M_ad"=M_ad, "M_skin"=M_skin, 
            "M_excreted"=M_excreted, "C_water"=C_water))
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
    
    #dM_i/dt: ug/d
    
    dM_gill <- kx*(C_water - C_gill*UF/Pbw) + F_card*(C_ven-C_gill)
    dM_exreted <- kx*C_gill*UF/Pbw
    
    # Water concentration: ug/ml
    dC_water = 0
    
    # Venous blood
    dM_ven <- (F_li + F_rp + F_git + F_go)*M_li/(V_li*P_li) + F_br*M_br/(V_br*P_br) +
      F_ad*M_ad/(V_ad*P_ad) + (F_ki + (1-a_fpp)*F_pp + (1-a_fs)*F_skin)*M_ki/(V_ki*P_ki) +
      a_fpp*F_pp*M_pp/(V_pp*P_pp) + a_fs*F_skin*M_skin/(V_skin*P_skin) -
      F_card*C_ven #+ kx*(C_water - C_ven*UF/Pbw)
    
    # Arterial blood
    dM_art <- F_card*(C_gill - C_art) #F_card*(C_ven - C_art)
    
    # GIT lumen
    dM_lumen <- -(ku + ke_feces) * M_lumen
    
    # GIT
    dM_git <- ku*M_lumen + F_git*(C_art - M_git/(V_git*P_git))
    
    # Liver
    dM_li <- F_li*(C_art - M_li/(V_li*P_li)) + F_rp*(M_rp/(V_rp*P_rp) - M_li/(V_li*P_li)) +
      F_git*(M_git/(V_git*P_git) - M_li/(V_li*P_li)) + F_go*(M_go/(V_go*P_go) - M_li/(V_li*P_li))
    
    # Kidney
    dM_ki <- F_ki*(C_art - M_ki/(V_ki*P_ki)) + (1-a_fpp)*F_pp*(M_pp/(V_pp*P_pp) - M_ki/(V_ki*P_ki)) +
      (1 - a_fs)*F_skin*(M_skin/(V_skin*P_skin) - M_ki/(V_ki*P_ki))
    
    # Brain
    dM_br <- F_br*(C_art - M_br/(V_br*P_br)) # ug/day
    
    # Gonads
    dM_go <- F_go*(C_art - M_go/(V_go*P_go)) # ug/day
    
    # Richly perfused tissues
    dM_rp <- F_rp*(C_art - M_rp/(V_rp*P_rp)) # ug/day
    
    # Poorly perfused tissues
    dM_pp <- F_pp*(C_art - M_pp/(V_pp*P_pp)) # ug/day
    
    # Adipose
    dM_ad <- F_ad*(C_art - M_ad/(V_ad*P_ad)) # ug/day
    
    # Skin
    dM_skin <- F_skin*(C_art - M_skin/(V_skin*P_skin)) # ug/day
    
    return(list(c("dM_gill"=dM_gill, "dM_ven"=dM_ven, "dM_art"=dM_art, "dM_lumen"=dM_lumen, "dM_git"=dM_git,
                  "dM_li"=dM_li, "dM_ki"=dM_ki, "dM_br"=dM_br, "dM_go"=dM_go, "dM_rp"=dM_rp,
                  "dM_pp"=dM_pp, "dM_ad"=dM_ad, "dM_skin"=dM_skin, "dM_exreted"=dM_exreted,
                  "dC_water"=dC_water),
                "C_gill"=C_gill,
                "C_art"=C_art, "C_ven"=C_ven, "C_git"=C_git, "C_li"=C_li,
                "C_ki"=C_ki, "C_br"=C_br, "C_go"=C_go, "C_rp"=C_rp, 
                "C_pp"=C_pp, "C_ad"=C_ad, "C_skin"=C_skin
    ))
    
  })
}

#============#
# User Input #
#============#

setwd('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/Zebrafish_PFAS_Experiment_Design')
BW = 1 #g
T_water = 20 + 273 #K 
sex = "Male"
pH_w = 7.2
C_water = 0
params <- create.params(BW, T_water, sex, pH_w, C_water)
inits <- create.inits(params)
sample_time <- seq(0,10,0.5)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
# check mass balance 
rowSums(solution[,2:15])
