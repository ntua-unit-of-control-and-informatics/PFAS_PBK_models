### PBK model for zebrafish female exposure to single PFAS ###
# Paper of the model: Golosovskai et al. 2024
# doi: https://doi.org/10.1016/j.scitotenv.2023.168738
# V. Minadakis

##############################
### Variable -(Notation)- Units
# Amounts:        -(A_x) - microg 
# Volumes:      -(V_x) - mL
# Time:            -(t)   - d
# Flows:           -(F_x) - mL/d
# Concentrations: -(C_x) - microg/mL
# Masses:         -(BW)  - g
# Temperature: -(TC_c)  - Celsius
# Temperature: -(TC_k)  - Kelvin
# Density of each tissue is considered equal to 1
###############################

#===============================================
### create.params function

create.params <- function(user_input){
  with(as.list(user_input),{# Physiological parameters
    Bw_Fcard_ref= 0.5      # Body weight of reference for F_card from Pery 2014 and Hamilton 2014
    F_card_ref   = (((11.1+12)/2)*10^(-3))*(24*60)/0.5 #ml/d/g converted from Pery et al (11.1 ul/min) /BW_Fcard_ref
    BW_VO2_ref = 0.4       # # Reference BW for O2 extraction rate value Vergauwen et al
    # Environmental condition
    TC_c         = 26
    TA           = 3000     # Arrhenius temperature  in Kelvin for zebrafish
    TR_Fcard     = 26 + 273.15  # Reference T for cardiac output (Kelvin) Average of Pery et al 2014 (26) and Hamilton et al 2014 (28)
    TR_VO2       = 27 + 273.15  # Reference T for O2 consumption rate (Kelvin) Vergauwen et al
    V_water      = V_water    # Volume of aquarium (mL) large to account for flow-through system
    
    # Effective respiratory volume & cardiac output
    VO2_ref      = 9.84 # mg O2/d/g --> calculated from ventilation volume proposed by P?ry 2014 (Qw= 0.55 mL/min/g)
    OEE          = 0.71 # Oxygen extraction efficiency of 71% proposed by Erickson 1990
    Sat          = 0.90 # dissolved oxygen saturation of 90% proposed by Erickson 1990
    art_ven_frac = 1/3 # fraction of arterial blood
    
    sc_blood  = 0.0222 # volume scaling factor : fraction of BW (%)
    sc_gon    = 0.08
    sc_brain  = 0.01
    sc_liv    = 0.0205
    sc_fat    = 0.022  # fathead minnow data
    sc_skin   = 0.10   # rainbow trout data
    sc_git    = 0.099  # rainbow trout data
    sc_kidney = 0.0021
    sc_rp     = 0.0266
    sc_pp     = 0.6176   # (1 -  sc_blood - sc_gon - sc_brain - sc_liv - sc_fat -  sc_skin - sc_git - sc_kidney -sc_rp)
    # 
    gon_frac    = 0.01252    # Fraction of arterial blood flow
    brain_frac  = 0.0324   
    liv_frac    = 0.01942  
    fat_frac    = 0.0135   
    skin_frac   = 0.05639  
    git_frac    = 0.174    
    kidney_frac = 0.02233  
    rp_frac     = 0.13255  
    pp_frac     = 0.53689   #(1 - gon_frac - brain_frac - liv_frac - fat_frac - skin_frac - git_frac -  kidney_frac - rp_frac)
    
    a_Fpp = 0.4 # Fraction of PPT blood going to venous
    a_Fs  = 0.1 # Fraction of skin blood going to venous
    
    Bt_M = 0.0008 # amount of proteins in zf blood (M)
    
    # Egg #
    Egggrowth	=	0.027	
    spawnrate = 1.5 # d
    V_one_egg = 0.000212 #V_egg
    
    #Oral absorption
    Ku_1       = 0   # Diffusion coefficient
    frac_abs_1 = 0   # Absorption fraction.  /!\ Parametrisation with frac_abs implies Ke_feces null !
    
    Cl_bile_1  = 0 # ml/d
    Ke_urine_1 = 0 # 1/d
    K_BG_1     = 0 # 1/d
    Ke_feces_1 = 0.83  # 1/d estiamted from Nichols  et al. 2004 0.83
    
    urine_rate     = 0.057947686 # V_burst = 1.2 mL.kg-1 every 29.82 minutes proposed by Curtis 1991 --> 1.2e-03 mL.g BW-1?
    urination_interval = 30/60/24      # 30min proposed by Curtis 1991 --> 30/60/24hr
    
    # Substance specific params
    # Chemical parameters
    substance = substance
    if(substance == "PFOA"){
      MW_1 = 414.07 # g/mol PFOA
      pKa =  0.5
      i=1
      logKow = 4.35
      delta_logKow = 3.1
      # Dmw = 10^(3.51-2.0)
      # Dpw = 10^(3.43-2.0)
      Pegggon_1 = 0.127 # our data
      Pbb_1 = 0.2
      Plivb_1 = 0.6
      Pgonb_1 = 0.4
      Pgitb_1 = 0.64
      Pppb_1 = 0.08
      Prpb_1 = 0.26
      Pfatb_1 = 0.5
      Pskinb_1 = 0.5
      Pkidb_1 = 0.64
      Kd = 0.83 * 1e-03
      
      K_1 = 1/Kd # 1/Kd with Kd in M  PFOA maybe multiply by the number of active sites
      
    }else if(substance == "PFHpA"){
      MW_1 = 364.06 # g/mol 
      pKa =  0.47
      i=1
      logKow = 3.94
      delta_logKow = 3.1
      # Dmw = 10^(3.51-2.0)
      # Dpw = 10^(3.43-2.0)
      Pegggon_1 = 0.17 # our data
      Pbb_1 = 0.08
      Plivb_1 = 1.4
      Pgonb_1 = 0.13
      Pgitb_1 = 0.93
      Pppb_1 = 0.08
      Prpb_1 = 0.23
      Pfatb_1 = 0.5
      Pskinb_1 = 0.5
      Pkidb_1 = 0.93
      Kd = 0.44 * 1e-03
      
      K_1 = 1/Kd # 1/Kd with Kd in M  PFHpA maybe multiply by the number of active sites
      
    }else if(substance == "PFNA"){
      MW_1 = 464.08 # g/mol 
      pKa =  0.52
      i=1
      logKow = 4.97
      delta_logKow = 3.1
      # Dmw = 10^(3.51-2.0)
      # Dpw = 10^(3.43-2.0)
      Pegggon_1 = 0.13 # our data
      Pbb_1 = 0.07
      Plivb_1 = 0.26
      Pgonb_1 = 0.13
      Pgitb_1 = 0.93
      Pppb_1 = 0.05
      Prpb_1 = 0.25
      Pfatb_1 = 0.3
      Pskinb_1 = 0.5
      Pkidb_1 = 0.31
      Kd = 0.58 * 1e-03
      
      K_1 = 1/Kd # 1/Kd with Kd in M  PFHpA maybe multiply by the number of active sites
      
    }else if (substance == "PFHxS"){
      MW_1 = 400.12 # g/mol 
      pKa =   -3.34
      i=1
      logKow = 3.41
      delta_logKow = 3.1
      # Dmw = 10^(3.51-2.0)
      # Dpw = 10^(3.43-2.0)
      Pegggon_1 = 0.19 # our data
      Pbb_1 = 0.06
      Plivb_1 = 0.46
      Pgonb_1 = 0.62
      Pgitb_1 = 0.64
      Pppb_1 = 0.06
      Prpb_1 = 0.26
      Pfatb_1 = 0.5
      Pskinb_1 = 0.5
      Pkidb_1 = 0.64
      Kd = 0.71 * 1e-03
      
      K_1 = 1/Kd # 1/Kd with Kd in M  PFHpA maybe multiply by the number of active sites
      
    }else if(substance == "PFBS"){
      MW_1 = 300.1 # g/mol 
      pKa =   -3.57
      i=1
      logKow = 3.9
      delta_logKow = 3.1
      # Dmw = 10^(3.51-2.0)
      # Dpw = 10^(3.43-2.0)
      Pegggon_1 = 0.342 # our data
      Pbb_1 = 0.09
      Plivb_1 = 0.9
      Pgonb_1 = 0.9
      Pgitb_1 = 0.7
      Pppb_1 = 0.09
      Prpb_1 = 0.27
      Pfatb_1 = 0.4
      Pskinb_1 = 0.5
      Pkidb_1 = 0.7
      Kd = 0.71 * 1e-03
      
      K_1 = 1/Kd # 1/Kd with Kd in M  PFHpA maybe multiply by the number of active sites
    }
    # Exposure quantity (microg)
    WaterExposure_1  = WaterExposure_1 # microg/mL 
    
    return(list(  "Bw_Fcard_ref"=Bw_Fcard_ref, "F_card_ref" =F_card_ref,
                  "BW_VO2_ref" = BW_VO2_ref, "TC_c" = TC_c,
                  "TA" = TA, "TR_Fcard" = TR_Fcard, "TR_VO2" = TR_VO2,
                  "V_water" = V_water, "VO2_ref" = VO2_ref,
                  "OEE" = OEE, "Sat" = Sat, "art_ven_frac" = art_ven_frac,
                  "sc_blood" = sc_blood, "sc_gon" = sc_gon,
                  "sc_brain" = sc_brain, "sc_liv" = sc_liv,
                  "sc_fat" = sc_fat, "sc_skin" = sc_skin,
                  "sc_git" = sc_git, "sc_kidney" = sc_kidney,
                  "sc_rp" = sc_rp, "sc_pp" = sc_pp,
                  "gon_frac" = gon_frac, "brain_frac" = brain_frac,
                  "liv_frac" = liv_frac, "fat_frac" = fat_frac,
                  "skin_frac" = skin_frac, "git_frac" = git_frac,
                  "kidney_frac" = kidney_frac, "rp_frac" = rp_frac,
                  "pp_frac" = pp_frac, 
                  "a_Fpp" = a_Fpp, "a_Fs"  = a_Fs,
                  "Bt_M" = Bt_M,
                  "Egggrowth"	=	Egggrowth,
                  "spawnrate" = spawnrate,
                  "V_one_egg" = V_one_egg,
                  "Ku_1" = Ku_1, "frac_abs_1" = frac_abs_1,
                  "Cl_bile_1" = Cl_bile_1, "Ke_urine_1" = Ke_urine_1,
                  "K_BG_1" = K_BG_1, "Ke_feces_1" = Ke_feces_1,
                  "urine_rate" = urine_rate, "urination_interval" = urination_interval,
                  #"MW_1" = MW_1,
                  "pKa" =  pKa, "i"=i,  "logKow" = logKow,
                  "delta_logKow" = delta_logKow, "Pegggon_1" = Pegggon_1,
                  "Pbb_1" = Pbb_1, "Plivb_1" = Plivb_1, "Pgonb_1" = Pgonb_1,
                  "Pgitb_1" = Pgitb_1, "Pppb_1" = Pppb_1, "Prpb_1" = Prpb_1,
                  "Pfatb_1" = Pfatb_1, "Pskinb_1" = Pskinb_1, "Pkidb_1" = Pkidb_1,
                  "Kd" = Kd, "K_1" = K_1, "WaterExposure_1" = WaterExposure_1,
                  "BW" = BW, "start_point"= start_point, "end_point" = end_point,
                  "time_of_exposure"=time_of_exposure, "WaterExposure_1"=WaterExposure_1
                  
    ))})
  
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    V_egg = 0.000212; 
    V_urine = 0; A_admin_gill_1 = 0; A_lumen_1 = 0; A_excr_gill_1 = 0; 
    A_urine_1 = 0; A_urine_cum_1 = 0; A_feces_1 = 0; A_excr_water_1 = 0; 
    A_bile_1 = 0; A_blood_1 = 0; A_liv_1 = 0; A_egg_1 = 0; A_egg_cum_1 = 0; 
    A_gon_1 = 0; A_fat_1 = 0; A_git_1 = 0; A_brain_1 = 0; A_kidney_1 = 0; 
    A_skin_1 = 0; A_rp_1 = 0; A_pp_1 = 0; 
    A_water_1 = 0#WaterExposure_1 * V_water
    
    return(c("V_egg" = V_egg, 
             "V_urine" = V_urine, "A_admin_gill_1" = A_admin_gill_1, 
             "A_lumen_1" = A_lumen_1, "A_excr_gill_1" = A_excr_gill_1, 
             "A_urine_1" = A_urine_1, "A_urine_cum_1" = A_urine_cum_1, 
             "A_feces_1" = A_feces_1, "A_excr_water_1" = A_excr_water_1, 
             "A_bile_1" = A_bile_1, "A_blood_1" = A_blood_1, "A_liv_1" = A_liv_1,
             "A_egg_1" = A_egg_1, "A_egg_cum_1" = A_egg_cum_1, 
             "A_gon_1" = A_gon_1, "A_fat_1" = A_fat_1, "A_git_1" = A_git_1,
             "A_brain_1" = A_brain_1, "A_kidney_1" = A_kidney_1, 
             "A_skin_1" = A_skin_1, "A_rp_1" = A_rp_1, "A_pp_1" = A_pp_1, 
             "A_water_1" = A_water_1))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    events_end_exposure_w <- list(data = rbind(data.frame(var = c("A_water_1"),
                                                          time  = time_of_exposure,
                                                          value = WaterExposure_1*V_water,
                                                          method = c("add"))
    ))
    
    #return(events_end_exposure_w)
    
    events_spawn <-          list(data = rbind(data.frame(var = c("A_egg_1"),
                                                          time  = seq(spawnrate,end_point,spawnrate),
                                                          value = 0,
                                                          method = c("replace")),
                                               data.frame(var = c("V_egg"),
                                                          time  = seq(spawnrate,end_point,spawnrate),
                                                          value = 2.12E-04,
                                                          method = c("replace"))
    ))
  
    return(list(data = rbind(events_end_exposure_w[[1]], events_spawn[[1]])) )
    })
}

custom.func <- function(){
  return()
}

ode.func <- function(t,initial_v, parameters, custom.func) {
  
  with(as.list(c(initial_v, parameters)),{
    
    PS       = pi^(1/3)*(6*V_one_egg)^(2/3)*(V_egg/V_one_egg)
    
    #Volumes (ml)
    V_blood = sc_blood * BW
    V_liv    = sc_liv    * BW
    V_gon    = sc_gon    * BW
    V_fat    = sc_fat    * BW
    V_git    = sc_git    * BW
    V_brain  = sc_brain  * BW
    V_kidney = sc_kidney * BW
    V_skin   = sc_skin   * BW
    V_rp     = sc_rp     * BW
    V_pp     = BW-(V_blood+V_liv+V_gon+V_fat+V_git+V_brain+V_kidney+V_skin+V_rp)
    
    V_total = V_blood+V_liv+V_gon+V_fat+V_git+V_brain+V_kidney+V_skin+V_rp+V_pp
    V_bal = BW - V_total
    
    KT <- function(T, TR, TA){ exp ( (TA / TR) - (TA / T) ) }
    
    #Parameter equations
    TC_k = TC_c + 273.15 # (degree K)
    Fcard = F_card_ref * KT(T=TC_k , TR=TR_Fcard,TA=TA) * ((BW/Bw_Fcard_ref)^(-0.1))* BW #cardiac output adjusted for temperature and BW of the simulation (ml/d)
    ####################################################################################################
    
    ##Flows to organs (ml/d)
    Fliv    = liv_frac    * Fcard
    Fgon    = gon_frac    * Fcard
    Fgit    = git_frac    * Fcard
    Ffat	  = fat_frac    * Fcard
    Fbrain  = brain_frac  * Fcard
    Fkidney = kidney_frac * Fcard
    Fskin   = skin_frac   * Fcard
    Frp     = rp_frac     * Fcard
    Fpp     = (Fcard-(Fliv+Fgon+Fgit+Ffat+Fbrain+Fkidney+Fskin+Frp))
    
    Ftotal  = (Fliv+Fgon+Fgit+Ffat+Fbrain+Fkidney+Fskin+Frp+Fpp)
    Fbal    = Fcard- Ftotal
    Fegg    = Fgon
    
    VO2     = VO2_ref*KT(T=TC_k,TR=TR_VO2,TA=TA)*((BW/BW_VO2_ref)^(-0.1))*BW   # O2 consumption rate (mg/d) adjusted for T and BW of the simmulation
    Co2w    = ((-0.24 * TC_c + 14.04) * Sat)/(10^3)                            # (mg O2/mL)C of O2 in water at T in celsius
    Fwater  = VO2/(OEE*Co2w)                                                   #Effective respratory volume. Equation from Barber et al.
    
    #### Plasma protein binding parameters ####
    
    # Kd - dissociation constant, mM;
    # K - association constant, taken as 1/Kd in parameters vector
    # Bt_M - total concentration of protein, M
    
    Free_b_1 = 1/(1+(K_1*Bt_M)) # equation from https://pubs.acs.org/doi/full/10.1021/acs.chemrestox.1c00193 for albumin (Bt - total protein concentration, equals to maximum binding capacity). Linear. For one compound
    
    logKow_ion = logKow-delta_logKow
    
    fn_fish = 1/(1+10^(i*(7.4-pKa)))  # (i = 1 for acid and -1 for base)
    
    dow = fn_fish * 10^(logKow) + (1-fn_fish)*10^(logKow_ion)
    
    PBW_dep = 0.008 * 0.3 * dow + 0.007*2.0 * dow^(0.94) + 0.134 * 2.9 * dow^(0.63) + 0.851 # PBW for ionic substances Wang 2022
    
    kow = 10^(logKow)
    
    y_water = VO2/(OEE*Co2w)* ( 1/(1000^(0.25) * BW^(0.75)) )  # gill ventilation coefficient OEE=0.71 ml/d
    y_blood = Fcard * PBW_dep*( 1/(1000^(0.25) * BW^(0.75)) )   # blood perfusion coefficient
    
    kx = ((BW/1000)^(0.75)) /( 2.8*10^(-3) + 68/dow + 1/y_water + 1/y_blood )*1000
    
    kout = (1/(68*(dow-1)+1) * kx)
    
    ############# FRACTION UNBOUND APPROACH 1. INDEPENDENT ASSESSMENT OF COMPOUNDS #####################
    
    Pskinb_1 = (0.037*0.3*dow + 0.012*2.0*(dow^(0.94)) + 0.259*2.9*(dow^(0.63)) + 0.692)/PBW_dep
    Pkidb_1  = (0.059*0.3*dow + 0.023*2.0*(dow^(0.94)) + 0.215*2.9*(dow^(0.63)) + 0.703)/PBW_dep
    
    ################### CONCENTRATIONS #############################
    
    ## COMPOUND 1 ##
    
    C_blood_1 = A_blood_1/V_blood
    C_plasma_1 = (C_blood_1)*0.86 # 0.86 ratio_blood_plasma
    C_liv_1 = A_liv_1/V_liv
    C_fat_1 = A_fat_1/V_fat
    C_egg_1 = A_egg_1/V_egg
    C_gon_1 = A_gon_1/V_gon
    C_gon_egg_1 = (A_gon_1+A_egg_1)/(V_gon+V_egg)
    C_git_1 = A_git_1/V_git
    C_brain_1 = A_brain_1/V_brain
    C_kidney_1 = A_kidney_1/ V_kidney 
    C_skin_1 = A_skin_1/V_skin
    C_rp_1  = A_rp_1/V_rp
    C_pp_1  = A_pp_1/V_pp
    C_whole_body_1 = ((A_blood_1 +  A_liv_1 + A_bile_1 + A_gon_1 + A_brain_1 + A_fat_1 +A_skin_1 + A_kidney_1 + A_git_1 + A_pp_1 + A_rp_1 )
                      / (V_blood + V_liv + V_gon + V_brain + V_fat + V_skin + V_git + V_kidney + V_pp + V_rp))
    C_carcass_1 = (A_blood_1 +  A_bile_1 + A_fat_1 +A_skin_1 + A_kidney_1 + A_git_1 + A_pp_1 + A_rp_1) /
      (V_blood + V_fat + V_skin + V_git + V_kidney + V_pp + V_rp)
    # C_carcass_1 = (A_bile_1 + A_fat_1 +A_skin_1 + A_kidney_1 + A_git_1 + A_pp_1 + A_rp_1) /
    # (V_blood + V_fat + V_skin + V_git + V_kidney + V_pp + V_rp)
    
    ############## END OF CONCENTRATIONS #############
    
    dV_egg = Egggrowth
    
    dV_urine = urine_rate * BW
    
    ####################### CHEMICAL IN THE BODY ######################
    
    ### ABSORPTION ###
    
    dA_admin_gill_1 = kx*(A_water_1/V_water)
    
    dA_lumen_1 = A_bile_1*K_BG_1 - Ku_1*A_lumen_1 - Ke_feces_1*A_lumen_1
    
    ### ELIMINATION/ECXRETION ###
    dA_excr_gill_1  = kout*(A_blood_1* Free_b_1/V_blood)
    
    dA_urine_1 = Ke_urine_1 * A_kidney_1
    dA_urine_cum_1 = Ke_urine_1 *A_kidney_1
    
    dA_feces_1 = Ke_feces_1 * A_lumen_1
    
    dA_excr_water_1 = dA_excr_gill_1 + dA_urine_cum_1 + dA_feces_1 
    
    dA_bile_1 = Cl_bile_1 * (A_liv_1/V_liv) - A_bile_1*K_BG_1
    
    ### DISTRIBUTION ###
    dA_blood_1 = dA_admin_gill_1 - dA_excr_gill_1 - Free_b_1*((A_blood_1/V_blood) * (Fliv + Ffat + Fskin + Fgon + Fgit + Fbrain + Fkidney + Frp + Fpp)) +
      Free_b_1*(Ffat*((A_fat_1/V_fat)/Pfatb_1) + a_Fs*Fskin*((A_skin_1/V_skin)/Pskinb_1)+ Fbrain*((A_brain_1/V_brain)/Pbb_1)+
                  (Fliv+Fgon+Fgit+Frp)*((A_liv_1/V_liv)/Plivb_1) + (Fkidney+(1-a_Fs)*Fskin+(1-a_Fpp)*Fpp)*((A_kidney_1/V_kidney)/Pkidb_1) + a_Fpp*Fpp*((A_pp_1/V_pp)/Pppb_1))
    
    
    dA_liv_1    = Free_b_1*(Fliv*(A_blood_1/V_blood) + Fgon*((A_gon_1/V_gon)/Pgonb_1) + Frp*((A_rp_1/V_rp)/Prpb_1) + Fgit*((A_git_1/V_git)/Pgitb_1)) - 
      Free_b_1*(Fliv + Frp + Fgit + Fgon)*((A_liv_1/V_liv)/Plivb_1) - 
      Cl_bile_1*(A_liv_1/V_liv)
    
    dA_egg_1    = PS/V_egg* (((A_gon_1/V_gon)/Pgonb_1)-((A_egg_1/V_egg)/Pegggon_1))
    
    dA_egg_cum_1  = PS/V_egg* (((A_gon_1/V_gon)/Pgonb_1)-((A_egg_1/V_egg)/Pegggon_1))
    
    dA_gon_1      = Free_b_1*(Fgon*(A_blood_1/V_blood) - Fgon*((A_gon_1/V_gon)/Pgonb_1)) + PS/V_egg*(((A_egg_1/V_egg)/Pegggon_1) - ((A_gon_1/V_gon)/Pgonb_1))
    
    dA_fat_1    = Free_b_1*Ffat *((A_blood_1/V_blood) - ((A_fat_1/V_fat)/Pfatb_1))
    
    dA_git_1    = Ku_1*A_lumen_1 + Free_b_1*Fgit*((A_blood_1/V_blood) - (A_git_1/V_git)/Pgitb_1)
    
    dA_brain_1  = Free_b_1*Fbrain*(A_blood_1/V_blood - (A_brain_1/V_brain)/Pbb_1)
    
    dA_kidney_1 = Free_b_1*(Fkidney*(A_blood_1/V_blood) + (1-a_Fpp)*Fpp*((A_pp_1/V_pp)/Pppb_1) + (1-a_Fs)*Fskin*((A_skin_1/V_skin)/Pskinb_1)) - Free_b_1*(Fkidney + (1-a_Fpp)*Fpp + (1-a_Fs)*Fskin)*((A_kidney_1/V_kidney)/Pkidb_1) - dA_urine_1
    
    dA_skin_1   = Free_b_1*Fskin* (A_blood_1/V_blood - (A_skin_1/V_skin)/Pskinb_1)
    dA_rp_1     = Free_b_1*Frp* (A_blood_1/V_blood - (A_rp_1/V_rp)/Prpb_1)
    dA_pp_1     = Free_b_1*Fpp* (A_blood_1/V_blood - (A_pp_1/V_pp)/Pppb_1)
    
    # Chemical in aquarium water
    dA_water_1 = dA_excr_water_1 - dA_admin_gill_1
    
    ### MASS BALANCES ###
    
    A_body_1 = (A_blood_1 +  A_liv_1 + A_bile_1 + A_gon_1 + A_brain_1 + A_fat_1 +A_skin_1 + A_kidney_1 + A_git_1 + A_pp_1 + A_rp_1 )
    
    A_elim_1 = (A_excr_water_1 + A_egg_cum_1 + A_lumen_1)
    
    Mass_bal_1 = (A_admin_gill_1 - A_body_1 - A_elim_1)
    
    ########## END OF COMPOUND #############
    
    list(c(dV_egg, 
           dV_urine, 
           
           dA_admin_gill_1, 
           dA_lumen_1, 
           dA_excr_gill_1, 
           dA_urine_1, 
           dA_urine_cum_1, 
           dA_feces_1, 
           dA_excr_water_1, 
           dA_bile_1, 
           dA_blood_1,
           dA_liv_1, 
           dA_egg_1, 
           dA_egg_cum_1, 
           dA_gon_1, 
           dA_fat_1,
           dA_git_1, 
           dA_brain_1, 
           dA_kidney_1, 
           dA_skin_1, 
           dA_rp_1, 
           dA_pp_1, 
           dA_water_1
    ),
    "Fcard" = Fcard,
    "Fcard_PBW" = Fcard * PBW_dep,
    "Fwater" = Fwater,
    "kx" = kx,
    "dow" = dow,
    "logKow" =logKow,
    "Mass_balance_1" = Mass_bal_1,
    "Free_b_1" = Free_b_1,
    "C_whole_body_1" = C_whole_body_1,
    "C_blood_1" = C_blood_1,
    "C_plasma_1" = C_plasma_1,
    "C_liv_1" = C_liv_1,
    "C_fat_1" = C_fat_1,
    "C_egg_1" = C_egg_1,
    "C_gon_egg_1"= C_gon_egg_1,
    "C_gon_1" = C_gon_1,
    "C_git_1" = C_git_1,
    "C_brain_1" = C_brain_1,
    "C_kidney_1" = C_kidney_1,
    "C_skin_1" = C_skin_1,
    "C_rp_1" =  C_rp_1,
    "C_pp_1" = C_pp_1,
    "C_carcass_1" = C_carcass_1,
    "Plivb_1" = Plivb_1, 
    "Pgonb_1" =  Pgonb_1, 
    "Pegggon_1" = Pegggon_1,
    "Pbb_1" = Pbb_1, 
    "Pfatb_1" = Pfatb_1,
    "Pskinb_1" = Pskinb_1,
    "Pgitb_1" = Pgitb_1,
    "Pkidb_1" =  Pkidb_1,
    "Prpb_1" = Prpb_1,
    "Pppb_1"= Pppb_1
    )
  }
  )
}

### simulation parameters ###
user_input <- list(
  'BW'       = 0.8,
  'substance' = "PFNA",
  'WaterExposure_1' = 100, #8.91*10^(-3), #8.91*10^(-3), # microg/mL ,
  'time_of_exposure' = 0,
  'V_water' = 100, #1E12,
  start_point  = 0 ,      # days
  end_point = 10
  # period = NA,
  # time_final_dose = NA,
  # time_first_dose = 0,
  # time_end_exposure = 32
)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

sample_time <- seq(0,10,1)
### RUNNING THE SIMULATION ###
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
#print(tail(solution))


#====================
#7. Upload on Jaqpot 
#===================
# Subset of features to be displayed on the user interface
predicted.feats <- c("V_egg", "V_urine", "A_admin_gill_1", "A_lumen_1",
                     "A_excr_gill_1", "A_urine_1", "A_urine_cum_1",
                     "A_feces_1", "A_excr_water_1", "A_bile_1", "A_blood_1",
                     "A_liv_1", "A_egg_1", "A_egg_cum_1", "A_gon_1", "A_fat_1",
                     "A_git_1", "A_brain_1", "A_kidney_1", "A_skin_1", "A_rp_1",
                     "A_pp_1", "A_water_1",
#                      "Fcard",
#                      "Fcard_PBW",
#                      "Fwater",
#                      "kx",
#                      "dow",
#                      "logKow",
#                      "Mass_balance_1",
#                      "Free_b_1",
                     "C_whole_body_1", "C_blood_1", "C_plasma_1",
                     "C_liv_1", "C_fat_1", "C_egg_1", "C_gon_egg_1",
                     "C_gon_1", "C_git_1", "C_brain_1", "C_kidney_1",
                     "C_skin_1", "C_rp_1", "C_pp_1", "C_carcass_1")#,
                     # "Plivb_1" , 
                     # "Pgonb_1" , 
                     # "Pegggon_1" ,
                     # "Pbb_1" , 
                     # "Pfatb_1",
                     # "Pskinb_1",
                     # "Pgitb_1",
                     # "Pkidb_1",
                     # "Prpb_1" ,
                     # "Pppb_1")

# Log in Jaqpot server
jaqpotr::login.api()
# Deploy the model on the Jaqpot server to create a web service
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func,
                     method = "bdf",url = "https://api.jaqpot.org/jaqpot/")
