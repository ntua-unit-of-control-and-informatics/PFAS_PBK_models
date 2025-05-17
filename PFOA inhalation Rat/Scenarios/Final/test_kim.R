
create_variable_params <- function(BW,sex,  estimated_params, fixed_params){
  
  # BW in kg
  # Cheng and Ng 2017 Table S1
  # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
  #assuming the density of tissue is 1 g/mL.
  # Estimated parameters
  if(sex == "M"){
    RAFOatp_k <- estimated_params[1]
  }else{
    RAFOatp_k <-estimated_params[2]
  }
  RAFOat3 <- estimated_params[3]
  CL_int <- estimated_params[4] #uL/min/million hepatocytes
  RAFOatp_l <- estimated_params[5]
  RAFUrat <- RAFOatp_k
  RAFOat1 <- RAFOat3
  RAFOatp2_l <- RAFOatp_l
  RAFOatp_lu_ap <-estimated_params[6]
  RAFOatp_lu_bas <- RAFOatp_lu_ap
  RAFNtcp <- RAFOatp_l
  RAFOatp2_Int <- 6.799068e-07
  
  
  wall_width <-  0.5e-6*100 #[m] -->[cm], capillary wall width. From Ashoor et al. (2018) [10.1016/j.rpor.2018.09.007]. Similar value from Sosula (1974) [10.1016/0026-2862(74)90011-9]
  basement_membrane <-  0.5e-6*100 #assumption
  dif <- 5.46e-6 #[cm^2/s], from  Gauthier et al. (2024) [10.1021/acsestwater.4c00631]
  Pgap_continuous <- dif/ (wall_width+basement_membrane) #cm/s
  Pgap_fenestrated <- Pgap_continuous
  # In organs with sinusoid there is no basement membrane at gaps,
  # so we use basemenent membrane thickness = 0
  Pgap_sinusoidal <- dif/ wall_width #cm/s
  Papp_RYU <- 1.46e-6 # cm/s, at pH = 7.4 from Ryu et al. (2024) [https://doi.org/10.1016/j.chemosphere.2024.142390]
  Peff_monolayer <- Papp_RYU*3600 #cm/h
  
  #fraction of gaps in capilalry surface
  f_kidney <- 0.5 #Bulger et al. (1983) [10.1172/jci110950] 
  f_liver <- 0.08 #  Simon-Santamaria et al. (2010)[10.1093/gerona/glq108] (Antwi et al. (2023) give a range 2-20% [10.1371/journal.pone.0293526] )
  f_spleen <- 0.3 #assumption
  f_intestine <- 0.2 #assumption
  f_non_fenestrated <- 0.01
  
  VmK_api <- 0
  VmK_baso <- 0
  KmK_baso <- 1e20
  KmK_api <-   5e4
  KLfabp <-1/8e-6#mol/L, Sheng et al. (2018) [10.1007/s00204-017-2055-1]
  Ka <- 5.8e5# from Ryu  et al. (2024)#mol/L
  CLfeces_unscaled <- estimated_params[8]#in L/h/BW^(-0.25), scaling similar to Loccisano et al. (2012)
  CLfeces <- CLfeces_unscaled*BW^(-0.25) 
  
  f_alb_avail <-  estimated_params[9]
  f_fabp_avail  <- 2
  
  f_a2u_avail <- 0
  
  koff_alb <-  2#estimated_params[9]
  koff_fabp <-     0.017#estimated_params[10]#
  koff_a2u <- 1#estimated_params[11]
  #In order to scale transporter Vmax, we need to have the tissue weight to estimate
  # tissue protein
  PVIN <- 2.24e-2 #Brown et al. 1997, p 416, Table 5: 1.4+0.84
  VIN <- PVIN * BW #intestine volume kg=L
  
  #Liver
  PVL <- 3.66e-2 #Brown et al. 1997
  VL <- PVL * BW #liver volume kg=L
  
  PVLu <- 0.48e-2 #Brown et al. 1997, p 418-19 mean of values for male rat
  VLu <- PVLu * BW
  
  #Kidney
  PVK <- 7.3e-3 #Brown et al. 1997
  VK <- PVK * BW #kidney volume kg=L 
  
  # These are explained thoroughly in a later section
  f_tubular <- 0.8
  f_PTC_prot_to_tub_prot <- 0.6939
  
  MW = 414.07 #g/mol, PFOA molecular weight
  muscle_protein <- 158.45 #mg/g muscle (protein data from Cheek et al.,1971 
  #and muscle mass from Caster et al.,1956) ***DO WE BELIEVE THAT THIS IS PER GRAM OF TOTAL ORGAN OR TISSUE?
  intestine_protein <- muscle_protein
  intestine_protein_total <- intestine_protein*(1000*VIN)
  
  #Kidney
  kidney_protein_per_rat <- 1000*(0.218+0.225+0.212)/3#mg of protein per rat  (Addis 1936)
  rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
  rat_kidney_weight_addis <- rat_weight_addis*0.0073 # kidney fraction to BW, Brown (1997)
  kidney_protein_per_gram <- kidney_protein_per_rat/rat_kidney_weight_addis #mg of protein/g kidney
  
  kidney_cells = 1.47e07 #cells/g https://doi.org/10.1038/s41598-024-53270-2
  kidney_cells_total <- kidney_cells* (1000*VK)
  kidney_protein_total <- kidney_protein_per_gram* (1000*VK) #mg
  PTC_protein <- f_tubular*f_PTC_prot_to_tub_prot*kidney_protein_total
  
  #Oatp kidney
  VmK_Oatp_in_vitro <- 9.3 #nmol/mg protein/min (Weaver et al. 2010)
  VmK_Oatp_scaled <- 60*VmK_Oatp_in_vitro*MW*PTC_protein/1000  #physiologically scaled to in vivo, ug/h
  VmK_Oatp <- VmK_Oatp_scaled*RAFOatp_k #in vivo value, in  ug/h
  KmK_Oatp=  126.4*MW# [umol/L] * g/mol  --> ug/L, from Weaver et al. (2010)
  
  #oat1 kidney
  VmK_Oat1_in_vitro= 2.6 #nmol/mg protein/min (Weaver et al. 2010)
  VmK_Oat1_scaled = 60*VmK_Oat1_in_vitro*MW*PTC_protein/1000 #physiologically scaled to in vivo, ug/h
  VmK_Oat1= VmK_Oat1_scaled*RAFOat1 #in vivo value, in   ug/h
  KmK_Oat1= 43.2 * MW #umol/L (Weaver et al. 2010) --> ug/L
  
  #oat3 kidney
  VmK_Oat3_in_vitro= 3.8 #nmol/mg protein/min  (Weaver et al. 2010)
  VmK_Oat3_scaled = 60*VmK_Oat3_in_vitro*MW*PTC_protein/1000 #physiologically scaled to in vivo, ug/h
  VmK_Oat3 = VmK_Oat3_scaled*RAFOat3 #in vivo value, in   ug/h
  KmK_Oat3= 65.7 * MW #umol/L (Weaver et al. 2010) --> ug/L
  
  #Urat1 kidney
  VmK_Urat_in_vitro= 1520e-3 #nmol/mg protein/min  (Lin et al. 2023)
  VmK_Urat_scaled = 60*VmK_Urat_in_vitro*MW*PTC_protein/1000 #physiologically scaled to in vivo, ug/h
  VmK_Urat = VmK_Urat_scaled*RAFUrat #in vivo value, in   ug/h
  KmK_Urat = 820.04 * MW #umol/L (Lin et al. 2023) --> ug/L
  
  #Liver
  liver_protein_per_rat <- 1000*(1.52+1.53+1.52)/3#mg of protein per rat  (Addis 1936)
  rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
  rat_liver_weight_addis <- rat_weight_addis*0.0366 # liver fraction to BW, Brown (1997)
  liver_protein_per_gram <- liver_protein_per_rat/rat_liver_weight_addis #mg or protein/g liver
  #Total liver cells:2e09  https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=110895)
  HEPGL <- 104  #million hepatocytes/gram of rat liver (Fattah et al., 2016. https://doi.org/10.1124/dmd.115.066381 )
  #HEPGL_human = 117*1e6 #hepatocytes per g of liver (Sohlenius-Sternbeck et al. 2006) 
  # Scaled hepatobiliary clearance
  CL_hepatobiliary <- CL_int*1e-6*HEPGL*VL*60 #L/h
  
  #oatp1-liver
  VmL_Oatp_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
  VmL_Oatp_scaled = 60*VmL_Oatp_in_vitro*MW*liver_protein_per_gram*(VL*1000)/1000   #physiologically scaled to in vivo, ug/h
  VmL_Oatp = VmL_Oatp_scaled*RAFOatp_l #in vivo value, in  ug/h
  KmL_Oatp = KmK_Oatp #same as kidney
  
  #oatp2b1-liver
  VmL_Oatp2_in_vitro= 1493e-3 #nmol/mg protein/min  (Lin et al. 2023)
  #physiologically scaled to in vivo
  VmL_Oatp2_scaled = 60*VmL_Oatp2_in_vitro*MW*liver_protein_per_gram*(VL*1000)/1000  #ug/h
  VmL_Oatp2 = VmL_Oatp2_scaled*RAFOatp2_l #in vivo value, in  ug/h
  KmL_Oatp2 = 148.68*MW #umol/L (Lin et al. 2023) --> ug/L
  
  #Ntcp liver
  VmL_Ntcp_in_vitro= 3#nmol/mg protein/min   Ruggiero et al. 2021
  #physiologically scaled to in vivo
  VmL_Ntcp_scaled = 60*VmL_Ntcp_in_vitro*MW*liver_protein_per_gram*(VL*1000)/1000 # ug/h 
  VmL_Ntcp = VmL_Ntcp_scaled*RAFNtcp #in vivo value, in  ug/h
  KmL_Ntcp= 20 * MW #umol/L, Ruggiero et al. 2021 --> ug/L
  
  #Lung
  lung_protein_per_gram <- 134 # 134 mg/mL tissue --> 134 mg/g tissue, Figure 2, https://doi.org/10.1007/s00580-021-03242-z 
  
  #oatp-lung-ap (from ALF to tissue)
  VmLu_Oatp_ap_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
  VmLu_Oatp_ap_scaled = 60*VmLu_Oatp_ap_in_vitro*MW*lung_protein_per_gram*(VLu*1000)/1000   #physiologically scaled to in vivo, ug/h
  VmLu_Oatp_ap = VmLu_Oatp_ap_scaled*RAFOatp_lu_ap #in vivo value, in  ug/h
  KmLu_Oatp_ap = KmK_Oatp #same as kidney
  
  #oatp-lung-bas (from IS to tissue)
  VmLu_Oatp_bas_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
  VmLu_Oatp_bas_scaled = 60*VmLu_Oatp_bas_in_vitro*MW*lung_protein_per_gram*(VLu*1000)/1000   #physiologically scaled to in vivo, ug/h
  VmLu_Oatp_bas = VmLu_Oatp_bas_scaled*RAFOatp_lu_bas #in vivo value, in  ug/h
  KmLu_Oatp_bas = KmK_Oatp #same as kidney
  
  #Intestine
  #oatp2b1-intestine
  VmIn_Oatp2_in_vitro= 456.63e-3 #nmol/mg protein/min  (Kimura et al., 2017) 
  #assuming that the mediated transport is performed only by this transporter
  VmIn_Oatp2_scaled = 60*VmIn_Oatp2_in_vitro*MW*intestine_protein_total/1000   #physiologically scaled to in vivo, ug/h
  VmIn_Oatp2 = VmIn_Oatp2_scaled*RAFOatp2_Int #in vivo value, in  ug/h, same RAF as in liver
  KmIn_Oatp2 = 8.3*MW #umol/L (Kimura et al., 2017) --> ug/L
  
  Mr_albumin <- 66500#g/mol
  CalbB_init  <- f_alb_avail*mean(c(593,551, 591, 509, 535))*1e-6 #mol/L, from Rose & Klemcke (2015) [PMID: 26424242]
  #Albumin concentration in blood and interstitial fluid compartments(mol/m^3 = 1e-6* nmol/g)
  #CalbB_init <- f_alb_avail*486*1e-06 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  
  CalbKB_init <- CalbB_init
  CalbLB_init <- CalbB_init
  CalbSTB_init <- CalbB_init
  CalbINB_init <- CalbB_init
  CalbMB_init <- CalbB_init
  CalbAB_init <- CalbB_init
  CalbRB_init <- CalbB_init
  CalbBoB_init <- CalbB_init
  CalbLuB_init <- CalbB_init
  CalbSPB_init <- CalbB_init
  CalbGoB_init <- CalbB_init
  CalbHB_init <- CalbB_init
  CalbBrB_init <- CalbB_init
  CalbSKB_init <- CalbB_init
  
  #Interstitial/plasma concentration ratio (IPR)
  #values from Kawai et al., 1994, Table C-I
  IPR_K = 0.5
  IPR_L = 0.5
  IPR_ST = 0.5
  IPR_IN = 0.9
  IPR_M = 0.6
  IPR_A = 0.5
  IPR_Lu = 0.5
  IPR_Sp = 0.5
  IPR_H = 0.5
  IPR_SK = 1
  IPR_Br = 0.5
  IPR_Go = 0.5 #assumption
  IPR_Bo = 0.5 #assumption
  IPR_R = (IPR_K+IPR_L+IPR_ST+IPR_IN+IPR_M+IPR_A+IPR_Lu+IPR_Sp+IPR_H+IPR_SK+IPR_Br+IPR_Go+IPR_Bo)/13 #average IPR of all the included organs (kg=L)
  
  Hct <- 0.41
  
  #Convert blood concentration to plasma concentration first to derive the interstitial concentration
  CalbKF_init  <- (CalbKB_init/(1-Hct)) * IPR_K
  CalbLF_init <- (CalbLB_init/(1-Hct))* IPR_L 
  CalbSTF_init <- (CalbSTB_init/(1-Hct))* IPR_ST
  CalbINF_init <- (CalbINB_init/(1-Hct))* IPR_IN
  CalbMF_init <- (CalbMB_init/(1-Hct))* IPR_M
  CalbAF_init <- (CalbAB_init/(1-Hct))* IPR_A
  CalbRF_init <- (CalbRB_init/(1-Hct))* IPR_R
  CalbBoF_init <- (CalbBoB_init/(1-Hct))* IPR_Bo
  CalbLuF_init <- (CalbLuB_init/(1-Hct))* IPR_Lu
  CalbSPF_init <- (CalbSPB_init/(1-Hct))* IPR_Sp
  CalbGoF_init <- (CalbGoB_init/(1-Hct))* IPR_Go
  CalbHF_init <- (CalbHB_init/(1-Hct))* IPR_H
  CalbBrF_init <- (CalbBrB_init/(1-Hct))* IPR_Br
  CalbSKF_init <- (CalbSKB_init/(1-Hct))* IPR_SK
  CalbLuAF_init <- (10/100/(1-Hct)) * CalbB_init #based on Woods et al. 2015 statement https://doi.org/10.1016/j.jconrel.2015.05.269
  
  
  #Alpha2mu-globulin concentration in kidney tissue (mol/L)
  if (sex == "M"){
    a2u_globulin_k = 8.77*kidney_protein_total*1e-3/VK #mg/L, 8.77 mg/g kidney protein from https://doi.org/10.1016/0300-483X(86)90197-6 
    Ca2uKT_init <- f_a2u_avail*(a2u_globulin_k*1e-3/15.5e3) #[mol/L]
    
    #Ca2uKT_init <- 321.51*1e-3 #[umol/L]*1e-3 -->(mol/m3), from Cheng et al. (2017)
    
  }else if(sex == "F"){
    Ca2uKT_init <- 0 #mol/L
  }
  
  
  #LFABP concentration in kidney and liver tissue (mol/m^3)
  L_FABP_L = 28.2e-3*liver_protein_per_rat/VL #mg/L, 28.2 ug/mg cytosolic protein from https://doi.org/10.1016/S0021-9258(18)34463-6
  #cytosolic protein is 96.3% of the total liver protein, https://doi.org/10.18632/aging.101009
  CFabpLT_init = f_fabp_avail*(L_FABP_L*1e-3/14e3) #[mol/L]
  
  
  #LFABP concentration in kidney and liver tissue (mol/m^3)
  CFabpKT_init <- f_fabp_avail*2.65*1e-6  #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  
  #======Table S2=======#
  #Equilibrium association constant (m^3/mol= 10^-3*M-1) for albumin(Ka), LFABP(KL_fabp),
  #and alpha2mu-globulin(Ka2u). See SI section S2-2 for details
  
  #Ka <-  24.18 #3.1*7.8 m3/mol multiplying by number of binding sites (Cheng et al. 2021)
  #Ka <-  1e05*1e-3 #[L/mol]*1e-3--->m3/mol
  Ka2u <- 5*1e02 #[L/mol]*1e-3--->m3/mol, value from Cheng et al. (2017)
  
  
  kon_alb <- Ka * koff_alb #1/M/s
  kon_a2u <- Ka2u * koff_a2u#1/M/s
  kon_fabp <- KLfabp * koff_fabp #1/M/s
  
  # Following the calculations  of Lin et al. (2023) for Caco-2 cells
  ClINFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
  Awell = 9 #cm^2 (for a 35 mm culture dish)
  Swell = 1.12 #cm^2
  well_protein = 0.346 #mg protein
  protein_per_well = (well_protein * Awell)/Swell #mg protein/well
  RAF_papp <- 1 
  Peff_LIN = RAF_papp*(ClINFT_unscaled*60*1e-06*1e3*protein_per_well)/(Awell *2) #cm/h,at  pH = 6.0
  #  Lin et al. (2023) used data from kimura et al. (2017) [http://dx.doi.org/10.1016/j.toxlet.2017.05.012]
  # which is more appropriate for 
  # For endothelial and cellular permeability we use the Ryu et al. (2024) value
  #Papp = Peff_RYU
  k_gut_in = ( (2*Peff_monolayer/100) * fixed_params$AINL)*1000 #L/h
  k_gut_out = ( (2*Peff_monolayer/100) * fixed_params$AINL)*1000 #L/h
  
  
  #passive diffusion rates, in L/h
  kLFLT = ((2*Peff_monolayer/100) * fixed_params$AcL)*1000 #m^3/h * 1000 --> L/h
  #kLTLbile = ((Papp/100) * AcLBilec)*1000 #m^3/h * 1000 --> L/h
  kMFMT = ((2*Peff_monolayer/100) * fixed_params$AcM)*1000 #m^3/h * 1000 --> L/h
  kSTFSTT = ((2*Peff_monolayer/100) * fixed_params$AcST)*1000 #m^3/h * 1000 --> L/h 
  kINFINT = ((2*Peff_monolayer/100) * fixed_params$AcIN)*1000 #m^3/h * 1000 --> L/h 
  kAFAT = ((2*Peff_monolayer/100) * fixed_params$AcA)*1000 #m^3/h * 1000 --> L/h 
  kLuTLuF = ((2*Peff_monolayer/100) * fixed_params$AcLu)*1000 #m^3/h * 1000 --> L/h
  kLuTLuAF = ((2*Peff_monolayer/100) * fixed_params$AcALF)*1000 #m^3/h * 1000 --> L/h
  kSPFSPT = ((2*Peff_monolayer/100) * fixed_params$AcSP)*1000 #m^3/h * 1000 --> L/h 
  kHFHT = ((2*Peff_monolayer/100) * fixed_params$AcH)*1000 #m^3/h * 1000 --> L/h 
  kBrFBrT = ((2*Peff_monolayer/100) * fixed_params$AcBr)*1000 #m^3/h * 1000 --> L/h 
  kGoFGoT = ((2*Peff_monolayer/100) * fixed_params$AcGo)*1000 #m^3/h * 1000 --> L/h 
  kSKFSKT = ((2*Peff_monolayer/100) * fixed_params$AcSK)*1000 #m^3/h * 1000 --> L/h
  kBoFBoT = ((2*Peff_monolayer/100) * fixed_params$AcBo)*1000 #m^3/h * 1000 --> L/h
  kRFRT = ((2*Peff_monolayer/100) * fixed_params$AcR)*1000 #m^3/h*1000 --> L/h 
  
  #Diffusion rates in L/h between renal tubule filtrate and tubule cells
  kPtcTu <- ((2*Peff_monolayer/100) * fixed_params$APT) *1000 #diffusion between proximal tubule cells and tubule filtrate
  kDalcTu <- ((2*Peff_monolayer/100) * fixed_params$ADAL) *1000 #diffusion between descending/ascending cells and tubule filtrate
  kDtcTu <- ((2*Peff_monolayer/100) * fixed_params$ADT) *1000 #diffusion between distal tubule cells and tubule filtrate
  kCdcTu <- ((2*Peff_monolayer/100) * fixed_params$ACD) *1000 #diffusion between collecting duct cells and tubule filtrate
  
  #Diffusion rates in L/h between  tubule cells and interstitial space
  kDtcF <- ((2*Peff_monolayer/100) * fixed_params$AcK_DTC) *1000
  kPtcF <- ((2*Peff_monolayer/100) * fixed_params$AcK_PTC) *1000
  kDalcF <- ((2*Peff_monolayer) * fixed_params$AcK_DALC) *1000 #diffusion between proximal tubule cells and interstitial space
  kCdcF <- ((2*Peff_monolayer/100) * fixed_params$AcK_CDC) *1000 #diffusion between descending/ascending cells and interstitial space
  kKTrestF  <- ((2*Peff_monolayer/100) * fixed_params$AcKTrest) *1000 #diffusion between rest of kidney cells and interstitial space
  
  # All effective permeabilities are in cm/h and are multiplied by 10, because they are divided by 100 to be converted to meters
  # and later in the code they mutiply SA which is in m^2 becoming m^3, which is multiplied by 1000 to become Liters
  # Transendothelial diffusion
  Ptrans_diff_K <- Peff_monolayer*10 * (1-f_kidney) #mm/h 
  Ptrans_diff_L <- Peff_monolayer*10* (1-f_liver) #mm/h
  Ptrans_diff_ST <- Peff_monolayer*10* (1-f_intestine) #mm/h
  Ptrans_diff_IN <- Peff_monolayer*10* (1-f_non_fenestrated)  #mm/h
  Ptrans_diff_A <- Peff_monolayer*10 * (1-f_non_fenestrated)#mm/h
  Ptrans_diff_M <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  Ptrans_diff_R <- Peff_monolayer*10* (1-f_non_fenestrated)#mm/h
  Ptrans_diff_Lu <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  Ptrans_diff_SP <- Peff_monolayer*10* (1-f_spleen)  #mm/h
  Ptrans_diff_H <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  Ptrans_diff_Br <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  Ptrans_diff_Go <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  Ptrans_diff_SK <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  Ptrans_diff_Bo <- Peff_monolayer*10* (1-f_non_fenestrated) #mm/h
  
  
  
  #Estimation of permeability through capillary fenestra and discontinuities:
  PparaKi <- Pgap_fenestrated*3600*10*f_kidney#mm/h
  PparaLi <- Pgap_sinusoidal*3600*10 * f_liver #mm/h
  PparaSt <- Pgap_fenestrated*3600*10 * f_non_fenestrated #mm/h
  PparaIn <- Pgap_fenestrated*3600*10 * f_intestine #mm/h
  PparaMu <- Pgap_continuous*3600*10 * f_non_fenestrated #mm/h
  PparaAd <- Pgap_continuous*3600*10 * f_non_fenestrated #mm/h
  PparaRe <- Pgap_continuous*3600*10 * f_non_fenestrated #mm/h
  PparaLu <- Pgap_continuous*3600*10 * f_non_fenestrated #mm/h
  PparaSp <- Pgap_sinusoidal*3600*10 * f_spleen  #mm/h
  PparaHt <- Pgap_continuous*3600*10 * f_non_fenestrated #mm/h
  PparaBr <- 0
  PparaGo <- Pgap_fenestrated*3600*10 * f_non_fenestrated #mm/h
  PparaSk <- Pgap_continuous*3600*10 * f_non_fenestrated #mm/h
  PparaBo <- Pgap_continuous*3600*10  * f_non_fenestrated#mm/h
  
  return(list(
    
    'CalbB_init'= CalbB_init, 'CalbKF_init'=CalbKF_init, 'CalbLF_init'=CalbLF_init,
    'CalbMF_init'=CalbMF_init, 'CalbAF_init'=CalbAF_init, 'CalbRF_init'=CalbRF_init,
    'CalbBoF_init'=CalbBoF_init, 'CalbLuF_init' =CalbLuF_init,
    'CalbLuAF_init'=CalbLuAF_init, 'CalbSPF_init' =CalbSPF_init,
    'CalbGoF_init' =CalbGoF_init, 'CalbHF_init' =CalbHF_init,
    'CalbBrF_init' =CalbBrF_init, 'CalbSTF_init' =CalbSTF_init,
    'CalbINF_init' =CalbINF_init, 'CalbSKF_init' =CalbSKF_init, 
    
    'CalbKB_init'=CalbKB_init,'CalbLB_init'=CalbLB_init,'CalbSTB_init'=CalbSTB_init,
    'CalbINB_init'=CalbINB_init, 'CalbMB_init'=CalbMB_init,'CalbAB_init'=CalbAB_init,
    'CalbRB_init'=CalbRB_init,'CalbBoB_init'=CalbBoB_init,
    'CalbLuB_init'=CalbLuB_init, 'CalbSPB_init'=CalbSPB_init,'CalbGoB_init'=CalbGoB_init,
    'CalbHB_init'=CalbHB_init,'CalbBrB_init'=CalbBrB_init,'CalbSKB_init'=CalbSKB_init,
    
    'Ca2uKT_init'=Ca2uKT_init,'CFabpKT_init'=CFabpKT_init,'CFabpLT_init'=CFabpLT_init, 
    
    'Ka'=Ka, 'Ka2u'=Ka2u, 'KLfabp'=KLfabp,
    
    "koff_alb" = koff_alb, "koff_a2u" = koff_a2u, "koff_fabp" = koff_fabp,
    "kon_alb" = kon_alb, "kon_a2u" = kon_a2u, "kon_fabp" = kon_fabp,
    
    "CL_hepatobiliary" = CL_hepatobiliary, 'CLfeces' = CLfeces,
    'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'= KmL_Oatp, 'VmL_Ntcp'= VmL_Ntcp,
    'VmL_Oatp2'=VmL_Oatp2, 'KmL_Oatp2'= KmL_Oatp2, 
    'VmIn_Oatp2'=VmIn_Oatp2, 'KmIn_Oatp2'= KmIn_Oatp2,
    'KmL_Ntcp'= KmL_Ntcp,'VmK_Oatp'= VmK_Oatp, 
    'KmK_Oatp'=KmK_Oatp, 'VmLu_Oatp_ap'= VmLu_Oatp_ap,
    'KmLu_Oatp_ap'=KmLu_Oatp_ap, 'VmLu_Oatp_bas'= VmLu_Oatp_bas,
    'KmLu_Oatp_bas'=KmLu_Oatp_bas,
    
    
    'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 
    'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
    'VmK_Urat'=VmK_Urat, 'KmK_Urat'=KmK_Urat, 
    
    'KmK_baso' = KmK_baso, 'KmK_api' = KmK_api,
    'VmK_baso' = VmK_baso,'VmK_api' = VmK_api,
    
    'k_gut_in' = k_gut_in, 'k_gut_out' = k_gut_out,
    'kKTrestF'=kKTrestF, 'kCdcF' = kCdcF, 'kDalcF' = kDalcF, 'kPtcF' = kPtcF, 'kDtcF' = kDtcF,
    'kPtcTu'=kPtcTu, 'kDalcTu' = kDalcTu, 'kDtcTu' = kDtcTu, 'kCdcTu' = kCdcTu, 
    'kLFLT'=kLFLT,  'kAFAT'=kAFAT, 
    'kRFRT'=kRFRT,
    'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF, 'kLuTLuAF'=kLuTLuAF, 'kSPFSPT' =kSPFSPT,
    'kSTFSTT' =kSTFSTT, 'kINFINT' =kINFINT, 'kHFHT' =kHFHT,
    'kBrFBrT' =kBrFBrT, 'kGoFGoT' =kGoFGoT,
    'kSKFSKT' =kSKFSKT, 'kBoFBoT'=kBoFBoT,
    
    'Ptrans_diff_K'=Ptrans_diff_K, 'Ptrans_diff_L'=Ptrans_diff_L, 
    'Ptrans_diff_A'=Ptrans_diff_A, 'Ptrans_diff_M'=Ptrans_diff_M, 'Ptrans_diff_R'=Ptrans_diff_R, 'Ptrans_diff_Lu' = Ptrans_diff_Lu,
    'Ptrans_diff_SP'=Ptrans_diff_SP, 'Ptrans_diff_H'=Ptrans_diff_H, 'Ptrans_diff_Br'=Ptrans_diff_Br, 'Ptrans_diff_ST' = Ptrans_diff_ST,
    'Ptrans_diff_IN'=Ptrans_diff_IN, 'Ptrans_diff_Go'=Ptrans_diff_Go,
    'Ptrans_diff_SK' = Ptrans_diff_SK,  'Ptrans_diff_Bo' = Ptrans_diff_Bo,
    
    "PparaKi" = PparaKi,"PparaLi" = PparaLi,"PparaSt" = PparaSt,"PparaIn" = PparaIn,
    "PparaMu" = PparaMu,"PparaAd" = PparaAd,"PparaRe" = PparaRe,"PparaLu" = PparaLu,
    "PparaSp" = PparaSp,"PparaHt" = PparaHt,"PparaBr" = PparaBr,"PparaGo" = PparaGo,
    "PparaSk" = PparaSk,"PparaBo" = PparaBo
    
    
  ))
}  

# Set up simulations for the 9th case, i.e. Kim (2016) ORAL male blood
BW <- 0.25  #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[9]])
params <- c(fixed_params[[9]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_OR_Mblood <-  solution[, c("time", "Cplasma")]

######################################################################################
# Set up simulations for the 10th case, i.e. Kim (2016) IV male blood
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[10]])
params <- c(fixed_params[[10]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs

sample_time=c(0, 5/60, seq(1,288,1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_IV_Mblood <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 25th case, i.e. Kim (2016) ORAL female blood
BW <- 0.25  #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[25]])
params <- c(fixed_params[[25]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= c(seq(0, 1.2, 0.2), seq(1.5,24,0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_OR_Fblood <-  solution[, c("time", "Cplasma")]

########################################################################################
# Set up simulations for the 26th case, i.e. Kim (2016) IV female blood
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[26]])
params <- c(fixed_params[[26]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0, 25, 0.5)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]


preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] <- preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] /1000
preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] <- preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] /1000
preds_kim_IV_Mblood <- preds_kim_IV_Mblood[2:dim(preds_kim_IV_Mblood)[1],]
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000
preds_kim_IV_Fblood <- preds_kim_IV_Fblood[2:dim(preds_kim_IV_Fblood)[1],]

experiment9 <- list()
# Convert Kim ORAL male blood from long to wide format using reshape
experiment9[[1]] <- reshape(kim_OR_Mblood[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9[[1]]) <- c("Time",unique(kim_OR_Mblood$Tissue))


# Convert Kim IV male blood from long to wide format using reshape
experiment9[[2]] <- reshape(kim_IV_Mblood[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9[[2]]) <- c("Time",unique(kim_IV_Mblood$Tissue))




#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment9[[3]] <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9[[3]]) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment9[[4]]<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9[[4]]) <- c("Time",unique(kim_IV_Fblood$Tissue))



# Put the experiments in a list
experiments <- list(experiment9 = experiment9)


colnames(preds_kim_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_kim_IV_Mblood) <- c ("Time", "Plasma")
colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")



# Create a list containing the corresponding predictions
simulations <- list(predictions9 = list(preds_kim_OR_Mblood = preds_kim_OR_Mblood,  
                                        preds_kim_IV_Mblood = preds_kim_IV_Mblood,
                                        preds_kim_OR_Fblood = preds_kim_OR_Fblood, 
                                          preds_kim_IV_Fblood = preds_kim_IV_Fblood))



plot_names <- list(
                   predictions9 = list(preds_kim_OR_Mblood = "Oral",  
                                       preds_kim_IV_Mblood = "IV",
                                       preds_kim_OR_Fblood = "Oral", 
                                         preds_kim_IV_Fblood = "IV") )


library(ggplot2)
library(patchwork)

create.plots <- function(predictions, observations, compartment, plot_name) {
  cls <- c("Prediction" = "#0072B2", "Observation" = "#D55E00")
  
  if(plot_name %in% c( "ELF",  "Brain",  "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                       "Heart", "Brain", "Gonads", "Stomach", "Intestine") ){
    y_name <-  expression("PFOA concentration (" * mu * "g/g tissue)")
  }else if(plot_name %in% c( "Urine Female 1 mg/kg", "Feces Female 1 mg/kg",  "Urine Female 5 mg/kg", 
                             "Feces Female 5 mg/kg", "Urine Female 25 mg/kg", "Feces Female 25 mg/kg",
                             "Urine Male 1 mg/kg","Feces Male 1 mg/kg",  "Urine Male 5 mg/kg", 
                             "Feces Male 5 mg/kg","Urine Male 25 mg/kg", "Feces Male 25 mg/kg")){
    y_name <- expression("PFOA Mass (" * mu * "g)")
  }else{
    y_name <- expression("PFOA concentration (" * mu * "g/mL)")
  }
  ggplot(data = predictions) +
    geom_line(aes_string(x = "Time", y = rlang::expr(!!compartment), color = '"Prediction"'), 
              size = 1.2, alpha = 0.9) +
    geom_point(data = observations, 
               aes_string(x = "Time", y = rlang::expr(!!compartment), color = '"Observation"'), 
               size = 3.5, shape = 21, stroke = 1, fill = "white") +
    labs(title = plot_name, 
         y = y_name,
         x = "Time (hours)",
         color = "Type") +
    scale_color_manual(values = cls) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
}

# Iterate over experiments
for (i in seq_along(experiments)) {
  all_plots <- list()
  
  for (j in seq_along(experiments[[i]])) {
    observations <- experiments[[i]][[j]]
    predictions <- simulations[[i]][[j]]
    compartments <- names(predictions)[-1]  # Skip "Time"
    
    if(length(compartments) > 1){
      compartment_plots <- lapply(compartments, function(compartment) {
        create.plots(
          predictions = predictions,
          observations = observations,
          compartment = compartment,
          plot_name = compartment
        )
      })
    }else{
      compartment_plots <- lapply(compartments, function(compartment) {
        create.plots(
          predictions = predictions,
          observations = observations,
          compartment = compartment,
          plot_name = plot_names[[i]][[j]]
        )
      })
    }
    
    all_plots <- c(all_plots, compartment_plots)  # Flatten here
  }
  
  # Combine with patchwork
  if (length(all_plots) == 1) {
    final_plot <- all_plots[[1]] + theme(legend.position = "bottom")
  } else {
    final_plot <- wrap_plots(all_plots, ncol = 3) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
  }
  
  print(final_plot)
}