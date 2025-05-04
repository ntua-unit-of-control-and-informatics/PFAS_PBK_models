
create_variable_params <- function(BW,sex,  estimated_params, fixed_params){
  
  # BW in kg
  # Cheng and Ng 2017 Table S1
  # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
  #assuming the density of tissue is 1 g/mL.
  # Estimated parameters
  if(sex == "M"){
    RAFOatp_k <- 3.416684e+01
    VmK_api <- 0
    f_fabp_avail  <-1.825424e-01
    
  }else{
    RAFOatp_k <-1.681949e-06
    VmK_api <- 0#estimated_params[15]
    f_fabp_avail  <-4.371305e-02
    
  }
  RAFOat3 <- 8.640175e+03
  CL_int <- 9.217010e+01 #uL/min/million hepatocytes
  RAFOatp_l <- 2.815961e+00
  RAFUrat <- RAFOatp_k
  RAFOat1 <- RAFOat3
  RAFOatp2_l <- RAFOatp_l
  RAFNtcp <- RAFOatp_l
  RAFOatp2_Int <- 6.799068e-07
  
  RAF_papp <- 1
  Papp_RYU <- RAF_papp*1.46e-6*3600 # cm/h, at pH = 7.4 from Ryu et al. (2024) [https://doi.org/10.1016/j.chemosphere.2024.142390]
  #Aqueous boundary layer permeability
  Paq <- 100*1e-6*3600 #cm/h, from Korzekwa et al. (2022) [doi:10.1111/cts.13314]
  #Paracellular permeability
  Ppara <- 0.2*1e-6*3600 #cm/h, from Korzekwa et al. (2022) [doi:10.1111/cts.13314]
  Peff_monolayer <- (Paq*Papp_RYU - Paq*Ppara)/(Paq+Ppara-Papp_RYU)
  
  Pgap_continuous <- 2.774701e-02 *3600#, cm/h, permeability by diffusion through tight junctions in continuous capillary types
  Pgap_fenestrated <- Pgap_continuous#estimated_params[9] # cm/s,permeability by diffusion through fenestra in fenestrated capillary types
  Pgap_sinusoidal <-Pgap_continuous#estimated_params[10] # cm/s,permeability by diffusion through discontinuities in sinusoidal capillary types
  # pH_isf <- 6; pH_iw <- 7.4; pKaa <- 3.8; pKab <- 1
  # Y = 1 + 10^(pKab - pH_isf) + 10^(pH_isf - pKaa) #1/Y is the neutral fraction of PFOA in ISF
  # X = 1 + 10^(pKab - pH_iw) + 10^(pH_iw - pKaa) #1/X is the neutral fraction of PFOA in cytosol
  
  f_alb_avail<- 3.765961e+00
  
  koff_alb <-  2.796183691
  koff_fabp <-   0.004720769#
  koff_a2u <- koff_alb
  
  VmK_baso <- 0
  KmK_baso <- 1e20
  KmK_api <-   5e4
  KLfabp <- 5.261153e+06#(1.2e5+4e4+1.9e4) #[L/mol] , value from Cheng et al. (2017)
  Ka <- 1.507899e+04 # 6e5 from Rue et al. (2024)#mol/L
  CLfeces_unscaled <-9.551628e-05#in L/h/BW^(-0.25), scaling similar to Loccisano et al. (2012)
  CLfeces <- CLfeces_unscaled*BW^(-0.25) 
  

  kabsUA <- 1.423399e-05 #L/h/m^2
  kCLEal <- 1.148107e-02 #L/h/m^2
  kCLEua <- 1.165535e-01 #L/h/m^2
  k_desorption <- estimated_params[2]
  RAFOatp_lu_ap <- estimated_params[1]
  RAFOatp_lu_bas <- RAFOatp_lu_ap
  
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
  
  #oatp-lung-ap (from ELF to tissue)
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
  
  CalbKF_init  <- CalbKB_init* IPR_K
  CalbLF_init <- CalbLB_init* IPR_L 
  CalbSTF_init <- CalbSTB_init* IPR_ST
  CalbINF_init <- CalbINB_init* IPR_IN
  CalbMF_init <- CalbMB_init* IPR_M
  CalbAF_init <- CalbAB_init* IPR_A
  CalbRF_init <- CalbRB_init* IPR_R
  CalbBoF_init <- CalbBoB_init* IPR_Bo
  CalbLuF_init <- CalbLuB_init* IPR_Lu
  CalbSPF_init <- CalbSPB_init* IPR_Sp
  CalbGoF_init <- CalbGoB_init* IPR_Go
  CalbHF_init <- CalbHB_init* IPR_H
  CalbBrF_init <- CalbBrB_init* IPR_Br
  CalbSKF_init <- CalbSKB_init* IPR_SK
  CalbLuELF_init <- 10/100 * CalbB_init #based on Woods et al. 2015 statement https://doi.org/10.1016/j.jconrel.2015.05.269
  
  
  #Alpha2mu-globulin concentration in kidney tissue (mol/L)
  if (sex == "M"){
    a2u_globulin_k = 8.77*kidney_protein_total*1e-3/VK #mg/L, 8.77 mg/g kidney protein from https://doi.org/10.1016/0300-483X(86)90197-6 
    Ca2uKT_init <- f_alb_avail*(a2u_globulin_k*1e-3/15.5e3) #[mol/L]
    
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
  
  kUAB <- kabsUA * fixed_params$Nasal_SA #absorption rate from upper airways to blood, in L/h
  CLEal <- kCLEal * fixed_params$Alveolar_SA #clearance rate from alveolar lining fluid to stomach, in L/h
  CLEua <- kCLEua * fixed_params$Nasal_SA #clearance rate rate from upper airways to stomach, in L/h
  
  
  #passive diffusion rates, in L/h
  kLFLT = ((2*Peff_monolayer/100) * fixed_params$AcL)*1000 #m^3/h * 1000 --> L/h
  #kLTLbile = ((Papp/100) * AcLBilec)*1000 #m^3/h * 1000 --> L/h
  kMFMT = ((2*Peff_monolayer/100) * fixed_params$AcM)*1000 #m^3/h * 1000 --> L/h
  kSTFSTT = ((2*Peff_monolayer/100) * fixed_params$AcST)*1000 #m^3/h * 1000 --> L/h 
  kINFINT = ((2*Peff_monolayer/100) * fixed_params$AcIN)*1000 #m^3/h * 1000 --> L/h 
  kAFAT = ((2*Peff_monolayer/100) * fixed_params$AcA)*1000 #m^3/h * 1000 --> L/h 
  kLuTLuF = ((2*Peff_monolayer/100) * fixed_params$AcLu)*1000 #m^3/h * 1000 --> L/h
  kLuTLuELF = ((2*Peff_monolayer/100) * fixed_params$AcELF)*1000 #m^3/h * 1000 --> L/h
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
  Ptrans_diff_K <- Peff_monolayer*10 #mm/h
  Ptrans_diff_L <- Peff_monolayer*10 #mm/h
  Ptrans_diff_ST <- Peff_monolayer*10 #mm/h
  Ptrans_diff_IN <- Peff_monolayer*10 #mm/h
  Ptrans_diff_A <- Peff_monolayer*10 #mm/h
  Ptrans_diff_M <- Peff_monolayer*10 #mm/h
  Ptrans_diff_R <- Peff_monolayer*10#mm/h
  Ptrans_diff_Lu <- Peff_monolayer*10 #mm/h
  Ptrans_diff_SP <- Peff_monolayer*10 #mm/h
  Ptrans_diff_H <- Peff_monolayer*10 #mm/h
  Ptrans_diff_Br <- Peff_monolayer*10 #mm/h
  Ptrans_diff_Go <- Peff_monolayer*10 #mm/h
  Ptrans_diff_SK <- Peff_monolayer*10 #mm/h
  Ptrans_diff_Bo <- Peff_monolayer*10 #mm/h
  
  #Estimation of permeability through capillary fenestra and discontinuities:
  PparaKi <- Pgap_fenestrated*10 #mm/h
  PparaLi <- Pgap_sinusoidal*10 #mm/h
  PparaSt <- Pgap_fenestrated*10 #mm/h
  PparaIn <- Pgap_fenestrated*10 #mm/h
  PparaMu <- Pgap_continuous*10 #mm/h
  PparaAd <- Pgap_continuous*10 #mm/h
  PparaRe <- Pgap_continuous*10 #mm/h
  PparaLu <- Pgap_continuous*10 #mm/h
  PparaSp <- Pgap_sinusoidal*10 #mm/h
  PparaHt <- Pgap_continuous*10 #mm/h
  PparaBr <- 0
  PparaGo <- Pgap_fenestrated*10 #mm/h
  PparaSk <- Pgap_continuous*10 #mm/h
  PparaBo <- Pgap_continuous*10 #mm/h
  
  return(list(
    
    'CalbB_init'= CalbB_init, 'CalbKF_init'=CalbKF_init, 'CalbLF_init'=CalbLF_init,
    'CalbMF_init'=CalbMF_init, 'CalbAF_init'=CalbAF_init, 'CalbRF_init'=CalbRF_init,
    'CalbBoF_init'=CalbBoF_init, 'CalbLuF_init' =CalbLuF_init,
    'CalbLuELF_init'=CalbLuELF_init, 'CalbSPF_init' =CalbSPF_init,
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
    
    'kUAB'= kUAB, 'CLEal'=CLEal, 'CLEua'=CLEua,
    'k_gut_in' = k_gut_in, 'k_gut_out' = k_gut_out,
    'kKTrestF'=kKTrestF, 'kCdcF' = kCdcF, 'kDalcF' = kDalcF, 'kPtcF' = kPtcF, 'kDtcF' = kDtcF,
    'kPtcTu'=kPtcTu, 'kDalcTu' = kDalcTu, 'kDtcTu' = kDtcTu, 'kCdcTu' = kCdcTu, 
    'kLFLT'=kLFLT,  'kAFAT'=kAFAT, 
    'kRFRT'=kRFRT,
    'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF, 'kLuTLuELF'=kLuTLuELF, 'kSPFSPT' =kSPFSPT,
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
    "PparaSk" = PparaSk,"PparaBo" = PparaBo,
    
    "k_desorption" = k_desorption
    
    
  ))
}  

create_fixed_params <- function(user.input){
  with(as.list(user.input),{
    
    #permeabilities correction factor
    kabs_st <- 0 #m/h
    #units conversion from Cheng 2017R, time-> h, PFOA mass->ng, tissues mass-> g
    Hct <- 0.41 #hematocrit for rats, https://doi.org/10.1080/13685538.2017.1350156 mean value for both males and females
    #======Table S1=======#    
    
    #Blood
    PVB <- 54e-3 #13.5 mL/244 g=0.055 mL/g~55e-3 mL/g (kg=L), Davies et al. 1993, for BW = 0.25 kg
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 31.2e-3 
    Vplasma <- PVplasma * BW #plasma volume kg=L
    VVen <- BW*11.3/250 	#volume of venous plasma (L); from doi:10.1007/bf02353860
    VArt <- BW*5.6/250	#volume of arterial plasma (L); from doi:10.1007/bf02353860
    
    #Kidney
    PVK <- 7.3e-3 #Brown et al. 1997
    VK <- PVK * BW #kidney volume kg=L 
    PVKB <- 0.16 #Brown et al. 1997
    VKB <- PVKB * PVK * BW #kidney blood volume kg=L
    PVKF <- 0.13 # Wolgast et al. (1981) [https//doi.org/10.1152/ajprenal.1981.241.2.F105]
    VKF <- PVKF * PVK * BW #kidney interstitial fluid volume kg=L
    
    # Here we assume that the length of each segment of the renal tubule is 
    # analogous to the body weight of the kidney in order to obtain
    # volumes and surface areas as a function of BW
    Vki_ref <- 2.9/1000 # Gilmer et al. (2018) handled data from Sperber, 1944 (Thesis: “Studies on the Mammalian Kidney")
    #Total Length of renal tubules, from Gilmer et al. (2018) [https://doi.org/10.1152/ajprenal.00219.2018.]
    LPT_slp <- (VK/Vki_ref)*(9886*1e-6)*26980*2 #m, Proximal tubule (short-looped nephrons)
    LPT_llp <- (VK/Vki_ref)*(11077*1e-6 )*11020*2#m, Proximal tubule (long-looped nephrons)
    LPT <- LPT_slp+LPT_llp
    LTDL_slp <- (VK/Vki_ref)*(1500*1e-6)*26980*2 #m, Thin descending limb (short-looped nephrons)
    LTDL_llp <- (VK/Vki_ref)*(6200*1e-6)*11020*2 #m, Thin descending limb (long-looped nephrons)
    LThinAL <-(VK/Vki_ref)* (4700*1e-6)*11020*2 #m, Thin ascending limb (long-looped nephrons)
    LThickAL_cort <- (VK/Vki_ref)*(1450*1e-6)*38000*2 #m, Cortical thick ascending limb
    LThickAL_med <- (VK/Vki_ref)*(2100*1e-6)*38000*2 #m, Medullary thick ascending limb 
    LDT_sup <- (VK/Vki_ref)*(1452*1e-6)*26980*2 #m, Distal tubule* (superficial)
    LDT_deep <- (VK/Vki_ref)*(1650*1e-6)*11020*2 #m, Distal tubule* (deep) 
    LDT <- LDT_sup + LDT_deep
    Lduct_cort <- (VK/Vki_ref)*(2900*1e-6)*6000*2#m, Cortical collecting duct 
    Lduct_med <-  (VK/Vki_ref)*(2100*1e-6)*6000*2 #m, Outer medullary collecting duct 
    
    # Renal tubule radius per compartment
    RPT_slp <- (22.9/2)*1e-6 #m, Proximal tubule (short-looped nephrons)
    RPT_llp <-(23.1/2)*1e-6#m, Proximal tubule (long-looped nephrons)
    RTDL_slp <- (15/2)*1e-6 #m, Thin descending limb (short-looped nephrons)
    RTDL_llp <- (15/2)*1e-6 #m, Thin descending limb (long-looped nephrons)
    RThinAL <- (15/2)*1e-6 #m, #Thin ascending limb (long-looped nephrons)
    RThickAL_cort <- (25.4/2)*1e-6 #m, Cortical thick ascending limb
    RThickAL_med <- (29/2)*1e-6 #m, Medullary thick ascending limb 
    RDT_sup <- (39/2)*1e-6 #m, Distal tubule* (superficial)
    RDT_deep <- (43/2)*1e-6 #m, Distal tubule* (deep) 
    Rduct_cort <- (24/2)*1e-6 #m, Cortical collecting duct 
    Rduct_med <-  (24/2)*1e-6 #m, Outer medullary collecting duct 
    
    #volumes of filtrate compartments,[m^3]*1e3 --> L
    VPT <- 2*pi*((RPT_slp^2)*LPT_slp + (RPT_llp^2)*LPT_llp)*1e3 # Proximal tubule (short- and long-looped)
    VTDL <- 2*pi*( (RTDL_slp^2)*LTDL_slp+(RTDL_llp^2)*LTDL_llp)*1e3 #Thin descending limb (short- and long-looped)
    VThinAL <-  2*pi*((RThinAL^2)*LThinAL)*1e3 #Thin ascending limb (long-looped)
    VThickAL <- 2*pi*((RThickAL_cort^2)*LThickAL_cort+(RThickAL_med^2)*LThickAL_med)*1e3 #Thick ascending limb (Cortical and Medullary)
    VDAL <- VTDL + VThinAL+ VThickAL # Loop of Henle
    VDT <- 2*pi*((RDT_sup^2)*LDT_sup+(RDT_deep^2)*LDT_deep)*1e3#Distal tubule (superficial+deep)
    VCD <- 2*pi*((Rduct_cort^2)*Lduct_cort+(Rduct_med^2)*Lduct_med)*1e3 #collecting duct (Cortical+Outer)
    VFil <-  VPT+VDAL+VDT+VCD #L
    # We hypothesize that when kidney is weighted the renal tubule content remains inside 
    VKT <- VK - VKF - VFil#kidney tissue volume kg=L
    
    #Wang et al. (2024) [https://doi.org/10.1038/s41598-024-53270-2] state that "More than 80% of 
    # renal cortical cells are tubular epithelial cells". By assuming that this percentage holds for
    # the medulla region, we have that the total fraction of tubular cells is 0.8
    f_tubular <- 0.8
    # Clark et al. (2019) [https://doi.org/10.1016/j.kint.2018.11.028] 
    # state that "Proximal tubule cells account for roughly 52% of the estimated
    # 206 million tubule epithelial cells per kidney. However, they account for approximately 69% of 
    # total tubule protein mass by virtue of their large size compared with other renal tubule cells".
    # Thus, here we assume that the volume is analogous to the protein mass:
    f_PTC_prot_to_tub_prot <- 0.6939
    f_DALC_prot_to_tub_prot <- (0.0339 + 0.1156)
    f_DTC_prot_to_tub_prot <- 0.08
    f_CDC_prot_to_tub_prot <- 0.0766
    
    VPTC <- f_tubular*f_PTC_prot_to_tub_prot*VKT #for comparison, for 300g male rat we have 0.38mL and Worley and Fisher have 0.34 mL
    VDALC <- f_tubular*f_DALC_prot_to_tub_prot*VKT #for comparison, for 300g male rat we have 0.38mL and Worley and Fisher have 0.34 mL
    VDTC <- f_tubular*f_DTC_prot_to_tub_prot*VKT #for comparison, for 300g male rat we have 0.38mL and Worley and Fisher have 0.34 mL
    VCDC <- f_tubular*f_CDC_prot_to_tub_prot*VKT #for comparison, for 300g male rat we have 0.38mL and Worley and Fisher have 0.34 mL
    VKTrest <- (1-f_tubular)*VKT
    
    VBladder <- 0.001 #https://doi.org/10.1152/physrev.00038.2003 (CHECK)
    
    #Liver
    PVL <- 3.66e-2 #Brown et al. 1997
    VL <- PVL * BW #liver volume kg=L
    PVLB <- 0.21  #Brown et al. 1997
    VLB <- PVLB * PVL * BW #liver blood volume kg=L
    PVLF <- 0.16  #pkSim
    VLF <- PVLF * PVL* BW #liver interstitial fluid volume kg=L
    VLT <- VL - VLF #liver tissue volume kg=L
    PVLbile <- 47e-3/200 #mL/g BW,  https://doi.org/10.1016/S0002-9440(10)64679-2
    VLbile <- PVLbile * BW #L
    
    #Intestine (small and large)
    PVIN <- 2.24e-2 #Brown et al. 1997, p 416, Table 5: 1.4+0.84
    VIN <- PVIN * BW #intestine volume kg=L
    # In the following we add the plasma and blood cell volumes of the small and large intestine from Shah and betts, (2012)
    PVINB <- 0.02
    VINB <- PVINB * PVIN * BW #intestine  blood volume kg=L
    PVINF <- (0.867+0.5)/280 #mL/g BW 
    VINF <- PVINF * BW #intestine interstitial fluid volume kg=L
    VINT <- VIN - VINF #intestine tissue volume kg=L
    
    #Stomach
    PVST <- 0.46e-2 #Brown et al. 1997, p 416, Table 5
    VST <- PVST * BW #stomach volume kg=L
    PVSTB <- 0.032 #from pkSim
    VSTB <- PVSTB * PVST * BW 
    PVSTF <-  0.10 # from pkSim
    VSTF <- PVSTF * PVST * BW 
    VSTT <- VST - VSTF #stomach tissue volume kg=L
    
    #Stomach and intestine lumen
    PVSTL <- 3.4/175 #mL/g BW, Connell et al., 2008, https://doi.org/10.1211/jpp.60.1.0008
    VSTL <- PVSTL * BW #stomach lumen volume kg=L
    PVINL <- (0.894+0.792+0.678+0.598+0.442)/230 # mL/g BW, Funai et al., 2023 https://doi.org/10.1038/s41598-023-44742-y --> Figure 3C
    VINL <- PVINL * BW #intestine lumen volume kg=L
    
    #Muscle
    PVM <- 40.43e-2 #Brown et al. 1997
    VM <- PVM * BW #muscle volume kg=L
    PVMB <- 0.03 #pkSim
    VMB <- PVMB * PVM * BW #muscle blood volume kg=L
    PVMF <- 0.12 #pkSim
    VMF <- PVMF * PVM * BW #muscle interstitial fluid volume kg=L
    VMT <- VM - VMF #muscle tissue volume kg=L
    
    #Adipose
    PVA <- 7e-2 #Brown et al. 1997
    VA <- PVA * BW #adipose volume kg=L
    PVAB <- 0.01 #pkSim
    VAB <- PVAB * PVA * BW #% adipose blood volume kg=L
    PVAF <- 0.14 #pkSim
    VAF <- PVAF * PVA * BW #adipose interstitial fluid volume kg=L
    VAT <- VA - VAF #adipose tissue volume kg=L
    
    #Upper airways
    PVUA <- 257e-3/288 #mm^3/g, for a 16-wk old male 288 g, Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    VUA <- PVUA * BW  #total volume of nasal cavity
    
    #Lung
    PVLu <- 0.48e-2 #Brown et al. 1997, p 418-19 mean of values for male rat
    VLu <- PVLu * BW
    PVLuB <- 0.63 #pkSim
    VLuB <- PVLuB * PVLu * BW #Brown et al. 1997, p 459 --> capillary blood occupied 9% of the lung volume
    PVLuF <- 0.263/280 #0.263 ml, Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VLuF <- PVLuF * BW #lung interstitial fluid volume
    PVLuELF <- 0.4/275 #0.4 mL Leslie et al, 1989 https://doi.org/10.1164/ajrccm/139.2.360 --> Watkins & Rannels 1979 https://doi.org/10.1152/jappl.1979.47.2.325  
    VLuELF <- PVLuELF * BW #lung alveolar lining fluid volume kg=LL
    VLuT <- VLu - VLuF - VLuELF #lung tissue volume kg=L
    
    #Spleen
    PVSP <- 0.2e-2  #Brown et al. 1997, p 416, Table 5
    VSP <- PVSP * BW
    PVSPB <- 0.28 #pkSim
    VSPB <- PVSPB * PVSP * BW #volume of the blood of spleen kg=L
    PVSPF <- 0.554/280 #ml/g BW-> Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VSPF <- PVSPF * BW #spleen interstitial fluid volume kg=L
    VSPT <- VSP - VSPF #spleen tissue volume kg=L
    
    #Heart
    PVH <- 0.33e-2  #Brown et al. 1997, p 416, Table 5
    VH <- PVH * BW
    PVHB <- 0.26 #pkSim
    VHB <- PVHB * PVH * BW #volume of the blood of heart kg=L
    PVHF <- 0.1 #pkSim
    VHF <- PVHF * PVH * BW #heart interstitial fluid volume kg=L
    VHT <- VH - VHF #heart tissue volume kg=L
    
    #Brain
    PVBr <- 0.57e-2  #Brown et al. 1997, p 416, Table 5
    VBr <- PVBr * BW
    PVBrB <- 0.04 #pkSim
    VBrB <- PVBrB * PVBr * BW #volume of the blood of brain kg=L
    PVBrF <- 17.5/100 * PVBr
    VBrF <- PVBrF * BW #https://doi.org/10.1016/j.pneurobio.2015.12.007 --> The IS occupies 15% to 20% of the total brain volume, brain IF volume kg=L 
    VBrT <- VBr - VBrF #brain tissue volume kg=L
    
    #gonads
    PVGo <- 0.25e-2/0.230 #pKsim, L/kg
    VGo <- PVGo * BW
    PVGoB <- 0.14 #pKsim
    VGoB <-PVGoB * PVGo * BW #volume of the blood of gonads kg=L
    PVGoF <- 0.07 #pKsim
    VGoF <- PVGoF * PVGo * BW #gonads interstitial fluid volume kg=L
    VGoT <- VGo - VGoF #gonads tissue volume kg=L
    
    #Skin
    PVSK <- 19.03e-2 #Brown et al. 1997, p 416, Table 5
    VSK <- PVSK * BW
    PVSKB <- 0.02 #pkSim
    VSKB <-PVSKB * PVSK * BW #volume of the blood of skin kg=L
    PVSKF <- 0.4  #https://doi.org/10.1111/j.1748-1716.1981.tb06901.x 40 mL/100 g tissue, BW = 200-250 g
    VSKF <- PVSKF * PVSK * BW #skin interstitial fluid volume kg=L
    VSKT <- VSK - VSKF #skin tissue volume kg=L
    
    #Bones
    PVBo <-1.59e-2/0.230 #pkSim
    VBo <- PVBo * BW
    PVBoB <- 0.04 #pkSim
    VBoB <-PVBoB * PVBo * BW #volume of the blood of bones kg=L
    PVBoF <- 0.1 #pkSim
    VBoF <- PVBoF * PVBo * BW #bones interstitial fluid volume kg=L
    VBoT <- VBo - VBoF #bones tissue volume kg=L
    
    #RoB
    PVR <- 1 - PVB - PVK - PVL - PVM - PVA - PVSP - PVH - PVBr - PVGo - PVST - PVIN - PVSK - PVBo
    VR <- PVR * BW #volume of the rest of the body kg=LL
    PVRB <-(PVKB+PVLB+PVLuB+PVMB+PVAB+PVSPB+PVHB+PVBrB+PVGoB+PVSTB+PVINB+PVSKB+PVBoB)/13 #average VF of all the included organs (kg=L)
    VRB <- PVRB * PVR * BW #volume of the blood of RoB kg=L
    PVRF <-(PVKF+PVLF+PVLuF+PVMF+PVAF+PVSPF+PVHF+PVBrF+PVGoF+PVSTF+PVINF+PVSKF+PVBoF)/13 #average VF of all the included organs (kg=L)
    VRF <- PVRF * PVR * BW #RoB of the blood of rest of body kg=L
    VRT <- VR - VRF #tissue volume of the rest of body kg=L
    
    #respiratory surface areas
    Nasal_SA <-  19*1e-4 #m2 Oller and Oberdörster (2010) [https://doi.org/10.1016/j.yrtph.2010.02.006]
    Alveolar_SA <- 0.40 #m2 Oller and Oberdörster (2010) [https://doi.org/10.1016/j.yrtph.2010.02.006]
    
    
    ##Capillary surface area for each tissue (Ai) as percentage of body weight (m^2/kg),
    #values from pkSim "Endothelial Surface area", Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    
    #Peritubular capillary density in From Gazzard et al. (2024) [https://doi.org/10.1002/ar.25576] 
    d_peritubular <- 0.024 #um^2/um^3,
    VK_gazzard <- 0.00363#L
    Vcortex <- 2*1295* VK/VK_gazzard #mm^3, NEED TO RECHECK IF IT REFERS TO ONE OF BOTH KIDNEYS
    A_peritubular <- (d_peritubular*Vcortex*1e3/1e6) #m^2
    A_peritubular_PTC <- A_peritubular * VPTC/ (VPTC + VDALC)
    A_peritubular_DTC <- A_peritubular * VDALC/ (VPTC + VDALC)
    
    
    BW_PKSim <- 0.23 #kg
    SA_scaling_factor <- BW/BW_PKSim  
    
    PAL <- 1136e-4 #m^2
    AL <- PAL *  SA_scaling_factor #liver surface area (m^2)
    
    PAST <- 33.77e-4 #m^2
    AST <- PAST * SA_scaling_factor #stomach surface area (m^2)
    PASTL<- 33.77e-4 #m^2
    ASTL<- PASTL * SA_scaling_factor #stomach lumen surface area (m^2)
    
    PAIN <- (74.82e-4+20.4e-4) #m^2 #large and small intestine
    AIN <- PAIN * SA_scaling_factor #intestine surface area (m^2)
    
    # Gut surface areas from PkSim
    Duodenum <- 5.50 #cm^2
    Upper_jejunum <- 27.33 #cm^2
    Lower_jejunum <- 27.33 #cm^2
    Upper_ileum <- 0.82 #cm^2
    Lower_ileum <- 0.59 #cm^2
    Cecum <- 1.63 #cm^2
    colon_ascendens <- 6.36 #cm^2
    colon_transversum <- 3.22 #cm^2
    colon_descendens <-2.27 #cm^2
    colon_sigmoid <-2.14 #cm^2
    
    SA_pksim <- Duodenum+Upper_jejunum+Lower_jejunum+Upper_ileum+Lower_ileum+Cecum+     
      colon_ascendens+colon_transversum+colon_descendens+colon_sigmoid
    
    #Enhancement factors
    EF_Duodenum <- 13.86
    EF_Upper_jejunum <-10.68
    EF_Lower_jejunum <- 8.31
    EF_Upper_ileum <- 8.11
    EF_Lower_ileum <- 11.22
    EF_Cecum <- 1.47
    EF_colon_ascendens <- 1.7 
    EF_colon_transversum <- 1.92
    EF_colon_descendens <- 1.95
    EF_colon_sigmoid <- 1.95
    
    
    
    SA_pksim_effective <- (Duodenum*EF_Duodenum+Upper_jejunum*EF_Upper_jejunum+Lower_jejunum*EF_Lower_jejunum+
                             Upper_ileum*EF_Upper_ileum+Lower_ileum*EF_Lower_ileum+Cecum*EF_Cecum+
                             colon_ascendens*EF_colon_ascendens+colon_transversum*EF_colon_transversum+
                             colon_descendens*EF_colon_descendens+colon_sigmoid*EF_colon_sigmoid)
    
    
    PAINL <- SA_pksim_effective * 1e-4 #m^2, from reference body weight
    AINL <- PAINL*SA_scaling_factor #m^2
    
    PAM <- 3042e-4 #m2
    AM <- PAM * SA_scaling_factor #muscle surface area (m^2)
    PAA <- 95.93e-4/0.23 #m2
    AA <- PAA * SA_scaling_factor #adipose surface area (m^2)
    
    PAR <- 100e-4#m^2, assumption
    AR <- PAR * SA_scaling_factor#surface area of rest of body (m^2)
    
    PLu <- 600.5e-4 #m^2 
    ALu <- PLu * SA_scaling_factor #lung surface area (m^2)
    PSP <- 162.3e-4 #m^2
    ASP <- PSP * SA_scaling_factor #spleen surface area (m^2)
    PH <-  201.1e-4 #m^2 
    AH <- PH * SA_scaling_factor #heart surface area (m^2)
    PBr <- 60.34e-4 #m^2
    ABr <- PBr * SA_scaling_factor #brain surface area (m^2)
    PGo <- 335.8e-4#m^2
    AGo <- PGo * SA_scaling_factor #gonads surface area (m^2)
    PSK <- 729.1e-4#m^2
    ASK <- PSK * SA_scaling_factor #skin surface area (m^2)
    PBo <- 621.4e-4#m^2
    ABo <- PBo * SA_scaling_factor #skin surface area (m^2)
    
    ####################################
    #----------------------------------#
    #             Flow Rates           #
    #----------------------------------#
    ####################################
    
    
    ####################################
    #            Blood flow rates      #
    ####################################
    
    #(QBi, in L/h) to different tissues (i=L, K, G, A, M, R)
    #as a percentage of cardiac output (Qcardiac L/h), which itself is a function
    #of body weight (BW)
    
    Qcardiac <- 0.235 * (BW^0.75) *60 #L/min->*60-> L/h
    PQBK <- 14.1/100 #Brown et al. 1997, p 438, Table 23
    QBK <- PQBK * Qcardiac #L/h
    PQBL <- 2.1/100 
    QBL <- PQBL * Qcardiac #L/h Brown et al. 1997, p 438, Table 23
    PQBST <- 0.16/100 #https://doi.org/10.3390/pharmaceutics6010097 Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBST = (8/1000)/Qcard *100
    QBST <- PQBST * Qcardiac #L/h
    PQBIN <- 9/100 #https://doi.org/10.3390/pharmaceutics6010097 Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBIN = (451/1000)/Qcard *100
    QBIN <- PQBIN * Qcardiac #L/h
    PQBM <- 27.8/100 #Brown et al. 1997, p 438, Table 23
    QBM <- PQBM * Qcardiac #L/h
    PQBA <- 7/100 #Brown et al. 1997, p 438, Table 23
    QBA <- PQBA * Qcardiac #L/h
    PQBLu <- 1 
    QBLu <- PQBLu * Qcardiac #L/h
    PQBSP <- 0.75/100 #https://doi.org/10.3390/pharmaceutics6010097 Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBSP = (37.5/1000)/Qcard *100
    QBSP <- PQBSP * Qcardiac #L/h
    PQBH <- 5.1/100 #Brown et al. 1997, p 438, Table 23
    QBH <- PQBH * Qcardiac #L/h
    PQBBr <- 2.0/100 #Brown et al. 1997, p 438, Table 23
    QBBr <- PQBBr * Qcardiac #L/h 
    PQBGo <- 0.28/100 #https://doi.org/10.1152/ajpregu.1987.253.2.R228 Qcard=0.235*(0.335^0.75)*60 (L/h) and PQBT = (0.295*60/1000)/Qcard *100
    QBGo <- PQBGo * Qcardiac #L/h
    PQBSK <- 1.35/100 #https://doi.org/10.1111/1523-1747.ep12277181 Qcard=0.235*(0.2^0.75)*60 (L/h) and PQBSK = (0.95*60/1000)/Qcard *100
    QBSK <- PQBSK * Qcardiac #L/h
    PBBo <- 12.2/100 #Brown et al. 1997, p 438, Table 23
    QBBo <- PBBo * Qcardiac #L/h
    
    # Total blood outflow from liver
    QBLtot <- QBL+QBSP+QBIN+QBST
    
    PQBR = 1 - PQBK - PQBL - PQBST - PQBIN - PQBM - PQBA - PQBH - PQBSK - PQBSP - PQBGo - PQBBr - PBBo
    QBR <- PQBR * Qcardiac #L/h
    
    
    ##################################
    #     Other fluids flow rates    #
    ##################################
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration rate (GFR), in L/h
    
    #PQbile <- 90/1000/24 # [mL/d/kg BW]/1000/24 --> L/h/kg BW
    #Qbile <- PQbile * BW #L/h
    if (sex == "M"){
      PQbile = 0.206/2 #L/kg liver/h #source: https://doi.org/10.1038/s41598-019-46150-7
      Qbile = PQbile* VL  #L/h
    }else if (sex == "F"){
      # Females have 44% more bile flow, source: doi:10.1042/cs0550253
      PQbile = 0.206/2 #L/kg liver/h #
      Qbile = 1.44* PQbile* VL  #L/h
    }
    Qfeces <- (8.18/0.21)*BW/24/1000 #g/kg BW, based on Cui et al.(2010)
    feces_density <- 1.29 #g/cm^3 --> g/mL from Lupton 1986, Fig 1. Fiber free control diet, https://doi.org/10.1093/jn/116.1.164
    
    if (sex == "M"){
      PQGFR <- 62.1  #L/h/kg   Corley et al., 2005 https://doi.org/10.1093/toxsci/kfi119, GFRC --> scaled fraction of kidney weight
      QGFR <- PQGFR * VK #L/h
      Qurine <- (60/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, https://doi.org/10.1002/nau.1006
    }else if(sex == "F"){
      PQGFR <- 41.04  #L/h/kg  Corley et al., 2005 https://doi.org/10.1093/toxsci/kfi119, GFRC --> scaled fraction of kidney weight
      QGFR <- PQGFR * VK #L/h
      Qurine <- (85/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, https://doi.org/10.1002/nau.1006  
    }
    QGE<- 1.41*BW^(-0.25) #gastric emptying time (1/(h*BW^0.25)); from Yang, 2013 [https://doi.org/10.1016/j.taap.2013.03.022]
    
    #flows of filtrate compartments from Gilmer et al. (2018) 
    # [https://doi.org/10.1152/ajprenal.00219.2018], nL/min ---> L/h 
    QPT_ref <- (40.9*26980+39.7*11020)*2*1e-9*60 #L/h, proximal tubule flow (short +long)
    QTDL_ref <- (26.4*26980+26.3*11020)*2*1e-9*60 #L/h, Thin descending limb (short +long)
    QTAL_ref <- 8 * 38000*2*1e-9*60 #L/h, Thin+thick (short +long) ascending 
    QDT_ref <- (6*26980 + 6*11020)*2*1e-9*60  #  Distal tubule (short + long)
    QCD_ref <- 14.4* 6000*2*1e-9*60  #  collecting duct (Cortical + medullary)
    
    # The values from Gilmer et al. (2018) are BW invariant. So since the GFR
    # we use is scaled to BW and GFR = QPT, we scale these values based on the
    # current GFR.
    Q_scaling_factor = QGFR/QPT_ref
    
    QPT <- Q_scaling_factor * QPT_ref #L/h, proximal tubule flow (short +long)
    QDAL <- Q_scaling_factor * QTDL_ref #L/h, Thin descending limb (short +long)
    QTAL <- Q_scaling_factor * QTAL_ref #L/h, Thin+thick (short +long) ascending 
    QDT <- Q_scaling_factor * QDT_ref  #L/h, Distal tubule (short + long)
    QCD <- Q_scaling_factor * QCD_ref  #L/h, collecting duct (Cortical + medullary)
    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    #Surface areas Interstitial - Intracellular (m^2), from PKSim 
    BW_ref <- 0.23
    AcK_total= 437.16*BW/BW_ref
    AcK_PTC <- AcK_total*VPTC/VFil # surface area of proximal tubule cells
    AcK_DALC <- AcK_total*VDALC/VFil # surface area of decending/ascending limb cells (loop of Henle)
    AcK_DTC <- AcK_total*VDTC/VFil # surface area of distal tubule cells
    AcK_CDC <- AcK_total*VCDC/VFil # surface area of collecting duct cells 
    AcKTrest <- AcK_total* VKTrest/VFil
    AcL= 84.45*BW/BW_ref
    AcST= 1007.31*BW/BW_ref
    AcIN= (400.94+152.39) *BW/BW_ref # small+large intestine
    AcM= 8.2*BW/BW_ref
    AcA= 3.87*BW/BW_ref
    AcLu= 0.05*BW/BW_ref
    AcSP= 564.05*BW/BW_ref
    AcH= 5.60*BW/BW_ref
    AcBr= 6.12e-4*BW/BW_ref
    AcGo= 2.01*BW/BW_ref
    AcSK= 0.11*BW/BW_ref
    AcBo= 6.52*BW/BW_ref
    # We don't have data for the surface area of IS-IC for the rest of the body, thus 
    # we naively assume an average:
    AcR= mean(c(AcK_total,AcL,AcST,AcIN,AcM,AcA,AcLu,AcSP,AcH,AcBr,AcGo,AcSK,AcBo))
    
    n <- 5 #enlargement factor of the apical membrane of tubule cells
    # Surface areas of the different subcompartments of kidney filtrate, m^2
    APT <-  2*pi*(RPT_slp*LPT_slp + RPT_llp*LPT_llp)*n # Proximal tubule (short- and long-looped)
    ATDL <- 2*pi*(RTDL_slp*LTDL_slp + RTDL_llp*LTDL_llp)*n #Thin descending limb (short- and long-looped)
    AThinAL <-  2*pi*RThinAL*LThinAL*n #Thin ascending limb (long-looped)
    AThickAL <- 2*pi*(RThickAL_cort*LThickAL_cort + RThickAL_med*LThickAL_med)*n #Thick ascending limb (Cortical and Medullary)
    ADAL <- ATDL + AThinAL + AThickAL
    ADT <- 2*pi*(RDT_sup*LDT_sup + RDT_deep*LDT_deep)*n #Distal tubule (superficial+deep)
    ACD <- 2*pi*(Rduct_cort*Lduct_cort + Rduct_med*Lduct_med)*n #collecting duct (Cortical+Outer)
    AFil <- APT + ADAL + ADT + ACD
    
    #Alveolar cells surface area (Type I and II), m^2
    AcELF = ((78.8*2*5320*1e-6) + (125*2*123*1e-6))*BW/0.29  #Stone et al., 1992, BW_ref = 0.29, values for each lung , https://doi.org/10.1165/ajrcmb/6.2.235
    
    # #canalicular surface area, m^2
    # rat_hep_surf_area = 22.95 * 1e2 # 22.95*1e6 cm2 --> m2,  https://doi.org/10.1074/jbc.271.12.6702
    # AcLBilec = 0.01 * 22.95 * 1e2 # m2 , canalicular membrane 1% of the surface area of the hepatocyte,https://www.ncbi.nlm.nih.gov/books/NBK470209/
    
    #Stomach
    # For identifiability reasons we assume that absorption takes place only through the intestines
    kabST <- (kabs_st* ASTL)*1000 #L/h
    
    MW = 414.07 #g/mol, PFOA molecular weight
    
    return(list(
      
      'VB'=VB, 'Vplasma'=Vplasma, 'VK'=VK, 'VKB'=VKB, 
      'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,'VBladder' = VBladder, 
      'VPT' = VPT, 'VDAL' = VDAL, 'VDT' = VDT, 'VCD' = VCD,
      'VPTC' = VPTC, 'VDALC' = VDALC, 'VDTC' = VDTC, 'VCDC' = VCDC,
      'VKTrest' = VKTrest,
      
      'VL'=VL, 'VLB'=VLB, 'VLF'=VLF, 'VLT'=VLT, 'VLbile'=VLbile,
      'VM'=VM, 'VMB'=VMB, 'VMF'=VMF, 'VMT'=VMT, 'VA'=VA, 'VAB'=VAB, 
      'VAF'=VAF, 'VAT'=VAT, 'VR'=VR, 'VRB'=VRB, 
      'VRF'=VRF, 'VRT'=VRT, 'VVen' = VVen,
      'VArt' = VArt, 'VLu'=VLu, 'VLuB'=VLuB, 'VLuF'=VLuF,
      'VLuELF'=VLuELF, 'VLuT'=VLuT,'VUA'=VUA, 
      'VSP'=VSP, 'VSPB'=VSPB, 'VSPF'=VSPF, 'VSPT'=VSPT,
      'VH'=VH, 'VHB'=VHB, 'VHF'=VHF, 'VHT'=VHT,
      'VBr'=VBr, 'VBrB'=VBrB, 'VBrF'=VBrF, 'VBrT'=VBrT,
      'VGo'=VGo, 'VGoB'=VGoB, 'VGoF'=VGoF, 'VGoT'=VGoT,
      'VIN'=VIN, 'VINB'=VINB, 'VINF'=VINF, 'VINT'=VINT,
      'VST'=VST, 'VSTB'=VSTB, 'VSTF'=VSTF, 'VSTT'=VSTT,
      'VSTL'=VSTL, 'VINL'=VINL,
      'VSK'=VSK,'VSKB'=VSKB, 'VSKF'=VSKF, 'VSKT'=VSKT,
      'VBo'=VBo,'VBoB'=VBoB, 'VBoF'=VBoF, 'VBoT'=VBoT,
      
      "Alveolar_SA" = Alveolar_SA, "Nasal_SA" = Nasal_SA,
      
      'A_peritubular_PTC' = A_peritubular_PTC, 
      'A_peritubular_DTC' = A_peritubular_DTC,
      'AL'=AL, 'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
      'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
      'AIN'=AIN, 'AGo'=AGo,
      'ASK'= ASK, 'ABo'=ABo,
      
      'AINL' = AINL, 'AcL' = AcL, 'AcM' = AcM, 'AcST' = AcST, 
      'AcIN' = AcIN, 'AcA' = AcA, 'AcLu' = AcLu, 'AcELF' = AcELF, 
      'AcSP' = AcSP, 'AcH' = AcH, 'AcBr' = AcBr, 'AcGo' = AcGo, 
      'AcSK' = AcSK, 'AcBo' = AcBo, 'AcR' = AcR, 'APT' = APT, 
      'ADAL' = ADAL, 'ADT' = ADT, 'ACD' = ACD, 'AcK_DTC' = AcK_DTC,
      'AcK_PTC' = AcK_PTC, 'AcK_DALC' = AcK_DALC,  
      'AcK_CDC' = AcK_CDC, 'AcKTrest' = AcKTrest,
      
      'Qcardiac'=Qcardiac, 'QBK'=QBK, 
      'QBL'=QBL, 'QBLtot'=QBLtot,
      'QBM'=QBM, 'QBA'=QBA,
      'QBR'=QBR, 'QBLu'=QBLu, 'Qfeces'=Qfeces,  'feces_density'=feces_density,
      'Qbile'=Qbile, 'QGFR'=QGFR,'Qurine'=Qurine, 
      'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST,
      'QBIN'=QBIN, 'QGE'=QGE,
      'QBGo'=QBGo,
      'QBSK'=QBSK, 'QBBo'=QBBo, 'Hct' = Hct,
      
      
      "QPT" = QPT, "QDAL" = QDAL, "QTAL" = QTAL, "QDT" = QDT, "QCD" = QCD,
      
      'f_tubular' =  f_tubular,  'f_PTC_prot_to_tub_prot' = f_PTC_prot_to_tub_prot, 
      'f_DALC_prot_to_tub_prot' = f_DALC_prot_to_tub_prot, 
      'f_DTC_prot_to_tub_prot' = f_DTC_prot_to_tub_prot , 
      'f_CDC_prot_to_tub_prot' = f_CDC_prot_to_tub_prot,
      'kabST'=kabST, 
      
      
      
      "admin.time" = admin.time, "admin.dose" = admin.dose,
      "admin.type" = admin.type, "MW"=MW, "sex" = sex, "BW" = BW,
      "depfr_head" = depfr_head, "depfr_AF" = depfr_AF
    ))
    
  })
}   


estimate_BFn_TVn <- function(sex, BW){
  
  #Nose-only Breathing frequency
  if (sex == "M"){
    PBFn <- 168*60/0.219  #1/min -->  1/h/kg, nose-only inhalation, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    BFn <- PBFn * BW  # 1/h
  }else if(sex == "F"){
    PBFn <- 177*60/0.145  #1/min -->  1/h/kg, nose-only inhalation, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    BFn <- PBFn * BW  # 1/h
  }
  
  #Nose-only Tidal Volume
  if (sex == "M"){
    PTVn <- 1.6/219   #mL/g --> mL=g, unitless,nose-only inhalation, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    TVn <- PTVn * BW  #L
  }else if(sex == "F"){
    PTVn <- 1.2/145   #mL/g --> mL=g, unitless,nose-only inhalation, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    TVn <- PTVn * BW  #L
  }
  return(c("BFn"=BFn, "TVn"=TVn))
}


estimate_BFi_TVi <- function(sex, BW){
  
  #Inhalation Breathing frequency
  if (sex == "M"){
    PBFi <- 86*60/0.387  #1/min -->  1/h/kg, inhalation, Moss et al., 2006, https://doi.org/10.1203/01.pdr.0000203104.45807.23
    #(95 in the study of Walker et al., 1997), https://doi.org/10.1016/S0034-5687(96)02520-0
    BFi <- PBFi * BW  # 1/h
  }else if(sex == "F"){
    PBFi <- 79*60/0.223  #1/min -->  1/h/kg, inhalation, Mauderly et al., 1986, https://doi.org/10.1002/jat.2550060106
    BFi <- PBFi * BW  # 1/h
  }
  
  #Inhalation Tidal Volume
  if (sex == "M"){
    PTVi <- 7.2/1e3   #mL/kg --> mL=g, unitless,inhalation, Walker et al., 1997, https://doi.org/10.1016/S0034-5687(96)02520-0
    TVi <- PTVi * BW  #L
  }else if(sex == "F"){
    PTVi <- 1.61/230   #mL/g --> mL=g, unitless,inhalation, Deng et al., 2023,  https://doi.org/10.1111/os.13630
    TVi <- PTVi * BW  #L
  }
  return(c("BFi"=BFi, "TVi"=TVi))
}


ode.func <- function(time, inits, params){
  t_24 <- floor(time/24)
  BW_init <- unlist(params["BW"])
  if (params["sex"] == "M"){
    growth_rate = 5.9/1000 #g/d
    BW_new <- BW_init + growth_rate*t_24
    if(BW_new>0.597){
      BW_new<-0.597
    }
  }else{
    growth_rate = 3.5/1000 #g/d
    BW_new <- BW_init + growth_rate*t_24
    if(BW_new>0.365){
      BW_new<-0.365
    }
  }
  BW_scaled_params <- c ('VB', 'Vplasma', 'VK', 'VKB', 'VKF',  'VFil','VPT' , 'VDAL' , 'VDT' , 'VCD' ,'VPTC' , 'VDALC' , 'VDTC' , 'VCDC' ,'VL', 'VLB', 'VLF',  'VLbile',
                         'VM', 'VMB', 'VMF', 'VA', 'VAB', 'VAF', 'VAT', 'VR', 'VRB',  'VRF',  'VVen' ,'VArt', 'VLu', 'VLuB', 'VLuF',
                         'VLuELF','VSP', 'VSPB', 'VSPF','VH', 'VHB', 'VHF','VBr', 'VBrB', 'VBrF','VGo', 'VGoB', 'VGoF',
                         'VIN', 'VINB', 'VINF', 'VST', 'VSTB', 'VSTF','VSTL', 'VINL','VSK','VSKB', 'VSKF',
                         'VBo','VBoB', 'VBoF','VLT', 'VKTrest','VINT', 'VSTT','VMT', 'VAT', 
                         'VLuT', 'VSPT','VHT', 'VBrT','VGoT', 'VSKT','VBoT', 'VRT',
                         'VKT', 'A_peritubular_PTC', 'A_peritubular_DTC','AL', 'AM', 'AA', 'AR', 'ALu','ASP', 'AH', 'ABr', 'AST',
                         'AIN', 'AGo','ASK', 'ABo','AINL', 'AcL' , 'AcM' , 'AcST' ,'AcIN', 'AcA' , 'AcLu' , 'AcELF' , 
                         'AcSP', 'AcH' , 'AcBr' , 'AcGo','AcSK', 'AcBo' , 'AcR' , 'APT' , 'ADAL' , 'ADT', 'ACD' , 'AcK_DALC','AcK_CDC' , 'AcKTrest',
                         'Qfeces','Qbile', 'QGFR','Qurine','kabST',"QPT" , "QDAL", "QTAL" , "QDT", "QCD", 'CLfeces', "CL_hepatobiliary", 
                         'VmL_Oatp', 'VmL_Ntcp','VmL_Oatp2',  'VmIn_Oatp2', 'VmK_Oatp','VmLu_Oatp_ap', 'VmLu_Oatp_bas',
                         'VmK_Oat1', 'VmK_Oat3','VmK_Urat','kKTrestF', 'kCdcF' , 'kDalcF' , 'kPtcF' , 'kDtcF' ,
                         'kPtcTu', 'kDalcTu' , 'kDtcTu' , 'kCdcTu' , 'kLFLT',  'kAFAT', 'kRFRT','kMFMT', 'kLuTLuF', 'kLuTLuELF', 'kSPFSPT' ,
                         'kSTFSTT' , 'kINFINT' , 'kHFHT' ,'kBrFBrT' , 'kGoFGoT' ,'kSKFSKT' , 'kBoFBoT', "k_gut_in", 'k_gut_out')
  Q_scaled_params <- c('QBK', 'QBL', 'QBLtot','QBM', 'QBA',
                       'QBR', 'QBLu','QBSP', 'QBH', 'QBBr', 'QBST','QBIN', 'QGE','QBGo','QBSK', 'QBBo', "PparaKi" ,"PparaLi" ,"PparaSt" ,"PparaIn" ,
                       "PparaMu" ,"PparaAd" ,"PparaRe" ,"PparaLu" ,"PparaSp","PparaHt" ,"PparaBr" ,"PparaGo","PparaSk" ,"PparaBo")
  Qcardiac_init <-  unlist(params["Qcardiac"])
  params["Qcardiac"] <- unlist(params["Qcardiac"])* (BW_new/BW_init)^0.75
  params["QGE"] <-   unlist(params["QGE"])* (BW_new/BW_init)^(-0.25) 
  params["CLfeces"] <-   unlist(params["CLfeces"])* (BW_new/BW_init)^(-0.25) 
  
  params[BW_scaled_params] <- lapply(params[BW_scaled_params], function(x) x * unname(BW_new/BW_init))
  params[Q_scaled_params] <- lapply(params[Q_scaled_params], function(x) x * unname(unlist(params["Qcardiac"])/Qcardiac_init))
  
  with(as.list(c(inits, params)),{
    
    #====================PFOA mass balance at each tissue or fluid compartment==============================     
    
    # Concentrations in ug/L
    
    # Blood concentration
    MVen <-  MVenf + MVenb
    MArt <- MArtf +MArtb
    
    CVen <- MVen/VVen
    CVenb <- MVenb/VVen
    CVenf <- MVenf/VVen
    
    CArt <- MArt/VArt
    CArtf <-MArtf/VArt
    CArtb <-MArtb/VArt
    
    # Kidney 
    
    MKB <- MKBf + MKBb
    MKF <- MKFf + MKFb
    CKB <- MKB/VKB # blood concentration
    CKBf <- MKBf/VKB
    CKBb <- MKBb/VKB
    CKF <- MKF/VKF  #interstitial fluid concentration
    CKFf <- MKFf/VKF
    CKFb <- MKFb/VKF
    
    CKTrestf <- MKTrestf/VKTrest
    CKTrestb <- MKTrestb/VKTrest
    CPTCf <- MPTCf/VPTC
    CPTCb <- MPTCb/VPTC
    CDALCf <- MDALCf/VDALC
    CDALCb <- MDALCb/VDALC
    CDTCf <- MDTCf/VDTC
    CDTCb <- MDTCb/VDTC
    CCDCf <- MCDCf/VCDC
    CCDCb <- MCDCb/VCDC
    
    MKTf <- MKTrestf + MPTCf + MDALCf + MDTCf + MCDCf
    MKTb <- MKTrestb + MPTCb + MDALCb + MDTCb + MCDCb
    MKT <- MKTf + MKTb
    Mfil <- MPT+MDAL+MDT+MCD
    CKTb <- MKTf/VKT
    CKTf <- MKTb/VKT
    CKT <- MKT/VKT # tissue concentration
    #renal filtrate
    CPT <- MPT/VPT #proximal tubule
    CDAL <- MDAL/ VDAL # tubule
    CDT <- MDT/ VDT #distal tubule
    CCD <- MCD/ VCD #proximal tubule
    CBladder <- MBladder/VBladder
    
    #Liver
    MLB <- MLBf + MLBb
    MLF <- MLFf + MLFb
    MLT <- MLTf + MLTb
    CLB <- MLB/VLB # blood concentration
    CLBf <- MLBf/VLB
    CLBb <- MLBb/VLB
    CLF <- MLF/VLF  #interstitial fluid concentration
    CLFf <- MLFf/VLF
    CLFb <- MLFb/VLF
    CLT <- MLT/VLT # tissue concentration
    CLTf <- MLTf/VLT
    CLTb <- MLTb/VLT
    Cbile <- MLbile/VLbile #Bile  canaliculi
    
    #Stomach
    
    MSTB <- MSTBf + MSTBb
    MSTF <- MSTFf + MSTFb
    MSTT <- MSTTf
    CSTB <- MSTB/VSTB # blood concentration
    CSTBf <- MSTBf/VSTB
    CSTBb <- MSTBb/VSTB
    CSTF <- MSTF/VSTF  #interstitial fluid concentration
    CSTFf <- MSTFf/VSTF
    CSTFb <- MSTFb/VSTF
    CSTTf <- MSTTf/VSTT # tissue concentration
    
    #Intestine
    
    MINB <- MINBf + MINBb
    MINF <- MINFf + MINFb
    MINT <- MINTf 
    CINB <- MINB/VINB # blood concentration
    CINBf <- MINBf/VINB
    CINBb <- MINBb/VINB
    CINF <- MINF/VINF  #Interstitial fluid concentration
    CINFf <- MINFf/VINF
    CINFb <- MINFb/VINF
    CINTf <- MINTf/VINT # tissue concentration
    
    #Stomach and Intestine lumens
    CSTL = MSTL/VSTL # Stomach Lumen concentration
    CINL = MINL/VINL # Intestine Lumen concentration
    
    #Muscle
    
    MMB <- MMBf + MMBb
    MMF <- MMFf + MMFb
    MMT <- MMTf 
    CMB <- MMB/VMB # blood concentration
    CMBf <- MMBf/VMB
    CMBb <- MMBb/VMB
    CMF <- MMF/VMF  #Interstitial fluid concentration
    CMFf <- MMFf/VMF
    CMFb <- MMFb/VMF
    CMTf <-  MMTf/VMT # tissue concentration
    
    #Adipose
    
    MAB <- MABf + MABb
    MAF <- MAFf + MAFb
    MAT <- MATf 
    CAB <- MAB/VAB # blood concentration
    CABf <- MABf/VAB
    CABb <- MABb/VAB
    CAF <- MAF/VAF  #Interstitial fluid concentration
    CAFf <- MAFf/VAF
    CAFb <- MAFb/VAF
    CATf <- MATf/VAT # tissue concentration
    
    #Rest-of-the-body
    
    MRB <- MRBf + MRBb
    MRF <- MRFf + MRFb
    MRT <- MRTf 
    CRB <- MRB/VRB # blood concentration
    CRBf <- MRBf/VRB
    CRBb <- MRBb/VRB
    CRF <- MRF/VRF  #Interstitial fluid concentration
    CRFf <- MRFf/VRF
    CRFb <- MRFb/VRF
    CRTf <- MRTf/VRT # tissue concentration
    
    #Lung
    #Upper airways
    CUA = MUA/VUA # upper airways concentration
    
    MLuB <- MLuBf + MLuBb
    MLuF <- MLuFf + MLuFb
    MLuT <- MLuTf 
    MLuELF <- MLuELFf + MLuELFb + MLuELFdust
    CLuB <- MLuB/VLuB # blood concentration
    CLuBf <- MLuBf/VLuB
    CLuBb <- MLuBb/VLuB
    CLuF <- MLuF/VLuF  #Interstitial fluid concentration
    CLuFf <- MLuFf/VLuF
    CLuFb <- MLuFb/VLuF
    CLuTf <- MLuTf/VLuT #tissue concentration
    CLuELF <- MLuELF/VLuELF #alveolar lining fluid concentration
    CLuELFf <- MLuELFf/VLuELF
    CLuELFb <- MLuELFb/VLuELF
    CLuELFdust <- MLuELFdust/VLuELF
    
    #Spleen
    MSPB <- MSPBf + MSPBb
    MSPF <- MSPFf + MSPFb
    MSPT <- MSPTf 
    CSPB <- MSPB/VSPB # blood concentration
    CSPBf <- MSPBf/VSPB
    CSPBb <- MSPBb/VSPB
    CSPF <- MSPF/VSPF  #Interstitial fluid concentration
    CSPFf <- MSPFf/VSPF
    CSPFb <- MSPFb/VSPF
    CSPTf <-  MSPTf/VSPT # tissue concentration
    
    #Heart
    MHB <- MHBf + MHBb
    MHF <- MHFf + MHFb
    MHT <- MHTf 
    CHB <- MHB/VHB # blood concentration
    CHBf <- MHBf/VHB
    CHBb <- MHBb/VHB
    CHF <- MHF/VHF  #Interstitial fluid concentration
    CHFf <- MHFf/VHF
    CHFb <- MHFb/VHF
    CHTf <- MHTf/VHT # tissue concentration
    
    #Brain
    
    MBrB <- MBrBf + MBrBb
    MBrF <- MBrFf + MBrFb
    MBrT <- MBrTf 
    CBrB <- MBrB/VBrB # blood concentration
    CBrBf <- MBrBf/VBrB
    CBrBb <- MBrBb/VBrB
    CBrF <- MBrF/VBrF  #Interstitial fluid concentration
    CBrFf <- MBrFf/VBrF
    CBrFb <- MBrFb/VBrF
    CBrTf <-  MBrTf/VBrT # tissue concentration
    
    #gonads
    
    MGoB <- MGoBf + MGoBb
    MGoF <- MGoFf + MGoFb
    MGoT <- MGoTf 
    CGoB <- MGoB/VGoB # blood concentration
    CGoBf <- MGoBf/VGoB
    CGoBb <- MGoBb/VGoB
    CGoF <- MGoF/VGoF  #Interstitial fluid concentration
    CGoFf <- MGoFf/VGoF
    CGoFb <- MGoFb/VGoF
    CGoTf <-  MGoTf/VGoT # tissue concentration
    
    #Skin
    
    MSKB <- MSKBf + MSKBb
    MSKF <- MSKFf + MSKFb
    MSKT <- MSKTf 
    CSKB <- MSKB/VSKB # blood concentration
    CSKBf <- MSKBf/VSKB
    CSKBb <- MSKBb/VSKB
    CSKF <- MSKF/VSKF  #Interstitial fluid concentration
    CSKFf <- MSKFf/VSKF
    CSKFb <- MSKFb/VSKF
    CSKTf <- MSKTf/VSKT # tissue concentration
    
    #Bones
    
    MBoB <- MBoBf + MBoBb
    MBoF <- MBoFf + MBoFb
    MBoT <- MBoTf 
    CBoB <- MBoB/VBoB # blood concentration
    CBoBf <- MBoBf/VBoB
    CBoBb <- MBoBb/VBoB
    CBoF <- MBoF/VBoF  #Interstitial fluid concentration
    CBoFf <- MBoFf/VBoF
    CBoFb <- MBoFb/VBoF
    CBoTf <- MBoTf/VBoT # tissue concentration
    
    #Calculation of free and bound PFOA in venous blood
    dCalbVenf <- koff_alb*CVenb/MW/1e6 - kon_alb*CalbVenf*CVenf/MW/1e6
    
    #Calculation of free and bound PFOA in arterial blood
    dCalbArtf <- koff_alb*CArtb/MW/1e6 - kon_alb*CalbArtf*CArtf/MW/1e6
    
    #--------------------------------------------------------------
    #Calculation of free concentrations in organ blood
    #--------------------------------------------------------------
    
    #Calculation of free and bound PFOA in kidney blood
    dCalbKBf <- koff_alb*CKBb/MW/1e6 - kon_alb*CalbKBf*CKBf/MW/1e6
    
    #Calculation of free and bound PFOA in liver blood
    dCalbLBf <- koff_alb*CLBb/MW/1e6 - kon_alb*CalbLBf*CLBf/MW/1e6
    
    #Calculation of free and bound PFOA in stomach blood
    dCalbSTBf <- koff_alb*CSTBb/MW/1e6 - kon_alb*CalbSTBf*CSTBf/MW/1e6
    
    #Calculation of free and bound PFOA in intestine blood
    dCalbINBf <- koff_alb*CINBb/MW/1e6 - kon_alb*CalbINBf*CINBf/MW/1e6
    
    #Calculation of free and bound PFOA in muscle blood
    dCalbMBf <- koff_alb*CMBb/MW/1e6 - kon_alb*CalbMBf*CMBf/MW/1e6
    
    #Calculation of free and bound PFOA in adipose blood
    dCalbABf <- koff_alb*CABb/MW/1e6 - kon_alb*CalbABf*CABf/MW/1e6
    
    #Calculation of free and bound PFOA in Rest-of-the-body blood
    dCalbRBf <- koff_alb*CRBb/MW/1e6 - kon_alb*CalbRBf*CRBf/MW/1e6
    
    #Calculation of free and bound PFOA in lungs blood
    dCalbLuBf <- koff_alb*CLuBb/MW/1e6 - kon_alb*CalbLuBf*CLuBf/MW/1e6 
    
    #Calculation of free and bound PFOA in spleen blood
    dCalbSPBf <- koff_alb*CSPBb/MW/1e6  - kon_alb*CalbSPBf*CSPBf/MW/1e6 
    
    #Calculation of free and bound PFOA in heart blood
    dCalbHBf <- koff_alb*CHBb/MW/1e6  - kon_alb*CalbHBf*CHBf/MW/1e6 
    
    #Calculation of free and bound PFOA in brain blood
    dCalbBrBf <- koff_alb*CBrBb/MW/1e6  - kon_alb*CalbBrBf*CBrBf/MW/1e6 
    
    #Calculation of free and bound PFOA in gonad blood
    dCalbGoBf <- koff_alb*CGoBb/MW/1e6  - kon_alb*CalbGoBf*CGoBf/MW/1e6 
    
    #Calculation of free and bound PFOA in skin blood
    dCalbSKBf <- koff_alb*CSKBb/MW/1e6  - kon_alb*CalbSKBf*CSKBf/MW/1e6 
    
    #Calculation of free and bound PFOA in bone blood
    dCalbBoBf <- koff_alb*CBoBb/MW/1e6  - kon_alb*CalbBoBf*CBoBf/MW/1e6
    
    #--------------------------------------------------------------
    #Calculation of free concentrations in organ interstitial fluid
    #--------------------------------------------------------------
    
    #Calculation of free and bound PFOA in kidney interstitial fluid
    dCalbKFf <- koff_alb*CKFb/MW/1e6 - kon_alb*CalbKFf*CKFf/MW/1e6
    
    #Calculation of free and bound PFOA in liver interstitial fluid
    dCalbLFf <- koff_alb*CLFb/MW/1e6 - kon_alb*CalbLFf*CLFf/MW/1e6
    
    #Calculation of free and bound PFOA in stomach interstitial fluid
    dCalbSTFf <- koff_alb*CSTFb/MW/1e6 - kon_alb*CalbSTFf*CSTFf/MW/1e6
    
    #Calculation of free and bound PFOA in intestine interstitial fluid
    dCalbINFf <- koff_alb*CINFb/MW/1e6 - kon_alb*CalbINFf*CINFf/MW/1e6
    
    #Calculation of free and bound PFOA in muscle interstitial fluid
    dCalbMFf <- koff_alb*CMFb/MW/1e6 - kon_alb*CalbMFf*CMFf/MW/1e6
    
    #Calculation of free and bound PFOA in adipose interstitial fluid
    dCalbAFf <- koff_alb*CAFb/MW/1e6 - kon_alb*CalbAFf*CAFf/MW/1e6
    
    #Calculation of free and bound PFOA in Rest-of-the-body interstitial fluid
    dCalbRFf <- koff_alb*CRFb/MW/1e6 - kon_alb*CalbRFf*CRFf/MW/1e6
    
    #Calculation of free and bound PFOA in lungs interstitial fluid
    dCalbLuFf <- koff_alb*CLuFb/MW/1e6 - kon_alb*CalbLuFf*CLuFf/MW/1e6
    
    #Calculation of free and bound PFOA in spleen interstitial fluid
    dCalbSPFf <- koff_alb*CSPFb/MW/1e6 - kon_alb*CalbSPFf*CSPFf/MW/1e6
    
    #Calculation of free and bound PFOA in heart interstitial fluid
    dCalbHFf <- koff_alb*CHFb/MW/1e6 - kon_alb*CalbHFf*CHFf/MW/1e6
    
    #Calculation of free and bound PFOA in brain interstitial fluid
    dCalBrFf <- koff_alb*CBrFb/MW/1e6 - kon_alb*CalBrFf*CBrFf/MW/1e6
    
    #Calculation of free and bound PFOA in gonad interstitial fluid
    dCalbGoFf <- koff_alb*CGoFb/MW/1e6 - kon_alb*CalbGoFf*CGoFf/MW/1e6
    
    #Calculation of free and bound PFOA in skin interstitial fluid
    dCalbSKFf <- koff_alb*CSKFb/MW/1e6 - kon_alb*CalbSKFf*CSKFf/MW/1e6
    
    #Calculation of free and bound PFOA in bone interstitial fluid
    dCalbBoFf <- koff_alb*CBoFb/MW/1e6 - kon_alb*CalbBoFf*CBoFf/MW/1e6
    
    #-------------------------------------------------------------------
    #Calculation of free concentrations in organ where we have tissue binding
    #-------------------------------------------------------------------
    #Calculation of free and bound PFOA in kidney Tissue
    dCa2uPTCf <- koff_a2u*CPTCb/MW/1e6 - kon_a2u*Ca2uPTCf*CPTCf/MW/1e6
    dCa2uDALCf <- koff_a2u*CDALCb/MW/1e6 - kon_a2u*Ca2uDALCf*CDALCf/MW/1e6
    dCa2uDTCf <- koff_a2u*CDTCb/MW/1e6 - kon_a2u*Ca2uDTCf*CDTCf/MW/1e6
    dCa2uCDCf <- koff_a2u*CCDCb/MW/1e6 - kon_a2u*Ca2uCDCf*CCDCf/MW/1e6
    dCa2uKTrestf <- koff_a2u*CKTrestb/MW/1e6 - kon_a2u*Ca2uKTrestf*CKTrestf/MW/1e6
    
    dCFabpPTCf <- koff_fabp*CPTCb/MW/1e6 - kon_fabp*CFabpPTCf*CPTCf/MW/1e6
    dCFabpDALCf <- koff_fabp*CDALCb/MW/1e6 - kon_fabp*CFabpDALCf*CDALCf/MW/1e6
    dCFabpDTCf <- koff_fabp*CDTCb/MW/1e6 - kon_fabp*CFabpDTCf*CDTCf/MW/1e6
    dCFabpCDCf <- koff_fabp*CCDCb/MW/1e6 - kon_fabp*CFabpCDCf*CCDCf/MW/1e6
    dCFabpKTrestf <- koff_fabp*CKTrestb/MW/1e6 - kon_fabp*CFabpKTrestf*CKTrestf/MW/1e6
    
    #Calculation of free and bound PFOA in liver tissue
    dCFabpLTf <- koff_fabp*CLTb/MW/1e6 - kon_fabp*CFabpLTf*CLTf/MW/1e6
    
    #Calculation of free and bound PFOA in alveolar lining fluid
    dCalbLuELFf = koff_alb*CLuELFb/MW/1e6 - kon_alb*CalbLuELFf*CLuELFf/MW/1e6
    
    # Bound PFOA
    #Blood
    dMVenb <-  kon_alb*CalbVenf*CVenf*VVen -  koff_alb*CVenb*VVen 
    dMArtb <- kon_alb*CalbArtf*CArtf*VArt -  koff_alb*CArtb*VArt 
    dMKBb <- kon_alb*CalbKBf*CKBf*VKB - koff_alb*CKBb*VKB 
    dMLBb <- kon_alb*CalbLBf*CLBf*VLB - koff_alb*CLBb*VLB 
    dMSTBb <- kon_alb*CalbSTBf*CSTBf*VSTB - koff_alb*CSTBb*VSTB
    dMINBb <- kon_alb*CalbINBf*CINBf*VINB - koff_alb*CINBb*VINB 
    dMMBb <- kon_alb*CalbMBf*CMBf*VMB - koff_alb*CMBb*VMB 
    dMABb <- kon_alb*CalbABf*CABf*VAB - koff_alb*CABb*VAB 
    dMRBb <- kon_alb*CalbRBf*CRBf*VRB - koff_alb*CRBb*VRB 
    dMLuBb <- kon_alb*CalbLuBf*CLuBf*VLuB - koff_alb*CLuBb*VLuB 
    dMSPBb <- kon_alb*CalbSPBf*CSPBf*VSPB - koff_alb*CSPBb*VSPB 
    dMHBb <- kon_alb*CalbHBf*CHBf*VHB - koff_alb*CHBb*VHB
    dMBrBb <- kon_alb*CalbBrBf*CBrBf*VBrB - koff_alb*CBrBb*VBrB
    dMGoBb <- kon_alb*CalbGoBf*CGoBf*VGoB - koff_alb*CGoBb*VGoB 
    dMSKBb <- kon_alb*CalbSKBf*CSKBf*VSKB - koff_alb*CSKBb*VSKB 
    dMBoBb <- kon_alb*CalbBoBf*CBoBf*VBoB - koff_alb*CBoBb*VBoB 
    
    #Interstitial fluid
    dMKFb <- kon_alb*CalbKFf*CKFf*VKF - koff_alb*CKFb*VKF 
    dMLFb <- kon_alb*CalbLFf*CLFf*VLF - koff_alb*CLFb*VLF 
    dMSTFb <- kon_alb*CalbSTFf*CSTFf*VSTF - koff_alb*CSTFb*VSTF 
    dMINFb <- kon_alb*CalbINFf*CINFf*VINF - koff_alb*CINFb*VINF
    dMMFb <- kon_alb*CalbMFf*CMFf*VMF - koff_alb*CMFb*VMF 
    dMAFb <- kon_alb*CalbAFf*CAFf*VAF - koff_alb*CAFb*VAF 
    dMRFb <- kon_alb*CalbRFf*CRFf*VRF - koff_alb*CRFb*VRF 
    dMLuFb <- kon_alb*CalbLuFf*CLuFf*VLuF - koff_alb*CLuFb*VLuF 
    dMSPFb <- kon_alb*CalbSPFf*CSPFf*VSPF - koff_alb*CSPFb*VSPF 
    dMHFb <- kon_alb*CalbHFf*CHFf*VHF - koff_alb*CHFb*VHF 
    dMBrFb <- kon_alb*CalBrFf*CBrFf*VBrF - koff_alb*CBrFb*VBrF 
    dMGoFb <- kon_alb*CalbGoFf*CGoFf*VGoF - koff_alb*CGoFb*VGoF 
    dMSKFb <- kon_alb*CalbSKFf*CSKFf*VSKF - koff_alb*CSKFb*VSKF 
    dMBoFb <- kon_alb*CalbBoFf*CBoFf*VBoF - koff_alb*CBoFb*VBoF 
    
    #Tissue
    dMPTCb <- kon_a2u*Ca2uPTCf*CPTCf*VPTC + kon_fabp*CFabpPTCf*CPTCf*VPTC -
      koff_fabp*CPTCb*VPTC - koff_a2u*CPTCb*VPTC 
    
    dMDALCb <- kon_a2u*Ca2uDALCf*CDALCf*VDALC + kon_fabp*CFabpDALCf*CDALCf*VDALC -
      koff_fabp*CDALCb*VDALC - koff_a2u*CDALCb*VDALC 
    
    dMDTCb <- kon_a2u*Ca2uDTCf*CDTCf*VDTC + kon_fabp*CFabpDTCf*CDTCf*VDTC -
      koff_fabp*CDTCb*VDTC - koff_a2u*CDTCb*VDTC 
    
    dMCDCb <- kon_a2u*Ca2uCDCf*CCDCf*VCDC + kon_fabp*CFabpCDCf*CCDCf*VCDC -
      koff_fabp*CCDCb*VCDC - koff_a2u*CCDCb*VCDC
    
    dMKTrestb <- kon_a2u*Ca2uKTrestf*CKTrestf*VKTrest + kon_fabp*CFabpKTrestf*CKTrestf*VKTrest -
      koff_fabp*CKTrestb*VKTrest - koff_a2u*CKTrestb*VKTrest 
    
    dMLTb <- kon_fabp*CFabpLTf*CLTf*VLT - koff_fabp*CLTb*VLT 
    
    #Alveolar lining fluid
    dMLuELFb <-  kon_alb*CalbLuELFf*CLuELFf*VLuELF -  koff_alb*CLuELFb*VLuELF 
    
    #====================================================================================================================
    
    #Arterial Blood
    dMArtf = QBLu*CLuBf - CArtf*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+
                                   QBST+QBIN+QBGo+QBSK+QBBo) - QGFR*CArtf +
      koff_alb*CArtb*VArt - kon_alb*CalbArtf*CArtf*VArt 
    
    #Venous Blood
    dMVenf =  kUAB*CUA - CVenf*QBLu + QBK*CKBf + QBLtot*CLBf + QBM*CMBf + QBA*CABf + QBR*CRBf+
      QBH*CHBf + QBBr*CBrBf+ QBGo*CGoBf + QBSK*CSKBf + QBBo*CBoBf+
      koff_alb*CVenb*VVen - kon_alb*CalbVenf*CVenf*VVen 
    
    #Kidney
    #blood subcompartment
    dMKBf = QBK*CArtf - QBK*CKBf- (Ptrans_diff_K+PparaKi)*(A_peritubular_PTC+A_peritubular_DTC)*(CKBf-CKFf) + 
      (VmK_baso*CPTCf/(KmK_baso+CPTCf))+
      koff_alb*CKBb*VKB - kon_alb*CalbKBf*CKBf*VKB
    
    #interstitial fluid subcompartment
    dMKFf =  (Ptrans_diff_K+PparaKi)*(A_peritubular_PTC+A_peritubular_DTC)*(CKBf-CKFf)- 
      kPtcF*(CKFf-CPTCf) - kDalcF*(CKFf-CDALCf) -
      kDtcF*(CKFf-CDTCf) - kCdcF*(CKFf-CCDCf)   -  kKTrestF*(CKFf-CKTrestf) + 
      koff_alb*CKFb*VKF - kon_alb*CalbKFf*CKFf*VKF-
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) - (VmK_Oat3*CKFf/(KmK_Oat3+CKFf))
    
    #proximal tubule  cells subcompartment
    dMPTCf = kPtcF*(CKFf-CPTCf) - kPtcTu*(CPTCf - CPT)  +
      (VmK_Oatp*CPT/(KmK_Oatp+CPT)) + (VmK_Urat*CPT/(KmK_Urat+CPT))+
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) + (VmK_Oat3*CKFf/(KmK_Oat3+CKFf)) - 
      (VmK_baso*CPTCf/(KmK_baso+CPTCf)) -(VmK_api*CPTCf/(KmK_api+CPTCf))-
      (kon_a2u*Ca2uPTCf*CPTCf*VPTC + kon_fabp*CFabpPTCf*CPTCf*VPTC -
         koff_fabp*CPTCb*VPTC - koff_a2u*CPTCb*VPTC) 
    
    #Tubule cells in Loop of Henle 
    dMDALCf =  kDalcF*(CKFf-CDALCf) - kDalcTu*(CDALCf - CDAL)- 
      (kon_a2u*Ca2uDALCf*CDALCf*VDALC + kon_fabp*CFabpDALCf*CDALCf*VDALC -
         koff_fabp*CDALCb*VDALC - koff_a2u*CDALCb*VDALC )
    
    #Distal convoluted tubule cells 
    dMDTCf =   kDtcF*(CKFf-CDTCf)- kDtcTu*(CDTCf - CDT)-
      (kon_a2u*Ca2uDTCf*CDTCf*VDTC + kon_fabp*CFabpDTCf*CDTCf*VDTC -
         koff_fabp*CDTCb*VDTC - koff_a2u*CDTCb*VDTC )
    
    #Collecting duct cells 
    dMCDCf =  kCdcF*(CKFf-CCDCf) - kCdcTu*(CCDCf - CCD)- 
      (kon_a2u*Ca2uCDCf*CCDCf*VCDC + kon_fabp*CFabpCDCf*CCDCf*VCDC -
         koff_fabp*CCDCb*VCDC - koff_a2u*CCDCb*VCDC)
    
    #Rest of kidney tissue subcompartment
    dMKTrestf =  kKTrestF*(CKFf-CKTrestf)- 
      (kon_a2u*Ca2uKTrestf*CKTrestf*VKTrest + kon_fabp*CFabpKTrestf*CKTrestf*VKTrest -
         koff_fabp*CKTrestb*VKTrest - koff_a2u*CKTrestb*VKTrest)
    
    #Proximal convoluted tubule
    dMPT =  QGFR*CArtf + kPtcTu*(CPTCf - CPT) - (VmK_Oatp*CPT/(KmK_Oatp+CPT)) - 
      (VmK_Urat*CPT/(KmK_Urat+CPT)) + (VmK_api*CPTCf/(KmK_api+CPTCf))- QDAL*CPT
    
    #Descending limb, Ascending limb (Loop of Henle )
    dMDAL =  QDAL*CPT + kDalcTu*(CDALCf - CDAL)  -  QDT*CDAL
    
    # Distal convoluted tubule 
    dMDT =   QDT*CDAL + kDtcTu*(CDTCf - CDT) -  QCD*CDT
    
    #Collecting duct
    dMCD = QCD*CDT + kCdcTu*(CCDCf - CCD) - Qurine*CCD
    
    # Bladder
    dMBladder = Qurine*CCD - Qurine*CBladder
    
    
    #Liver
    #blood subcompartment
    dMLBf = QBL*CArtf + QBSP*CSPBf + QBIN*CINBf + QBST*CSTBf - 
      QBLtot*CLBf - (Ptrans_diff_L+PparaLi)*AL*(CLBf-CLFf)  +
      koff_alb*CLBb*VLB -kon_alb*CalbLBf*CLBf*VLB 
    #interstitial fluid subcompartment 
    dMLFf = (Ptrans_diff_L+PparaLi)*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) - 
      (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) - (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf)) -
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_alb*CLFb*VLF -kon_alb*CalbLFf*CLFf*VLF 
    #Liver tissue subcompartment
    dMLTf = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) + (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf))+
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_fabp*CLTb*VLT-
      kon_fabp*CFabpLTf*CLTf*VLT - CL_hepatobiliary*CLTf
    #Bile  canaliculi subcompartment
    dMLbile = CL_hepatobiliary*CLTf - Cbile*Qbile
    
    
    
    #Stomach
    #blood subcompartment
    dMSTBf = QBST*CArtf - QBST*CSTBf - (Ptrans_diff_ST+PparaSt)*AST*(CSTBf-CSTFf)  +
      koff_alb*CSTBb*VSTB-kon_alb*CalbSTBf*CSTBf*VSTB 
    #interstitial fluid subcompartment 
    dMSTFf = (Ptrans_diff_ST+PparaSt)*AST*(CSTBf-CSTFf) - kSTFSTT*(CSTFf-CSTTf) +
      koff_alb*CSTFb*VSTF - kon_alb*CalbSTFf*CSTFf*VSTF 
    #Stomach tissue subcompartment
    dMSTTf = kSTFSTT*(CSTFf-CSTTf) + kabST*CSTL 
    #Stomach lumen
    dMSTL = - QGE*CSTL -kabST*CSTL  + CLEal*CLuELFf + CLEua*CUA 
    
    
    #Intestine
    #blood subcompartment
    dMINBf = QBIN*CArtf - QBIN*CINBf - (Ptrans_diff_IN+PparaIn)*AIN*(CINBf-CINFf) +
      koff_alb*CINBb*VINB - kon_alb*CalbINBf*CINBf*VINB
    #interstitial fluid subcompartment 
    dMINFf = (Ptrans_diff_IN+PparaIn)*AIN*(CINBf-CINFf) - kINFINT*(CINFf-CINTf) +
      koff_alb*CINFb*VINF - kon_alb*CalbINFf*CINFf*VINF
    #Intestine tissue subcompartment
    dMINTf = kINFINT*(CINFf-CINTf) + k_gut_in*CINL - k_gut_out*CINTf + (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL)) 
    #Intestine lumen
    dMINL = QGE*CSTL - (CLfeces*CINL) - k_gut_in*CINL + k_gut_out*CINTf + Cbile*Qbile - 
      (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))
    
    
    #Muscle
    #blood subcompartment
    dMMBf = QBM*CArtf - QBM*CMBf - (Ptrans_diff_M+PparaMu)*AM*(CMBf-CMFf)  +
      koff_alb*CMBb*VMB - kon_alb*CalbMBf*CMBf*VMB 
    #interstitial fluid subcompartment 
    dMMFf =(Ptrans_diff_M+PparaMu)*AM*(CMBf-CMFf) - kMFMT*(CMFf- CMTf) +
      koff_alb*CMFb*VMF - kon_alb*CalbMFf*CMFf*VMF
    #Muscle tissue subcompartment 
    dMMTf = kMFMT*(CMFf- CMTf) 
    
    
    #Adipose
    #blood subcompartment
    dMABf = QBA*CArtf - QBA*CABf - (Ptrans_diff_A+PparaAd)*AA*(CABf-CAFf)  +
      koff_alb*CABb*VAB - kon_alb*CalbABf*CABf*VAB
    #interstitial fluid subcompartment 
    dMAFf = (Ptrans_diff_A+PparaAd)*AA*(CABf-CAFf)  - kAFAT*(CAFf-CATf) +
      koff_alb*CAFb*VAF - kon_alb*CalbAFf*CAFf*VAF 
    #Adipose tissue subcompartment 
    dMATf =  kAFAT*(CAFf-CATf)
    
    
    #Rest of body
    #blood subcompartment
    dMRBf = QBR*CArtf - QBR*CRBf - (Ptrans_diff_R+PparaRe)*AR*(CRBf-CRFf)  +
      koff_alb*CRBb*VRB - kon_alb*CalbRBf*CRBf*VRB
    #interstitial fluid subcompartment 
    dMRFf = (Ptrans_diff_R+PparaRe)*AR*(CRBf-CRFf) - kRFRT*(CRFf -CRTf) +
      koff_alb*CRFb*VRF - kon_alb*CalbRFf*CRFf*VRF 
    #Rest of body tissue subcompartment 
    dMRTf = kRFRT*(CRFf -CRTf) 
    
    #Upper airways
    dMUA = - CLEua*CUA  - kUAB * CUA 
    
    #Lung 
    #blood subcompartment
    dMLuBf = CVenf*QBLu - QBLu*CLuBf - (Ptrans_diff_Lu+PparaLu)*ALu*(CLuBf-CLuFf) +
      koff_alb*CLuBb*VLuB - kon_alb*CalbLuBf*CLuBf*VLuB 
    #interstitial fluid subcompartment
    dMLuFf = (Ptrans_diff_Lu+PparaLu)*ALu*(CLuBf-CLuFf)+ kLuTLuF*(CLuTf-CLuFf) + 
      koff_alb*CLuFb*VLuF - kon_alb*CalbLuFf*CLuFf*VLuF - (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) 
    #Lung tissue
    dMLuTf =  - kLuTLuF*(CLuTf-CLuFf) -  kLuTLuELF*(CLuTf-CLuELFf) + (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) +
      (VmLu_Oatp_ap*CLuELFf/(KmLu_Oatp_ap+CLuELFf)) 
    #Alveolar lining fluid
    dMLuELFf =  k_desorption * CLuELFdust + kLuTLuELF*(CLuTf-CLuELFf) + 
      koff_alb*CLuELFb*VLuELF - kon_alb*CalbLuELFf*CLuELFf*VLuELF -
      (VmLu_Oatp_ap*CLuELFf/(KmLu_Oatp_ap+CLuELFf)) - CLEal*CLuELFf
    dMLuELFdust <- -k_desorption * CLuELFdust
    
    #Spleen
    #blood subcompartment
    dMSPBf = QBSP*CArtf - QBSP*CSPBf - (Ptrans_diff_SP+PparaSp)*ASP*(CSPBf-CSPFf) + 
      koff_alb*CSPBb*VSPB - kon_alb*CalbSPBf*CSPBf*VSPB
    #interstitial fluid subcompartment 
    dMSPFf = (Ptrans_diff_SP+PparaSp)*ASP*(CSPBf-CSPFf) - kSPFSPT*(CSPFf -CSPTf) +
      koff_alb*CSPFb*VSPF - kon_alb*CalbSPFf*CSPFf*VSPF 
    #Spleen tissue subcompartment 
    dMSPTf = kSPFSPT*(CSPFf -CSPTf)
    
    
    #Heart
    #blood subcompartment
    dMHBf = QBH*CArtf - QBH*CHBf - (Ptrans_diff_H+PparaHt)*AH*(CHBf-CHFf) + 
      koff_alb*CHBb*VHB - kon_alb*CalbHBf*CHBf*VHB 
    #interstitial fluid subcompartment 
    dMHFf = (Ptrans_diff_H+PparaHt)*AH*(CHBf-CHFf) - kHFHT*(CHFf -CHTf) + 
      koff_alb*CHFb*VHF - kon_alb*CalbHFf*CHFf*VHF 
    #Heart tissue subcompartment 
    dMHTf = kHFHT*(CHFf -CHTf) 
    
    
    #Brain
    #blood subcompartment
    dMBrBf = QBBr*CArtf - QBBr*CBrBf - (Ptrans_diff_Br+PparaBr)*ABr*(CBrBf-CBrFf)  + 
      koff_alb*CBrBb*VBrB - kon_alb*CalbBrBf*CBrBf*VBrB 
    #interstitial fluid subcompartment 
    dMBrFf = (Ptrans_diff_Br+PparaBr)*ABr*(CBrBf-CBrFf) - kBrFBrT*(CBrFf -CBrTf) +
      koff_alb*CBrFb*VBrF - kon_alb*CalBrFf*CBrFf*VBrF 
    #Brain tissue subcompartment 
    dMBrTf = kBrFBrT*(CBrFf -CBrTf) 
    
    
    #Gonads
    #blood subcompartment
    dMGoBf = QBGo*CArtf - QBGo*CGoBf - (Ptrans_diff_Go+PparaGo)*AGo*(CGoBf-CGoFf) +
      koff_alb*CGoBb*VGoB - kon_alb*CalbGoBf*CGoBf*VGoB 
    #interstitial fluid subcompartment 
    dMGoFf = (Ptrans_diff_Go+PparaGo)*AGo*(CGoBf-CGoFf) - kGoFGoT*(CGoFf -CGoTf) +
      koff_alb*CGoFb*VGoF - kon_alb*CalbGoFf*CGoFf*VGoF 
    #gonads tissue subcompartment 
    dMGoTf = kGoFGoT*(CGoFf -CGoTf) 
    
    
    #Skin
    #blood subcompartment
    dMSKBf = QBSK*CArtf - QBSK*CSKBf - (Ptrans_diff_SK+PparaSk)*ASK*(CSKBf-CSKFf) +
      koff_alb*CSKBb*VSKB - kon_alb*CalbSKBf*CSKBf*VSKB 
    #interstitial fluid subcompartment
    dMSKFf =(Ptrans_diff_SK+PparaSk)*ASK*(CSKBf-CSKFf) - kSKFSKT*(CSKFf -CSKTf) +
      koff_alb*CSKFb*VSKF - kon_alb*CalbSKFf*CSKFf*VSKF
    #Skin tissue subcompartment
    dMSKTf = kSKFSKT*(CSKFf -CSKTf) 
    
    
    #Bones
    #blood subcompartment
    dMBoBf = QBBo*CArtf - QBBo*CBoBf - (Ptrans_diff_Bo+PparaBo)*ABo*(CBoBf-CBoFf) +
      koff_alb*CBoBb*VBoB - kon_alb*CalbBoBf*CBoBf*VBoB 
    #interstitial fluid subcompartment
    dMBoFf = (Ptrans_diff_Bo+PparaBo)*ABo*(CBoBf-CBoFf) - kBoFBoT*(CBoFf -CBoTf) +
      koff_alb*CBoFb*VBoF -  kon_alb*CalbBoFf*CBoFf*VBoF 
    #Bones tissue subcompartment
    dMBoTf = kBoFBoT*(CBoFf -CBoTf) 
    
    #Excreta#
    dMfeces <- CLfeces*CINL
    dMurine <- Qurine*CBladder
    dVurine = Qurine
    dVfeces = Qfeces
    
    #Concentration calculation in each compartment 
    
    Cblood <- (MVen +MArt)/ (VVen+VArt)
    Mblood <- MVen +MArt
    Cplasma <- Cblood/(1-Hct)
    
    Mkidney <- MKB + MKF+ MKT + Mfil
    Ckidney <- Mkidney/VK    
    
    Cliver <- (MLB + MLF+ MLT + MLbile )/(VLB+VLF+VLT+VLbile)
    Mliver <- MLB + MLF+ MLT + MLbile
    
    Cstomach <-  (MSTB + MSTF+ MSTT + MSTL)/(VSTB+VSTF+VSTT+VSTL)
    Cintestine <-  (MINB + MINF+ MINT+MINL)/(VINB+VINF+VINT+VINL)
    Cmuscle <-  (MMB + MMF+ MMT)/(VMB+VMF+VMT)
    Cadipose <-  (MAB + MAF+ MAT)/(VAB+VAF+VAT)
    
    CUpperair <- MUA/VUA
    CLuELF <- (MLuELFf+MLuELFb+MLuELFdust)/VLuELF
    Clungs <-  (MLuB + MLuF+ MLuT + MLuELF +MLuELFdust )/(VLuB+VLuF+VLuT+VLuELF)
    Clungtissue <- (MLuB + MLuF+ MLuT)/(VLuB+VLuF+VLuT)
    Crest <-  (MRB + MRF+ MRT)/(VRB+VRF+VRT)
    Ccarcass <- (MMB+MMF+MMT+MAB+MAF+MAT+MRB+MRF+MRT+MBoB+MBoF+MBoT+MSKB+MSKF+MSKT)/(VM+VA+VR+VBo+VSK)
    Cfeces <- Mfeces/(Vfeces*feces_density)
    Curine <- Murine/Vurine
    Cspleen <-  (MSPB + MSPF+ MSPT)/(VSPB+VSPF+VSPT)
    Cheart <-  (MHB + MHF+ MHT)/(VHB+VHF+VHT)
    
    Cbrain <-  (MBrB + MBrF+ MBrT)/(VBrB+VBrF+VBrT)
    Mbrain <- MBrB + MBrF+ MBrT
    
    Cgonads <-  (MGoB + MGoF+ MGoT)/(VGoB+VGoF+VGoT)
    Cskin <-  (MSKB + MSKF+ MSKT)/(VSKB+VSKF+VSKT)
    Cbones <-  (MBoB + MBoF+ MBoT)/(VBoB+VBoF+VBoT)
    
    
    #Concentration calculation in each compartment 
    
    
    list(c( 'dCalbVenf' = dCalbVenf, 'dCalbArtf' = dCalbArtf, 
            'dCalbKBf' = dCalbKBf, 'dCalbLBf' = dCalbLBf, 'dCalbSTBf' = dCalbSTBf, 
            'dCalbINBf' = dCalbINBf, 'dCalbMBf' = dCalbMBf, 'dCalbABf' = dCalbABf,
            'dCalbRBf' = dCalbRBf, 'dCalbLuBf' = dCalbLuBf, 'dCalbSPBf' = dCalbSPBf,
            'dCalbHBf' = dCalbHBf, 'dCalbBrBf' = dCalbBrBf, 'dCalbGoBf' = dCalbGoBf, 
            'dCalbSKBf' = dCalbSKBf, 'dCalbBoBf'=dCalbBoBf, 'dCalbKFf' = dCalbKFf,
            'dCalbLFf' = dCalbLFf, 'dCalbSTFf' = dCalbSTFf,'dCalbINFf' = dCalbINFf,
            'dCalbMFf' = dCalbMFf, 'dCalbAFf' = dCalbAFf,'dCalbRFf' = dCalbRFf,
            'dCalbLuFf' = dCalbLuFf, 'dCalbSPFf' = dCalbSPFf, 'dCalbHFf' = dCalbHFf,
            'dCalBrFf' = dCalBrFf, 'dCalbGoFf' = dCalbGoFf, 'dCalbSKFf' = dCalbSKFf,
            'dCalbBoFf' = dCalbBoFf, 'dCa2uPTCf' = dCa2uPTCf, 'dCa2uDALCf' = dCa2uDALCf,
            'dCa2uDTCf' = dCa2uDTCf, 'dCa2uCDCf' = dCa2uCDCf, 'dCa2uKTrestf' = dCa2uKTrestf,
            'dCFabpPTCf' = dCFabpPTCf, 'dCFabpDALCf' = dCFabpDALCf, 'dCFabpDTCf' = dCFabpDTCf,
            'dCFabpCDCf' = dCFabpCDCf, 'dCFabpKTrestf' = dCFabpKTrestf,
            'dCFabpLTf' = dCFabpLTf, 'dCalbLuELFf' = dCalbLuELFf,
            
            
            'dMVenb' = dMVenb, 'dMArtb' = dMArtb, 'dMKBb' = dMKBb, 
            'dMLBb' = dMLBb,'dMSTBb' = dMSTBb,'dMINBb' = dMINBb,'dMMBb' = dMMBb,
            'dMABb' = dMABb, 'dMRBb' = dMRBb,'dMLuBb' = dMLuBb, 'dMSPBb' = dMSPBb, 
            'dMHBb' = dMHBb,  'dMBrBb' = dMBrBb,  'dMGoBb' = dMGoBb, 
            'dMSKBb' = dMSKBb, 'dMBoBb' = dMBoBb,'dMKFb' = dMKFb, 'dMLFb' = dMLFb, 
            'dMSTFb' = dMSTFb,  'dMINFb' = dMINFb, 'dMMFb' = dMMFb, 
            'dMAFb' = dMAFb, 'dMRFb' = dMRFb, 'dMLuFb' = dMLuFb, 
            'dMSPFb' = dMSPFb,  'dMHFb' = dMHFb, 'dMBrFb' = dMBrFb, 
            'dMGoFb' = dMGoFb, 'dMSKFb' = dMSKFb, 'dMBoFb' = dMBoFb,
            'dMPTCb' = dMPTCb, 'dMDALCb' = dMDALCb, 'dMDTCb' = dMDTCb,
            'dMCDCb' = dMCDCb, 'dMKTrestb' = dMKTrestb,
            'dMLTb' = dMLTb, 'dMLuELFb'=dMLuELFb,
            
            
            'dMArtf'=dMArtf, 'dMVenf'=dMVenf, 'dMKBf'=dMKBf, 
            'dMKFf'=dMKFf,  'dMPTCf' = dMPTCf, 'dMDALCf' = dMDALCf,  
            'dMDTCf' = dMDTCf, 'dMCDCf' = dMCDCf,
            'dMKTrestf' = dMKTrestf, 'dMPT' = dMPT,  'dMDAL' = dMDAL,
            'dMDT' =  dMDT, 'dMCD' = dMCD,
            
            'dMBladder' = dMBladder, 'dMLBf'=dMLBf, 
            'dMLFf'=dMLFf, 'dMLTf'=dMLTf, 'dMLbile'=dMLbile,
            
            'dMSTBf'=dMSTBf, 'dMSTFf'=dMSTFf, 'dMSTTf'=dMSTTf, 'dMSTL'=dMSTL,
            'dMINBf'=dMINBf, 'dMINFf'=dMINFf, 'dMINTf'=dMINTf,'dMINL'=dMINL,
            
            'dMMBf'=dMMBf, 'dMMFf'=dMMFf, 'dMMTf'=dMMTf,
            'dMABf'=dMABf, 'dMAFf'=dMAFf, 'dMATf'=dMATf, 
            'dMRBf'=dMRBf, 'dMRFf'=dMRFf,'dMRTf'=dMRTf, 'dMUA'=dMUA,
            'dMLuBf'=dMLuBf, 'dMLuFf'=dMLuFf,'dMLuTf'=dMLuTf,'dMLuELFf' = dMLuELFf,
            "dMLuELFdust" = dMLuELFdust,
            
            'dMSPBf'=dMSPBf, 'dMSPFf'=dMSPFf, 'dMSPTf'=dMSPTf,
            'dMHBf'=dMHBf, 'dMHFf'=dMHFf, 'dMHTf'=dMHTf,
            'dMBrBf'=dMBrBf, 'dMBrFf'=dMBrFf, 'dMBrTf'=dMBrTf,
            'dMGoBf'=dMGoBf, 'dMGoFf'=dMGoFf, 'dMGoTf'=dMGoTf,
            'dMSKBf'=dMSKBf, 'dMSKFf'=dMSKFf, 'dMSKTf'=dMSKTf,
            'dMBoBf'=dMBoBf, 'dMBoFf'=dMBoFf, 'dMBoTf'=dMBoTf,
            'dMfeces'=dMfeces,'dMurine'=dMurine,'dVfeces'=dVfeces,'dVurine'=dVurine 
    ),
    
    
    'CVen'=CVen, 'CVenb'=CVenb, 'CVenf'=CVenf, 'CArt'=CArt, 'CArtf'=CArtf, 'CArtb'=CArtb,
    'CKB'=CKB, 'CKBf'=CKBf, 'CKBb'=CKBb, 'CKF'=CKF, 'CKFf'=CKFf, 'CKFb'=CKFb, 'CKT'=CKT,
    'CKTf'=CKTf, 'CKTb'=CKTb,  'CBladder' = CBladder,'CKTrestf' = CKTrestf,
    'CKTrestb' =CKTrestb, 'CPTCf' = CPTCf,
    'CPTCb' = CPTCb, 'CDALCf' = CDALCf, 'CDALCb' = CDALCb,
    'CDTCf' = CDTCf, 'CDTCb' = CDTCb, 'CCDCf' = CCDCf,
    'CCDCb' = CCDCb, 'CPT' = CPT, 'CDAL' = CDAL, 'CDT' = CDT, 'CCD' = CCD,
    'CLB'=CLB, 'CLBf'=CLBf, 'CLBb'=CLBb, 
    'CLF'=CLF, 'CLFf'=CLFf, 'CLFb'=CLFb, 'CLT'=CLT, 'CLTf'=CLTf, 'CLTb'=CLTb, 
    'CSTB'=CSTB, 'CSTBf'=CSTBf, 'CSTBb'=CSTBb, 'CSTF'=CSTF, 'CSTFf'=CSTFf, 'CSTFb'=CSTFb,
    'CSTTf'=CSTTf, 'CINB'=CINB, 'CINBf'=CINBf, 'CINBb'=CINBb, 'CINF'=CINF, 'CINFf'=CINFf,
    'CINFb'=CINFb, 'CINTf'=CINTf, 'CSTL'=CSTL, 'CINL'=CINL, 'CMB'=CMB, 'CMBf'=CMBf,
    'CMBb'=CMBb, 'CMF'=CMF, 'CMFf'=CMFf, 'CMFb'=CMFb, 'CMTf'=CMTf, 'CAB'=CAB, 'CABf'=CABf,
    'CABb'=CABb, 'CAF'=CAF, 'CAFf'=CAFf, 'CAFb'=CAFb, 'CATf'=CATf, 'CRB'=CRB, 'CRBf'=CRBf,
    'CRBb'=CRBb, 'CRF'=CRF, 'CRFf'=CRFf, 'CRFb'=CRFb, 'CRTf'=CRTf, 'CUA'=CUA, 'CLuB'=CLuB,
    'CLuBf'=CLuBf, 'CLuBb'=CLuBb, 'CLuF'=CLuF, 'CLuFf'=CLuFf, 'CLuFb'=CLuFb, 'CLuTf'=CLuTf,
    'CLuELF'=CLuELF, 'CLuELFf'=CLuELFf, 'CLuELFb'=CLuELFb, 'CSPB'=CSPB, 'CSPBf'=CSPBf, 
    'CSPBb'=CSPBb, 'CSPF'=CSPF, 'CSPFf'=CSPFf, 'CSPFb'=CSPFb, 'CSPTf'=CSPTf,  'CHB'=CHB,
    'CHBf'=CHBf, 'CHBb'=CHBb, 'CHF'=CHF, 'CHFf'=CHFf, 'CHFb'=CHFb, 'CHTf'=CHTf, 
    'CBrB'=CBrB, 'CBrBf'=CBrBf, 'CBrBb'=CBrBb, 'CBrF'=CBrF, 'CBrFf'=CBrFf, 'CBrFb'=CBrFb,
    'CBrTf'=CBrTf, 'CGoB'=CGoB, 'CGoBf'=CGoBf, 'CGoBb'=CGoBb, 'CGoF'=CGoF, 'CGoFf'=CGoFf, 
    'CGoFb'=CGoFb, 'CGoTf'=CGoTf, 'CSKB'=CSKB, 'CSKBf'=CSKBf, 'CSKBb'=CSKBb, 'CSKF'=CSKF,
    'CSKFf'=CSKFf, 'CSKFb'=CSKFb, 'CSKTf'=CSKTf, 'CBoB'=CBoB, 'CBoBf'=CBoBf, 'CBoBb'=CBoBb,
    'CBoF'=CBoF, 'CBoFf'=CBoFf, 'CBoFb'=CBoFb, 'CBoTf'=CBoTf,     
    
    'Cblood'=Cblood, 'Mblood'=Mblood, 'Cplasma'=Cplasma, 
    'Ckidney'=Ckidney, 'Mkidney'=Mkidney, 'Cliver'=Cliver, 'Mliver'=Mliver, 
    'Cstomach'=Cstomach, 'Cintestine'=Cintestine, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose,
    'CUpperair'=CUpperair,  "Clungtissue" = Clungtissue, 'Clungs'=Clungs,
    'Crest'=Crest, 'Ccarcass'=Ccarcass, 'Cfeces'=Cfeces,
    'Curine'=Curine, 'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain, 
    'Mbrain'=Mbrain, 'Cgonads'=Cgonads, 'Cskin'=Cskin, 'Cbones'=Cbones, 
    "Cbile" = Cbile
    
    )
    
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    CalbVenf<- CalbB_init; MVenf<- 0; MVenb<- 0;
    CalbArtf<- CalbB_init; MArtf<- 0; MArtb<- 0; 
    CalbKBf<- CalbB_init; MKBf<- 0; MKBb<- 0;
    CalbLBf<- CalbB_init; MLBf<- 0; MLBb<- 0;
    CalbSTBf<- CalbB_init; MSTBf<- 0; MSTBb<- 0;
    CalbINBf<- CalbB_init; MINBf<- 0; MINBb<- 0;
    CalbMBf<- CalbB_init; MMBf<- 0; MMBb<- 0; 
    CalbABf<- CalbB_init; MABf<- 0; MABb<- 0;
    CalbRBf<- CalbB_init; MRBf<- 0; MRBb<- 0; 
    CalbLuBf<- CalbB_init; MLuBf<- 0; MLuBb<- 0;
    CalbSPBf<- CalbB_init; MSPBf<- 0; MSPBb<- 0;
    CalbHBf<- CalbB_init; MHBf<- 0; MHBb<- 0; 
    CalbBrBf<- CalbB_init; MBrBf<- 0; MBrBb<- 0;
    CalbGoBf<- CalbB_init; MGoBf<- 0; MGoBb<- 0; 
    CalbSKBf<- CalbB_init; MSKBf<- 0; MSKBb<- 0;
    CalbBoBf<- CalbB_init; MBoBf<- 0; MBoBb<- 0; 
    
    CalbKFf<- CalbKF_init; MKFf<- 0; MKFb<- 0;
    CalbLFf<- CalbLF_init; MLFf<- 0; MLFb<- 0; 
    CalbSTFf<- CalbSTF_init; MSTFf<- 0; MSTFb<- 0;
    CalbINFf<- CalbINF_init; MINFf<- 0; MINFb<- 0; 
    CalbMFf<- CalbMF_init; MMFf<- 0; MMFb<- 0;
    CalbAFf<- CalbAF_init; MAFf<- 0; MAFb<- 0; 
    CalbRFf<- CalbRF_init; MRFf<- 0; MRFb<- 0; 
    CalbLuFf<- CalbLuF_init; MLuFf<- 0; MLuFb<- 0;
    CalbSPFf<- CalbSPF_init; MSPFf<- 0; MSPFb<- 0; 
    CalbHFf<- CalbHF_init; MHFf<- 0; MHFb<- 0;
    CalBrFf<- CalbBrF_init; MBrFf<- 0; MBrFb<- 0;
    CalbGoFf<- CalbGoF_init; MGoFf<- 0; MGoFb<- 0; 
    CalbSKFf<- CalbSKF_init; MSKFf<- 0; MSKFb<- 0;
    CalbBoFf<- CalbBoF_init; MBoFf<- 0; MBoFb<- 0;
    Ca2uPTCf = Ca2uKT_init; Ca2uDALCf =  Ca2uKT_init;
    Ca2uDTCf =  Ca2uKT_init; Ca2uCDCf =  Ca2uKT_init; Ca2uKTrestf = Ca2uKT_init;
    CFabpPTCf =  CFabpKT_init; CFabpDALCf = CFabpKT_init; CFabpDTCf = CFabpKT_init;
    CFabpCDCf = CFabpKT_init; CFabpKTrestf =CFabpKT_init;
    
    MPTCf <- 0;MDALCf <- 0;
    MDTCf <- 0; MCDCf <- 0;
    MKTrestf <- 0; MPT <- 0;  MDAL<- 0;
    MDT <- 0; MCD <- 0;
    
    CFabpLTf<- CFabpLT_init; CalbLuELFf<- CalbLuELF_init;
    MLTf<- 0; MLTb<- 0;  MPTCb = 0; MDALCb = 0; MDTCb = 0
    MCDCb = 0; MKTrestb = 0; MLbile <-0; MSTTf <- 0;  MINTf <- 0; MLuTf <- 0; 
    MLuELFf <- 0; MUA <- 0; MLuELFb<- 0; MLuELFdust <- 0; MSPTf <- 0; MHTf <- 0;  MBrTf <- 0;
    MGoTf <- 0;  MSKTf <- 0; MBoTf <- 0; MMTf <- 0; MATf <- 0; MRTf <- 0;
    
    MPT <- 0; MBladder <- 0; Murine <-0;MSTL <-0;  MINL <-0;
    Mfeces <-0;  Vurine <-0; Vfeces <-0
    
    return(c('CalbVenf' = CalbVenf, 'CalbArtf' = CalbArtf, 
             'CalbKBf' = CalbKBf, 'CalbLBf' = CalbLBf, 'CalbSTBf' = CalbSTBf, 
             'CalbINBf' = CalbINBf, 'CalbMBf' = CalbMBf, 'CalbABf' = CalbABf,
             'CalbRBf' = CalbRBf, 'CalbLuBf' = CalbLuBf, 'CalbSPBf' = CalbSPBf,
             'CalbHBf' = CalbHBf, 'CalbBrBf' = CalbBrBf, 'CalbGoBf' = CalbGoBf, 
             'CalbSKBf' = CalbSKBf, 'CalbBoBf'=CalbBoBf, 'CalbKFf' = CalbKFf,
             'CalbLFf' = CalbLFf, 'CalbSTFf' = CalbSTFf,'CalbINFf' = CalbINFf,
             'CalbMFf' = CalbMFf, 'CalbAFf' = CalbAFf,'CalbRFf' = CalbRFf,
             'CalbLuFf' = CalbLuFf, 'CalbSPFf' = CalbSPFf, 'CalbHFf' = CalbHFf,
             'CalBrFf' = CalBrFf, 'CalbGoFf' = CalbGoFf, 'CalbSKFf' = CalbSKFf,
             'CalbBoFf' = CalbBoFf, 
             'Ca2uPTCf' = Ca2uPTCf, 'Ca2uDALCf' =  Ca2uDALCf,
             'Ca2uDTCf' =  Ca2uDTCf, 'Ca2uCDCf' =  Ca2uCDCf, 'Ca2uKTrestf' = Ca2uKTrestf,
             'CFabpPTCf' =  CFabpPTCf, 'CFabpDALCf' = CFabpDALCf, 'CFabpDTCf' = CFabpDTCf,
             'CFabpCDCf' = CFabpCDCf, 'CFabpKTrestf' =CFabpKTrestf,
             'CFabpLTf' = CFabpLTf, 'CalbLuELFf' = CalbLuELFf,
             
             
             'MVenb' = MVenb, 'MArtb' = MArtb, 'MKBb' = MKBb, 
             'MLBb' = MLBb,'MSTBb' = MSTBb,'MINBb' = MINBb,'MMBb' = MMBb,
             'MABb' = MABb, 'MRBb' = MRBb,'MLuBb' = MLuBb, 'MSPBb' = MSPBb, 
             'MHBb' = MHBb,  'MBrBb' = MBrBb,  'MGoBb' = MGoBb, 
             'MSKBb' = MSKBb, 'MBoBb' = MBoBb, 'MKFb' = MKFb, 'MLFb' = MLFb, 
             'MSTFb' = MSTFb,  'MINFb' = MINFb, 'MMFb' = MMFb, 
             'MAFb' = MAFb, 'MRFb' = MRFb, 'MLuFb' = MLuFb, 
             'MSPFb' = MSPFb,  'MHFb' = MHFb, 'MBrFb' = MBrFb, 
             'MGoFb' = MGoFb, 'MSKFb' = MSKFb, 'MBoFb' = MBoFb,
             'MPTCb' = MPTCb, 'MDALCb' = MDALCb, 'MDTCb' = MDTCb,
             'MCDCb' = MCDCb, 'MKTrestb' = MKTrestb, 'MLTb' = MLTb, 
             
             'MLuELFb'=MLuELFb,
             
             'MArtf'=MArtf, 'MVenf'=MVenf, 'MKBf'=MKBf, 
             'MKFf'=MKFf,  'MPTCf' = MPTCf, 'MDALCf' = MDALCf,  
             'MDTCf' = MDTCf, 'MCDCf' = MCDCf,
             'MKTrestf' = MKTrestf, 'MPT' = MPT,  'MDAL' = MDAL,
             'MDT' =  MDT, 'MCD' = MCD,
             
             'MBladder' = MBladder,  'MLBf'=MLBf, 
             'MLFf'=MLFf, 'MLTf'=MLTf, 'MLbile'= MLbile,
             
             'MSTBf'=MSTBf, 'MSTFf'=MSTFf, 'MSTTf'=MSTTf, 'MSTL'=MSTL,
             'MINBf'=MINBf, 'MINFf'=MINFf, 'MINTf'=MINTf,'MINL'=MINL,
             
             'MMBf'=MMBf, 'MMFf'=MMFf, 'MMTf'=MMTf,
             'MABf'=MABf, 'MAFf'=MAFf, 'MATf'=MATf, 
             'MRBf'=MRBf, 'MRFf'=MRFf,'MRTf'=MRTf,'MUA'=MUA,
             'MLuBf'=MLuBf, 'MLuFf'=MLuFf,'MLuTf'=MLuTf,'MLuELFf' = MLuELFf,
             "MLuELFdust" = MLuELFdust,
             
             
             'MSPBf'=MSPBf, 'MSPFf'=MSPFf, 'MSPTf'=MSPTf,
             'MHBf'=MHBf, 'MHFf'=MHFf, 'MHTf'=MHTf,
             'MBrBf'=MBrBf, 'MBrFf'=MBrFf, 'MBrTf'=MBrTf,
             'MGoBf'=MGoBf, 'MGoFf'=MGoFf, 'MGoTf'=MGoTf,
             'MSKBf'=MSKBf, 'MSKFf'=MSKFf, 'MSKTf'=MSKTf,
             'MBoBf'=MBoBf, 'MBoFf'=MBoFf, 'MBoTf'=MBoTf,
             'Mfeces'=Mfeces,'Murine'=Murine,'Vfeces'=Vfeces,'Vurine'=Vurine  
             
    ))
    
    
  })
}


create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "iv"){
        events <- list(data = rbind(data.frame(var = c("MVenf"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("MSTL"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "inh"){
        
        events <- list(data = rbind(data.frame(var = c("MLuELFdust"),  time = admin.time, 
                                               
                                               value = admin.dose*depfr_AF, method = c("add")) ))
      }else if (admin.type == "nasal"){
        
        events <- list(data = rbind(data.frame(var = c("MUA"),  time = admin.time, 
                                               
                                               value = c(admin.dose*depfr_head), method = c("add")),
                                    
                                    data.frame(var = c("MLuELFf"),  time = admin.time, 
                                               
                                               value = c(admin.dose*depfr_AF), method = c("add")) ))
        
      }
    }
    return(events)
  })
}

create_all_fixed_params <- function(){
  params <- list()
  
  # Set up simulations for the 1st case, i.e. kudo (2007) high dose, tissues
  BW <- 0.29  # body weight (kg)
  admin.dose_per_g <- 16.56 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in hours
  admin.type <- "iv"
  sex <- "M" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[1]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 2nd case, i.e. kudo (2007) low dose, tissues
  BW <- 0.29  # body weight (kg)
  admin.dose_per_g <- 0.041 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 # time when doses are administered, in hours
  admin.type <- "iv"
  sex <- "M" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[2]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 3rd case, i.e. kim (2016) IV male tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "iv"
  sex <- "M" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[3]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 4th case, i.e. kim (2016) ORAL male tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[4]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "iv"
  sex <- "F" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[5]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 6th case, i.e. kim (2016) IV female tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "F" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[6]] <- create_fixed_params(user_input)
  
  
  # Set up simulations for the 7th case, i.e. Gustafsson (2022) oral male tissues
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  admin.dose_per_g <- 0.364 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" =0)
  params[[7]] <- create_fixed_params(user_input)
  
  
  # Set up simulations for the 8th case, i.e. Gustafsson (2022) Inhalation male blood
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  duration <- 0.345 #hours, 22.5 min
  k = 9#partition of administration packages
  admin.dose <- rep((512.5*0.62*0.335)/k, length.out = k) #ug PFOA, for 22.5 min inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "inh"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" = 1)
  params[[8]] <- create_fixed_params(user_input)
  
  
  
  # Set up simulations for the 9th case, i.e. Hinderliter Inhalation male single low
  BW <- (0.311+0.195)/2  #kg, https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  BFn = 166*60# 1/h
  TVn = 1.71*1e-3# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 1.2 # administered dose in mg/m^3
  depfr_head <- 0.3073
  depfr_AF <- (0.1537+0.0281)
  k = 1*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF)
  params[[9]] <- create_fixed_params(user_input)
  
  
  # Set up simulations for the 10th case, i.e. Hinderliter Inhalation male single medium
  BW <- (0.311+0.195)/2  #kg, https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  BFn = 166*60# 1/h
  TVn = 1.71*1e-3# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 9.8 # administered dose in mg/m^3
  depfr_head <- 0.3272
  depfr_AF <- (0.1289+0.0276)
  k = 1*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF)
  params[[10]] <- create_fixed_params(user_input) 
  
  
  # Set up simulations for the 11th case, i.e. Hinderliter Inhalation male single high
  BW <- (0.311+0.195)/2 #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  BFn = 166*60# 1/h
  TVn = 1.71*1e-3# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 27 # administered dose in mg/m^3
  depfr_head <- 0.3812
  depfr_AF <- (0.1694+0.0259)
  k = 1*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF ) 
  params[[11]] <- create_fixed_params(user_input)
  
  
  # Set up simulations for the 12th case, i.e. Hinderliter Inhalation female single low
  BW <- (0.197+0.145)/2  #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  BFn = 166*60# 1/h
  TVn = 1.05*1e-3# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 1.2 # administered dose in mg/m^3
  depfr_head <- 0.2991
  depfr_AF <- (0.1305+0.0243)
  k = 1*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  params[[12]] <- create_fixed_params(user_input)
  
  
  # Set up simulations for the 13th case, i.e. Hinderliter Inhalation female single medium
  BW <- (0.197+0.145)/2  #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  BFn = 166*60# 1/h
  TVn = 1.05*1e-3# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 9.8 # administered dose in mg/m^3
  depfr_head <- 0.3376
  depfr_AF <- (0.1074+0.0228)
  k = 1*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  params[[13]] <- create_fixed_params(user_input) 
  
  
  # Set up simulations for the 14th case, i.e. Hinderliter Inhalation female single high
  BW <- (0.197+0.145)/2  #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  BFn = 166*60# 1/h
  TVn = 1.05*1e-3# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 27 # administered dose in mg/m^3
  depfr_head <- 0.3404
  depfr_AF <- (0.1510+0.0242)
  k = 1*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  params[[14]] <- create_fixed_params(user_input) 
  
  
  
  return(params)
  
}
##########
#########
########
#####
###
##
#
obj.func <- function(x, dataset, fixed_params){
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- x
  ##########################
  #-------------------------
  # Kudo high
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. kudo (2007) high dose, tissues
  BW <- 0.29  # body weight (kg)
  sex <- "M" 
  variable_params <- create_variable_params(BW, sex,estimated_params, fixed_params[[1]])
  params <- c(fixed_params[[1]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time=seq(0,2,0.1)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Kudo et al. 2007 high dose
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kudo_high <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Clungs")]/1000) #convert ug/kg to ug/g
  
  obs_kudo_high <- list(exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  score[1] <- AAFE(predictions = preds_kudo_high, observations = obs_kudo_high)
  
  ##########################
  #-------------------------
  # Kudo low
  #-------------------------
  ##########################
  # Set up simulations for the 2nd case, i.e. kudo (2007) low dose, tissues
  BW <- 0.29  # body weight (kg)
  sex <- "M" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[2]])
  params <- c(fixed_params[[2]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,2,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Kudo et al. 2007 low dose
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kudo_low <- solution[solution$time %in% unique(exp_data$time), c("Clungs")]
  
  preds_kudo_low<- as.data.frame(preds_kudo_low /1000) #convert ug/kg to ug/g
  
  
  obs_kudo_low <- list(exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  score[2] <- AAFE(predictions = preds_kudo_low, observations = obs_kudo_low)
  
  ##########################
  #-------------------------
  # Kim IV male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 3rd case, i.e. kim (2016) IV male tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[3]])
  params <- c(fixed_params[[3]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,288,1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of kim (2016) IV male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_IV_Mtissues <- solution[solution$time %in% unique(exp_data$time), c("Clungs")]
  
  preds_kim_IV_Mtissues<- as.data.frame(preds_kim_IV_Mtissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_IV_Mtissues <- list(exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  score[3] <- AAFE(predictions = preds_kim_IV_Mtissues, observations = obs_kim_IV_Mtissues)
  
  ##########################
  #-------------------------
  # Kim ORAL male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 4th case, i.e. kim (2016) ORAL male tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  sex <- "M" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[4]])
  params <- c(fixed_params[[4]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,288,1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of kim (2016) ORAL male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_OR_Mtissues <- solution[solution$time %in% unique(exp_data$time), c("Clungs")]
  
  preds_kim_OR_Mtissues<- as.data.frame(preds_kim_OR_Mtissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_OR_Mtissues <- list(exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  score[4] <- AAFE(predictions = preds_kim_OR_Mtissues, observations = obs_kim_OR_Mtissues)
  
  ##########################
  #-------------------------
  # Kim IV female tissues
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  sex <- "F" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[5]])
  params <- c(fixed_params[[5]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,24,1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of kim (2016) ORAL female tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_IV_Ftissues <- solution[solution$time %in% unique(exp_data$time), c("Clungs")]
  
  preds_kim_IV_Ftissues<- as.data.frame(preds_kim_IV_Ftissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_IV_Ftissues <- list(exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  score[5] <- AAFE(predictions = preds_kim_IV_Ftissues, observations = obs_kim_IV_Ftissues)
  
  ##########################
  #-------------------------
  # Kim ORAL female tissues
  #-------------------------
  ##########################
  # Set up simulations for the 6th case, i.e. kim (2016) IV female tissues
  BW <- 0.25 #kg, from Kim et al. 2018
  sex <- "F" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[6]])
  params <- c(fixed_params[[6]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,24,1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of kim (2016) ORAL female tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_OR_Ftissues <- solution[solution$time %in% unique(exp_data$time), c("Clungs")]
  
  preds_kim_OR_Ftissues<- as.data.frame(preds_kim_OR_Ftissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_OR_Ftissues <- list(exp_data[exp_data$Tissue == "Lung", "concentration"])
  
  score[6] <- AAFE(predictions = preds_kim_OR_Ftissues, observations = obs_kim_OR_Ftissues)
  
  
  ##########################
  #-------------------------
  # Gustafsson Oral male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 28st case, i.e. Gustafsson (2022) oral male tissues
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[7]])
  params <- c(fixed_params[[7]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time= seq(0,48,1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  #======================================df7=========================================================
  exp_data <- dataset$df7 # retrieve data of Gustafsson (2022) oral male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CLuELF","Clungtissue")
  
  preds_gus_OR_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_OR_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_gus_OR_Mtissues <- list( exp_data[exp_data$Tissue == "ELF", "concentration"],
                               exp_data[exp_data$Tissue == "Lung", "concentration"]) 
  
  score[7] <- AAFE(predictions = preds_gus_OR_Mtissues, observations = obs_gus_OR_Mtissues)
  
  
  ##########################
  #-------------------------
  # Gustafsson Inhalation male blood
  #-------------------------
  ##########################
  # Set up simulations for the 8th case, i.e. Gustafsson (2022) Inhalation male blood
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[8]])
  params <- c(fixed_params[[8]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(union(events$data$time, seq(0,48,1)))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df8=========================================================
  
  exp_data <- dataset$df8 # retrieve data of Gustafsson (2022)
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_gus_INH_Mblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_INH_Mblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_gus_INH_Mblood <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  score[8] <- AAFE(predictions = preds_gus_INH_Mblood, observations = obs_gus_INH_Mblood)
  
  
  ##########################
  #-------------------------
  # Gustafsson Inhalation male tissues
  #-------------------------
  ##########################
  
  #======================================df9=========================================================
  
  exp_data <- dataset$df9 # retrieve data of Gustafsson (2022) Inhalation male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CLuELF","Cliver","Clungtissue", "Ckidney")
  
  preds_gus_INH_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, "time"]
    
    preds_gus_INH_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_gus_INH_Mtissues <- list(exp_data[exp_data$Tissue == "ELF", "concentration"],
                               exp_data[exp_data$Tissue == "Liver", "concentration"],
                               exp_data[exp_data$Tissue == "Lung", "concentration"], 
                               exp_data[exp_data$Tissue == "Kidney", "concentration"]) 
  
  
  
  score[9] <- AAFE(predictions = preds_gus_INH_Mtissues, observations = obs_gus_INH_Mtissues)
  
  column_names <- c("CLuELF")
  preds_gus_INH_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, "time"]
    
    preds_gus_INH_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_gus_INH_Mtissues <- list(exp_data[exp_data$Tissue == "ELF", "concentration"]) 
  score[16] <- AAFE(predictions = preds_gus_INH_Mtissues, observations = obs_gus_INH_Mtissues)
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single low
  #-------------------------
  ##########################
  # Set up simulations for the 10th case, i.e. Hinderliter Inhalation male single low
  BW <- (0.311+0.195)/2  #kg, https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex = "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[9]])
  params <- c(fixed_params[[9]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(union(events$data$time, seq(0,30,1)))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df10=========================================================
  
  exp_data <- dataset$df10 # retrieve data of Hinderliter Inhalation male single low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_hind_INH_Mblood_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Mblood_low [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_hind_INH_Mblood_low <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[10] <- AAFE(predictions = preds_hind_INH_Mblood_low, observations = obs_hind_INH_Mblood_low)
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single medium
  #-------------------------
  ##########################
  # Set up simulations for the 11th case, i.e. Hinderliter Inhalation male single medium
  BW <- (0.311+0.195)/2  #kg, https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[10]])
  params <- c(fixed_params[[10]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(union(events$data$time, seq(0,30,1)))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df11=========================================================
  
  exp_data <- dataset$df11 # retrieve data of Hinderliter Inhalation male single medium
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_hind_INH_Mblood_medium <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Mblood_medium [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_hind_INH_Mblood_medium <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[11] <- AAFE(predictions = preds_hind_INH_Mblood_medium, observations = obs_hind_INH_Mblood_medium)
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single high
  #-------------------------
  ##########################
  # Set up simulations for the 12th case, i.e. Hinderliter Inhalation male single high
  BW <- (0.311+0.195)/2  #kg, https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[11]])
  params <- c(fixed_params[[11]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(c(0.5, union(events$data$time, seq(0,30,1))))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df12=========================================================
  
  exp_data <- dataset$df12 # retrieve data of Hinderliter Inhalation male single high
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_hind_INH_Mblood_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Mblood_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_hind_INH_Mblood_high <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[12] <- AAFE(predictions = preds_hind_INH_Mblood_high, observations = obs_hind_INH_Mblood_high)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single low
  #-------------------------
  ##########################
  # Set up simulations for the 13th case, i.e. Hinderliter Inhalation female single low
  BW <- (0.197+0.145)/2  #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[12]])
  params <- c(fixed_params[[12]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(union(events$data$time, seq(0,9,1)))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df13=========================================================
  
  exp_data <- dataset$df13 # retrieve data of Hinderliter Inhalation male single low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_hind_INH_Fblood_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Fblood_low [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_hind_INH_Fblood_low <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[13] <- AAFE(predictions = preds_hind_INH_Fblood_low, observations = obs_hind_INH_Fblood_low)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single medium
  #-------------------------
  ##########################
  # Set up simulations for the 14th case, i.e. Hinderliter Inhalation female single low
  BW <- (0.197+0.145)/2  #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[13]])
  params <- c(fixed_params[[13]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(c(0.5, union(events$data$time, seq(0,12,1))))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df14=========================================================
  
  exp_data <- dataset$df14 # retrieve data of Hinderliter Inhalation male single medium
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_hind_INH_Fblood_medium <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Fblood_medium [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_hind_INH_Fblood_medium <- list (exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[14] <- AAFE(predictions = preds_hind_INH_Fblood_medium, observations = obs_hind_INH_Fblood_medium)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single high
  #-------------------------
  ##########################
  # Set up simulations for the 15th case, i.e. Hinderliter Inhalation female single high
  BW <- (0.197+0.145)/2  #kg,  https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[14]])
  params <- c(fixed_params[[14]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= sort(c(0.5,union(events$data$time, seq(0,12,1))))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-02, atol = 1e-02))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df15=========================================================
  
  exp_data <- dataset$df15 # retrieve data of Hinderliter Inhalation male single high
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_hind_INH_Fblood_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_hind_INH_Fblood_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_hind_INH_Fblood_high <- list (exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
  
  score[15] <- AAFE(predictions = preds_hind_INH_Fblood_high, observations = obs_hind_INH_Fblood_high)
  
  
  
  ########################################################################################
  # Estimate final score
  
  final_score <- mean(score, na.rm = TRUE)
  return(final_score)
  
}

################################################################################

setwd('/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat')
#setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")

MW <- 414.07 #g/mol
source("Goodness-of-fit-metrics.R")

# Read data
kudo_high_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_high_kudo_2007.xlsx")
kudo_low_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_low_kudo_2007.xlsx")
kim_IV_Mtissues <- openxlsx::read.xlsx("Data/PFOA_male_tissues_IV_kim_2016.xlsx")
kim_OR_Mtissues <- openxlsx::read.xlsx("Data/PFOA_male_tissues_ORAL_kim_2016.xlsx")
kim_IV_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_IV_kim_2016.xlsx")
kim_OR_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_ORAL_kim_2016.xlsx")
gus_OR_Mblood <- openxlsx::read.xlsx("Data/Gustafsson 2022_PFOA_Plasma Male rats_Oral.xlsx")
gus_OR_Mtissues <- openxlsx::read.xlsx("Data/Gustafsson 2022_PFOA_Tissues Male rats_Oral.xlsx")
gus_OR_Mtissues$Tissue[1] <- "ELF"
gus_INH_Mblood <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Plasma Male rats_Inhalation.xlsx")
gus_INH_Mtissues <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Tissues Male rats_Inhalation.xlsx")
gus_INH_Mtissues$Tissue[1] <- "ELF"
hind_INH_Mblood_low <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_Low_dose.xlsx")
hind_INH_Mblood_medium <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_Medium_dose.xlsx")
hind_INH_Mblood_high <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_High_dose.xlsx")
hind_INH_Fblood_low <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_Low_dose.xlsx")
hind_INH_Fblood_medium <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_Medium_dose.xlsx")
hind_INH_Fblood_high <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_High_dose.xlsx")

#setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Inhalation_based_on_scenario15_PT/Training/AAFE/Inhalation_scenario24e")

setwd('/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Final_Inhalation_model/Training/AAFE/inhalation_permeability/fine_tune_raf')

dataset <- list("df1" = kudo_high_dose, "df2" = kudo_low_dose, "df3" = kim_IV_Mtissues, "df4" = kim_OR_Mtissues,
                "df5" = kim_IV_Ftissues, "df6" = kim_OR_Ftissues, "df7" = gus_OR_Mtissues, "df8" = gus_INH_Mblood,
                "df9" = gus_INH_Mtissues, "df10" = hind_INH_Mblood_low, "df11" = hind_INH_Mblood_medium,
                "df12" = hind_INH_Mblood_high,"df13" = hind_INH_Fblood_low,
                "df14" = hind_INH_Fblood_medium, "df15" = hind_INH_Fblood_high)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #""
              "xtol_rel" = 1e-07,
              "ftol_rel" = 1e-07,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0, 
              "maxeval" = 1000, 
              "print_level" = 1)

ClINFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
Awell = 9 #cm^2 (for a 35 mm culture dish)
Swell = 1.12 #cm^2
well_protein = 0.346 #mg protein
protein_per_well = (well_protein * Awell)/Swell #mg protein/well
Papp_Kimura = (ClINFT_unscaled*60*1e-06*1e3*protein_per_well)/Awell  #cm/h,at  pH = 6.0
Papp_RYU = 1.46e-6*3600 # cm/h, at pH = 7.4 from Ryu et al. (2024) [https://doi.org/10.1016/j.chemosphere.2024.142390]


# Create initial conditions (zero initialisation)
#Parameter names:
#  kabsUA, kCLEal, kCLEua, RAFlu

N_pars <- 2# Number of parameters to be fitted
fit <- c(1, 4.213369e-05)

lb	= c(1e-2, 2e-05)
ub = c(10, 6e-05)

fixed_params <- create_all_fixed_params()
# Run the optimization algorithm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             lb	= lb,
                             ub = ub,
                             opts = opts,
                             dataset = dataset,
                             fixed_params = fixed_params)

#estimated_params <- exp(optimizer$solution)
estimated_params <- (optimizer$solution)
save.image("fine_tune_raf.RData")



# Set up simulations for the 1st case, i.e. kudo (2007) high dose, tissues
BW <- 0.29  # body weight (kg)
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[1]])
params <- c(fixed_params[[1]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time=seq(0,2,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kudo_high <-  solution[, c("time","Clungs")]
###############################################################################
# Set up simulations for the 2nd case, i.e. kudo (2007) low dose, tissues
BW <- 0.29  # body weight (kg)
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[2]])
params <- c(fixed_params[[2]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,2,0.1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kudo_low <- solution[, c("time","Clungs")]
############################################################################
# Set up simulations for the 3rd case, i.e. kim (2016) IV male tissues
BW <- 0.25 #kg, from Kim et al. 2018
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[3]])
params <- c(fixed_params[[3]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_IV_Mtissues <-  solution[, c("time","Clungs")]

#########################################################################
# Set up simulations for the 4th case, i.e. kim (2016) ORAL male tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[4]])
params <- c(fixed_params[[4]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_OR_Mtissues <-  solution[, c("time","Clungs")]

###################################################################
# Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[5]])
params <- c(fixed_params[[5]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_IV_Ftissues <-  solution[, c("time","Clungs")]

#################################################################################
# Set up simulations for the 6th case, i.e. kim (2016) IV female tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[6]])
params <- c(fixed_params[[6]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_kim_OR_Ftissues <-  solution[, c("time","Clungs")]


##########################################################################################
# Set up simulations for the 7th case, i.e. Gustafsson (2022) oral male tissues
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[7]])
params <- c(fixed_params[[7]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0,48,0.2)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_gus_OR_Mtissues <-  solution[, c("time", "CLuELF","Cliver", "Clungtissue", "Ckidney")]

##########################################################################################
# Set up simulations for the 8th case, i.e. Gustafsson Inhalation male blood
BW <- 0.5125  #kg, from Gustafsson et al., 2022
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[8]])
params <- c(fixed_params[[8]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

sample_time=seq(0,48,0.2)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_gus_INH_Mblood <-  solution[, c("time", "Cplasma")]
preds_gus_INH_Mtissues <-  solution[, c("time", "CLuELF","Cliver", "Clungtissue", "Ckidney")]


#################################################################################
# Set up simulations for the 10th case, i.e. Hinderliter Inhalation male single low
BW <- (0.311+0.195)/2  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[9]])
params <- c(fixed_params[[9]], variable_params)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_hind_INH_Mblood_low <-  solution[, c("time", "Cplasma")]


#################################################################################
# Set up simulations for the 11th case, i.e. Hinderliter Inhalation male single medium
BW <- (0.311+0.195)/2  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[10]])
params <- c(fixed_params[[10]], variable_params)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_hind_INH_Mblood_medium <-  solution[, c("time", "Cplasma")]


#################################################################################
# Set up simulations for the 12th case, i.e. Hinderliter Inhalation male single high
BW <- (0.311+0.195)/2  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
depfr_head <- 0.2822
depfr_AF <- (0.1148+0.0177)
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[11]])
params <- c(fixed_params[[11]], variable_params)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_hind_INH_Mblood_high <-  solution[, c("time", "Cplasma")]

#################################################################################
# Set up simulations for the 13th case, i.e. Hinderliter Inhalation female single low
BW <- (0.197+0.145)/2  #kg
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[12]])
params <- c(fixed_params[[12]], variable_params)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,9,0.04)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_hind_INH_Fblood_low <-  solution[, c("time", "Cplasma")]

# Set up simulations for the 14th case, i.e. Hinderliter Inhalation female single medium
BW <- (0.197+0.145)/2  #
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[13]])
params <- c(fixed_params[[13]], variable_params)
inits <- create.inits(params)
events <- create.events(params)



sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-02, atol = 1e-02))

preds_hind_INH_Fblood_medium <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 15th case, i.e. Hinderliter Inhalation female single high
BW <- (0.197+0.145)/2  #kg,
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[14]])
params <- c(fixed_params[[14]], variable_params)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-3, atol = 1e-3))

preds_hind_INH_Fblood_high <-  solution[, c("time", "Cplasma")]


##########################################################################################################
#########################################################################################################

#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets
preds_kudo_high[,2:dim(preds_kudo_high)[2]] <- preds_kudo_high[,2:dim(preds_kudo_high)[2]] /1000 
preds_kudo_low[,2:dim(preds_kudo_low)[2]] <- preds_kudo_low[,2:dim(preds_kudo_low)[2]] /1000 
preds_kim_IV_Mtissues[,2:dim(preds_kim_IV_Mtissues)[2]] <- preds_kim_IV_Mtissues[,2:dim(preds_kim_IV_Mtissues)[2]] /1000
preds_kim_OR_Mtissues[,2:dim(preds_kim_OR_Mtissues)[2]] <- preds_kim_OR_Mtissues[,2:dim(preds_kim_OR_Mtissues)[2]] /1000
preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] <- preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] /1000
preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] <- preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] /1000
preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] <- preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] /1000
preds_gus_INH_Mblood[,2:dim(preds_gus_INH_Mblood)[2]] <- preds_gus_INH_Mblood[,2:dim(preds_gus_INH_Mblood)[2]] /1000
preds_gus_INH_Mtissues[,2:dim(preds_gus_INH_Mtissues)[2]] <- preds_gus_INH_Mtissues[,2:dim(preds_gus_INH_Mtissues)[2]] /1000
preds_hind_INH_Mblood_low[,2:dim(preds_hind_INH_Mblood_low)[2]] <- preds_hind_INH_Mblood_low[,2:dim(preds_hind_INH_Mblood_low)[2]] /1000
preds_hind_INH_Mblood_medium[,2:dim(preds_hind_INH_Mblood_medium)[2]] <- preds_hind_INH_Mblood_medium[,2:dim(preds_hind_INH_Mblood_medium)[2]] /1000
preds_hind_INH_Mblood_high[,2:dim(preds_hind_INH_Mblood_high)[2]] <- preds_hind_INH_Mblood_high[,2:dim(preds_hind_INH_Mblood_high)[2]] /1000
preds_hind_INH_Fblood_low[,2:dim(preds_hind_INH_Fblood_low)[2]] <- preds_hind_INH_Fblood_low[,2:dim(preds_hind_INH_Fblood_low)[2]] /1000
preds_hind_INH_Fblood_medium[,2:dim(preds_hind_INH_Fblood_medium)[2]] <- preds_hind_INH_Fblood_medium[,2:dim(preds_hind_INH_Fblood_medium)[2]] /1000
preds_hind_INH_Fblood_high[,2:dim(preds_hind_INH_Fblood_high)[2]] <- preds_hind_INH_Fblood_high[,2:dim(preds_hind_INH_Fblood_high)[2]] /1000


# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("PFOA concentration (" * mu* "g/g tissue)" ),
         x = "Time (hours)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("", values=cls,
                       guide = guide_legend(override.aes =
                                              list(shape = c(16,NA),
                                                   linetype = c(0,1))))+
    theme_light() + 
    theme(legend.position=c(1,1), 
          legend.justification=c(0, 1), 
          legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          axis.title=element_text(size=14),
          legend.text = element_text(size=14)
    )
  
}


# Convert Kudo High dose from long to wide format using reshape
experiment_inh_1 <- reshape(kudo_high_dose[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_1) <- c("Time",kudo_high_dose$Tissue )

# Convert Kudo Low dose from long to wide format using reshape
experiment_inh_2 <- reshape(kudo_low_dose[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_2) <- c("Time",kudo_low_dose$Tissue )

# Convert Kim IV Male tissues from long to wide format using reshape
experiment_inh_3 <- reshape(kim_IV_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_3) <- c("Time",kim_IV_Mtissues$Tissue )

# Convert Kim ORAL Male tissues from long to wide format using reshape
experiment_inh_4 <- reshape(kim_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_4) <- c("Time",kim_OR_Mtissues$Tissue )

# Convert Kim IV female tissues from long to wide format using reshape
experiment_inh_5 <- reshape(kim_IV_Ftissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_5) <- c("Time",kim_IV_Ftissues$Tissue )

# Convert Kim ORAL female tissues from long to wide format using reshape
experiment_inh_6 <- reshape(kim_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_6) <- c("Time",kim_OR_Ftissues$Tissue )


# Convert Gustafsson Oral male tissues from long to wide format using reshape
experiment_inh_7 <- reshape(gus_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_7) <- c("Time",unique(gus_OR_Mtissues$Tissue))

# Convert Gustafsson Inhalation male blood from long to wide format using reshape
experiment_inh_8 <- reshape(gus_INH_Mblood[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_8) <- c("Time",unique(gus_INH_Mblood$Tissue))


# Convert Gustafsson Inhalation male tissues from long to wide format using reshape
experiment_inh_9 <- reshape(gus_INH_Mtissues[c("Tissue" ,"Time_hours", 
                                               "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_9) <- c("Time",unique(gus_INH_Mtissues$Tissue))


# Convert Hinderliter Inhalation male single low from long to wide format using reshape
experiment_inh_10 <- reshape(hind_INH_Mblood_low[c("Tissue" ,"Time_hours", 
                                                   "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_10) <- c("Time",unique(hind_INH_Mblood_low$Tissue))


# Convert Hinderliter Inhalation male single medium from long to wide format using reshape
experiment_inh_11 <- reshape(hind_INH_Mblood_medium[c("Tissue" ,"Time_hours", 
                                                      "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_11) <- c("Time",unique(hind_INH_Mblood_medium$Tissue))


# Convert Hinderliter Inhalation male single high from long to wide format using reshape
experiment_inh_12 <- reshape(hind_INH_Mblood_high[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_12) <- c("Time",unique(hind_INH_Mblood_high$Tissue))


# Convert Hinderliter Inhalation female single low from long to wide format using reshape
experiment_inh_13 <- reshape(hind_INH_Fblood_low[c("Tissue" ,"Time_hours", 
                                                   "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_13) <- c("Time",unique(hind_INH_Fblood_low$Tissue))


# Convert Hinderliter Inhalation female single medium from long to wide format using reshape
experiment_inh_14 <- reshape(hind_INH_Fblood_medium[c("Tissue" ,"Time_hours", 
                                                      "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_14) <- c("Time",unique(hind_INH_Fblood_medium$Tissue))


# Convert Hinderliter Inhalation female single high from long to wide format using reshape
experiment_inh_15 <- reshape(hind_INH_Fblood_high[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_15) <- c("Time",unique(hind_INH_Fblood_high$Tissue))

# Put the experiments in a list
experiments <- list(experiment_inh_1 = experiment_inh_1, experiment_inh_2 = experiment_inh_2,
                    experiment_inh_3 = experiment_inh_3, experiment_inh_4 = experiment_inh_4,
                    experiment_inh_5 = experiment_inh_5, experiment_inh_6 = experiment_inh_6,
                    experiment_inh_7 = experiment_inh_7,
                    experiment_inh_8 = experiment_inh_8,experiment_inh_9 = experiment_inh_9, 
                    experiment_inh_10 = experiment_inh_10, experiment_inh_11 = experiment_inh_11,  
                    experiment_inh_12 = experiment_inh_12, experiment_inh_13 = experiment_inh_13, 
                    experiment_inh_14 = experiment_inh_14, experiment_inh_15 = experiment_inh_15)


# Rename predictions so that they share the same name as the names of the experimental data dataframe
colnames(preds_kudo_high) <- c( "Time", "Lung")
colnames(preds_kudo_low) <-  colnames(preds_kudo_high) 
colnames(preds_kim_IV_Mtissues) <- c( "Time","Lung")
colnames(preds_kim_OR_Mtissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_IV_Ftissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_OR_Ftissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_gus_OR_Mtissues) <- c ("Time", "ELF", "Liver", "Lung", "Kidney")
colnames(preds_gus_INH_Mblood) <- c ("Time", "Plasma")
colnames(preds_gus_INH_Mtissues) <- c ("Time", "ELF", "Liver", "Lung", "Kidney")
colnames(preds_hind_INH_Mblood_low) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Mblood_medium) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Mblood_high) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Fblood_low) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Fblood_medium) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Fblood_high) <- c ("Time", "Plasma")

# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_kudo_high,  predictions2 = preds_kudo_low, 
                    predictions3 = preds_kim_IV_Mtissues, predictions4 = preds_kim_OR_Mtissues,
                    predictions5 = preds_kim_IV_Ftissues, predictions6 = preds_kim_OR_Ftissues,
                    predictions7 = preds_gus_OR_Mtissues, predictions8 = preds_gus_INH_Mblood, 
                    predictions9 = preds_gus_INH_Mtissues, predictions10 = preds_hind_INH_Mblood_low,
                    predictions11 = preds_hind_INH_Mblood_medium, predictions12 = preds_hind_INH_Mblood_high, 
                    predictions13 = preds_hind_INH_Fblood_low,predictions14 = preds_hind_INH_Fblood_medium,
                    predictions15 = preds_hind_INH_Fblood_high)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment )
  })
  if(length(compartments) == 1){
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                               common.legend = TRUE, legend = "right"))
    
  }else{
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                               common.legend = TRUE, legend = "right"))
  }
  
  
  plot.margin=unit(c(0,0,0,0), "pt")
  
  
  # Save the plot with dynamically adjusted dimensions
  ggsave(paste0("experiment", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}




