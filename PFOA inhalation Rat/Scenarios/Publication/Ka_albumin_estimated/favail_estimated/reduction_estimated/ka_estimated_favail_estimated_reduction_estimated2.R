# Check notes for description of scenario
#female increase in BW: 3.5 g/d
#Male increase in BW: 5.9 g.d

library(deSolve)

create_variable_params <- function(BW,sex,  estimated_params, fixed_params){
  
  if(sex == "M"){
    RAFOatp_k <- estimated_params[1]
    f_fabp_avail  <- estimated_params[11]
    
  }else{
    RAFOatp_k <-estimated_params[2]
    f_fabp_avail  <- estimated_params[12]
    
  }
  RAFOat3 <- estimated_params[3]
  
  CL_int <- estimated_params[4] #uL/min/million hepatocytes
  HEPGL <- 104  #million hepatocytes/gram of rat liver (Fattah et al., 2016. [doi: 10.1124/dmd.115.066381] )
  # Scaled hepatobiliary clearance
  CL_hepatobiliary <- CL_int*1e-6*HEPGL*(fixed_params$MKi_drained*1000) *60 #L/h
  
  RAFOatp_l <- estimated_params[5]
  RAFUrat <- RAFOatp_k
  RAFOat1 <- RAFOat3
  RAFOatp2_l <- RAFOatp_l
  RAFOatp_lu_ap <-estimated_params[6]
  RAFOatp_lu_bas <- RAFOatp_lu_ap
  RAFNtcp <- RAFOatp_l
  RAFOatp2_Int <- 6.799068e-07
  
  # The parameters below were part of tests involving hypothetical transporters in the kidneys.
  # In the final model these transporters are switched off.
  VmK_api <- 0
  VmK_baso <- 0
  KmK_baso <- 1e20
  KmK_api <-   5e4
  
  KLfabp <-  1.2e5#[L/mol]. From Sheng et al. (2018) [doi:10.1007/s00204-017-2055-1]
  Ka <- estimated_params[9]#[L/mol]. From Rue et al. (2024)
  Ka2u <- 1e03 #[L/mol], from Han et al. (2004) [doi:10.1081/DCT-200039725]
  
  CLfeces_unscaled <- estimated_params[7]#in L/h/BW^(-0.25), scaling similar to Loccisano et al. (2012)
  CLfeces <- CLfeces_unscaled*BW^(-0.25) 
  
  f_alb_avail<-  estimated_params[10]
  
  f_a2u_avail <- 1
  
  koff_alb <-  1e-2*3600#[1/h]
  koff_fabp <-  1e-2*3600#[1/h]
  koff_a2u <- koff_alb
  
  reduction_factor <- estimated_params[8]
  
  f_tubular <- 0.8
  f_PTC_prot_to_tub_prot <- 0.6939
  
  # For the section below, we assume that all measurement were made relative to rinsed organs that 
  # contained only capillary blood weight. Consequently, we multiply by the mass of the drained organ
  
  MW = 414.07 #g/mol, PFOA molecular weight
  intestine_protein <-109 #mg/g intestine. from van de Kerkhof (2007) [https://research.rug.nl/en/publications/drug-metabolism-in-human-and-rat-intestine-an-in-vitro-approach]
  intestine_protein_total <- intestine_protein*(1000*fixed_params$MIn_drained) #mg
  
  #Kidney
  kidney_protein_per_rat <- 1000*(0.218+0.225+0.212)/3 #mg of protein per rat  (Addis 1936)
  rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
  rat_kidney_weight_addis <- rat_weight_addis*0.0073 # kidney fraction to BW, Brown (1997)
  kidney_protein_per_gram <- kidney_protein_per_rat/rat_kidney_weight_addis #mg of protein/g kidney
  kidney_protein_total <- kidney_protein_per_gram* (1000*fixed_params$MKi_drained) #mg
  PTC_protein <- f_tubular*f_PTC_prot_to_tub_prot*kidney_protein_total #mg
  
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
  liver_protein_total <- liver_protein_per_gram* (1000*fixed_params$MLi_drained) #mg
  
  
  #oatp1-liver
  VmL_Oatp_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
  VmL_Oatp_scaled = 60*VmL_Oatp_in_vitro*MW*liver_protein_total/1000   #physiologically scaled to in vivo, ug/h
  VmL_Oatp = VmL_Oatp_scaled*RAFOatp_l #in vivo value, in  ug/h
  KmL_Oatp = KmK_Oatp #same as kidney
  
  #oatp2b1-liver
  VmL_Oatp2_in_vitro= 1493e-3 #nmol/mg protein/min  (Lin et al. 2023)
  #physiologically scaled to in vivo
  VmL_Oatp2_scaled = 60*VmL_Oatp2_in_vitro*MW*liver_protein_total/1000  #ug/h
  VmL_Oatp2 = VmL_Oatp2_scaled*RAFOatp2_l #in vivo value, in  ug/h
  KmL_Oatp2 = 148.68*MW #umol/L (Lin et al. 2023) --> ug/L
  
  #Ntcp liver
  VmL_Ntcp_in_vitro= 3#nmol/mg protein/min   Ruggiero et al. 2021
  #physiologically scaled to in vivo
  VmL_Ntcp_scaled = 60*VmL_Ntcp_in_vitro*MW*liver_protein_total/1000 # ug/h 
  VmL_Ntcp = VmL_Ntcp_scaled*RAFNtcp #in vivo value, in  ug/h
  KmL_Ntcp= 20 * MW #umol/L, Ruggiero et al. 2021 --> ug/L
  
  #Lung
  lung_protein_per_gram <- 134 #  mg/g tissue, Figure 2, [doi: 10.1007/s00580-021-03242-z] 
  lung_protein_total <- lung_protein_per_gram*(fixed_params$MLu_drained*1000) #mg
  
  #oatp-lung-ap (from ALF to tissue)
  VmLu_Oatp_ap_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
  VmLu_Oatp_ap_scaled = 60*VmLu_Oatp_ap_in_vitro*MW*lung_protein_total/1000   #physiologically scaled to in vivo, ug/h
  VmLu_Oatp_ap = VmLu_Oatp_ap_scaled*RAFOatp_lu_ap #in vivo value, in  ug/h
  KmLu_Oatp_ap = KmK_Oatp #same as kidney
  
  #oatp-lung-bas (from IS to tissue)
  VmLu_Oatp_bas_in_vitro= 9.3 #nmol/mg protein/min  (Weaver et al. 2010)
  VmLu_Oatp_bas_scaled = 60*VmLu_Oatp_bas_in_vitro*MW*lung_protein_total/1000   #physiologically scaled to in vivo, ug/h
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
  Mr_a2u <- 15500 #g/mol, from  Hai and Kizilbash (2013) [doi:10.6026/97320630009145]
  Mr_fabp <- 12000 #g/mol, from Ockner and  Manning (1976) [doi:10.1172/JCI108510] 
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
  #values from Kawai et al.(1994) 
  IPR_K = 0.5
  IPR_L = 0.5
  IPR_ST = 0.9 #same as gut
  IPR_IN = 0.9 #same as gut
  IPR_M = 0.6
  IPR_A = 0.5
  IPR_Lu = 0.5
  IPR_Sp = 0.5
  IPR_H = 0.5
  IPR_SK = 1.0
  IPR_Br = 0.5
  IPR_Go = 0.5 #assumption
  IPR_Bo = 1.0 
  IPR_R = 0.5 # assumption
  
  Hct <- 0.46 # From Davies and Morris (1993) [doi.10.1023/A:1018943613122].
  # Worley et al. (2015) used the same value [doi:10.1016/j.taap.2015.10.017]. 
  # Ehresman et al. (2007) found no binding to RBCs and a concentration of plasma:blood 2:1 [doi:10.1016/j.envres.2006.06.008]
  # Filho et al. (2017) provide age-based values [10.1080/13685538.2017.1350156]
  
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
  CalbLuAF_init <- (10/100)* (CalbB_init/(1-Hct))  #based on Woods et al. 2015 statement [doi: 10.1016/j.jconrel.2015.05.269]
 
   #Total cytosolic protein in liver tissue
  cytosolic_protein_L <-105*(1000*fixed_params$MLi_drained) #[mg], from Krähenbühl and Krähenbühl (2023) [doi: 10.3390/ijms24054365]
  cytosolic_protein_K <- 53.3*(1000*fixed_params$Mcortex) #[mg], assuming same value as humans. From Scotcher et al. (2017) [doi: 10.1124/dmd.117.075242]
  if (sex == "M"){
    #Alpha2mu-globulin concentration in male kidney tissue
    a2u_globulin_k = 22.7*(kidney_protein_total/1000)/fixed_params$VKi_drained #mg/L. From Stonard et al. (1986) [doi: 10.1016/0300-483X(86)90197-6]
    Ca2uKT_init <- f_a2u_avail*(a2u_globulin_k*1e-3/Mr_a2u) #[mol/L]
    
    #LFABP concentration in male liver tissue
    L_FABP_L = 30.1*1e-3*cytosolic_protein_L/fixed_params$VLi_drained #[mg/L]. From  Ockner and  Manning (1976) [doi: 10.1016/S0021-9258(18)34463-6]
    CFabpLT_init = f_fabp_avail*(L_FABP_L*1e-3/Mr_fabp) #[mol/L]
    
    #LFABP concentration in male kidney tissue
    L_FABP_K = 0.078*1e-9*cytosolic_protein_K/fixed_params$VKi_drained #[mol/L]. From  Maatman et al. (1992) [doi: 10.1042/bj2880285]
    CFabpKT_init <- f_fabp_avail*L_FABP_K  #[mol/L]
    
  }else if(sex == "F"){
    #Alpha2mu-globulin concentration in female kidney tissue
    Ca2uKT_init <- 0 #mol/L
    
    #LFABP concentration in female liver tissue
    L_FABP_L = 51.7*1e-3*cytosolic_protein_L/fixed_params$VLi_drained #[mg/L]. From  Ockner and  Manning (1976) [doi: 10.1016/S0021-9258(18)34463-6]
    CFabpLT_init = f_fabp_avail*(L_FABP_L*1e-3/Mr_fabp) #[mol/L]
    
    #LFABP concentration in female kidney tissue
    L_FABP_K = 0.051*1e-9*cytosolic_protein_K/fixed_params$VKi_drained #[mol/L]. From  Maatman et al. (1992) [doi: 10.1042/bj2880285]
    CFabpKT_init <- f_fabp_avail*L_FABP_K  #[mol/L]
    
  }
  
  kon_alb <- Ka * koff_alb #1/h/M
  kon_a2u <- Ka2u * koff_a2u#1/h/M
  kon_fabp <- KLfabp * koff_fabp #1/h/M
  
  # Following the calculations  of Lin et al. (2023) for Caco-2 cells
  ClINFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017 [doi: 10.1016/j.toxlet.2017.05.012]
  Awell = 9 #cm^2 (for a 35 mm culture dish)
  Swell = 1.12 #cm^2
  well_protein = 0.346 #mg protein
  protein_per_well = (well_protein * Awell)/Swell #mg protein/well
  RAF_papp <- 1 
  Peff_KIMURA = RAF_papp*(ClINFT_unscaled*60*1e-06*1e3*protein_per_well)/(Awell *2) #cm/h,at  pH = 6.0
  Papp_RYU <- 1.46e-6*3600 # cm/h, at pH = 7.4 from Ryu et al. (2024) [doi: 10.1016/j.chemosphere.2024.142390]
  # The permeability rate of Kimura et al. (2017)is 30 times greater than the one from Ryu et al. (2024)
  # For endothelial and cellular permeability we use the Ryu et al. (2024) value
  Peff_monolayer <- Papp_RYU #cm/h
  
  k_gut_in = ( (2*Peff_monolayer/100) * fixed_params$AINL)*1000 #L/h
  k_gut_out = ( (2*Peff_monolayer/100) * fixed_params$AINL)*1000 #L/h
  
  #Stomach
  # For identifiability reasons we assume that absorption takes place only through the intestines
  kabST <- 0#(kabs_st* fixed_params$ASTL)*1000 #L/h
 
  
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
  
  # Physiologic upper limit of pore size from Sarin et al. (2010) [10.1186/2040-2384-2-14]
  pore_diameters <- c(
    Ki = 9,
    Li = 135,     # use 135 as average for rat liver
    St = 9,       # midpoint of 6-12 nm for fenestrated with diaphragm
    In = 9,
    Mu = 5,
    Ad = 5,
    Re = 5,
    Lu = 5,
    Sp = 5000,    # 5 microns = 5000 nm
    Ht = 5,
    Br = 1,       # tight junction
    Go = 9,
    Sk = 9,
    Bo = 5
  )
  
  dif <- 5.46e-6 #[cm^2/s], diffusion coefficient of PFOA in water, from  Gauthier et al. (2024) [doi: 10.1021/acsestwater.4c00631]
  kboltzman <- 1.38e-23 #J/kelvin
  Temp <- 37+273 #kelvin 
  dyn_visc <- 6.9e-4 #Pa*s, Dynamic viscosity of water at 37 o C
  # from Stokes–Einstein:
  R_H  <-  kboltzman * Temp/ (6*pi*dif*1e-4*dyn_visc)*1e9 #nm, hydrodynamic radius of PFOA 

  #Renkin equation (Renkin, (1954), [PMC2147404])
  #Diffusion reduction due to steric hindrance at the entrance to the pores and frictional resistance within the pores 
  lambda <- R_H/(pore_diameters/2)
  lambda["Br"] <- 1 #Upper limit
  renkin_reduction <- (1-lambda)^2 * (1 - 2.104*lambda + 2.09*lambda^3 - 0.95*lambda^5)
  wall_width <-  0.5e-6*100 #[m] -->[cm], capillary wall width. From Ashoor et al. (2018) [10.1016/j.rpor.2018.09.007]. Similar value from Sosula (1974) [10.1016/0026-2862(74)90011-9]
  basement_membrane <-  0.5e-6*100 #assumption
  
  Pgap <-reduction_factor*dif*renkin_reduction/ (wall_width+basement_membrane) #cm/s
  # In organs with sinusoid there is no basement membrane at gaps,
  # so we use basemenent membrane thickness = 0
  Pgap_sinusoidal <-  reduction_factor*dif*renkin_reduction/ wall_width #cm/s
  
  #fraction of gaps in capillary surface
  f_kidney <- 0.3 #Bulger et al. (1983) [10.1172/jci110950] Mou et al. (2024) [10.3390/ijms25169107]
  f_liver <- 0.08 #  Simon-Santamaria et al. (2010)[10.1093/gerona/glq108] (Antwi et al. (2023) give a range 2-20% [10.1371/journal.pone.0293526] )
  f_spleen <- 0.08 #assumption
  f_intestine <- 0.095 #Simionescu et al. (1974) [doi: 10.1083/jcb.60.1.128]
  f_non_fenestrated <- 0.0048 # Clough & Michel (1988) [doi: 10.1113/jphysiol.1988.sp017348]
  
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
  PparaKi <- Pgap["Ki"]*3600*10*f_kidney#mm/h
  PparaLi <- Pgap_sinusoidal["Li"]*3600*10 * f_liver #mm/h
  PparaSt <- Pgap["St"]*3600*10 * f_non_fenestrated #mm/h
  PparaIn <- Pgap["In"]*3600*10 * f_intestine #mm/h
  PparaMu <- Pgap["Mu"]*3600*10 * f_non_fenestrated #mm/h
  PparaAd <- Pgap["Ad"]*3600*10 * f_non_fenestrated #mm/h
  PparaRe <- Pgap["Re"]*3600*10 * f_non_fenestrated #mm/h
  PparaLu <- Pgap["Lu"]*3600*10 * f_non_fenestrated #mm/h
  PparaSp <- Pgap_sinusoidal["Sp"]*3600*10 * f_spleen  #mm/h
  PparaHt <- Pgap["Ht"]*3600*10 * f_non_fenestrated #mm/h
  PparaBr <- 0
  PparaGo <- Pgap["Go"]*3600*10 * f_non_fenestrated #mm/h
  PparaSk <- Pgap["Sk"]*3600*10 * f_non_fenestrated #mm/h
  PparaBo <- Pgap["Bo"]*3600*10  * f_non_fenestrated#mm/h
  
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
    
    'k_gut_in' = k_gut_in, 'k_gut_out' = k_gut_out,'kabST'=kabST,
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
    "PparaSk" = PparaSk,"PparaBo" = PparaBo,
    
    "MW" = MW, "Hct" = Hct
    
  ))
}  

create_fixed_params <- function(user.input){
  with(as.list(user.input),{
    
    # Total Volumes, Residual blood, Capillary blood, Interstitial fluid and Cellular volumes of organs
    # All volume values are in [L], densities are in kg/L
    
    # The fraction of organ mass relative to body weight in Brown (1997) refers to drained organ mass.
    # We assume that this means that the blood from large vessels was removed, but capillary blood in the organ still remains.
    
    # We distinguish  blood into total residual blood, which is the blood remaining in both large vessels and capillaries 
    # in nonbled animals and capillary residual blood, which is the blood remaining in capillaries in bled animals. We take 
    # the residual blood fraction from Brown when possible (average of multiple studies in most cases). 
    
    # The capillary residual blood was drawn from Bernareggi and  Rowland (1991) [doi: 10.1007/BF01062191], which was measured after bleeding the animals. 
    # Thus, we assume that the reported values were calculated relative to bled tissue that contains residual blood only in the capillaries and not in the large vessels.
   
    # We assume that the fraction of interstitial space from Kawai et al. (1994) is relative to the parenchyma organ volume.
    # We base this on Tsuji et al. , which is one of the references used. In the latter reference, the 
    # ISF is calculated on a blood-free organ weight basis.

    #Kidney
    PVKi <- 7.3e-3 #fraction of drained kidney mass relative to body weight. from Brown et al. (1997)
    d_kidney <- 1.05 #density of kidneys [kg/L]
    VKi_drained <- PVKi * BW / d_kidney #volume of  drained kidneys
    MKi_drained <- VKi_drained * d_kidney
    fcap_ki <- 0.046 #blood remaining in rat organs after bleeding, in [ml/g drained tissue].from Bernareggi and  Rowland (1991)
    VKiB_cap <- fcap_ki * MKi_drained  #volume of kidney blood in capillaries 
    PVKiB <- 0.16 #fraction of residual blood in kidneys. from Brown et al. (1997)
    VKi_tot <- (VKi_drained-VKiB_cap)/(1-PVKiB) #total kidney volume
    MK <- (VKi_drained - VKiB_cap) *d_kidney #total kidney parenchyma mass
    VKiB <- PVKiB * VKi_tot #residual blood volume  in kidneys
    fcap_to_blood_ki<-  VKiB_cap/VKiB #fraction of capillary to total organ blood
    PVKF <- 0.13 # fraction of interstitial fluid in kidneys.From Larson et al. (1984) [doi: 10.1111/j.1748-1716.1984.tb00137.x]
    VKiF <- PVKF *(VKi_drained - VKiB_cap) #kidney interstitial fluid volume 
    VKiB_lv <-  VKiB - VKiB_cap #volume of kidney blood in large vessels (lv)
      
    # Here we assume that the length of each segment of the renal tubule is 
    # analogous to the body weight of the kidney in order to obtain
    # volumes and surface areas as a linear function of BW
    VKi_ref <- (2*2.9/ d_kidney)/1000 # [L]. Gilmer et al. (2018) handled data from Sperber, 1944 (Thesis: “Studies on the Mammalian Kidney").
    # In Sperber, p.303, it says that the kidneys of an old female white rat were used. We assume that the 
    # weight provided refers to a single kidney and not both.
    
    #Total Length of renal tubules, from Gilmer et al. (2018) [doi:10.1152/ajprenal.00219.2018]
    LPT_slp <- (VKi_tot/VKi_ref)*(9886*1e-6)*26980*2 #[m], Proximal tubule (short-looped nephrons)
    LPT_llp <- (VKi_tot/VKi_ref)*(11077*1e-6 )*11020*2#[m], Proximal tubule (long-looped nephrons)
    LPT <- LPT_slp+LPT_llp
    LTDL_slp <- (VKi_tot/VKi_ref)*(1500*1e-6)*26980*2 #[m], Thin descending limb (short-looped nephrons)
    LTDL_llp <- (VKi_tot/VKi_ref)*(6200*1e-6)*11020*2 #[m], Thin descending limb (long-looped nephrons)
    LThinAL <-(VKi_tot/VKi_ref)* (4700*1e-6)*11020*2 #[m], Thin ascending limb (long-looped nephrons)
    LThickAL_cort <- (VKi_tot/VKi_ref)*(1450*1e-6)*38000*2 #[m], Cortical thick ascending limb
    LThickAL_med <- (VKi_tot/VKi_ref)*(2100*1e-6)*38000*2 #[m], Medullary thick ascending limb 
    LDT_sup <- (VKi_tot/VKi_ref)*(1452*1e-6)*26980*2 #[m], Distal tubule* (superficial)
    LDT_deep <- (VKi_tot/VKi_ref)*(1650*1e-6)*11020*2 #[m], Distal tubule* (deep) 
    LDT <- LDT_sup + LDT_deep
    Lduct_cort <- (VKi_tot/VKi_ref)*(2900*1e-6)*6000*2 #[m], Cortical collecting duct 
    Lduct_med <-  (VKi_tot/VKi_ref)*(2100*1e-6)*6000*2 #[m], Outer medullary collecting duct 
    
    # Renal tubule radius per compartment
    RPT_slp <- (22.9/2)*1e-6 #[m], Proximal tubule (short-looped nephrons)
    RPT_llp <-(23.1/2)*1e-6#[m], Proximal tubule (long-looped nephrons)
    RTDL_slp <- (15/2)*1e-6 #[m], Thin descending limb (short-looped nephrons)
    RTDL_llp <- (15/2)*1e-6 #[m], Thin descending limb (long-looped nephrons)
    RThinAL <- (15/2)*1e-6 #[m], #Thin ascending limb (long-looped nephrons)
    RThickAL_cort <- (25.4/2)*1e-6 #[m], Cortical thick ascending limb
    RThickAL_med <- (29/2)*1e-6 #[m], Medullary thick ascending limb 
    RDT_sup <- (39/2)*1e-6 #[m], Distal tubule* (superficial)
    RDT_deep <- (43/2)*1e-6 #[m], Distal tubule* (deep) 
    Rduct_cort <- (24/2)*1e-6 #[m], Cortical collecting duct 
    Rduct_med <-  (24/2)*1e-6 #[m], Outer medullary collecting duct 
    
    #volumes of filtrate compartments,[m^3]*1e3 --> [L]
    VPT <- 2*pi*((RPT_slp^2)*LPT_slp + (RPT_llp^2)*LPT_llp)*1e3 # Proximal tubule (short- and long-looped)
    VTDL <- 2*pi*( (RTDL_slp^2)*LTDL_slp+(RTDL_llp^2)*LTDL_llp)*1e3 #Thin descending limb (short- and long-looped)
    VThinAL <-  2*pi*((RThinAL^2)*LThinAL)*1e3 #Thin ascending limb (long-looped)
    VThickAL <- 2*pi*((RThickAL_cort^2)*LThickAL_cort+(RThickAL_med^2)*LThickAL_med)*1e3 #Thick ascending limb (Cortical and Medullary)
    VDAL <- VTDL + VThinAL+ VThickAL # Loop of Henle
    VDT <- 2*pi*((RDT_sup^2)*LDT_sup+(RDT_deep^2)*LDT_deep)*1e3#Distal tubule (superficial+deep)
    VCD <- 2*pi*((Rduct_cort^2)*Lduct_cort+(Rduct_med^2)*Lduct_med)*1e3 #collecting duct (Cortical+Outer)
    VFil <-  VPT+VDAL+VDT+VCD #L
    # We hypothesize that when kidney is weighted the renal tubule content remains inside 
    VKT <- VKi_drained - VKiF - VFil - VKiB_cap#kidney tissue volume kg=L
    
    #Wang et al. (2024) [doi:10.1038/s41598-024-53270-2] state that "More than 80% of 
    # renal cortical cells are tubular epithelial cells". By assuming that this percentage holds for
    # the medulla region as well, we estimate that the total fraction of tubular cells to total kidney cells is 0.8
    f_tubular <- 0.8
    
    # Clark et al. (2019) [doi: 10.1016/j.kint.2018.11.028] 
    # state that "Proximal tubule cells account for roughly 52% of the estimated
    # 206 million tubule epithelial cells per kidney. However, they account for approximately 69% of 
    # total tubule protein mass by virtue of their large size compared with other renal tubule cells".
    # Thus, here we assume that the volume is analogous to the protein mass, and from Clark et al. (2019), we get:
    f_PTC_prot_to_tub_prot <- 0.6939
    f_DALC_prot_to_tub_prot <- (0.0339 + 0.1156)
    f_DTC_prot_to_tub_prot <- 0.08
    f_CDC_prot_to_tub_prot <- 0.0766
    
    VPTC <- f_tubular*f_PTC_prot_to_tub_prot*VKT #for comparison, for 300g male rat we get 0.76mL and Worley and Fisher(2015) estimate 0.30 mL
    VDALC <- f_tubular*f_DALC_prot_to_tub_prot*VKT
    VDTC <- f_tubular*f_DTC_prot_to_tub_prot*VKT 
    VCDC <- f_tubular*f_CDC_prot_to_tub_prot*VKT 
    VKTrest <- (1-f_tubular)*VKT
    
    VBladder <- 0.001 #[L], from Andersson & Arner (2004) [doi:10.1152/physrev.00038.2003]
    
    #Liver
    PVLi <- 3.66e-2 # fraction of drained liver mass relative to body weight. From Brown et al. (1997)
    d_liver <- 1.05 #  Assumption
    VLi_drained <- PVLi * BW/d_liver
    MLi_drained <- VLi_drained * d_liver
    fcap_li <- 0.057 # from Bernareggi and  Rowland (1991)
    VLiB_cap <- fcap_li * MLi_drained 
    PVLiB <- 0.21  #from Brown et al. 1997
    VLi_tot <- (VLi_drained-VLiB_cap)/(1-PVLiB)
    MLi <- (VLi_drained - VLiB_cap) * d_liver
    VLiB <- PVLiB * VLi_tot
    fcap_to_blood_li<-  VLiB_cap/VLiB
    PVLiF <- 0.163  # Kawai et al. (1994)
    VLiF <- PVLiF * (VLi_drained - VLiB_cap)
    PVBile <- 47e-3/200 # fraction of biliary system of the normal rat liver relative to body weight. FromT  Masyuk et al. (2001) [doi: 10.1016/S0002-9440(10)64679-2]
    VBile <- PVBile * BW #L
    VLiT <- VLi_drained - VLiF - VLiB_cap - VBile
    VLiB_lv <-  VLiB -VLiB_cap
    
 
    # Intestine (small and large)
    PVIn <- 2.24e-2 # from Brown et al. (1997)
    d_int <- 1.04 # from Brown et al. (1997)
    VIn_drained <- PVIn * BW / d_int
    MIn_drained <- VIn_drained * d_int
    fcap_in <- 0.010 # from Bernareggi and  Rowland (1991)
    VInB_cap <- fcap_in * MIn_drained
    PVInB <- 0.024   # from  Kawai et al. (1994)
    VIn_tot <- (VIn_drained - VInB_cap) / (1 - PVInB)
    MIn <- (VIn_drained - VInB_cap) * d_int
    VInB <- PVInB * VIn_tot
    fcap_to_blood_in<-  VInB_cap/VInB 
    PVInF <- 0.094 #from Kawai et al. (1994)
    VInF <- PVInF * (VIn_drained - VInB_cap)
    VInT <- VIn_drained - VInF - VInB_cap
    VInB_lv <- VInB - VInB_cap
    
 
    # Stomach
    PVSt <- 0.46e-2 #from Brown et al. (1997)
    d_stomach <- 1.05 #from Brown et al. (1997)
    VSt_drained <- PVSt * BW / d_stomach
    MSt_drained <- VSt_drained * d_stomach
    fcap_st <- 0.011 # from Bernareggi and  Rowland (1991)
    VStB_cap <- fcap_st * MSt_drained
    PVStB <- 0.032 #from  Kawai et al. (1994)
    VSt_tot <- (VSt_drained - VStB_cap) / (1 - PVStB)
    MSt <-  (VSt_drained - VStB_cap)* d_stomach
    VStB <- PVStB * VSt_tot
    fcap_to_blood_st<-  VStB_cap/VStB 
    PVStF <- 0.10 #from  Kawai et al. (1994)
    VStF <- PVStF * (VSt_drained - VStB_cap)
    VStT <- VSt_drained - VStF - VStB_cap
    VStB_lv <- VStB - VStB_cap
    
    #Stomach and intestine lumen
    # The values below are from Conell et al. (2008) [doi: 10.1211/jpp.60.1.0008] and represent measured average water content of fed and fasted state from Figure 4.
    PVStL <- 0.53/175 #mL/g BW
    VStL <- PVStL * BW 
    PVInL <- (5.4- 0.53)/175 # mL/g BW
    VInL <- PVInL * BW 
    
    # Muscle
    PVMu <- 40.43e-2 #from Brown et al. (1997)
    d_muscle <- 1.04 #from Brown et al. (1997)
    VMu_drained <- PVMu * BW / d_muscle
    fcap_mu <- 0.004 # from Bernareggi and  Rowland (1991)
    VMuB_cap <- fcap_mu * VMu_drained* d_muscle
    PVMuB <- 0.04 #from  Brown et al. (1997)
    VMu_tot <- (VMu_drained - VMuB_cap) / (1 - PVMuB)
    MMu <- (VMu_drained - VMuB_cap)  * d_muscle
    VMuB <- PVMuB * VMu_tot
    fcap_to_blood_mu <-  VMuB_cap/VMuB 
    PVMuF <- 0.120 #from  Kawai et al. (1994)
    VMuF <- PVMuF * (VMu_drained - VMuB_cap)
    VMuT <- VMu_drained - VMuF - VMuB_cap
    VMuB_lv <- VMuB - VMuB_cap
    
    # Adipose
    PVAd <- mean(c(5.5,7))*1e-2 #from Brown et al. (1997)
    d_adipose <- 0.92    # from Brown et al. (1997)
    VAd_drained <- PVAd * BW / d_adipose
    fcap_ad <- 0.005 #from Bernareggi and  Rowland (1991)
    VAdB_cap <- fcap_ad * VAd_drained * d_adipose
    PVAdB <- 0.01 #from  Kawai et al. (1994)
    VAd_tot <- (VAd_drained - VAdB_cap) / (1 - PVAdB)
    MAd <- (VAd_drained - VAdB_cap) * d_adipose
    VAdB <- PVAdB * VAd_tot
    fcap_to_blood_ad<-  VAdB_cap/VAdB 
    PVAdF <- 0.135 #from  Kawai et al. (1994)
    VAdF <- PVAdF * (VAd_drained - VAdB_cap)
    VAdT <- VAd_drained - VAdF - VAdB_cap
    VAdB_lv <- VAdB - VAdB_cap
    
    # Lung
    PVLu <- 0.50e-2 #from Brown et al. (1997)
    d_lung <- 1.05 #from Brown et al. (1997)    
    VLu_drained <- PVLu * BW / d_lung
    MLu_drained <- VLu_drained * d_lung
    fcap_lu <- 0.175 #from  Bernareggi and  Rowland (1991)
    VLuB_cap <- fcap_lu *MLu_drained
    PVLuB <- 0.36  #from  Brown et al. (1997)
    VLu_tot <- (VLu_drained - VLuB_cap) / (1 - PVLuB)
    MLu <- (VLu_drained - VLuB_cap) * d_lung
    VLuB <- PVLuB * VLu_tot
    fcap_to_blood_lu<-  VLuB_cap/VLuB 
    PVLuF <- 0.188 #from  Kawai et al. (1994)
    VLuF <- PVLuF * (VLu_drained - VLuB_cap)
    VLuAF <- ((0.014+0.018)/2)/((335 + 290)/2) * BW # from Boger et al, 2023 [doi: 10.1016/j.xphs.2023.01.001] 
    VLuT <- VLu_drained - VLuF - VLuAF - VLuB_cap
    VLuB_lv <- VLuB - VLuB_cap
    
    #  For the spleen assume that all of the blood is capillary blood and thus we use the value of total blood from Brown et al. (1997),
    # as capillary blood, because the value of residual blood after bleeding from  Bernareggi and  Rowland (1991) was higher.
    PVSp <- 0.2e-2 # from Brown et al. (1997)
    d_spleen <- 1.05  # from Brown et al. (1997)
    VSp_drained <- PVSp * BW / d_spleen
    PVSpB <- 0.22 #from Brown et al. (1997)
    fcap_sp <- PVSpB
    VSpB_cap <- fcap_sp * VSp_drained * d_spleen
    VSp_tot <- VSp_drained 
    MSp <- (VSp_drained-VSpB_cap)*d_spleen
    VSpB <- VSpB_cap
    fcap_to_blood_sp<-  VSpB_cap/VSpB 
    PVSpF <- 0.150 #from  Kawai et al. (1994)
    VSpF <- PVSpF * (VSp_drained - VSpB_cap)
    VSpT <- VSp_drained - VSpF - VSpB_cap
    VSpB_lv <- 0
    
    # Heart
    PVHt <- 0.33e-2  #Brown et al. (1997)
    d_heart <- 1.03  #Brown et al. (1997)
    VHt_drained <- PVHt * BW / d_heart
    fcap_ht <- 0.061  #from  Bernareggi and  Rowland (1991)
    VHtB_cap <- fcap_ht * VHt_drained * d_heart
    PVHtB <- 0.26  #Brown et al. (1997)
    VHt_tot <- (VHt_drained - VHtB_cap) / (1 - PVHtB)
    MHt <- (VHt_drained - VHtB_cap)*d_heart
    VHtB <- PVHtB * VHt_tot
    fcap_to_blood_ht <-  VHtB_cap/VHtB 
    PVHtF <- 0.10 #from  Kawai et al. (1994)
    VHtF <- PVHtF * (VHt_drained - VHtB_cap)
    VHtT <- VHt_drained - VHtF - VHtB_cap
    VHtB_lv <- VHtB - VHtB_cap
    
    # Brain
    PVBr <- 0.57e-2  #Brown et al. (1997)
    d_brain <- 1.035
    VBr_drained <- PVBr * BW / d_brain
    fcap_br <- 0.014  #from  Bernareggi and  Rowland (1991)
    VBrB_cap <- fcap_br * VBr_drained * d_brain
    PVBrB <- 0.03  #Brown et al. (1997)
    VBr_tot <- (VBr_drained - VBrB_cap) / (1 - PVBrB)
    MBr <- (VBr_drained - VBrB_cap) * d_brain
    VBrB <- PVBrB * VBr_tot
    fcap_to_blood_br <-  VBrB_cap/VBrB 
    PVBrF <- 0.004 #from  Kawai et al. (1994).
    VBrF <- PVBrF * (VBr_drained - VBrB_cap)
    VBrT <- VBr_drained - VBrF - VBrB_cap
    VBrB_lv <- VBrB - VBrB_cap
    
    # Gonads
    PVGo <- 0.25e-2/0.230 #from PK-Sim
    d_gonads <- 1.05 #assumption
    VGo_drained <- PVGo * BW / d_gonads
    fcap_go <- 0.007  # value for testis, from  Bernareggi and  Rowland (1991)
    VGoB_cap <- fcap_go * VGo_drained * d_gonads
    PVGoB <- 0.14 #from PK-Sim
    VGo_tot <- (VGo_drained - VGoB_cap) / (1 - PVGoB)
    MGo <- (VGo_drained - VGoB_cap) * d_gonads
    VGoB <- PVGoB * VGo_tot
    fcap_to_blood_go <-  VGoB_cap/VGoB 
    PVGoF <- 0.07 #from PK-Sim
    VGoF <- PVGoF * (VGo_drained - VGoB_cap)
    VGoT <- VGo_drained - VGoF - VGoB_cap
    VGoB_lv <- VGoB - VGoB_cap
    
    # Skin
    PVSk <- 19.03e-2  #Brown et al. (1997)
    d_skin <- 1.12 #density of dermis. From Brown et al. (1997)
    VSk_drained <- PVSk * BW / d_skin
    fcap_sk <- 0.002 #from  Bernareggi and  Rowland (1991)
    VSkB_cap <- fcap_sk * VSk_drained * d_skin
    PVSkB <- 0.02  #Brown et al. (1997)
    VSk_tot <- (VSk_drained - VSkB_cap) / (1 - PVSkB)
    MSk <-  (VSk_drained - VSkB_cap)*d_skin
    VSkB <- PVSkB * VSk_tot
    fcap_to_blood_sk <-  VSkB_cap/VSkB 
    PVSkF <- 0.3 # from  Kawai et al. (1994). Wiig and Ree (1981) report a similar value (0.4) #[doi: 10.1111/j.1748-1716.1981.tb06901.x] 40 mL/100 g tissue, BW = 200-250 g
    VSkF <- PVSkF * (VSk_drained - VSkB_cap)
    VSkT <- VSk_drained - VSkF - VSkB_cap
    VSkB_lv <- VSkB - VSkB_cap
    
    # Bones
    PVBo <- mean(c(5,7))*1e-2  #from Brown et al. (1997)
    d_bone <- 1.96  #average of trabecular and cortical. From Brown et al. (1997)
    VBo_drained <- PVBo * BW / d_bone
    fcap_bo <- 0.019   #from  Bernareggi and  Rowland (1991)
    VBoB_cap <- fcap_bo * VBo_drained * d_bone
    PVBoB <- 0.04  #Brown et al. (1997)
    VBo_tot <- (VBo_drained - VBoB_cap) / (1 - PVBoB)
    MBo <- (VBo_drained - VBoB_cap) * d_bone
    VBoB <- PVBoB * VBo_tot
    fcap_to_blood_bo <-  VBoB_cap/VBoB 
    PVBoF <- 0.10 # from  Kawai et al. (1994).
    VBoF <- PVBoF * (VBo_drained - VBoB_cap)
    VBoT <- VBo_drained - VBoF - VBoB_cap
    VBoB_lv <- VBoB - VBoB_cap
    
  
    Mtissue_tot <- MK +  MLi + MIn + MSt +  MMu + MAd + MLu + MSp + MHt + MBr + MGo + MSk + MBo #[kg] total tissue mass without the RoB (and with no blood)
    Mlumen <- 1*VInL + 1*VStL #kg
    VBlood <- (0.06 * BW *1000+ 0.77)/1000    # Total blood. From Lee and Blaufox (1981) [PMID: 3965655]
    
    #Rest-of-the-body (re) !!!!! 
    MR_tot <- BW - Mtissue_tot - 1*VBlood - Mlumen
    d_re <- 1.05 # assumption
    VRe_tot <- MR_tot/d_re
    PVReB <- (PVKiB+PVLiB+PVLuB+PVMuB+PVAdB+PVSpB+PVHtB+PVBrB+PVGoB+PVStB+PVInB+PVSkB+PVBoB)/13 #average fraction of all the included organs
    VReB <- PVReB *VRe_tot
    fcap_to_blood_re <- (fcap_to_blood_ki + fcap_to_blood_li + fcap_to_blood_in + fcap_to_blood_st + fcap_to_blood_mu + fcap_to_blood_ad + fcap_to_blood_lu + 
                          fcap_to_blood_sp + fcap_to_blood_ht + fcap_to_blood_br + fcap_to_blood_go + fcap_to_blood_sk + fcap_to_blood_bo)/13
    VReB_cap <- fcap_to_blood_re*VReB
    VRe_drained <- VRe_tot - VReB + VReB_cap
    PVReF <-(PVKF+PVLiF+PVLuF+PVMuF+PVAdF+PVSpF+PVHtF+PVBrF+PVGoF+PVStF+PVInF+PVSkF+PVBoF)/13 #average interstitial fraction of all the included organs
    VReF <- PVReF * (VRe_drained -VReB_cap)
    VReT <- VRe_drained - VReF - VReB_cap 
    VReB_lv <- VReB - VReB_cap 
   
    Vcap_tot <-  VKiB_cap + VLiB_cap + VInB_cap + VStB_cap + VMuB_cap + VAdB_cap + VLuB_cap + VSpB_cap + VHtB_cap + VBrB_cap + 
                  VGoB_cap + VSkB_cap + VBoB_cap +  VReB_cap
    VB_central <- VBlood - Vcap_tot # blood not in capillaries
    f_ven <- 3/4 #from Brown et al. (1997)
    f_art <- 1/4 #from Brown et al. (1997)
    VVen <- f_ven *VB_central 	#volume of venous blood. From Kawai et al. (1994) [doi:10.1007/bf02353860]
    VArt <- f_art *VB_central 	#volume of arterial blood.  From Kawai et al. (1994)
    #-------------------------------------------------------------------------------------------------
    #----------------------------       Surface Areas       ------------------------------------------
    #-------------------------------------------------------------------------------------------------
    
    # For scaling capillary surface area we use simple linear scaling 
    # while for interstitial-intracellular surface area we scale with body weight to the 3/4,
    # following PK-Sim
    BW_ref <- 0.23
    linear_scaling_factor <- BW/BW_ref
    nonlinear_scaling_factor <- (BW/BW_ref)^0.75
    ##Capillary surface area for each tissue (Ai)(m^2),
    
    #Peritubular capillary density in From Gazzard et al. (2024) [doi: 10.1002/ar.25576] 
    d_peritubular <- 0.024 #um^2/um^3,
    VK_gazzard <- 0.00363#L
    Vcortex <- 2*1295* VKi_tot/VK_gazzard #mm^3,
    A_peritubular <- (d_peritubular*Vcortex*1e3/1e6) #m^2
    
    # The following values from PK-Sim "Endothelial Surface area", Niederalt et al., 2018, [doi: 10.1007/s10928-017-9559-4]
    PAL <- 1136e-4 #m^2
    AL <- PAL *  linear_scaling_factor #liver surface area (m^2)
    PAST <- 33.77e-4 #m^2
    AST <- PAST * linear_scaling_factor #stomach surface area (m^2)
    PAIN <- (74.82e-4+20.6e-4) #m^2 #large and small intestine
    AIN <- PAIN * linear_scaling_factor #intestine surface area (m^2)
    PAM <- 3043e-4 #m2
    AM <- PAM * linear_scaling_factor #muscle surface area (m^2)
    PAA <- 95.93e-4/0.23 #m2
    AA <- PAA * linear_scaling_factor #adipose surface area (m^2)
    PAR <- 100e-4#m^2, assumption
    AR <- PAR * linear_scaling_factor #surface area of rest of body (m^2)
    PLu <- 600.5e-4 #m^2 
    ALu <- PLu * linear_scaling_factor #lung surface area (m^2)
    PSP <- 162.3e-4 #m^2
    ASP <- PSP * linear_scaling_factor #spleen surface area (m^2)
    PH <-  201.1e-4 #m^2 
    AH <- PH * linear_scaling_factor #heart surface area (m^2)
    PBr <- 60.34e-4 #m^2
    ABr <- PBr * linear_scaling_factor #brain surface area (m^2)
    PGo <- 335.8e-4#m^2
    AGo <- PGo * linear_scaling_factor #gonads surface area (m^2)
    PSK <- 729.1e-4#m^2
    ASK <- PSK * linear_scaling_factor #skin surface area (m^2)
    PBo <- 621.4e-4#m^2
    ABo <- PBo * linear_scaling_factor #bones surface area (m^2)
    
    # Gut surface areas from PK-Sim
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

    SA_PKSim <- Duodenum+Upper_jejunum+Lower_jejunum+Upper_ileum+Lower_ileum+Cecum+     
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
    
    
    SA_PKSim_effective <- (Duodenum*EF_Duodenum+Upper_jejunum*EF_Upper_jejunum+Lower_jejunum*EF_Lower_jejunum+
                             Upper_ileum*EF_Upper_ileum+Lower_ileum*EF_Lower_ileum+Cecum*EF_Cecum+
                             colon_ascendens*EF_colon_ascendens+colon_transversum*EF_colon_transversum+
                             colon_descendens*EF_colon_descendens+colon_sigmoid*EF_colon_sigmoid)
    
    
    AINL <- SA_PKSim_effective * 1e-4 #m^2, no body weight scaling
    
    #Surface areas Interstitial - Intracellular (m^2), from PK-Sim 
    BW_ref <- 0.23/1e4
    AcK_total= 437.16*nonlinear_scaling_factor
    # for renal tubule, we assume that the total surface area is split proportionally to the volume of each segment.
    AcK_PTC <- AcK_total*VPTC/VFil # surface area of proximal tubule cells
    AcK_DALC <- AcK_total*VDALC/VFil # surface area of decending/ascending limb cells (loop of Henle)
    AcK_DTC <- AcK_total*VDTC/VFil # surface area of distal tubule cells
    AcK_CDC <- AcK_total*VCDC/VFil # surface area of collecting duct cells 
    AcKTrest <- AcK_total* VKTrest/VFil
    AcL= 84.45*nonlinear_scaling_factor
    AcST= 1007.31*nonlinear_scaling_factor
    AcIN= (400.94+152.39) *nonlinear_scaling_factor # small+large intestine
    AcM= 8.2*nonlinear_scaling_factor
    AcA= 3.87*nonlinear_scaling_factor
    AcLu= 0.05*nonlinear_scaling_factor
    AcSP= 564.05*nonlinear_scaling_factor
    AcH= 5.60*nonlinear_scaling_factor
    AcBr= 6.12e-4*nonlinear_scaling_factor
    AcGo= 2.01*nonlinear_scaling_factor
    AcSK= 0.11*nonlinear_scaling_factor
    AcBo= 6.52*nonlinear_scaling_factor
    # We don't have data for the surface area of IS-IC for the rest of the body, thus 
    # we naively assume an average:
    AcR= median(c(AcK_total,AcL,AcST,AcIN,AcM,AcA,AcLu,AcSP,AcH,AcBr,AcGo,AcSK,AcBo))
    
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
    
    #Alveolar cells surface area, m^2
    AcALF = 0.4* (BW/0.363)  #Stone et al., 1992 [doi:10.1165/ajrcmb/6.2.235]
    
    #----------------------------------------------------------------------------------------------
    #-------------------------             Flow Rates           -----------------------------------
    #----------------------------------------------------------------------------------------------

    ########
    #Blood flow rates 
    ###
    
    #(QBi, in L/h) to different tissues expressed using the corresponding percentage of cardiac output
    
    Qcardiac <- 0.235 * (BW^0.75) *60 #[L/h]
    PQBK <- 14.1/100 #Brown et al. 1997, p 438, Table 23
    QBK <- PQBK * Qcardiac #L/h
    PQBL <- 2.4/100 
    QBL <- PQBL * Qcardiac #L/h Brown et al. 1997, p 438, Table 23
    PQBST <- 0.16/100 #[doi: 10.3390/pharmaceutics6010097] Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBST = (8/1000)/Qcard 
    QBST <- PQBST * Qcardiac #L/h
    PQBIN <- 9.1/100 #[doi: 10.3390/pharmaceutics6010097] Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBIN = (451/1000)/Qcard 
    QBIN <- PQBIN * Qcardiac #L/h
    PQBM <- 27.8/100 #Brown et al. 1997, p 438, Table 23
    QBM <- PQBM * Qcardiac #L/h
    PQBA <- 7/100 #Brown et al. 1997, p 438, Table 23
    QBA <- PQBA * Qcardiac #L/h
    PQBLu <- 1 
    QBLu <- PQBLu * Qcardiac #L/h
    PQBSP <- 0.76/100 #[doi: 10.3390/pharmaceutics6010097] Qcard=0.235*(0.250^0.75)*60 (L/h) and PQBSP = (37.8/1000)/Qcard 
    QBSP <- PQBSP * Qcardiac #L/h
    PQBH <- 4.9/100 #Brown et al. 1997, p 438, Table 23
    QBH <- PQBH * Qcardiac #L/h
    PQBBr <- 2.0/100 #Brown et al. 1997, p 438, Table 23
    QBBr <- PQBBr * Qcardiac #L/h 
    PQBGo <- 1.04/100 #[doi: 10.1152/ajpregu.1987.253.2.R228]. Vgon =(0.25e-2/0.230)*335 [g]  Qgon = 60*(0.295/1000)*Vgon [L/h],  fraction:f_gon=  Qgon/Qcard
    QBGo <- PQBGo * Qcardiac #L/h
    PQBSK <- 5.8/100 # Brown et al. (1997)
    QBSK <- PQBSK * Qcardiac #L/h
    PBBo <- 12.2/100 #Brown et al. 1997, p 438, Table 23
    QBBo <- PBBo * Qcardiac #L/h
    
    # Total blood outflow from liver
    QBLtot <- QBL+QBSP+QBIN+QBST
    
    PQBR = 1 - PQBK - PQBL - PQBST - PQBIN - PQBM - PQBA - PQBH - PQBSK - PQBSP - PQBGo - PQBBr - PBBo
    QBR <- PQBR * Qcardiac #L/h
    
    
    #########
    #     Other fluids flow rates    
    ###
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration rate (GFR), in L/h
    
    #PQbile <- 90/1000/24 # [mL/d/kg BW]/1000/24 --> L/h/kg BW
    #Qbile <- PQbile * BW #L/h
    if (sex == "M"){
      PQbile = 0.206/2 #L/kg liver/h. The study uses male rats. From  Prasnicka et al. (2018) [doi: 10.1038/s41598-019-46150-7]
      Qbile = PQbile* VLi_drained  #L/h
    }else if (sex == "F"){
      # Females have 44% more bile flow. From Vonk et al. (1978) [doi:10.1042/cs0550253]
      PQbile = 0.206/2 #L/kg liver/h #
      Qbile = 1.44* PQbile* VLi_drained  #L/h
    }
    Qfeces <- (8.18/0.21)*BW/24/1000 #g/kg BW, based on Cui et al.(2010)
    feces_density <- 1.29 #g/cm^3 --> g/mL from Lupton 1986, Fig 1. Fiber free control diet, [doi: 10.1093/jn/116.1.164]
    
    if (sex == "M"){
      PQGFR <- 62.1  #L/h/kg kidney. From Corley et al., 2005 [doi: 10.1093/toxsci/kfi119]
      QGFR <- PQGFR * VKi_drained #L/h
      Qurine <- (60/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, [doi: 10.1002/nau.1006]
    }else if(sex == "F"){
      PQGFR <- 41.04 #L/h/kg kidney. From Corley et al., 2005 [doi: 10.1093/toxsci/kfi119]
      QGFR <- PQGFR * VKi_drained #L/h
      Qurine <- (85/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, [doi: 10.1002/nau.1006] 
    }
    QGE<- 1.4*BW^(-0.25) #gastric emptying time [1/(h/BW^(-0.25)]*BW^(-0.25)--> [1/h]. From Yang et al. (2013) [doi: 10.1016/j.taap.2013.03.022]
    
    #flows of filtrate compartments from Gilmer et al. (2018) [doi: 10.1152/ajprenal.00219.2018]
    # The values represent flow rates at the beginning of individual segments
    # Units: [nL/min]*2*1e-9*60 ---> L/h 
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
    QTDL <- Q_scaling_factor * QTDL_ref #L/h, Thin descending limb (short +long)
    QTAL <- Q_scaling_factor * QTAL_ref #L/h, Thin+thick (short +long) ascending 
    QDT <- Q_scaling_factor * QDT_ref  #L/h, Distal tubule (short + long)
    QCD <- Q_scaling_factor * QCD_ref  #L/h, collecting duct (Cortical + medullary)

    return(list(
      
       'VKiB'=VKiB, 'VKiB_cap'=VKiB_cap, "VKiB_lv" = VKiB_lv,
      'VKiF'=VKiF, 'VKT'=VKT, 'VFil'=VFil,'VBladder' = VBladder, 
      'VPT' = VPT, 'VDAL' = VDAL, 'VDT' = VDT, 'VCD' = VCD,
      'VPTC' = VPTC, 'VDALC' = VDALC, 'VDTC' = VDTC, 'VCDC' = VCDC,
      'VKTrest' = VKTrest,
      
      'VLiB'=VLiB, 'VLiB_cap'=VLiB_cap,'VLiB_lv'=VLiB_lv, 'VLiF'=VLiF, 'VLiT'=VLiT, 'VBile'=VBile,
      'VInB'=VInB, 'VInB_cap'=VInB_cap,'VInB_lv'=VInB_lv,'VInF'=VInF, 'VInT'=VInT,
      'VSt_drained'=VSt_drained, 'VStB_cap'=VStB_cap, 'VStB_lv'=VStB_lv,'VStB'=VStB,'VStF'=VStF, 'VStT'=VStT,
      'VStL'=VStL, 'VInL'=VInL,
      'VMu_drained'=VMu_drained, 'VMuB'=VMuB, 'VMuB_cap'=VMuB_cap, 'VMuB_lv'=VMuB_lv, 'VMuF'=VMuF, 'VMuT'=VMuT, 
      'VAd_drained'=VAd_drained, 'VAdB'=VAdB,'VAdB_cap'=VAdB_cap,'VAdB_lv'=VAdB_lv, 'VAdF'=VAdF, 'VAdT'=VAdT, 
      'VLu_drained'=VLu_drained, 'VLuB'=VLuB,  'VLuB_cap'=VLuB_cap, 'VLuB_lv'=VLuB_lv,'VLuF'=VLuF,'VLuAF'=VLuAF, 'VLuT'=VLuT,
      'VSp_drained'=VSp_drained, 'VSpB'=VSpB, 'VSpB_cap'=VSpB_cap,'VSpB_lv'=VSpB_lv,'VSpF'=VSpF, 'VSpT'=VSpT,
      'VHt_drained'=VHt_drained, 'VHtB'=VHtB, 'VHtB_cap'=VHtB_cap,'VHtB_lv'=VHtB_lv,'VHtF'=VHtF, 'VHtT'=VHtT,
      'VBr_drained'=VBr_drained, 'VBrB'=VBrB,'VBrB_cap'=VBrB_cap,'VBrB_lv'=VBrB_lv, 'VBrF'=VBrF, 'VBrT'=VBrT,
      'VGo_drained'=VGo_drained, 'VGoB'=VGoB, 'VGoB_cap'=VGoB_cap,'VGoB_lv'=VGoB_lv,'VGoF'=VGoF, 'VGoT'=VGoT,
      'VSk_drained'=VSk_drained,'VSkB'=VSkB, 'VSkB_cap'=VSkB_cap,'VSkB_lv'=VSkB_lv, 'VSkF'=VSkF, 'VSkT'=VSkT,
      'VBo_drained'=VBo_drained,'VBoB'=VBoB, 'VBoB_cap'=VBoB_cap, 'VBoB_lv'=VBoB_lv, 'VBoF'=VBoF, 'VBoT'=VBoT,
      'VRe_drained'=VRe_drained,'VReB'=VReB,'VReB_cap'=VReB_cap,'VReB_lv'=VReB_lv,'VReF'=VReF, 'VReT'=VReT, 
      
      "VKi_tot" = VKi_tot, "VLi_tot" = VLi_tot,  "VIn_tot" = VIn_tot, "VSt_tot" = VSt_tot, "VMu_tot" = VMu_tot, 
      "VAd_tot" = VAd_tot, "VLu_tot" = VLu_tot, "VSp_tot" = VSp_tot,
      "VHt_tot" = VHt_tot, "VBr_tot" = VBr_tot, "VGo_tot" = VGo_tot, "VSk_tot" = VSk_tot, "VBo_tot" = VBo_tot, "VRe_tot" = VRe_tot,
      
      "VBlood" = VBlood, 'VVen' = VVen, 'VArt' = VArt,
      
      "MLi_drained" = MLi_drained,  "MKi_drained" = MKi_drained,  "MIn_drained" = MIn_drained, "MLu_drained" = MLu_drained,
      "VLi_drained" = VLi_drained, "VKi_drained" = VKi_drained,
      'A_peritubular' = A_peritubular, 
      'AL'=AL, 'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
      'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
      'AIN'=AIN, 'AGo'=AGo,
      'ASK'= ASK, 'ABo'=ABo,
      
      'AINL' = AINL, 'AcL' = AcL, 'AcM' = AcM, 'AcST' = AcST, 
      'AcIN' = AcIN, 'AcA' = AcA, 'AcLu' = AcLu, 'AcALF' = AcALF, 
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
      'QBSK'=QBSK, 'QBBo'=QBBo,
      
      
      "QPT" = QPT, "QTDL" = QTDL, "QTAL" = QTAL, "QDT" = QDT, "QCD" = QCD,
      
      'f_tubular' =  f_tubular,  'f_PTC_prot_to_tub_prot' = f_PTC_prot_to_tub_prot, 
      'f_DALC_prot_to_tub_prot' = f_DALC_prot_to_tub_prot, 
      'f_DTC_prot_to_tub_prot' = f_DTC_prot_to_tub_prot , 
      'f_CDC_prot_to_tub_prot' = f_CDC_prot_to_tub_prot,
      
      
      
      "admin.time" = admin.time, "admin.dose" = admin.dose,
      "admin.type" = admin.type, "MW"=MW, "sex" = sex, "BW" = BW
      
    ))
    
  })
}   


ode.func <- function(time, inits, params){
  # t_24 <- floor(time/24)
  # BW_init <- unlist(params["BW"])
  # if (params["sex"] == "M"){
  #   growth_rate = 5.9/1000 #g/d
  #   BW_new <- BW_init + growth_rate*t_24
  #   if(BW_new>0.597){
  #     BW_new<-0.597
  #   }
  # }else{
  #   growth_rate = 3.5/1000 #g/d
  #   BW_new <- BW_init + growth_rate*t_24
  #   if(BW_new>0.365){
  #     BW_new<-0.365
  #   }
  # }
  # BW_linear_scaled_params <- c ('VBlood',  'VK', 'VKiB', 'VKiF',  'VFil','VPT' , 'VDAL' , 'VDT' , 'VCD' ,'VPTC' , 'VDALC' , 'VDTC' , 
  #                               'VCDC' ,'VLi', 'VLiB', 'VLiF',  'VBile',
  #                        'VMu_drained', 'VMuB', 'VMuF', 'VAd_drained', 'VAdB', 'VAdF', 'VAdT', 'VR', 'VReB',  'VReF',  'VVen' ,'VArt', 'VLu_drained', 'VLuB', 'VLuF',
  #                        'VLuAF','VSp_drained', 'VSpB', 'VSpF','VHt_drained', 'VHtB', 'VHtF','VBr_drained', 'VBrB', 'VBrF','VGo_drained', 'VGoB', 'VGoF',
  #                        'VIn', 'VInB', 'VInF', 'VSt_drained', 'VStB', 'VStF','VStL', 'VInL','VSk_drained','VSkB', 'VSkF',
  #                        'VBo_drained','VBoB', 'VBoF','VLiT', 'VKTrest','VInT', 'VStT','VMuT', 'VAdT', 
  #                        'VLuT', 'VSpT','VHtT', 'VBrT','VGoT', 'VSkT','VBoT', 'VReT',
  #                        'VKT', 'A_peritubular', 'AL', 'AM', 'AA', 'AR', 'ALu','ASP', 'AH', 'ABr', 'AST',
  #                        'AIN', 'AGo','ASK', 'ABo','AINL','APT' , 'ADAL' , 'ADT', 'ACD' ,
  #                        'Qfeces','Qbile', 'QGFR','Qurine','kabST',"QPT" , "QTDL", "QTAL" , "QDT", "QCD", 'CLfeces',
  #                        "CL_hepatobiliary", 
  #                        'VmL_Oatp', 'VmL_Ntcp','VmL_Oatp2',  'VmIn_Oatp2', 'VmK_Oatp','VmLu_Oatp_ap', 'VmLu_Oatp_bas',
  #                        'VmK_Oat1', 'VmK_Oat3','VmK_Urat','kPtcTu', 'kDalcTu' , 'kDtcTu' , 'kCdcTu' ,)
  # 
  #  BW_nonlinear_scaled_params <- c ('kKTrestF', 'kCdcF' , 'kDalcF' , 'kPtcF' , 'kDtcF' ,
  #                                   'kLFLT',  'kAFAT', 'kRFRT','kMFMT', 'kLuTLuF', 'kLuTLuAF', 'kSPFSPT' ,
  #                                   'kSTFSTT' , 'kINFINT' , 'kHFHT' ,'kBrFBrT' , 'kGoFGoT' ,'kSKFSKT' , 'kBoFBoT', 
  #                                   "k_gut_in", 'k_gut_out')
  # 
  #  Q_scaled_params <- c('QBK', 'QBL', 'QBLtot','QBM', 'QBA',
  #                      'QBR', 'QBLu','QBSP', 'QBH', 'QBBr', 'QBST','QBIN', 'QGE','QBGo','QBSK', 'QBBo', "PparaKi" ,"PparaLi" ,"PparaSt" ,"PparaIn" ,
  #                      "PparaMu" ,"PparaAd" ,"PparaRe" ,"PparaLu" ,"PparaSp","PparaHt" ,"PparaBr" ,"PparaGo","PparaSk" ,"PparaBo")
  # Qcardiac_init <-  unlist(params["Qcardiac"])
  # params["Qcardiac"] <- unlist(params["Qcardiac"])* (BW_new/BW_init)^0.75
  # params["QGE"] <-   unlist(params["QGE"])* (BW_new/BW_init)^(-0.25) 
  # params["CLfeces"] <-   unlist(params["CLfeces"])* (BW_new/BW_init)^(-0.25) 
  # 
  # params[BW_linear_scaled_params] <- lapply(params[BW_linear_scaled_params], function(x) x * unname(BW_new/BW_init))
  # params[BW_nonlinear_scaled_params] <- lapply(params[BW_nonlinear_scaled_params], function(x) x * unname(BW_new/BW_init)^0.75)
  # 
  # params[Q_scaled_params] <- lapply(params[Q_scaled_params], function(x) x * unname(unlist(params["Qcardiac"])/Qcardiac_init))
  # 
  with(as.list(c(inits, params)),{
    
    #====================PFOA mass balance at each tissue or fluid compartment==============================     
    
    # Concentrations in ug/L, mass in ug, time in hours, volumes in L, flows in L/h
    
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
    CKB <- MKB/VKiB_cap # blood concentration
    CKBf <- MKBf/VKiB_cap
    CKBb <- MKBb/VKiB_cap
    CKF <- MKF/VKiF  #interstitial fluid concentration
    CKFf <- MKFf/VKiF
    CKFb <- MKFb/VKiF
    
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
    CKTb <- MKTb/VKT
    CKTf <- MKTf/VKT
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
    CLB <- MLB/VLiB_cap # blood concentration
    CLBf <- MLBf/VLiB_cap
    CLBb <- MLBb/VLiB_cap
    CLF <- MLF/VLiF  #interstitial fluid concentration
    CLFf <- MLFf/VLiF
    CLFb <- MLFb/VLiF
    CLT <- MLT/VLiT # tissue concentration
    CLTf <- MLTf/VLiT
    CLTb <- MLTb/VLiT
    CBile <- MBile/VBile #Bile  canaliculi
    
    #Stomach
    
    MSTB <- MSTBf + MSTBb
    MSTF <- MSTFf + MSTFb
    MSTT <- MSTTf
    CSTB <- MSTB/VStB_cap # blood concentration
    CSTBf <- MSTBf/VStB_cap
    CSTBb <- MSTBb/VStB_cap
    CSTF <- MSTF/VStF  #interstitial fluid concentration
    CSTFf <- MSTFf/VStF
    CSTFb <- MSTFb/VStF
    CSTTf <- MSTTf/VStT # tissue concentration
    
    #Intestine
    
    MINB <- MINBf + MINBb
    MINF <- MINFf + MINFb
    MINT <- MINTf 
    CINB <- MINB/VInB_cap # blood concentration
    CINBf <- MINBf/VInB_cap
    CINBb <- MINBb/VInB_cap
    CINF <- MINF/VInF  #Interstitial fluid concentration
    CINFf <- MINFf/VInF
    CINFb <- MINFb/VInF
    CINTf <- MINTf/VInT # tissue concentration
    
    #Stomach and Intestine lumens
    CSTL = MSTL/VStL # Stomach Lumen concentration
    CINL = MINL/VInL # Intestine Lumen concentration
    
    #Muscle
    
    MMB <- MMBf + MMBb
    MMF <- MMFf + MMFb
    MMT <- MMTf 
    CMB <- MMB/VMuB_cap # blood concentration
    CMBf <- MMBf/VMuB_cap
    CMBb <- MMBb/VMuB_cap
    CMF <- MMF/VMuF  #Interstitial fluid concentration
    CMFf <- MMFf/VMuF
    CMFb <- MMFb/VMuF
    CMTf <-  MMTf/VMuT # tissue concentration
    
    #Adipose
    
    MAB <- MABf + MABb
    MAF <- MAFf + MAFb
    MAT <- MATf 
    CAB <- MAB/VAdB_cap # blood concentration
    CABf <- MABf/VAdB_cap
    CABb <- MABb/VAdB_cap
    CAF <- MAF/VAdF  #Interstitial fluid concentration
    CAFf <- MAFf/VAdF
    CAFb <- MAFb/VAdF
    CATf <- MATf/VAdT # tissue concentration
    
    #Rest-of-the-body
    
    MRB <- MRBf + MRBb
    MRF <- MRFf + MRFb
    MRT <- MRTf 
    CRB <- MRB/VReB_cap # blood concentration
    CRBf <- MRBf/VReB_cap
    CRBb <- MRBb/VReB_cap
    CRF <- MRF/VReF  #Interstitial fluid concentration
    CRFf <- MRFf/VReF
    CRFb <- MRFb/VReF
    CRTf <- MRTf/VReT # tissue concentration
    
    #Lung
    
    MLuB <- MLuBf + MLuBb
    MLuF <- MLuFf + MLuFb
    MLuT <- MLuTf 
    MLuAF <- MLuAFf + MLuAFb
    CLuB <- MLuB/VLuB_cap # blood concentration
    CLuBf <- MLuBf/VLuB_cap
    CLuBb <- MLuBb/VLuB_cap
    CLuF <- MLuF/VLuF  #Interstitial fluid concentration
    CLuFf <- MLuFf/VLuF
    CLuFb <- MLuFb/VLuF
    CLuTf <- MLuTf/VLuT #tissue concentration
    CLuAF <- MLuAF/VLuAF #alveolar lining fluid concentration
    CLuAFf <- MLuAFf/VLuAF
    CLuAFb <- MLuAFb/VLuAF
    
    #Spleen
    
    MSPB <- MSPBf + MSPBb
    MSPF <- MSPFf + MSPFb
    MSPT <- MSPTf 
    CSPB <- MSPB/VSpB_cap # blood concentration
    CSPBf <- MSPBf/VSpB_cap
    CSPBb <- MSPBb/VSpB_cap
    CSPF <- MSPF/VSpF  #Interstitial fluid concentration
    CSPFf <- MSPFf/VSpF
    CSPFb <- MSPFb/VSpF
    CSPTf <-  MSPTf/VSpT # tissue concentration
    
    #Heart
    MHB <- MHBf + MHBb
    MHF <- MHFf + MHFb
    MHT <- MHTf 
    CHB <- MHB/VHtB_cap # blood concentration
    CHBf <- MHBf/VHtB_cap
    CHBb <- MHBb/VHtB_cap
    CHF <- MHF/VHtF  #Interstitial fluid concentration
    CHFf <- MHFf/VHtF
    CHFb <- MHFb/VHtF
    CHTf <- MHTf/VHtT # tissue concentration
    
    #Brain
    
    MBrB <- MBrBf + MBrBb
    MBrF <- MBrFf + MBrFb
    MBrT <- MBrTf 
    CBrB <- MBrB/VBrB_cap # blood concentration
    CBrBf <- MBrBf/VBrB_cap
    CBrBb <- MBrBb/VBrB_cap
    CBrF <- MBrF/VBrF  #Interstitial fluid concentration
    CBrFf <- MBrFf/VBrF
    CBrFb <- MBrFb/VBrF
    CBrTf <-  MBrTf/VBrT # tissue concentration
    
    #gonads
    
    MGoB <- MGoBf + MGoBb
    MGoF <- MGoFf + MGoFb
    MGoT <- MGoTf 
    CGoB <- MGoB/VGoB_cap # blood concentration
    CGoBf <- MGoBf/VGoB_cap
    CGoBb <- MGoBb/VGoB_cap
    CGoF <- MGoF/VGoF  #Interstitial fluid concentration
    CGoFf <- MGoFf/VGoF
    CGoFb <- MGoFb/VGoF
    CGoTf <-  MGoTf/VGoT # tissue concentration
    
    #Skin
    
    MSKB <- MSKBf + MSKBb
    MSKF <- MSKFf + MSKFb
    MSKT <- MSKTf 
    CSKB <- MSKB/VSkB_cap # blood concentration
    CSKBf <- MSKBf/VSkB_cap
    CSKBb <- MSKBb/VSkB_cap
    CSKF <- MSKF/VSkF  #Interstitial fluid concentration
    CSKFf <- MSKFf/VSkF
    CSKFb <- MSKFb/VSkF
    CSKTf <- MSKTf/VSkT # tissue concentration
    
    #Bones
    
    MBoB <- MBoBf + MBoBb
    MBoF <- MBoFf + MBoFb
    MBoT <- MBoTf 
    CBoB <- MBoB/VBoB_cap # blood concentration
    CBoBf <- MBoBf/VBoB_cap
    CBoBb <- MBoBb/VBoB_cap
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
    dCalbBrFf <- koff_alb*CBrFb/MW/1e6 - kon_alb*CalbBrFf*CBrFf/MW/1e6
    
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
    dCalbLuAFf = koff_alb*CLuAFb/MW/1e6 - kon_alb*CalbLuAFf*CLuAFf/MW/1e6
    
    # Bound PFOA
    #Blood
    dMVenb <-  kon_alb*CalbVenf*CVenf*VVen -  koff_alb*CVenb*VVen 
    dMArtb <- kon_alb*CalbArtf*CArtf*VArt -  koff_alb*CArtb*VArt 
    dMKBb <- kon_alb*CalbKBf*CKBf*VKiB_cap - koff_alb*CKBb*VKiB_cap 
    dMLBb <- kon_alb*CalbLBf*CLBf*VLiB_cap - koff_alb*CLBb*VLiB_cap 
    dMSTBb <- kon_alb*CalbSTBf*CSTBf*VStB_cap - koff_alb*CSTBb*VStB_cap
    dMINBb <- kon_alb*CalbINBf*CINBf*VInB_cap - koff_alb*CINBb*VInB_cap 
    dMMBb <- kon_alb*CalbMBf*CMBf*VMuB_cap - koff_alb*CMBb*VMuB_cap 
    dMABb <- kon_alb*CalbABf*CABf*VAdB_cap - koff_alb*CABb*VAdB_cap 
    dMRBb <- kon_alb*CalbRBf*CRBf*VReB_cap - koff_alb*CRBb*VReB_cap 
    dMLuBb <- kon_alb*CalbLuBf*CLuBf*VLuB_cap - koff_alb*CLuBb*VLuB_cap 
    dMSPBb <- kon_alb*CalbSPBf*CSPBf*VSpB_cap - koff_alb*CSPBb*VSpB_cap 
    dMHBb <- kon_alb*CalbHBf*CHBf*VHtB_cap - koff_alb*CHBb*VHtB_cap
    dMBrBb <- kon_alb*CalbBrBf*CBrBf*VBrB_cap - koff_alb*CBrBb*VBrB_cap
    dMGoBb <- kon_alb*CalbGoBf*CGoBf*VGoB_cap - koff_alb*CGoBb*VGoB_cap 
    dMSKBb <- kon_alb*CalbSKBf*CSKBf*VSkB_cap - koff_alb*CSKBb*VSkB_cap 
    dMBoBb <- kon_alb*CalbBoBf*CBoBf*VBoB_cap - koff_alb*CBoBb*VBoB_cap 
    
    #Interstitial fluid
    dMKFb <- kon_alb*CalbKFf*CKFf*VKiF - koff_alb*CKFb*VKiF 
    dMLFb <- kon_alb*CalbLFf*CLFf*VLiF - koff_alb*CLFb*VLiF 
    dMSTFb <- kon_alb*CalbSTFf*CSTFf*VStF - koff_alb*CSTFb*VStF 
    dMINFb <- kon_alb*CalbINFf*CINFf*VInF - koff_alb*CINFb*VInF
    dMMFb <- kon_alb*CalbMFf*CMFf*VMuF - koff_alb*CMFb*VMuF 
    dMAFb <- kon_alb*CalbAFf*CAFf*VAdF - koff_alb*CAFb*VAdF 
    dMRFb <- kon_alb*CalbRFf*CRFf*VReF - koff_alb*CRFb*VReF 
    dMLuFb <- kon_alb*CalbLuFf*CLuFf*VLuF - koff_alb*CLuFb*VLuF 
    dMSPFb <- kon_alb*CalbSPFf*CSPFf*VSpF - koff_alb*CSPFb*VSpF 
    dMHFb <- kon_alb*CalbHFf*CHFf*VHtF - koff_alb*CHFb*VHtF 
    dMBrFb <- kon_alb*CalbBrFf*CBrFf*VBrF - koff_alb*CBrFb*VBrF 
    dMGoFb <- kon_alb*CalbGoFf*CGoFf*VGoF - koff_alb*CGoFb*VGoF 
    dMSKFb <- kon_alb*CalbSKFf*CSKFf*VSkF - koff_alb*CSKFb*VSkF 
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
    
    dMLTb <- kon_fabp*CFabpLTf*CLTf*VLiT - koff_fabp*CLTb*VLiT 
    
    #Alveolar lining fluid
    dMLuAFb <-  kon_alb*CalbLuAFf*CLuAFf*VLuAF -  koff_alb*CLuAFb*VLuAF 
    
    #====================================================================================================================
    
    #Arterial Blood
    dMArtf = QBLu*CLuBf - CArtf*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+
                                   QBST+QBIN+QBGo+QBSK+QBBo) - QGFR*CArtf +
      koff_alb*CArtb*VArt - kon_alb*CalbArtf*CArtf*VArt 
    
    #Venous Blood
    dMVenf = - CVenf*QBLu + QBK*CKBf + QBLtot*CLBf + QBM*CMBf + QBA*CABf + QBR*CRBf+
      QBH*CHBf + QBBr*CBrBf+ QBGo*CGoBf + QBSK*CSKBf + QBBo*CBoBf+
      koff_alb*CVenb*VVen - kon_alb*CalbVenf*CVenf*VVen 
    
    #Kidney
    #blood subcompartment
    dMKBf = QBK*CArtf - QBK*CKBf- (Ptrans_diff_K+PparaKi)*A_peritubular*(CKBf-CKFf) + 
      (VmK_baso*CPTCf/(KmK_baso+CPTCf))+
      koff_alb*CKBb*VKiB_cap - kon_alb*CalbKBf*CKBf*VKiB_cap
    
    #interstitial fluid subcompartment
    dMKFf =  (Ptrans_diff_K+PparaKi)*A_peritubular*(CKBf-CKFf)- 
      kPtcF*(CKFf-CPTCf) - kDalcF*(CKFf-CDALCf) -
      kDtcF*(CKFf-CDTCf) - kCdcF*(CKFf-CCDCf)   -  kKTrestF*(CKFf-CKTrestf) + 
      koff_alb*CKFb*VKiF - kon_alb*CalbKFf*CKFf*VKiF-
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
      (VmK_Urat*CPT/(KmK_Urat+CPT)) + (VmK_api*CPTCf/(KmK_api+CPTCf))- QTDL*CPT
    
    #Descending limb, Ascending limb (Loop of Henle )
    dMDAL =  QTDL*CPT + kDalcTu*(CDALCf - CDAL)  -  QDT*CDAL
    
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
      koff_alb*CLBb*VLiB_cap -kon_alb*CalbLBf*CLBf*VLiB_cap 
    #interstitial fluid subcompartment 
    dMLFf = (Ptrans_diff_L+PparaLi)*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) - 
      (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) - (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf)) -
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_alb*CLFb*VLiF -kon_alb*CalbLFf*CLFf*VLiF 
    #Liver tissue subcompartment
    dMLTf = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) + (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf))+
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_fabp*CLTb*VLiT-
      kon_fabp*CFabpLTf*CLTf*VLiT - CL_hepatobiliary*CLTf
    #Bile  canaliculi subcompartment
    dMLbile = CL_hepatobiliary*CLTf - CBile*Qbile
    
    #Stomach
    #blood subcompartment
    dMSTBf = QBST*CArtf - QBST*CSTBf - (Ptrans_diff_ST+PparaSt)*AST*(CSTBf-CSTFf)  +
      koff_alb*CSTBb*VStB_cap-kon_alb*CalbSTBf*CSTBf*VStB_cap 
    #interstitial fluid subcompartment 
    dMSTFf = (Ptrans_diff_ST+PparaSt)*AST*(CSTBf-CSTFf) - kSTFSTT*(CSTFf-CSTTf) +
      koff_alb*CSTFb*VStF - kon_alb*CalbSTFf*CSTFf*VStF 
    #Stomach tissue subcompartment
    dMSTTf = kSTFSTT*(CSTFf-CSTTf) + kabST*CSTL 
    #Stomach lumen
    dMSTL = - QGE*CSTL -kabST*CSTL 
    
    
    #Intestine
    #blood subcompartment
    dMINBf = QBIN*CArtf - QBIN*CINBf - (Ptrans_diff_IN+PparaIn)*AIN*(CINBf-CINFf) +
      koff_alb*CINBb*VInB_cap - kon_alb*CalbINBf*CINBf*VInB_cap
    #interstitial fluid subcompartment 
    dMINFf = (Ptrans_diff_IN+PparaIn)*AIN*(CINBf-CINFf) - kINFINT*(CINFf-CINTf) +
      koff_alb*CINFb*VInF - kon_alb*CalbINFf*CINFf*VInF
    #Intestine tissue subcompartment
    dMINTf = kINFINT*(CINFf-CINTf) + k_gut_in*CINL - k_gut_out*CINTf + (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL)) 
    #Intestine lumen
    dMINL = QGE*CSTL - (CLfeces*CINL) - k_gut_in*CINL + k_gut_out*CINTf + CBile*Qbile - 
      (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))
    
    
    #Muscle
    #blood subcompartment
    dMMBf = QBM*CArtf - QBM*CMBf - (Ptrans_diff_M+PparaMu)*AM*(CMBf-CMFf)  +
      koff_alb*CMBb*VMuB_cap - kon_alb*CalbMBf*CMBf*VMuB_cap 
    #interstitial fluid subcompartment 
    dMMFf =(Ptrans_diff_M+PparaMu)*AM*(CMBf-CMFf) - kMFMT*(CMFf- CMTf) +
      koff_alb*CMFb*VMuF - kon_alb*CalbMFf*CMFf*VMuF
    #Muscle tissue subcompartment 
    dMMTf = kMFMT*(CMFf- CMTf) 
    
    
    #Adipose
    #blood subcompartment
    dMABf = QBA*CArtf - QBA*CABf - (Ptrans_diff_A+PparaAd)*AA*(CABf-CAFf)  +
      koff_alb*CABb*VAdB_cap - kon_alb*CalbABf*CABf*VAdB_cap
    #interstitial fluid subcompartment 
    dMAFf = (Ptrans_diff_A+PparaAd)*AA*(CABf-CAFf)  - kAFAT*(CAFf-CATf) +
      koff_alb*CAFb*VAdF - kon_alb*CalbAFf*CAFf*VAdF 
    #Adipose tissue subcompartment 
    dMATf =  kAFAT*(CAFf-CATf)
    
    
    #Rest of body
    #blood subcompartment
    dMRBf = QBR*CArtf - QBR*CRBf - (Ptrans_diff_R+PparaRe)*AR*(CRBf-CRFf)  +
      koff_alb*CRBb*VReB_cap - kon_alb*CalbRBf*CRBf*VReB_cap
    #interstitial fluid subcompartment 
    dMRFf = (Ptrans_diff_R+PparaRe)*AR*(CRBf-CRFf) - kRFRT*(CRFf -CRTf) +
      koff_alb*CRFb*VReF - kon_alb*CalbRFf*CRFf*VReF 
    #Rest of body tissue subcompartment 
    dMRTf = kRFRT*(CRFf -CRTf) 
    
    
    #Lung 
    #blood subcompartment
    dMLuBf = CVenf*QBLu - QBLu*CLuBf - (Ptrans_diff_Lu+PparaLu)*ALu*(CLuBf-CLuFf) +
      koff_alb*CLuBb*VLuB_cap - kon_alb*CalbLuBf*CLuBf*VLuB_cap 
    #interstitial fluid subcompartment
    dMLuFf = (Ptrans_diff_Lu+PparaLu)*ALu*(CLuBf-CLuFf)+ kLuTLuF*(CLuTf-CLuFf) + 
      koff_alb*CLuFb*VLuF - kon_alb*CalbLuFf*CLuFf*VLuF - (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) 
    #Lung tissue
    dMLuTf =  - kLuTLuF*(CLuTf-CLuFf) -  kLuTLuAF*(CLuTf-CLuAFf) + (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) +
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf)) 
    #Alveolar lining fluid
    dMLuAFf =  kLuTLuAF*(CLuTf-CLuAFf) + koff_alb*CLuAFb*VLuAF - kon_alb*CalbLuAFf*CLuAFf*VLuAF -
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf)) 
    
    
    #Spleen
    #blood subcompartment
    dMSPBf = QBSP*CArtf - QBSP*CSPBf - (Ptrans_diff_SP+PparaSp)*ASP*(CSPBf-CSPFf) + 
      koff_alb*CSPBb*VSpB_cap - kon_alb*CalbSPBf*CSPBf*VSpB_cap
    #interstitial fluid subcompartment 
    dMSPFf = (Ptrans_diff_SP+PparaSp)*ASP*(CSPBf-CSPFf) - kSPFSPT*(CSPFf -CSPTf) +
      koff_alb*CSPFb*VSpF - kon_alb*CalbSPFf*CSPFf*VSpF 
    #Spleen tissue subcompartment 
    dMSPTf = kSPFSPT*(CSPFf -CSPTf)
    
    
    #Heart
    #blood subcompartment
    dMHBf = QBH*CArtf - QBH*CHBf - (Ptrans_diff_H+PparaHt)*AH*(CHBf-CHFf) + 
      koff_alb*CHBb*VHtB_cap - kon_alb*CalbHBf*CHBf*VHtB_cap 
    #interstitial fluid subcompartment 
    dMHFf = (Ptrans_diff_H+PparaHt)*AH*(CHBf-CHFf) - kHFHT*(CHFf -CHTf) + 
      koff_alb*CHFb*VHtF - kon_alb*CalbHFf*CHFf*VHtF 
    #Heart tissue subcompartment 
    dMHTf = kHFHT*(CHFf -CHTf) 
    
    
    #Brain
    #blood subcompartment
    dMBrBf = QBBr*CArtf - QBBr*CBrBf - (Ptrans_diff_Br+PparaBr)*ABr*(CBrBf-CBrFf)  + 
      koff_alb*CBrBb*VBrB_cap - kon_alb*CalbBrBf*CBrBf*VBrB_cap 
    #interstitial fluid subcompartment 
    dMBrFf = (Ptrans_diff_Br+PparaBr)*ABr*(CBrBf-CBrFf) - kBrFBrT*(CBrFf -CBrTf) +
      koff_alb*CBrFb*VBrF - kon_alb*CalbBrFf*CBrFf*VBrF 
    #Brain tissue subcompartment 
    dMBrTf = kBrFBrT*(CBrFf -CBrTf) 
    
    
    #Gonads
    #blood subcompartment
    dMGoBf = QBGo*CArtf - QBGo*CGoBf - (Ptrans_diff_Go+PparaGo)*AGo*(CGoBf-CGoFf) +
      koff_alb*CGoBb*VGoB_cap - kon_alb*CalbGoBf*CGoBf*VGoB_cap 
    #interstitial fluid subcompartment 
    dMGoFf = (Ptrans_diff_Go+PparaGo)*AGo*(CGoBf-CGoFf) - kGoFGoT*(CGoFf -CGoTf) +
      koff_alb*CGoFb*VGoF - kon_alb*CalbGoFf*CGoFf*VGoF 
    #gonads tissue subcompartment 
    dMGoTf = kGoFGoT*(CGoFf -CGoTf) 
    
    
    #Skin
    #blood subcompartment
    dMSKBf = QBSK*CArtf - QBSK*CSKBf - (Ptrans_diff_SK+PparaSk)*ASK*(CSKBf-CSKFf) +
      koff_alb*CSKBb*VSkB_cap - kon_alb*CalbSKBf*CSKBf*VSkB_cap 
    #interstitial fluid subcompartment
    dMSKFf =(Ptrans_diff_SK+PparaSk)*ASK*(CSKBf-CSKFf) - kSKFSKT*(CSKFf -CSKTf) +
      koff_alb*CSKFb*VSkF - kon_alb*CalbSKFf*CSKFf*VSkF
    #Skin tissue subcompartment
    dMSKTf = kSKFSKT*(CSKFf -CSKTf) 
    
    
    #Bones
    #blood subcompartment
    dMBoBf = QBBo*CArtf - QBBo*CBoBf - (Ptrans_diff_Bo+PparaBo)*ABo*(CBoBf-CBoFf) +
      koff_alb*CBoBb*VBoB_cap - kon_alb*CalbBoBf*CBoBf*VBoB_cap 
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
    Cplasma <- Cblood/(1-Hct)
    
    Mblood <- MVen +MArt
    Mkidney <- MKB + MKF+ MKT + Mfil + Cblood*VKiB_lv  
    Mkidney_bled <- MKB + MKF+ MKT + Mfil
    Mbrain <- MBrB + MBrF+ MBrT +  Cblood*VBrB_lv
    Mbrain_bled <- MBrB + MBrF+ MBrT
    Mliver <- MLB + MLF+ MLT + MBile +  Cblood*VLiB_lv
    Mliver_bled <- MLB + MLF+ MLT + MBile
    
   
    Ckidney <- (MKB + MKF+ MKT + Mfil+ Cblood*VKiB_lv )/VKi_tot
    Cliver <- (MLB + MLF+ MLT + MBile + Cblood*VLiB_lv)/(VLi_tot)
    Cintestine <-  (MINB + MINF+ MINT+MINL+ Cblood*VInB_lv)/VIn_tot
    Cstomach <-  (MSTB + MSTF+ MSTT + MSTL+ Cblood*VStB_lv)/VSt_tot
    Cmuscle <-  (MMB + MMF+ MMT+ Cblood*VMuB_lv)/VMu_tot
    Cadipose <-  (MAB + MAF+ MAT+ Cblood*VAdB_lv)/VAd_tot
    Clungs <-  (MLuB + MLuF+ MLuT + MLuAF+ Cblood*VLuB_lv)/VLu_tot
    Clungtissue <- (MLuB + MLuF+ MLuT+ Cblood*VLuB_lv)/(VLuB+VLuF+VLuT)
    Cspleen <-  (MSPB + MSPF+ MSPT+ Cblood*VSpB_lv)/VSp_tot
    Cheart <-  (MHB + MHF+ MHT+ Cblood*VHtB_lv)/VHt_tot
    Cbrain <-  (MBrB + MBrF+ MBrT+ Cblood*VBrB_lv)/VBr_tot
    Cgonads <-  (MGoB + MGoF+ MGoT+ Cblood*VGoB_lv)/VGo_tot
    Cskin <-  (MSKB + MSKF+ MSKT+ Cblood*VSkB_lv)/VSk_tot
    Cbones <-  (MBoB + MBoF+ MBoT+ Cblood*VBoB_lv)/VBo_tot
    Crest <-  (MRB + MRF+ MRT+ Cblood*VReB_lv)/VRe_tot
    Ccarcass <- (MMB+MMF+MMT+MAB+MAF+MAT+MRB+MRF+MRT+MBoB+MBoF+MBoT+MSKB+MSKF+MSKT+ 
                   Cblood*(VMuB_lv+VAdB_lv+VReB_lv+VBoB_lv+VSkB_lv))/(VMu_tot+VAd_tot+VRe_tot+VBo_tot+VSk_tot)
    
    Ckidney_bled <- (MKB + MKF+ MKT + Mfil)/VKi_tot
    Cliver_bled <- (MLB + MLF+ MLT + MBile )/(VLi_tot)
    Cintestine_bled <-  (MINB + MINF+ MINT+MINL)/VIn_tot
    Cstomach_bled <-  (MSTB + MSTF+ MSTT + MSTL)/VSt_tot
    Cmuscle_bled <-  (MMB + MMF+ MMT)/VMu_tot
    Cadipose_bled <-  (MAB + MAF+ MAT)/VAd_tot
    Clungs_bled <-  (MLuB + MLuF+ MLuT + MLuAF)/VLu_tot
    Clungtissue_bled <- (MLuB + MLuF+ MLuT)/(VLuB+VLuF+VLuT)
    Cspleen_bled <-  (MSPB + MSPF+ MSPT)/VSp_tot
    Cheart_bled <-  (MHB + MHF+ MHT)/VHt_tot
    Cbrain_bled <-  (MBrB + MBrF+ MBrT)/VBr_tot
    Cgonads_bled <-  (MGoB + MGoF+ MGoT)/VGo_tot
    Cskin_bled <-  (MSKB + MSKF+ MSKT)/VSk_tot
    Cbones_bled <-  (MBoB + MBoF+ MBoT)/VBo_tot
    Crest_bled <-  (MRB + MRF+ MRT)/VRe_tot
    Ccarcass_bled <- (MMB+MMF+MMT+MAB+MAF+MAT+MRB+MRF+MRT+MBoB+MBoF+MBoT+MSKB+MSKF+MSKT)/(VMu_tot+VAd_tot+VRe_tot+VBo_tot+VSk_tot)
    
    CalveolarLF <- MLuAF/VLuAF
    Cfeces <- Mfeces/(Vfeces*feces_density)
    Curine <- Murine/Vurine
    
    
    
    list(c( 'dCalbVenf' = dCalbVenf, 'dCalbArtf' = dCalbArtf, 
            'dCalbKBf' = dCalbKBf, 'dCalbLBf' = dCalbLBf, 'dCalbSTBf' = dCalbSTBf, 
            'dCalbINBf' = dCalbINBf, 'dCalbMBf' = dCalbMBf, 'dCalbABf' = dCalbABf,
            'dCalbRBf' = dCalbRBf, 'dCalbLuBf' = dCalbLuBf, 'dCalbSPBf' = dCalbSPBf,
            'dCalbHBf' = dCalbHBf, 'dCalbBrBf' = dCalbBrBf, 'dCalbGoBf' = dCalbGoBf, 
            'dCalbSKBf' = dCalbSKBf, 'dCalbBoBf'=dCalbBoBf, 'dCalbKFf' = dCalbKFf,
            'dCalbLFf' = dCalbLFf, 'dCalbSTFf' = dCalbSTFf,'dCalbINFf' = dCalbINFf,
            'dCalbMFf' = dCalbMFf, 'dCalbAFf' = dCalbAFf,'dCalbRFf' = dCalbRFf,
            'dCalbLuFf' = dCalbLuFf, 'dCalbSPFf' = dCalbSPFf, 'dCalbHFf' = dCalbHFf,
            'dCalbBrFf' = dCalbBrFf, 'dCalbGoFf' = dCalbGoFf, 'dCalbSKFf' = dCalbSKFf,
            'dCalbBoFf' = dCalbBoFf, 'dCa2uPTCf' = dCa2uPTCf, 'dCa2uDALCf' = dCa2uDALCf,
            'dCa2uDTCf' = dCa2uDTCf, 'dCa2uCDCf' = dCa2uCDCf, 'dCa2uKTrestf' = dCa2uKTrestf,
            'dCFabpPTCf' = dCFabpPTCf, 'dCFabpDALCf' = dCFabpDALCf, 'dCFabpDTCf' = dCFabpDTCf,
            'dCFabpCDCf' = dCFabpCDCf, 'dCFabpKTrestf' = dCFabpKTrestf,
            'dCFabpLTf' = dCFabpLTf, 'dCalbLuAFf' = dCalbLuAFf,
            
            
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
            'dMLTb' = dMLTb, 'dMLuAFb'=dMLuAFb,
            
            
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
            'dMRBf'=dMRBf, 'dMRFf'=dMRFf,'dMRTf'=dMRTf,
            'dMLuBf'=dMLuBf, 'dMLuFf'=dMLuFf,'dMLuTf'=dMLuTf,'dMLuAFf' = dMLuAFf,
            
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
    'CRBb'=CRBb, 'CRF'=CRF, 'CRFf'=CRFf, 'CRFb'=CRFb, 'CRTf'=CRTf, 'CLuB'=CLuB,
    'CLuBf'=CLuBf, 'CLuBb'=CLuBb, 'CLuF'=CLuF, 'CLuFf'=CLuFf, 'CLuFb'=CLuFb, 'CLuTf'=CLuTf,
    'CLuAF'=CLuAF, 'CLuAFf'=CLuAFf, 'CLuAFb'=CLuAFb, 'CSPB'=CSPB, 'CSPBf'=CSPBf, 
    'CSPBb'=CSPBb, 'CSPF'=CSPF, 'CSPFf'=CSPFf, 'CSPFb'=CSPFb, 'CSPTf'=CSPTf,  'CHB'=CHB,
    'CHBf'=CHBf, 'CHBb'=CHBb, 'CHF'=CHF, 'CHFf'=CHFf, 'CHFb'=CHFb, 'CHTf'=CHTf, 
    'CBrB'=CBrB, 'CBrBf'=CBrBf, 'CBrBb'=CBrBb, 'CBrF'=CBrF, 'CBrFf'=CBrFf, 'CBrFb'=CBrFb,
    'CBrTf'=CBrTf, 'CGoB'=CGoB, 'CGoBf'=CGoBf, 'CGoBb'=CGoBb, 'CGoF'=CGoF, 'CGoFf'=CGoFf, 
    'CGoFb'=CGoFb, 'CGoTf'=CGoTf, 'CSKB'=CSKB, 'CSKBf'=CSKBf, 'CSKBb'=CSKBb, 'CSKF'=CSKF,
    'CSKFf'=CSKFf, 'CSKFb'=CSKFb, 'CSKTf'=CSKTf, 'CBoB'=CBoB, 'CBoBf'=CBoBf, 'CBoBb'=CBoBb,
    'CBoF'=CBoF, 'CBoFf'=CBoFf, 'CBoFb'=CBoFb, 'CBoTf'=CBoTf,     
    
    'Cblood'=Cblood,  'Cplasma'=Cplasma, 'Mblood'=Mblood, 'Mkidney'=Mkidney, 'Mliver'=Mliver, 'Mbrain'=Mbrain,
    
    'Ckidney'=Ckidney, 'Cliver'=Cliver,
    'Cstomach'=Cstomach, 'Cintestine'=Cintestine, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose,
    'Clungs'=Clungs, 'Clungtissue'=Clungtissue, 'Crest'=Crest, 'Ccarcass'=Ccarcass,
    'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain, 
     'Cgonads'=Cgonads, 'Cskin'=Cskin, 'Cbones'=Cbones, 
    
    
    'Ckidney_bled'=Ckidney_bled, 'Cliver_bled'=Cliver_bled,
    'Cstomach_bled'=Cstomach_bled, 'Cintestine_bled'=Cintestine_bled, 'Cmuscle_bled'=Cmuscle_bled, 'Cadipose_bled'=Cadipose_bled,
    'Clungs_bled'=Clungs_bled, 'Clungtissue_bled'=Clungtissue_bled, 'Crest_bled'=Crest_bled, 'Ccarcass_bled'=Ccarcass_bled,
    'Cspleen_bled'=Cspleen_bled, 'Cheart_bled'=Cheart_bled, 'Cbrain_bled'=Cbrain_bled, 
    'Cgonads_bled'=Cgonads_bled, 'Cskin_bled'=Cskin_bled, 'Cbones_bled'=Cbones_bled, 
    
    
    "CBile" = CBile, 'CalveolarLF' = CalveolarLF, "Cfeces" = Cfeces,  "Curine" = Curine
     
    
    )
    
  })
}

#Initial condition for each compartment.

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
    CalbBrFf<- CalbBrF_init; MBrFf<- 0; MBrFb<- 0;
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
    
    CFabpLTf<- CFabpLT_init; CalbLuAFf<- CalbLuAF_init;
    MLTf<- 0; MLTb<- 0;  MPTCb = 0; MDALCb = 0; MDTCb = 0
    MCDCb = 0; MKTrestb = 0; MBile <-0; MSTTf <- 0;  MINTf <- 0; MLuTf <- 0; 
    MLuAFf <- 0; MLuAFb<- 0; MSPTf <- 0; MHTf <- 0;  MBrTf <- 0;
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
             'CalbBrFf' = CalbBrFf, 'CalbGoFf' = CalbGoFf, 'CalbSKFf' = CalbSKFf,
             'CalbBoFf' = CalbBoFf, 
             'Ca2uPTCf' = Ca2uPTCf, 'Ca2uDALCf' =  Ca2uDALCf,
             'Ca2uDTCf' =  Ca2uDTCf, 'Ca2uCDCf' =  Ca2uCDCf, 'Ca2uKTrestf' = Ca2uKTrestf,
             'CFabpPTCf' =  CFabpPTCf, 'CFabpDALCf' = CFabpDALCf, 'CFabpDTCf' = CFabpDTCf,
             'CFabpCDCf' = CFabpCDCf, 'CFabpKTrestf' =CFabpKTrestf,
             'CFabpLTf' = CFabpLTf, 'CalbLuAFf' = CalbLuAFf,
             
             
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
             
             'MLuAFb'=MLuAFb,
             
             'MArtf'=MArtf, 'MVenf'=MVenf, 'MKBf'=MKBf, 
             'MKFf'=MKFf,  'MPTCf' = MPTCf, 'MDALCf' = MDALCf,  
             'MDTCf' = MDTCf, 'MCDCf' = MCDCf,
             'MKTrestf' = MKTrestf, 'MPT' = MPT,  'MDAL' = MDAL,
             'MDT' =  MDT, 'MCD' = MCD,
             
             'MBladder' = MBladder,  'MLBf'=MLBf, 
             'MLFf'=MLFf, 'MLTf'=MLTf, 'MBile'= MBile,
             
             'MSTBf'=MSTBf, 'MSTFf'=MSTFf, 'MSTTf'=MSTTf, 'MSTL'=MSTL,
             'MINBf'=MINBf, 'MINFf'=MINFf, 'MINTf'=MINTf,'MINL'=MINL,
             
             'MMBf'=MMBf, 'MMFf'=MMFf, 'MMTf'=MMTf,
             'MABf'=MABf, 'MAFf'=MAFf, 'MATf'=MATf, 
             'MRBf'=MRBf, 'MRFf'=MRFf,'MRTf'=MRTf,
             'MLuBf'=MLuBf, 'MLuFf'=MLuFf,'MLuTf'=MLuTf,'MLuAFf' = MLuAFf,
             
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
                     "sex" = sex)
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
                     "sex" = sex)
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
                     "sex" = sex)
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
                     "sex" = sex)
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
                     "sex" = sex)
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
                     "sex" = sex)
  params[[6]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
  BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 12 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  
  params[[7]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
  BW <- 0.2  # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "F" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  
  params[[8]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 9th case, i.e. Kim (2016) ORAL male blood
  BW <- 0.25  #kg, from Kim et al. 2018
  admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  
  params[[9]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 10th case, i.e. Kim (2016) IV male blood
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
                     "sex" = sex)
  
  params[[10]] <- create_fixed_params(user_input)
  
  
  # Set up simulations for the 11th case, i.e.Kemper 2003 (Worley) ORAL female  LOW
  sex <- "F"
  BW <- 0.2 #kg
  admin.type <-"oral"
  admin.time <- 0
  #Female, oral 1mg/kg dose
  admin.dose <- 1 * BW * 1000 #ug
  params[[11]] <-   create_fixed_params(list('BW'=BW,
                                             "admin.dose"= admin.dose,
                                             "admin.time" = admin.time, 
                                             "admin.type" = admin.type,
                                             "sex" = sex))
  
  # Set up simulations for the 12th case, i.e.Kemper 2003 (Worley) ORAL female MEDIUM
  admin.dose <- 5 * BW * 1000 #ug
  params[[12]] <-   create_fixed_params(list('BW'=BW,
                                             "admin.dose"= admin.dose,
                                             "admin.time" = admin.time, 
                                             "admin.type" = admin.type,
                                             "sex" = sex))
  
  # Set up simulations for the 13th case, i.e.Kemper 2003 (Worley) ORAL female HIGH
  admin.dose <- 25 * BW * 1000 #ug
  params[[13]] <-   create_fixed_params(list('BW'=BW,
                                             "admin.dose"= admin.dose,
                                             "admin.time" = admin.time, 
                                             "admin.type" = admin.type,
                                             "sex" = sex))
  
  # Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL male LOW
  sex <- "M"
  BW <- 0.3 #kg
  admin.type <-"oral"
  admin.time <- 0
  #Male, oral 1mg/kg dose
  admin.dose <- 1 * BW * 1000 #ug
  params[[14]] <-   create_fixed_params(list('BW'=BW,
                                             "admin.dose"= admin.dose,
                                             "admin.time" = admin.time, 
                                             "admin.type" = admin.type,
                                             "sex" = sex))
  
  # Set up simulations for the 15th case, i.e.Kemper 2003 (Worley) ORAL male MEDIUM
  admin.dose <- 5 * BW * 1000 #ug
  params[[15]] <-   create_fixed_params(list('BW'=BW,
                                             "admin.dose"= admin.dose,
                                             "admin.time" = admin.time, 
                                             "admin.type" = admin.type,
                                             "sex" = sex))
  
  # Set up simulations for the 16th case, i.e.Kemper 2003 (Worley) ORAL male HIGH
  admin.dose <- 25 * BW * 1000 #ug
  params[[16]] <-   create_fixed_params(list('BW'=BW,
                                             "admin.dose"= admin.dose,
                                             "admin.time" = admin.time, 
                                             "admin.type" = admin.type,
                                             "sex" = sex))
  
  # Set up simulations for the 17th case, i.e. Dzierlenga 2021, IV male serum
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 6 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "iv"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[17]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 18st case, i.e. Dzierlenga 2021, ORAL male serum low
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 6 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[18]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 19nd case, i.e. Dzierlenga 2021, ORAL male serum medium
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 12 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[19]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 20d case, i.e. Dzierlenga 2021, ORAL male serum high
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 48 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[20]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 21th case, i.e. Dzierlenga 2021, IV female serum
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "iv"
  sex <- "F"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[21]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 22th case, i.e. Dzierlenga 2021, ORAL female serum low
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "F"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[22]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 23th case, i.e. Dzierlenga 2021, ORAL female serum medium
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "F"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[23]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 24th case, i.e. Dzierlenga 2021, ORAL female serum high
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  admin.dose_per_g <- 320 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "F"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[24]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 25th case, i.e. Kim (2016) ORAL female blood
  BW <- 0.25  #kg, from Kim et al. 2018
  admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
  admin.time <- 0 #time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "F" 
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "sex" = sex)
  params[[25]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 26th case, i.e. Kim (2016) IV male blood
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
                     "sex" = sex)
  params[[26]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 27th case, i.e. Gustafsson (2022) oral male blood
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
                     "sex" = sex)
  params[[27]] <- create_fixed_params(user_input)
  
  # Set up simulations for the 28st case, i.e. Gustafsson (2022) Inhalation male tissues
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
                     "sex" = sex )
  params[[28]] <- create_fixed_params(user_input)
  
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
  score <- rep(NA, N_data+9)
  N_points <- 0 # total number of points
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)
  tissue_multiplier <- 2
  last_point_multiplier <- 10
  dzirlenga_female_liver_multiplier <- 2
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Kudo et al. 2007 high dose
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kudo_high <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Cblood",
                                                                                        "Cliver","Ckidney", "Ccarcass",
                                                                                        "Clungs", "Cspleen", "Cheart",
                                                                                        "Cbrain", "Cgonads", "Cstomach", "Cintestine")]/1000) #convert ug/kg to ug/g
  
  obs_kudo_high <- list(exp_data[exp_data$Tissue == "Blood", "concentration"],
                        exp_data[exp_data$Tissue == "Liver", "concentration"],
                        exp_data[exp_data$Tissue == "Kidney", "concentration"],
                        exp_data[exp_data$Tissue == "Carcass", "concentration"],
                        exp_data[exp_data$Tissue == "Lung", "concentration"],
                        exp_data[exp_data$Tissue == "Spleen", "concentration"],
                        exp_data[exp_data$Tissue == "Heart", "concentration"],
                        exp_data[exp_data$Tissue == "Brain", "concentration"],
                        exp_data[exp_data$Tissue == "Gonads", "concentration"],
                        exp_data[exp_data$Tissue == "Stomach", "concentration"],
                        exp_data[exp_data$Tissue == "Intestine", "concentration"])
  N_kudo_high <- sum(unlist(lapply( obs_kudo_high,length)))
  score[1] <- N_kudo_high* tissue_multiplier*AAFE(predictions = preds_kudo_high, observations = obs_kudo_high)
  N_points <- N_points + N_kudo_high* tissue_multiplier
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Kudo et al. 2007 low dose
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kudo_low <- solution[solution$time %in% unique(exp_data$time), c("Cblood",
                                                                         "Cliver","Ckidney", "Ccarcass",
                                                                         "Clungs", "Cspleen", "Cheart",
                                                                         "Cbrain", "Cgonads", "Cstomach", "Cintestine")]
  preds_kudo_low<- as.data.frame(preds_kudo_low /1000) #convert ug/kg to ug/g
  
  
  obs_kudo_low <- list(exp_data[exp_data$Tissue == "Blood", "concentration"],
                       exp_data[exp_data$Tissue == "Liver", "concentration"],
                       exp_data[exp_data$Tissue == "Kidney", "concentration"],
                       exp_data[exp_data$Tissue == "Carcass", "concentration"],
                       exp_data[exp_data$Tissue == "Lung", "concentration"],
                       exp_data[exp_data$Tissue == "Spleen", "concentration"],
                       exp_data[exp_data$Tissue == "Heart", "concentration"],
                       exp_data[exp_data$Tissue == "Brain", "concentration"],
                       exp_data[exp_data$Tissue == "Gonads", "concentration"],
                       exp_data[exp_data$Tissue == "Stomach", "concentration"],
                       exp_data[exp_data$Tissue == "Intestine", "concentration"])
  
  N_kudo_low <- sum(unlist(lapply(obs_kudo_low, length)))  
  score[2] <- N_kudo_low* tissue_multiplier*AAFE(predictions = preds_kudo_low, observations = obs_kudo_low)  
  N_points <- N_points + N_kudo_low* tissue_multiplier
  
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of kim (2016) IV male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_IV_Mtissues <- solution[solution$time %in% unique(exp_data$time), c("Cliver","Ckidney", 
                                                                                "Clungs", "Cspleen", "Cheart"
  )]
  preds_kim_IV_Mtissues<- as.data.frame(preds_kim_IV_Mtissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_IV_Mtissues <- list(exp_data[exp_data$Tissue == "Liver", "concentration"],
                              exp_data[exp_data$Tissue == "Kidney", "concentration"],
                              exp_data[exp_data$Tissue == "Lung", "concentration"],
                              exp_data[exp_data$Tissue == "Spleen", "concentration"],
                              exp_data[exp_data$Tissue == "Heart", "concentration"])
  
  N_kim_IV_Mtissues <- sum(unlist(lapply(obs_kim_IV_Mtissues, length)))  
  score[3] <- N_kim_IV_Mtissues* tissue_multiplier*AAFE(predictions = preds_kim_IV_Mtissues, observations = obs_kim_IV_Mtissues)  
  N_points <- N_points + N_kim_IV_Mtissues*tissue_multiplier
  
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of kim (2016) ORAL male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_OR_Mtissues <- solution[solution$time %in% unique(exp_data$time), c("Cliver","Ckidney", 
                                                                                "Clungs", "Cspleen", "Cheart" )]
  
  preds_kim_OR_Mtissues<- as.data.frame(preds_kim_OR_Mtissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_OR_Mtissues <- list(exp_data[exp_data$Tissue == "Liver", "concentration"],
                              exp_data[exp_data$Tissue == "Kidney", "concentration"],
                              exp_data[exp_data$Tissue == "Lung", "concentration"],
                              exp_data[exp_data$Tissue == "Spleen", "concentration"],
                              exp_data[exp_data$Tissue == "Heart", "concentration"])
  
  N_kim_OR_Mtissues <- sum(unlist(lapply(obs_kim_OR_Mtissues, length))) 
   score[4] <- N_kim_OR_Mtissues*tissue_multiplier*AAFE(predictions = preds_kim_OR_Mtissues, observations = obs_kim_OR_Mtissues)  
   N_points <- N_points + N_kim_OR_Mtissues*tissue_multiplier
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of kim (2016) ORAL female tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_IV_Ftissues <- solution[solution$time %in% unique(exp_data$time), c("Cliver","Ckidney", 
                                                                                "Clungs", "Cspleen", "Cheart" )]
  
  preds_kim_IV_Ftissues<- as.data.frame(preds_kim_IV_Ftissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_IV_Ftissues <- list(exp_data[exp_data$Tissue == "Liver", "concentration"],
                              exp_data[exp_data$Tissue == "Kidney", "concentration"],
                              exp_data[exp_data$Tissue == "Lung", "concentration"],
                              exp_data[exp_data$Tissue == "Spleen", "concentration"],
                              exp_data[exp_data$Tissue == "Heart", "concentration"])
  
  N_kim_IV_Ftissues <- sum(unlist(lapply(obs_kim_IV_Ftissues, length)))  
  score[5] <- N_kim_IV_Ftissues*tissue_multiplier*AAFE(predictions = preds_kim_IV_Ftissues, observations = obs_kim_IV_Ftissues)  
  N_points <- N_points + N_kim_IV_Ftissues*tissue_multiplier
  
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of kim (2016) ORAL female tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_kim_OR_Ftissues <- solution[solution$time %in% unique(exp_data$time), c("Cliver","Ckidney", 
                                                                                "Clungs", "Cspleen", "Cheart" )]
  
  preds_kim_OR_Ftissues<- as.data.frame(preds_kim_OR_Ftissues /1000) #convert ug/kg to ug/g
  
  
  obs_kim_OR_Ftissues <- list(exp_data[exp_data$Tissue == "Liver", "concentration"],
                              exp_data[exp_data$Tissue == "Kidney", "concentration"],
                              exp_data[exp_data$Tissue == "Lung", "concentration"],
                              exp_data[exp_data$Tissue == "Spleen", "concentration"],
                              exp_data[exp_data$Tissue == "Heart", "concentration"])
  
  N_kim_OR_Ftissues <- sum(unlist(lapply(obs_kim_OR_Ftissues, length))) 
  score[6] <- N_kim_OR_Ftissues*tissue_multiplier*AAFE(predictions = preds_kim_OR_Ftissues, observations = obs_kim_OR_Ftissues)
  N_points <- N_points + N_kim_OR_Ftissues*tissue_multiplier
  
  ##########################
  #-------------------------
  # Dzierlenga ORAL male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
  BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[7]])
  params <- c(fixed_params[[7]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,864,1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df7=========================================================
  exp_data <- dataset$df7 # retrieve data of Dzierlenga (2021) ORAL male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cliver","Ckidney","Cbrain"  )
  
  preds_dzi_OR_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Mtissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_dzi_OR_Mtissues <- list( exp_data[exp_data$Tissue == "Liver", "concentration"],
                               exp_data[exp_data$Tissue == "Kidney", "concentration"],
                               exp_data[exp_data$Tissue == "Brain", "concentration"]) 
  
    N_dzi_OR_Mtissues <- sum(unlist(lapply(obs_dzi_OR_Mtissues, length)))  
    score[7] <- N_dzi_OR_Mtissues*AAFE(predictions = preds_dzi_OR_Mtissues, observations = obs_dzi_OR_Mtissues)  
    N_points <- N_points + N_dzi_OR_Mtissues*tissue_multiplier  
    
    N_dzi_OR_Mtissues2 <- length(preds_dzi_OR_Mtissues[[1]])
    score[39] <- N_dzi_OR_Mtissues2*AAFE(predictions = preds_dzi_OR_Mtissues[[1]], observations = obs_dzi_OR_Mtissues[[1]])
    N_points <- N_dzi_OR_Mtissues2 + N_points
    
    N_dzi_OR_Mtissues2 <- length(preds_dzi_OR_Mtissues[[2]])
    score[41] <- N_dzi_OR_Mtissues2*AAFE(predictions = preds_dzi_OR_Mtissues[[2]], observations = obs_dzi_OR_Mtissues[[2]])
    N_points <- N_dzi_OR_Mtissues2 + N_points
  ##########################
  #-------------------------
  # Dzierlenga ORAL female tissues
  #-------------------------
  ##########################
  # Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
  BW <- 0.2  # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[8]])
  params <- c(fixed_params[[8]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,24,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df8=========================================================
  
  exp_data <- dataset$df8 # retrieve data of Dzierlenga (2021) ORAL female tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cliver","Ckidney","Cbrain" )
  MW <-414.07
  
  preds_dzi_OR_Ftissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Ftissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_dzi_OR_Ftissues <- list( exp_data[exp_data$Tissue == "Liver", "concentration"],
                               exp_data[exp_data$Tissue == "Kidney", "concentration"],
                               exp_data[exp_data$Tissue == "Brain", "concentration"]) 
  
    N_dzi_OR_Ftissues <- sum(unlist(lapply(obs_dzi_OR_Ftissues, length)))  
    score[8] <- N_dzi_OR_Ftissues*AAFE(predictions = preds_dzi_OR_Ftissues, observations = obs_dzi_OR_Ftissues)  
    N_points <- N_points + N_dzi_OR_Ftissues*tissue_multiplier   
      
   
    N_dzi_OR_Ftissues2 <- length(preds_dzi_OR_Ftissues[[1]])
    score[40] <- dzirlenga_female_liver_multiplier*N_dzi_OR_Ftissues2*AAFE(predictions = preds_dzi_OR_Ftissues[[1]], observations = obs_dzi_OR_Ftissues[[1]])
    N_points <- dzirlenga_female_liver_multiplier*N_dzi_OR_Ftissues2 + N_points
    
    score[43] <-last_point_multiplier*AAFE(predictions = preds_dzi_OR_Ftissues[[1]][length(preds_dzi_OR_Ftissues[[1]])],
                                           observations = obs_dzi_OR_Ftissues[[1]][length(obs_dzi_OR_Ftissues[[1]])])
    N_points <- last_point_multiplier + N_points
  ##########################
  #-------------------------
  # Kim ORAL male blood
  #-------------------------
  ##########################
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df9=========================================================
  
  exp_data <- dataset$df9 # retrieve data of Kim (2016) ORAL male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_kim_OR_Mblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_kim_OR_Mblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_kim_OR_Mblood <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
    N_kim_OR_Mblood <- sum(unlist(lapply(obs_kim_OR_Mblood, length)))  
    score[9] <- N_kim_OR_Mblood*AAFE(predictions = preds_kim_OR_Mblood, observations = obs_kim_OR_Mblood)  
    N_points <- N_points + N_kim_OR_Mblood  
    
    score[35] <- last_point_multiplier*AAFE(predictions = list(preds_kim_OR_Mblood[[1]][length(preds_kim_OR_Mblood[[1]])]),                    
     observations = list(obs_kim_OR_Mblood[[1]][length(obs_kim_OR_Mblood[[1]])]))
     N_points <- N_points + last_point_multiplier
     
     score[42] <- last_point_multiplier*AAFE(predictions = list(preds_kim_OR_Mblood[[1]][2]),                    
                                             observations = list(obs_kim_OR_Mblood[[1]][2]))
     N_points <- N_points + last_point_multiplier
  ##########################
  #-------------------------
  # Kim IV male blood
  #-------------------------
  ##########################
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df10=========================================================
  
  exp_data <- dataset$df10 # retrieve data of Kim (2016) IV male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_kim_IV_Mblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_kim_IV_Mblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_kim_IV_Mblood <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  
    N_kim_IV_Mblood <- sum(unlist(lapply(obs_kim_IV_Mblood, length)))  
    score[10] <- N_kim_IV_Mblood*AAFE(predictions = preds_kim_IV_Mblood, observations = obs_kim_IV_Mblood)  
    N_points <- N_points + N_kim_IV_Mblood 
    
     score[36] <- last_point_multiplier*AAFE(predictions = list(preds_kim_IV_Mblood[[1]][length(preds_kim_IV_Mblood[[1]])]),
                         observations = list(obs_kim_IV_Mblood[[1]][length(preds_kim_IV_Mblood[[1]])]))
     N_points <- N_points + last_point_multiplier
  ##########################
  #-------------------------
  # Kemper 2003 (Worley)
  #-------------------------
  ##########################
  # Set up simulations for the 11th case, i.e.Kemper 2003 (Worley) ORAL female LOW
  
  sex <- "F"
  BW <- 0.2 #kg
  sample_time <- seq(0,192,1)
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[11]])
  params <- c(fixed_params[[11]], variable_params)
  events <- create.events(params)
  inits <- create.inits (params)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df11a=========================================================
  exp_data <- dataset$df11a # retrieve data of Kemper 2003  (from Worley) ORAL female urine LOW
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Murine")
  
  preds_Kemp_OR_Furine_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    
    preds_Kemp_OR_Furine_low [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Furine_low <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Furine_low <- (obs_Kemp_OR_Furine_low/100)*params$admin.dose/1000
  
  N_Kemp_OR_Furine_low <- length(obs_Kemp_OR_Furine_low)  
  score[11] <- N_Kemp_OR_Furine_low*AAFE(predictions = preds_Kemp_OR_Furine_low, observations = obs_Kemp_OR_Furine_low) 
   N_points <- N_points + N_Kemp_OR_Furine_low
  
  #======================================df11b=========================================================
  exp_data <- dataset$df11b # retrieve data of Kemper 2003  (from Loccisano) ORAL female feces LOW
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Mfeces")
  
  preds_Kemp_OR_Ffeces_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Ffeces_low [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Ffeces_low <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Ffeces_low <- (obs_Kemp_OR_Ffeces_low/100)*params$admin.dose/1000
  
  N_Kemp_OR_Ffeces_low <- length(obs_Kemp_OR_Ffeces_low)  
  score[12] <- N_Kemp_OR_Ffeces_low*AAFE(predictions = preds_Kemp_OR_Ffeces_low, observations = obs_Kemp_OR_Ffeces_low)  
  N_points <- N_points + N_Kemp_OR_Ffeces_low
  
  
  # Set up simulations for the 12th case, i.e.Kemper 2003 (Worley) ORAL female  MEDIUM
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[12]])
  params <- c(fixed_params[[12]], variable_params)
  events <- create.events(params)
  inits <- create.inits (params)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  
  #======================================df12a=========================================================
  exp_data <- dataset$df12a # retrieve data of Kemper 2003  (Worley) ORAL female urine MEDIUM
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Murine")
  
  preds_Kemp_OR_Furine_med <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Furine_med [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Furine_med <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Furine_med <- (obs_Kemp_OR_Furine_med/100)*params$admin.dose/1000
  
  N_Kemp_OR_Furine_med <- length(obs_Kemp_OR_Furine_med)  
  score[13] <- N_Kemp_OR_Furine_med*AAFE(predictions = preds_Kemp_OR_Furine_med, observations = obs_Kemp_OR_Furine_med)  
  N_points <- N_points + N_Kemp_OR_Furine_med
  
  
  #======================================df12b=========================================================
  exp_data <- dataset$df12b # retrieve data of Kemper 2003  (Loccisano) ORAL female feces MEDIUM
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Mfeces")
  
  preds_Kemp_OR_Ffeces_med <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Ffeces_med [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Ffeces_med <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Ffeces_med <- (obs_Kemp_OR_Ffeces_med/100)*params$admin.dose/1000
  
  N_Kemp_OR_Ffeces_med <- length(obs_Kemp_OR_Ffeces_med)  
  score[14] <- N_Kemp_OR_Ffeces_med*AAFE(predictions = preds_Kemp_OR_Ffeces_med, observations = obs_Kemp_OR_Ffeces_med) 
   N_points <- N_points + N_Kemp_OR_Ffeces_med
  
  ##########################################################################################
  # Set up simulations for the 13ath case, i.e.Kemper 2003 (Worley) ORAL female HIGH
  sample_time <- seq(0,672,1)
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[13]])
  params <- c(fixed_params[[13]], variable_params)
  events <- create.events(params)
  inits <- create.inits (params)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df13a=========================================================
  
  exp_data <- dataset$df13a # retrieve data of Kemper 2003  (from Worley) ORAL female urine HIGH
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Murine")
  
  preds_Kemp_OR_Furine_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Furine_high [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Furine_high <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Furine_high <- (obs_Kemp_OR_Furine_high/100)*params$admin.dose/1000
  
  N_Kemp_OR_Furine_high <- length(obs_Kemp_OR_Furine_high)  
  score[15] <- N_Kemp_OR_Furine_high*AAFE(predictions = preds_Kemp_OR_Furine_high, observations = obs_Kemp_OR_Furine_high) 
   N_points <- N_points + N_Kemp_OR_Furine_high
  
  #======================================df13b=========================================================
  exp_data <- dataset$df13b # retrieve data of Kemper 2003  (Loccisano) ORAL female feces HIGH
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Mfeces")
  
  preds_Kemp_OR_Ffeces_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Ffeces_high [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Ffeces_high <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Ffeces_high <- (obs_Kemp_OR_Ffeces_high/100)*params$admin.dose/1000
  
  N_Kemp_OR_Ffeces_high <- length(obs_Kemp_OR_Ffeces_high)  
  score[16] <- (N_Kemp_OR_Ffeces_high/2)*AAFE(predictions = preds_Kemp_OR_Ffeces_high, observations = obs_Kemp_OR_Ffeces_high) 
   N_points <- N_points + N_Kemp_OR_Ffeces_high/2
  
  ###################################################################################################
  # Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL male LOW
  
  sex <- "M"
  BW <- 0.3 #kg
  sample_time <- seq(0,673,1)
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[14]])
  params <- c(fixed_params[[14]], variable_params)
  inits <- create.inits (params)
  events <- create.events(params)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df14a=========================================================
  exp_data <- dataset$df14a # retrieve data of Kemper 2003  (Worley) ORAL male urine LOW
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Murine")
  
  preds_Kemp_OR_Murine_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Murine_low [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Murine_low <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Murine_low <- (obs_Kemp_OR_Murine_low/100)*params$admin.dose/1000
  
  N_Kemp_OR_Murine_low <- length(obs_Kemp_OR_Murine_low)  
  score[17] <- ( N_Kemp_OR_Murine_low/3)*AAFE(predictions = preds_Kemp_OR_Murine_low, observations = obs_Kemp_OR_Murine_low) 
   N_points <- N_points + N_Kemp_OR_Murine_low/3
  
  #======================================df14b=========================================================
  exp_data <- dataset$df14b # retrieve data of Kemper 2003  (Loccisano) ORAL male feces Low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Mfeces")
  
  preds_Kemp_OR_Mfeces_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Mfeces_low [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Mfeces_low <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Mfeces_low <- (obs_Kemp_OR_Mfeces_low/100)*params$admin.dose/1000
  
  N_Kemp_OR_Mfeces_low <- length(obs_Kemp_OR_Mfeces_low)  
  score[18] <- N_Kemp_OR_Mfeces_low*AAFE(predictions = preds_Kemp_OR_Mfeces_low, observations = obs_Kemp_OR_Mfeces_low) 
   N_points <- N_points + N_Kemp_OR_Mfeces_low
  
  ####################################################################################################
  # Set up simulations for the 15th case, i.e.Kemper 2003 (Worley) ORAL male MEDIUM
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[15]])
  params <- c(fixed_params[[15]], variable_params)
  inits <- create.inits (params)
  events <- create.events(params)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df15a========================================================= 
  
  exp_data <- dataset$df15a # retrieve data of Kemper 2003  (Worley) ORAL male urine MEDIUM
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Murine")
  
  preds_Kemp_OR_Murine_med <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Murine_med [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Murine_med <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Murine_med <- (obs_Kemp_OR_Murine_med/100)*params$admin.dose/1000
  
  N_Kemp_OR_Murine_med <- length(obs_Kemp_OR_Murine_med)  
  score[19] <- (N_Kemp_OR_Murine_med/3)*AAFE(predictions = preds_Kemp_OR_Murine_med, observations = obs_Kemp_OR_Murine_med) 
   N_points <- N_points + N_Kemp_OR_Murine_med/3
  
  
  #======================================df15b=========================================================
  exp_data <- dataset$df15b # retrieve data of Kemper 2003  (Loccisano) ORAL male feces MEDIUM
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Mfeces")
  
  preds_Kemp_OR_Mfeces_med <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Mfeces_med [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Mfeces_med <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Mfeces_med <- (obs_Kemp_OR_Mfeces_med/100)*params$admin.dose/1000
  
  N_Kemp_OR_Mfeces_med <- length(obs_Kemp_OR_Mfeces_med) 
   score[20] <- (N_Kemp_OR_Mfeces_med/2)*AAFE(predictions = preds_Kemp_OR_Mfeces_med, observations = obs_Kemp_OR_Mfeces_med)  
   N_points <- N_points + N_Kemp_OR_Mfeces_med/2
  
  ######################################################################################################
  # Set up simulations for the 16th case, i.e.Kemper 2003 (Worley) ORAL male  HIGH
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[16]])
  params <- c(fixed_params[[16]], variable_params)
  inits <- create.inits (params)
  events <- create.events(params)
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params, events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df16a=========================================================
  exp_data <- dataset$df16a # retrieve data of Kemper 2003  (Worley) ORAL male urine HIGH
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Murine")
  
  preds_Kemp_OR_Murine_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Murine_high [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Murine_high <- c(exp_data[exp_data$Tissue == "Urine", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Murine_high <- (obs_Kemp_OR_Murine_high/100)*params$admin.dose/1000
  
  N_Kemp_OR_Murine_high <- length(obs_Kemp_OR_Murine_high)  
  score[21] <- (N_Kemp_OR_Murine_high/3)*AAFE(predictions = preds_Kemp_OR_Murine_high, observations = obs_Kemp_OR_Murine_high) 
   N_points <- N_points + (N_Kemp_OR_Murine_high/3)
  
  #======================================df16b=========================================================
  exp_data <- dataset$df16b # retrieve data of Kemper 2003  (Loccisano) ORAL male feces
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Mfeces")
  
  preds_Kemp_OR_Mfeces_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    rounded_time <- round(exp_time)
    rounded_soltime <- round(solution$time)
    
    preds_Kemp_OR_Mfeces_high [[i]] <- solution[rounded_soltime %in% rounded_time, column_names[i]]/1000
  }
  
  obs_Kemp_OR_Mfeces_high <- c(exp_data[exp_data$Tissue == "Feces", "concentration"])
  # Estimate cumulative fecal mass
  obs_Kemp_OR_Mfeces_high <- (obs_Kemp_OR_Mfeces_high/100)*params$admin.dose/1000
  
  N_Kemp_OR_Mfeces_high <- length(obs_Kemp_OR_Mfeces_high)  
  score[22] <- (N_Kemp_OR_Mfeces_high/2)*AAFE(predictions = preds_Kemp_OR_Mfeces_high, observations = obs_Kemp_OR_Mfeces_high)  
  N_points <- N_points + N_Kemp_OR_Mfeces_high/2
  
  
  ##########################
  #-------------------------
  # Dzierlenga 2021 IV male serum
  #-------------------------
  ##########################
  
  # Set up simulations for the 17th case, i.e. Dzierlenga 2021, IV male serum
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[17]])
  params <- c(fixed_params[[17]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- c(0, 0.083, 0.25, 0.5, 1, 3, 6, seq(12, 1200, 4))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df17=========================================================
  
  exp_data <- dataset$df17 # retrieve data of Dzierlenga 2021, IV male serum
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_IV_Mserum <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_IV_Mserum [[i]] <- solution[solution$time %in% exp_time, column_names[i]] /1000
  }
  
  #we assume that clotting factors are negligible amount
  obs_dzi_IV_Mserum <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])
  
  N_dzi_IV_Mserum <- sum(unlist(lapply(obs_dzi_IV_Mserum, length)))  
  score[23] <- N_dzi_IV_Mserum*AAFE(predictions = preds_dzi_IV_Mserum, observations = obs_dzi_IV_Mserum)  
  N_points <- N_points + N_dzi_IV_Mserum
  
  ##########################
  #-------------------------
  # Dzierlenga 2021 ORAL male serum low
  #-------------------------
  ##########################
  
  # Set up simulations for the 18th case, i.e. Dzierlenga 2021, ORAL male serum low
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[18]])
  params <- c(fixed_params[[18]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df18=========================================================
  
  exp_data <- dataset$df18 # retrieve data of Dzierlenga 2021, ORAL male serum low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_OR_Mserum_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Mserum_low [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_OR_Mserum_low <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])
  
  N_dzi_OR_Mserum_low <- sum(unlist(lapply(obs_dzi_OR_Mserum_low, length)))  
  score[24] <- N_dzi_OR_Mserum_low*AAFE(predictions = preds_dzi_OR_Mserum_low, observations = obs_dzi_OR_Mserum_low)  
  N_points <- N_points + N_dzi_OR_Mserum_low
  
  ##########################
  #-------------------------
  # Dzierlenga 2021 ORAL male serum medium
  #-------------------------
  ##########################
  
  # Set up simulations for the 19th case, i.e. Dzierlenga 2021, ORAL male serum medium
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[19]])
  params <- c(fixed_params[[19]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df19=========================================================
  exp_data <- dataset$df19 # retrieve data of Dzierlenga 2021, ORAL male serum medium
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_OR_Mserum_medium <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Mserum_medium [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_OR_Mserum_medium <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])
  

  N_dzi_OR_Mserum_medium <- sum(unlist(lapply(obs_dzi_OR_Mserum_medium, length)))  
  score[25] <- N_dzi_OR_Mserum_medium* AAFE(predictions = preds_dzi_OR_Mserum_medium, observations = obs_dzi_OR_Mserum_medium)
  N_points <- N_points + N_dzi_OR_Mserum_medium
  
  ##########################
  #-------------------------
  # Dzierlenga 2021 ORAL male serum high
  #-------------------------
  ##########################
  # Set up simulations for the 20th case, i.e. Dzierlenga 2021, ORAL male serum high
  BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[20]])
  params <- c(fixed_params[[20]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df20=========================================================
  exp_data <- dataset$df20 # retrieve data of Dzierlenga 2021, ORAL male serum high
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_OR_Mserum_high <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Mserum_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_OR_Mserum_high <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])

  N_dzi_OR_Mserum_high <- sum(unlist(lapply(obs_dzi_OR_Mserum_high, length)))  
  score[26] <- N_dzi_OR_Mserum_high* AAFE(predictions = preds_dzi_OR_Mserum_high, observations = obs_dzi_OR_Mserum_high)
  N_points <- N_points + N_dzi_OR_Mserum_high
  ##########################
  #-------------------------
  # Dzierlenga 2021 IV female serum 
  #-------------------------
  ##########################
  
  # Set up simulations for the 21th case, i.e. Dzierlenga 2021, IV female serum
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[21]])
  params <- c(fixed_params[[21]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- c(0, 0.083, 0.25, seq(0.5, 192, 0.5))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df21=========================================================
  exp_data <- dataset$df21 # retrieve data of Dzierlenga 2021, IV female serum
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_IV_Fserum <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_IV_Fserum [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_IV_Fserum <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])

  N_dzi_IV_Fserum <- sum(unlist(lapply(obs_dzi_IV_Fserum, length)))
  score[27]       <- N_dzi_IV_Fserum * AAFE(preds_dzi_IV_Fserum, obs_dzi_IV_Fserum)
  N_points        <- N_points + N_dzi_IV_Fserum
  ##########################
  #-------------------------
  # Dzierlenga 2021 ORAL female serum low
  #-------------------------
  ##########################
  
  # Set up simulations for the 22th case, i.e. Dzierlenga 2021, ORAL female serum low
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[22]])
  params <- c(fixed_params[[22]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- seq(0, 96, 0.25)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df22=========================================================
  exp_data <- dataset$df22 # retrieve data of Dzierlenga 2021, ORAL female serum low
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_OR_Fserum_low <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Fserum_low [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_OR_Fserum_low <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])

  N_dzi_OR_Fserum_low <- sum(unlist(lapply(obs_dzi_OR_Fserum_low, length)))
  score[28]            <- N_dzi_OR_Fserum_low * AAFE(preds_dzi_OR_Fserum_low, obs_dzi_OR_Fserum_low)
  N_points             <- N_points + N_dzi_OR_Fserum_low
  ##########################
  #-------------------------
  # Dzierlenga 2021 ORAL female serum medium
  #-------------------------
  ##########################
  
  # Set up simulations for the 23th case, i.e. Dzierlenga 2021, ORAL female serum medium
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[23]])
  params <- c(fixed_params[[23]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- c(0, 0.25, seq(1, 192, 0.5))
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df23=========================================================
  exp_data <- dataset$df23 # retrieve data of Dzierlenga 2021, ORAL female serum medium
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_OR_Fserum_medium<- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Fserum_medium [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_OR_Fserum_medium <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])

  N_dzi_OR_Fserum_medium <- sum(unlist(lapply(obs_dzi_OR_Fserum_medium, length)))
  score[29]              <- N_dzi_OR_Fserum_medium * AAFE(preds_dzi_OR_Fserum_medium, obs_dzi_OR_Fserum_medium)
  N_points               <- N_points + N_dzi_OR_Fserum_medium
  ##########################
  #-------------------------
  # Dzierlenga 2021 ORAL female serum high
  #-------------------------
  ##########################
  
  # Set up simulations for the 24th case, i.e. Dzierlenga 2021, ORAL female serum high
  BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
  sex <- "F"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[24]])
  params <- c(fixed_params[[24]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time <- seq(0, 96, 0.25)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df24=========================================================
  exp_data <- dataset$df24 # retrieve data of Dzierlenga 2021, ORAL female serum high
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_dzi_OR_Fserum_high<- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_dzi_OR_Fserum_high [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000 
  }
  
  obs_dzi_OR_Fserum_high <- list(exp_data[exp_data$Tissue == "Serum", "concentration"])

  N_dzi_OR_Fserum_high <- sum(unlist(lapply(obs_dzi_OR_Fserum_high, length)))
  score[30]            <- N_dzi_OR_Fserum_high * AAFE(preds_dzi_OR_Fserum_high, obs_dzi_OR_Fserum_high)
  N_points             <- N_points + N_dzi_OR_Fserum_high
  ##########################
  #-------------------------
  # Kim ORAL female blood
  #-------------------------
  ##########################
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df25=========================================================
  exp_data <- dataset$df25 # retrieve data of Kim (2016) ORAL male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_kim_OR_Fblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_kim_OR_Fblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_kim_OR_Fblood <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])

  
  N_kim_OR_Fblood <- sum(unlist(lapply(obs_kim_OR_Fblood, length)))
  score[31]       <- N_kim_OR_Fblood * AAFE(preds_kim_OR_Fblood, obs_kim_OR_Fblood)
  N_points        <- N_points + N_kim_OR_Fblood
  
  score[37] <- last_point_multiplier * AAFE(predictions = list(preds_kim_OR_Fblood[[1]][length(preds_kim_OR_Fblood[[1]])]), 
                    observations = list(obs_kim_OR_Fblood[[1]][length(obs_kim_OR_Fblood[[1]])]))
  N_points        <- N_points + last_point_multiplier
  
  ##########################
  #-------------------------
  # Kim IV female blood
  #-------------------------
  ##########################
  # Set up simulations for the 26th case, i.e. Kim (2016) IV male blood
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
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df26=========================================================
  exp_data <- dataset$df26 # retrieve data of Kim (2016) IV male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_kim_IV_Fblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_kim_IV_Fblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_kim_IV_Fblood <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])

  
  N_kim_IV_Fblood <- sum(unlist(lapply(obs_kim_IV_Fblood, length)))
  score[32]       <- N_kim_IV_Fblood * AAFE(preds_kim_IV_Fblood, obs_kim_IV_Fblood)
  N_points        <- N_points + N_kim_IV_Fblood
  
  score[38] <- last_point_multiplier* AAFE(predictions = list(preds_kim_IV_Fblood[[1]][length(preds_kim_IV_Fblood[[1]])]), 
                    observations = list(obs_kim_IV_Fblood[[1]][length(obs_kim_IV_Fblood[[1]])]))
  N_points        <- N_points + last_point_multiplier
  
  ##########################
  #-------------------------
  # Gustafsson Oral male blood
  #-------------------------
  ##########################
  # Set up simulations for the 27th case, i.e. Gustafsson (2022) oral male blood
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M" 
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[27]])
  params <- c(fixed_params[[27]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df27=========================================================
  exp_data <- dataset$df27 # retrieve data of Gustafsson (2022) Oral male blood
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")
  
  preds_gus_OR_Mblood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_OR_Mblood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_gus_OR_Mblood <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])

  N_gus_OR_Mblood <- sum(unlist(lapply(obs_gus_OR_Mblood, length)))
  score[33]       <- N_gus_OR_Mblood * AAFE(preds_gus_OR_Mblood, obs_gus_OR_Mblood)
  N_points        <- N_points + N_gus_OR_Mblood
  ##########################
  #-------------------------
  # Gustafsson Oral male tissues
  #-------------------------
  ##########################
  # Set up simulations for the 28st case, i.e. Gustafsson (2022) Inhalation male tissues
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[28]])
  params <- c(fixed_params[[28]], variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df28=========================================================
  exp_data <- dataset$df28 # retrieve data of Gustafsson (2022) oral male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CalveolarLF","Cliver","Clungtissue", "Ckidney")
  
  preds_gus_OR_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_gus_OR_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  obs_gus_OR_Mtissues <- list( exp_data[exp_data$Tissue == "ALF", "concentration"],
                               exp_data[exp_data$Tissue == "Liver", "concentration"],
                               exp_data[exp_data$Tissue == "Lung", "concentration"], 
                               exp_data[exp_data$Tissue == "Kidney", "concentration"]) 
  

  N_gus_OR_Mtissues <- sum(unlist(lapply(obs_gus_OR_Mtissues, length)))
  score[34]         <- N_gus_OR_Mtissues*tissue_multiplier * AAFE(preds_gus_OR_Mtissues, obs_gus_OR_Mtissues)
  N_points          <- N_points + N_gus_OR_Mtissues*tissue_multiplier
  ########################################################################################

  # Estimate final score
  
  final_score <- sum(score)/N_points
  if(is.na(final_score)){
    final_score <- 100
  }
  return(final_score)
  
}

################################################################################


setwd("/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")


MW <- 414.07 #g/mol
source("Goodness-of-fit-metrics.R")

# Read data
kudo_high_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_high_kudo_2007.xlsx")
kudo_low_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_low_kudo_2007.xlsx")
kim_IV_Mtissues <- openxlsx::read.xlsx("Data/PFOA_male_tissues_IV_kim_2016.xlsx")
kim_OR_Mtissues <- openxlsx::read.xlsx("Data/PFOA_male_tissues_ORAL_kim_2016.xlsx")
kim_IV_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_IV_kim_2016.xlsx")
kim_OR_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_ORAL_kim_2016.xlsx")
dzi_OR_Mtissues <- openxlsx::read.xlsx("Data/Dzierlenga_tissue_male_ORAL_2021.xlsx")
dzi_OR_Mtissues$Concentration_microM <- dzi_OR_Mtissues$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Ftissues <- openxlsx::read.xlsx("Data/Dzierlenga_tissue_female_ORAL_2021.xlsx")
dzi_OR_Ftissues$Concentration_microM <- dzi_OR_Ftissues$Concentration_microM* MW/1000 #convert from uM to ug/g
kim_OR_Mblood <- openxlsx::read.xlsx("Data/PFOA_male_blood_ORAL_kim_2016.xlsx")
kim_IV_Mblood <- openxlsx::read.xlsx("Data/PFOA_male_blood_IV_kim_2016.xlsx")
Lup_OR_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_Lupton_2020.xlsx")
Kemp_OR_Ffeces_high <- openxlsx::read.xlsx("Data/PFOA_Feces_female_oral_25_mg_per_kg-Loc.xlsx")
Kemp_OR_Mfeces_high <- openxlsx::read.xlsx("Data/PFOA_Feces_male_oral_25_mg_per_kg-Loc.xlsx")
Kemp_OR_Furine_low <- openxlsx::read.xlsx("Data/PFOA_Urine_female_oral_1_mg_per_kg.xlsx")
Kemp_OR_Furine_med <- openxlsx::read.xlsx("Data/PFOA_Urine_female_oral_5_mg_per_kg.xlsx")
Kemp_OR_Furine_high <- openxlsx::read.xlsx("Data/PFOA_Urine_female_oral_25_mg_per_kg.xlsx")
Kemp_OR_Murine_low <- openxlsx::read.xlsx("Data/PFOA_Urine_male_oral_1_mg_per_kg.xlsx")
Kemp_OR_Murine_med <- openxlsx::read.xlsx("Data/PFOA_Urine_male_oral_5_mg_per_kg.xlsx")
Kemp_OR_Murine_high <- openxlsx::read.xlsx("Data/PFOA_Urine_male_oral_25_mg_per_kg.xlsx")
dzi_IV_Mserum <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_IV_2021.xlsx")
dzi_IV_Mserum$Concentration_microM <- dzi_IV_Mserum$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_OR_Mserum_low <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_ORAL_low_2021.xlsx")
dzi_OR_Mserum_low$Concentration_microM <- dzi_OR_Mserum_low$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Mserum_medium <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_ORAL_medium_2021.xlsx")
dzi_OR_Mserum_medium$Concentration_microM <- dzi_OR_Mserum_medium$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_OR_Mserum_high <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_ORAL_high_2021.xlsx")
dzi_OR_Mserum_high$Concentration_microM <- dzi_OR_Mserum_high$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_IV_Fserum <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_IV_2021.xlsx")
dzi_IV_Fserum$Concentration_microM <- dzi_IV_Fserum$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_OR_Fserum_low <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_ORAL_low_2021.xlsx")
dzi_OR_Fserum_low$Concentration_microM <- dzi_OR_Fserum_low$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Fserum_medium <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_ORAL_medium_2021.xlsx")
dzi_OR_Fserum_medium$Concentration_microM <- dzi_OR_Fserum_medium$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Fserum_high <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_ORAL_high_2021.xlsx")
dzi_OR_Fserum_high$Concentration_microM <- dzi_OR_Fserum_high$Concentration_microM * MW/1000 #convert from uM to ug/g
kim_OR_Fblood <- openxlsx::read.xlsx("Data/PFOA_female_blood_ORAL_kim_2016.xlsx")
kim_IV_Fblood <- openxlsx::read.xlsx("Data/PFOA_female_blood_IV_kim_2016.xlsx")
gus_OR_Mblood <- openxlsx::read.xlsx("Data/Gustafsson 2022_PFOA_Plasma Male rats_Oral.xlsx")
gus_OR_Mtissues <- openxlsx::read.xlsx("Data/Gustafsson 2022_PFOA_Tissues Male rats_Oral.xlsx")
Kemp_OR_Ffeces_med <- openxlsx::read.xlsx("Data/PFOA_Feces_female_oral_5_mg_per_kg-Loc.xlsx")
Kemp_OR_Mfeces_med <- openxlsx::read.xlsx("Data/PFOA_Feces_male_oral_5_mg_per_kg-Loc.xlsx")
Kemp_OR_Ffeces_low <- openxlsx::read.xlsx("Data/PFOA_Feces_female_oral_1_mg_per_kg-Loc.xlsx")
Kemp_OR_Mfeces_low <- openxlsx::read.xlsx("Data/PFOA_Feces_male_oral_1_mg_per_kg-Loc.xlsx")


setwd("/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Publication/Ka_albumin_estimated/favail_estimated/reduction_estimated")

dataset <- list("df1" = kudo_high_dose, "df2" = kudo_low_dose, "df3" = kim_IV_Mtissues, "df4" = kim_OR_Mtissues,
                "df5" = kim_IV_Ftissues, "df6" = kim_OR_Ftissues, "df7" = dzi_OR_Mtissues, "df8" = dzi_OR_Ftissues,
                "df9" = kim_OR_Mblood, "df10" = kim_IV_Mblood, 
                "df11a" = Kemp_OR_Furine_low, "df11b" = Kemp_OR_Ffeces_low,
                "df12a" = Kemp_OR_Furine_med, "df12b" = Kemp_OR_Ffeces_med,
                "df13a" = Kemp_OR_Furine_high,"df13b" = Kemp_OR_Ffeces_high,
                "df14a" = Kemp_OR_Murine_low,"df14b" = Kemp_OR_Mfeces_low,
                "df15a" = Kemp_OR_Murine_med,"df15b" = Kemp_OR_Mfeces_med,
                "df16a" = Kemp_OR_Murine_high,"df16b" = Kemp_OR_Mfeces_high,
                "df17" = dzi_IV_Mserum, "df18" = dzi_OR_Mserum_low, "df19" = dzi_OR_Mserum_medium,
                "df20" = dzi_OR_Mserum_high, "df21" = dzi_IV_Fserum, "df22" = dzi_OR_Fserum_low,
                "df23" = dzi_OR_Fserum_medium, "df24" = dzi_OR_Fserum_high, "df25" = kim_OR_Fblood, 
                "df26" = kim_IV_Fblood, "df27" = gus_OR_Mblood, "df28" = gus_OR_Mtissues)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_NEWUOA", #"NLOPT_LN_NEWUOA"   "NLOPT_LN_SBPLX"
              "xtol_rel" = 1e-06,
              "ftol_rel" = 1e-06,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0, 
              "maxeval" = 300, 
              "print_level" = 1)



N_pars <- 12# Number of parameters to be fitted
fit <- log(c(1.8e+01, 1e-03, 1e3, 0.02, 1, 1, 1e-4,  0.1, 1e5, 1, 1, 1))


fixed_params <- create_all_fixed_params()
# Run the optimization algorithm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             # lb	= lb,
                             # ub = ub,
                             opts = opts,
                             dataset = dataset,
                             fixed_params = fixed_params)

#estimated_params <- exp(optimizer$solution)
estimated_params <- exp(optimizer$solution)
save.image("ka_estimated_favail_estimated2_reduction_estimated.RData")


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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kudo_high <-  solution[, c("time","Cblood","Cliver","Ckidney", "Ccarcass","Clungs", 
                                 "Cspleen", "Cheart","Cbrain", "Cgonads", "Cstomach", 
                                 "Cintestine")]
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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kudo_low <- solution[, c("time","Cblood","Cliver","Ckidney", "Ccarcass","Clungs", 
                               "Cspleen", "Cheart","Cbrain", "Cgonads", "Cstomach", 
                               "Cintestine")]
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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kim_IV_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kim_OR_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kim_IV_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kim_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

###############################################################################################
# Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[7]])
params <- c(fixed_params[[7]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,864,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]

#############################################################################################
# Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
BW <- 0.2  # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[8]])
params <- c(fixed_params[[8]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,0.1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]

######################################################################################
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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kim_IV_Mblood <-  solution[, c("time", "Cplasma")]

##########################
#-------------------------
# Kemper 2003 (Worley)
#-------------------------
##########################
# Set up simulations for the 11th case, i.e.Kemper 2003 (Worley) ORAL female LOW

sex <- "F"
BW <- 0.2 #kg
sample_time <- seq(0,192,1)
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[11]])
params <- c(fixed_params[[11]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_Kemp_OR_Furine_low <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Ffeces_low <-  solution[, c("time", "Mfeces")]

##########################################################################################
# Set up simulations for the 12th case, i.e.Kemper 2003 (Worley) ORAL female  MEDIUM
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[12]])
params <- c(fixed_params[[12]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_Kemp_OR_Furine_med <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Ffeces_med <-  solution[, c("time", "Mfeces")]

##########################################################################################
# Set up simulations for the 13th case, i.e.Kemper 2003 (Worley) ORAL female HIGH
sample_time <- seq(0,672,1)
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[13]])
params <- c(fixed_params[[13]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_Kemp_OR_Furine_high <-  solution[1:169, c("time", "Murine")]
preds_Kemp_OR_Ffeces_high <-  solution[, c("time", "Mfeces")]

###################################################################################################
# Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL male LOW

sex <- "M"
BW <- 0.3 #kg
sample_time <- seq(0,673,1)
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[14]])
params <- c(fixed_params[[14]], variable_params)
inits <- create.inits (params)
events <- create.events(params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_Kemp_OR_Murine_low <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_low <-  solution[, c("time", "Mfeces")]

####################################################################################################
# Set up simulations for the 15th case, i.e.Kemper 2003 (Worley) ORAL male MEDIUM
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[15]])
params <- c(fixed_params[[15]], variable_params)
inits <- create.inits (params)
events <- create.events(params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_Kemp_OR_Murine_med <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_med <-  solution[, c("time", "Mfeces")]

######################################################################################################
# Set up simulations for the 16th case, i.e.Kemper 2003 (Worley) ORAL male  HIGH
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[16]])
params <- c(fixed_params[[16]], variable_params)
inits <- create.inits (params)
events <- create.events(params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_Kemp_OR_Murine_high <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_high <-  solution[, c("time", "Mfeces")]

#####################################################################################################
# Set up simulations for the 17th case, i.e. Dzierlenga 2021, IV male serum
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[17]])
params <- c(fixed_params[[17]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.083, 0.25, 0.5, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_IV_Mserum <-  solution[, c("time", "Cplasma")]

##############################################################################################
# Set up simulations for the 18th case, i.e. Dzierlenga 2021, ORAL male serum low
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[18]])
params <- c(fixed_params[[18]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Mserum_low <-  solution[, c("time", "Cplasma")]

########################################################################################################
# Set up simulations for the 19th case, i.e. Dzierlenga 2021, ORAL male serum medium
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[19]])
params <- c(fixed_params[[19]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Mserum_medium <-  solution[, c("time", "Cplasma")]

##############################################################################################
# Set up simulations for the 20th case, i.e. Dzierlenga 2021, ORAL male serum high
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[20]])
params <- c(fixed_params[[20]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Mserum_high <-  solution[, c("time", "Cplasma")]

###############################################################################################
# Set up simulations for the 21th case, i.e. Dzierlenga 2021, IV female serum
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[21]])
params <- c(fixed_params[[21]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.083, 0.25, seq(0.5, 192, 0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_IV_Fserum <-  solution[, c("time", "Cplasma")]

#########################################################################################
# Set up simulations for the 22th case, i.e. Dzierlenga 2021, ORAL female serum low
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[22]])
params <- c(fixed_params[[22]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0, 96, 0.25)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]

################################################################################################
# Set up simulations for the 23th case, i.e. Dzierlenga 2021, ORAL female serum medium
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[23]])
params <- c(fixed_params[[23]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, seq(1, 192, 0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Fserum_medium <-  solution[, c("time", "Cplasma")]

########################################################################################
# Set up simulations for the 24th case, i.e. Dzierlenga 2021, ORAL female serum high
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[24]])
params <- c(fixed_params[[24]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time <- seq(0, 96, 0.25)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_dzi_OR_Fserum_high <-  solution[, c("time", "Cplasma")]

##########################################################################################
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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

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
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]

###############################################################################################
# Set up simulations for the 27th case, i.e. Gustafsson (2022) oral male blood
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[27]])
params <- c(fixed_params[[27]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0,48,0.2)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_gus_OR_Mblood <-  solution[, c("time", "Cplasma")]

##########################################################################################
# Set up simulations for the 28st case, i.e. Gustafsson (2022) Inhalation male tissues
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params, fixed_params[[28]])
params <- c(fixed_params[[28]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0,48,0.2)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_gus_OR_Mtissues <-  solution[, c("time", "CalveolarLF","Cliver", "Clungtissue", "Ckidney")]
##########################################################################################################
#########################################################################################################

#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets
preds_kudo_high[,2:dim(preds_kudo_high)[2]] <- preds_kudo_high[,2:dim(preds_kudo_high)[2]] /1000 
preds_kudo_low[,2:dim(preds_kudo_low)[2]] <- preds_kudo_low[,2:dim(preds_kudo_low)[2]] /1000 
preds_kim_IV_Mtissues[,2:dim(preds_kim_IV_Mtissues)[2]] <- preds_kim_IV_Mtissues[,2:dim(preds_kim_IV_Mtissues)[2]] /1000
preds_kim_OR_Mtissues[,2:dim(preds_kim_OR_Mtissues)[2]] <- preds_kim_OR_Mtissues[,2:dim(preds_kim_OR_Mtissues)[2]] /1000
preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] <- preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] /1000
preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] <- preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] /1000
preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] <- preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000
preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] <- preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] /1000
preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] <- preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] /1000
preds_Kemp_OR_Ffeces_high[,2:dim(preds_Kemp_OR_Ffeces_high)[2]] <- preds_Kemp_OR_Ffeces_high[,2:dim(preds_Kemp_OR_Ffeces_high)[2]] /1000
preds_Kemp_OR_Mfeces_high[,2:dim(preds_Kemp_OR_Mfeces_high)[2]] <- preds_Kemp_OR_Mfeces_high[,2:dim(preds_Kemp_OR_Mfeces_high)[2]] /1000
preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] <- preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] /1000
preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] <- preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] /1000
preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] <- preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] /1000
preds_Kemp_OR_Murine_low[,2:dim(preds_Kemp_OR_Murine_low)[2]] <- preds_Kemp_OR_Murine_low[,2:dim(preds_Kemp_OR_Murine_low)[2]] /1000
preds_Kemp_OR_Murine_med[,2:dim(preds_Kemp_OR_Murine_med)[2]] <- preds_Kemp_OR_Murine_med[,2:dim(preds_Kemp_OR_Murine_med)[2]] /1000
preds_Kemp_OR_Murine_high[,2:dim(preds_Kemp_OR_Murine_high)[2]] <- preds_Kemp_OR_Murine_high[,2:dim(preds_Kemp_OR_Murine_high)[2]] /1000
preds_dzi_IV_Mserum[,2:dim(preds_dzi_IV_Mserum)[2]] <- preds_dzi_IV_Mserum[,2:dim(preds_dzi_IV_Mserum)[2]] /1000
preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] <- preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] /1000
preds_dzi_OR_Mserum_medium[,2:dim(preds_dzi_OR_Mserum_medium)[2]] <- preds_dzi_OR_Mserum_medium[,2:dim(preds_dzi_OR_Mserum_medium)[2]] /1000
preds_dzi_OR_Mserum_high[,2:dim(preds_dzi_OR_Mserum_high)[2]] <- preds_dzi_OR_Mserum_high[,2:dim(preds_dzi_OR_Mserum_high)[2]] /1000
preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] <- preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] /1000
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] <- preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] /1000
preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] <- preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] /1000
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000
preds_gus_OR_Mblood[,2:dim(preds_gus_OR_Mblood)[2]] <- preds_gus_OR_Mblood[,2:dim(preds_gus_OR_Mblood)[2]] /1000
preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] <- preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] /1000
preds_Kemp_OR_Ffeces_med[,2:dim(preds_Kemp_OR_Ffeces_med)[2]] <- preds_Kemp_OR_Ffeces_med[,2:dim(preds_Kemp_OR_Ffeces_med)[2]] /1000
preds_Kemp_OR_Mfeces_med[,2:dim(preds_Kemp_OR_Mfeces_med)[2]] <- preds_Kemp_OR_Mfeces_med[,2:dim(preds_Kemp_OR_Mfeces_med)[2]] /1000
preds_Kemp_OR_Ffeces_low[,2:dim(preds_Kemp_OR_Ffeces_low)[2]] <- preds_Kemp_OR_Ffeces_low[,2:dim(preds_Kemp_OR_Ffeces_low)[2]] /1000
preds_Kemp_OR_Mfeces_low[,2:dim(preds_Kemp_OR_Mfeces_low)[2]] <- preds_Kemp_OR_Mfeces_low[,2:dim(preds_Kemp_OR_Mfeces_low)[2]] /1000


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
experiment1 <- reshape(kudo_high_dose[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",kudo_high_dose$Tissue )

# Convert Kudo Low dose from long to wide format using reshape
experiment2 <- reshape(kudo_low_dose[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",kudo_low_dose$Tissue )

# Convert Kim IV Male tissues from long to wide format using reshape
experiment3 <- reshape(kim_IV_Mtissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment3) <- c("Time",kim_IV_Mtissues$Tissue )

# Convert Kim ORAL Male tissues from long to wide format using reshape
experiment4 <- reshape(kim_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment4) <- c("Time",kim_OR_Mtissues$Tissue )

# Convert Kim IV female tissues from long to wide format using reshape
experiment5 <- reshape(kim_IV_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",kim_IV_Ftissues$Tissue )

# Convert Kim ORAL female tissues from long to wide format using reshape
experiment6 <- reshape(kim_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",kim_OR_Ftissues$Tissue )

# Convert Dzierlenga ORAL male tissues from long to wide format using reshape
experiment7 <- reshape(dzi_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(dzi_OR_Mtissues$Tissue))


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))

# Convert Kim ORAL male blood from long to wide format using reshape
experiment9 <- reshape(kim_OR_Mblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9) <- c("Time",unique(kim_OR_Mblood$Tissue))


# Convert Kim IV male blood from long to wide format using reshape
experiment10 <- reshape(kim_IV_Mblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment10) <- c("Time",unique(kim_IV_Mblood$Tissue))


# Convert Kemper ORAL female urine low from long to wide format using reshape
experiment11a <- reshape(Kemp_OR_Furine_low [c("Tissue" ,"Time_h", 
                                               "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment11a) <- c("Time",unique(Kemp_OR_Furine_low$Tissue))
experiment11a$Urine = (experiment11a$Urine/100)*0.2*1


# Convert Kemper ORAL female feces low from long to wide format using reshape
experiment11b <- reshape(Kemp_OR_Ffeces_low[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment11b) <- c("Time",unique(Kemp_OR_Ffeces_low$Tissue))
experiment11b$Feces = (experiment11b$Feces/100)*0.2*1

# Convert Kemper ORAL female urine med from long to wide format using reshape
experiment12a <- reshape(Kemp_OR_Furine_med[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment12a) <- c("Time",unique(Kemp_OR_Furine_med$Tissue))
experiment12a$Urine = (experiment12a$Urine/100)*0.2*5


# Convert Kemper ORAL female feces medium from long to wide format using reshape
experiment12b <- reshape(Kemp_OR_Ffeces_med[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment12b) <- c("Time",unique(Kemp_OR_Ffeces_med$Tissue))
experiment12b$Feces = (experiment12b$Feces/100)*0.2*5

# Convert Kemper ORAL female urine high from long to wide format using reshape
experiment13a <- reshape(Kemp_OR_Furine_high[c("Tissue" ,"Time_h", 
                                               "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment13a) <- c("Time",unique(Kemp_OR_Furine_high$Tissue))
experiment13a$Urine = (experiment13a$Urine/100)*0.2*25

# Convert Kemper ORAL female feces from long to wide format using reshape
experiment13b <- reshape(Kemp_OR_Ffeces_high[c("Tissue" ,"Time_h", 
                                               "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment13b) <- c("Time",unique(Kemp_OR_Ffeces_high$Tissue))
experiment13b$Feces = (experiment13b$Feces/100)*0.2*25


# Convert Kemper ORAL male urine low from long to wide format using reshape
experiment14a <- reshape(Kemp_OR_Murine_low [c("Tissue" ,"Time_h", 
                                               "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment14a) <- c("Time",unique(Kemp_OR_Murine_low$Tissue))
experiment14a$Urine = (experiment14a$Urine/100)*0.3*1


# Convert Kemper ORAL male feces low from long to wide format using reshape
experiment14b <- reshape(Kemp_OR_Mfeces_low[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment14b) <- c("Time",unique(Kemp_OR_Mfeces_low$Tissue))
experiment14b$Feces = (experiment14b$Feces/100)*0.3*1

# Convert Kemper ORAL male urine med from long to wide format using reshape
experiment15a <- reshape(Kemp_OR_Murine_med[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment15a) <- c("Time",unique(Kemp_OR_Murine_med$Tissue))
experiment15a$Urine = (experiment15a$Urine/100)*0.3*5


# Convert Kemper ORAL male feces medium from long to wide format using reshape
experiment15b <- reshape(Kemp_OR_Mfeces_med[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment15b) <- c("Time",unique(Kemp_OR_Mfeces_med$Tissue))
experiment15b$Feces = (experiment15b$Feces/100)*0.3*5


# Convert Kemper ORAL male urine high from long to wide format using reshape
experiment16a <- reshape(Kemp_OR_Murine_high[c("Tissue" ,"Time_h", 
                                               "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment16a) <- c("Time",unique(Kemp_OR_Furine_high$Tissue))
experiment16a$Urine = (experiment16a$Urine/100)*0.3*25


# Convert Kemper ORAL male feces from long to wide format using reshape
experiment16b <- reshape(Kemp_OR_Mfeces_high[c("Tissue" ,"Time_h", 
                                               "Cum_dose_%")], 
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment16b) <- c("Time",unique(Kemp_OR_Mfeces_high$Tissue))
experiment16b$Feces = (experiment16b$Feces/100)*0.3*25

# Convert Dzierlenga 2021, IV male serum from long to wide format using reshape
experiment17 <- reshape(dzi_IV_Mserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment17) <- c("Time",unique(dzi_IV_Mserum$Tissue))

# Convert Dzierlenga 2021, ORAL male serum low from long to wide format using reshape
experiment18 <- reshape(dzi_OR_Mserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment18) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))

# Convert Dzierlenga 2021, ORAL male serum medium from long to wide format using reshape
experiment19 <- reshape(dzi_OR_Mserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment19) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL male serum high from long to wide format using reshape
experiment20 <- reshape(dzi_OR_Mserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment20) <- c("Time",unique(dzi_OR_Mserum_high$Tissue))

#Convert Dzierlenga 2021, IV female serum from long to wide format using reshape
experiment21 <- reshape(dzi_IV_Fserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment21) <- c("Time",unique(dzi_IV_Fserum$Tissue))

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment22 <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment22) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL female serum medium from long to wide format using reshape
experiment23 <- reshape(dzi_OR_Fserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment23) <- c("Time",unique(dzi_OR_Fserum_medium$Tissue))

#Convert Dzierlenga 2021, ORAL female serum high from long to wide format using reshape
experiment24 <- reshape(dzi_OR_Fserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment24) <- c("Time",unique(dzi_OR_Fserum_high$Tissue))

#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment25 <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment25) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment26<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment26) <- c("Time",unique(kim_IV_Fblood$Tissue))


# Convert Gustafsson Oral male blood from long to wide format using reshape
experiment27 <- reshape(gus_OR_Mblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment27) <- c("Time",unique(gus_OR_Mblood$Tissue))


# Convert Gustafsson Oral male tissues from long to wide format using reshape
experiment28 <- reshape(gus_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                          "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment28) <- c("Time",unique(gus_OR_Mtissues$Tissue))



# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment3, 
                    experiment4 = experiment4,experiment5 = experiment5, experiment6 = experiment6,
                    experiment7 = experiment7, experiment8 = experiment8,experiment9 = experiment9, 
                    experiment10 = experiment10, experiment11a = experiment11a, experiment11b = experiment11b, 
                    experiment12a = experiment12a, experiment12b = experiment12b,
                    experiment13a = experiment13a, experiment13b = experiment13b,
                    experiment14a = experiment14a, experiment14b = experiment14b,
                    experiment15a = experiment15a,experiment15b = experiment15b,
                    experiment16a = experiment16a,experiment16b = experiment16b,
                    experiment17 = experiment17, experiment18 = experiment18, 
                    experiment19 = experiment19, experiment20 = experiment20,experiment21 = experiment21, 
                    experiment22 = experiment22, experiment23 = experiment23, experiment24 = experiment24,
                    experiment25 = experiment25, experiment26 = experiment26, experiment27 = experiment27,
                    experiment28 = experiment28)


# Rename predictions so that they share the same name as the names of the experimental data dataframe
colnames(preds_kudo_high) <- c( "Time", "Blood", "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                                "Heart", "Brain", "Gonads", "Stomach", "Intestine")

colnames(preds_kudo_low) <-  colnames(preds_kudo_high) 
colnames(preds_kim_IV_Mtissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")

colnames(preds_kim_OR_Mtissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_IV_Ftissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_OR_Ftissues) <- colnames(preds_kim_IV_Mtissues)

colnames(preds_dzi_OR_Mtissues) <- c("Time","Liver","Kidney","Brain")
colnames(preds_dzi_OR_Ftissues) <- colnames(preds_dzi_OR_Mtissues)

colnames(preds_kim_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_kim_IV_Mblood) <- c ("Time", "Plasma")

colnames(preds_Kemp_OR_Ffeces_high) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_high) <- c ("Time", "Feces")

colnames(preds_Kemp_OR_Furine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_high) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_high) <- c ("Time", "Urine")

colnames(preds_dzi_IV_Mserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_high) <- c ("Time", "Serum")
colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")


colnames(preds_gus_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_gus_OR_Mtissues) <- c ("Time", "ALF", "Liver", "Lung", "Kidney")
colnames(preds_Kemp_OR_Ffeces_med) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_med) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Ffeces_low) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_low) <- c ("Time", "Feces")

# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_kudo_high,  predictions2 = preds_kudo_low, 
                    predictions3 = preds_kim_IV_Mtissues, 
                    predictions4 = preds_kim_OR_Mtissues, predictions5 = preds_kim_IV_Ftissues,
                    predictions6 = preds_kim_OR_Ftissues,
                    predictions7 = preds_dzi_OR_Mtissues, predictions8 = preds_dzi_OR_Ftissues, 
                    predictions9 = preds_kim_OR_Mblood,
                    predictions10 = preds_kim_IV_Mblood, 
                    predictions11a = preds_Kemp_OR_Furine_low, predictions11b = preds_Kemp_OR_Ffeces_low,
                    predictions12a = preds_Kemp_OR_Furine_med, predictions12b = preds_Kemp_OR_Ffeces_med,
                    predictions13a =preds_Kemp_OR_Furine_high, predictions13b = preds_Kemp_OR_Ffeces_high,
                    predictions14a = preds_Kemp_OR_Murine_low, predictions14b = preds_Kemp_OR_Mfeces_low, 
                    predictions15a = preds_Kemp_OR_Murine_med, predictions15b = preds_Kemp_OR_Mfeces_med,
                    predictions16a =preds_Kemp_OR_Murine_high, predictions16b = preds_Kemp_OR_Mfeces_high,
                    prediction17 =preds_dzi_IV_Mserum, predictions18 =preds_dzi_OR_Mserum_low,
                    predictions19 =preds_dzi_OR_Mserum_medium, predictions20 =preds_dzi_OR_Mserum_high, 
                    predictions21 =preds_dzi_IV_Fserum, predictions22 =preds_dzi_OR_Fserum_low,
                    predictions23 =preds_dzi_OR_Fserum_medium,
                    predictions24 =preds_dzi_OR_Fserum_high, predictions25 = preds_kim_OR_Fblood, 
                    predictions26 = preds_kim_IV_Fblood, predictions27 = preds_gus_OR_Mblood, 
                    predictions28 = preds_gus_OR_Mtissues)


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



