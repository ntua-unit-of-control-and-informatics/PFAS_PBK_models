# Check notes for description of scenario
#female increase in BW: 3.5 g/d
#Male increase in BW: 5.9 g.d

library(deSolve)

create_variable_params <- function(BW,sex,  estimated_params, fixed_params){
  
  # BW in kg
  # Cheng and Ng 2017 Table S1
  # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
  #assuming the density of tissue is 1 g/mL.
  # Estimated parameters
  RAF_api <- 1
  
  if(sex == "M"){
    RAFOatp_k <-2.088006e+00
  }else{
    RAFOatp_k <- 1E-7
  }
  RAFOat3 <- 100
  CL_int <- 1e+01 #uL/min/million hepatocytes
  RAFOatp_l <- 100
  RAFUrat <- RAFOatp_k
  RAFOat1 <- 0
  RAFOatp2_l <- RAFOatp_l
  RAFOatp_lu_ap <- 8.659643e+00 
  RAFOatp_lu_bas <- RAFOatp_lu_ap
  RAFNtcp <- RAFOatp_l
  RAFOatp2_Int <- 0
  
  Papp_gut <- 0.06
  Papp <- 0.05
  
  f_fabp_avail <- 1
  f_alb_avail <-1
  
  koff_alb <-  1e-0
  koff_fabp <-  koff_alb
  koff_a2u <- koff_alb
  
  
  VmK_baso <- 0
  KmK_baso <- 1e20
  KLfabp <- (1.2e5+4e4+1.9e4)  #[L/mol]*1e-3 , value from Cheng et al. (2017)
  Ka <-5e5 # 5.8e05 from Rue et al. (2024)#mol/L
  CLfeces_unscaled <- 4.091396e-04 #in L/h/BW^(-0.25), scaling similar to Loccisano et al. (2012)
  CLfeces <- CLfeces_unscaled*BW^(-0.25)  #in L/h
  
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
  
  KmK_api <-   KmK_Oatp
  VmK_api <- VmK_Oatp_scaled*RAF_api
  
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
  
  #Albumin concentration in interstitial fluid compartments(mol/m^3 = 1e-6* nmol/g)
  # CalbKF_init <- f_alb_avail*243*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  # CalbLF_init <- f_alb_avail*243*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  # CalbSTF_init <- f_alb_avail*146*1e-6 # [umol/L]*1e-6 -->(mol/L), same as Gut (assumption)
  # CalbINF_init <- f_alb_avail*146*1e-6 #[umol/L]*1e-6 -->(mol/L), same as Gut (assumption)
  # CalbMF_init <- f_alb_avail*146*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  # CalbAF_init <- f_alb_avail*.*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  # CalbRF_init <- f_alb_avail*73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  # CalbBoF_init <-f_alb_avail* 73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  # CalbLuF_init <- f_alb_avail*CalbINF_init #assumption 
  # 
  # CalbSPF_init <- f_alb_avail*243e-6 #[umol/L]*1e-6 -->(mol/L), same as liver (assumption)
  # CalbGoF_init <- f_alb_avail*41/Mr_albumin #mg/mL-->  (mol/L) from https://doi.org/10.1210/endo-116-5-1983 --> 41 mg/mL, MW=65 kg/mol
  # CalbHF_init <- f_alb_avail* 65/Mr_albumin##mg/mL--> (mol/L) https://doi.org/10.1007/s12291-010-0042-x --> 6.5 g/100 g tissue, MW=65 kg/mol 
  # CalbBrF_init <- f_alb_avail*8e-2/Mr_albumin  ##mg/mL--> (mol/L) https://doi.org/10.1016/0014-4886(90)90158-O --> 0.08 g/L, MW=65 kg/mol 
  # CalbSKF_init <- f_alb_avail*21/Mr_albumin ##mg/mL-->  (mol/L) https://doi.org/10.1111/j.1748-1716.1973.tb05464.x -->Table 2: 2.1 g/100 mL
  # 
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
  CalbLuAF_init <- 10/100 * CalbB_init #based on Woods et al. 2015 statement https://doi.org/10.1016/j.jconrel.2015.05.269
  
  
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
  Papp_LIN = RAF_papp*(ClINFT_unscaled*60*1e-06*1e3*protein_per_well)/(Awell *2) #cm/h,at  pH = 6.0
  Papp_RYU = RAF_papp*1.46e-6*3600 # cm/h, at pH = 7.4 from Ryu et al. (2024) [https://doi.org/10.1016/j.chemosphere.2024.142390]
  #  Lin et al. (2023) used data from kimura et al. (2017) [http://dx.doi.org/10.1016/j.toxlet.2017.05.012]
  # which is more appropriate for 
  # For endothelial and cellular permeability we use the Ryu et al. (2024) value
  #Papp = Papp_RYU
  k_gut_in = ( (Papp_gut/100) * fixed_params$AINL)*1000 #L/h
  k_gut_out = ( (Papp/100) * fixed_params$AINL)*1000 #L/h
  
  #passive diffusion rates, in L/h
  kLFLT = ((Papp/100) * fixed_params$AcL)*1000 #m^3/h * 1000 --> L/h
  #kLTLbile = ((Papp/100) * AcLBilec)*1000 #m^3/h * 1000 --> L/h
  kMFMT = ((Papp/100) * fixed_params$AcM)*1000 #m^3/h * 1000 --> L/h
  kSTFSTT = ((Papp/100) * fixed_params$AcST)*1000 #m^3/h * 1000 --> L/h 
  kINFINT = ((Papp/100) * fixed_params$AcIN)*1000 #m^3/h * 1000 --> L/h 
  kAFAT = ((Papp/100) * fixed_params$AcA)*1000 #m^3/h * 1000 --> L/h 
  kLuTLuF = ((Papp/100) * fixed_params$AcLu)*1000 #m^3/h * 1000 --> L/h
  kLuTLuAF = ((Papp/100) * fixed_params$AcALF)*1000 #m^3/h * 1000 --> L/h
  kSPFSPT = ((Papp/100) * fixed_params$AcSP)*1000 #m^3/h * 1000 --> L/h 
  kHFHT = ((Papp/100) * fixed_params$AcH)*1000 #m^3/h * 1000 --> L/h 
  kBrFBrT = ((Papp/100) * fixed_params$AcBr)*1000 #m^3/h * 1000 --> L/h 
  kGoFGoT = ((Papp/100) * fixed_params$AcGo)*1000 #m^3/h * 1000 --> L/h 
  kSKFSKT = ((Papp/100) * fixed_params$AcSK)*1000 #m^3/h * 1000 --> L/h
  kBoFBoT = ((Papp/100) * fixed_params$AcBo)*1000 #m^3/h * 1000 --> L/h
  kRFRT = ((Papp/100) * fixed_params$AcR)*1000 #m^3/h*1000 --> L/h 
  
  #Diffusion rates in L/h between renal tubule filtrate and tubule cells
  kPtcTu <- ((Papp/100) * fixed_params$APT) *1000 #diffusion between proximal tubule cells and tubule filtrate
  kDalcTu <- ((Papp/100) * fixed_params$ADAL) *1000 #diffusion between descending/ascending cells and tubule filtrate
  kDtcTu <- ((Papp/100) * fixed_params$ADT) *1000 #diffusion between distal tubule cells and tubule filtrate
  kCdcTu <- ((Papp/100) * fixed_params$ACD) *1000 #diffusion between collecting duct cells and tubule filtrate
  
  #Diffusion rates in L/h between  tubule cells and interstitial space
  kPtcF <- ((Papp/100) * fixed_params$AcK_PTC) *1000 #diffusion between proximal tubule cells and interstitial space
  kDtcF <- ((Papp/100) * fixed_params$AcK_DTC) *1000 #diffusion between descending/ascending cells and interstitial space
  
  kDalcF <- ((Papp/100) * fixed_params$AcK_DALC) *1000 #diffusion between proximal tubule cells and interstitial space
  kCdcF <- ((Papp/100) * fixed_params$AcK_CDC) *1000 #diffusion between descending/ascending cells and interstitial space
  kKTrestF  <- ((Papp/100) * fixed_params$AcKTrest) *1000 #diffusion between rest of kidney cells and interstitial space
  
  #Effective permeability (Peff, in mm/h) for blood (B), liver(L), kidney(K),
  #stomach(ST),intestine (IN), adipose(A), muscle(M), spleen (SP), heart (H), 
  #brain (Br), gonads (Go), rest of body(R). Similar to Lin et al. (2023) [https://doi.org/10.1021/acs.est.2c05642.]
  # we devide the effective permeability by 2 to account for the resistance at both
  # the apical and basolateral side of the endothelium
  
  PeffK <- Papp*10/2 #mm/h
  PeffL <- Papp*10/2 #mm/h
  PeffST <- Papp*10/2 #mm/h
  PeffIN <- Papp*10/2 #mm/h
  PeffA <- Papp*10/2 #mm/h
  PeffM <- Papp*10/2 #mm/h
  PeffR <- Papp*10/2 #mm/h
  PeffLu <- Papp*10/2 #mm/h
  PeffSP <- Papp*10/2 #mm/h
  PeffH <- Papp*10/2 #mm/h
  PeffBr <- Papp*10/2 #mm/h
  PeffGo <- Papp*10/2 #mm/h
  PeffSK <- Papp*10/2 #mm/h
  PeffBo <- Papp*10/2 #mm/h
  
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
    
    'Papp' = Papp, 'k_gut_in' = k_gut_in, 'k_gut_out' = k_gut_out,
    'kKTrestF'=kKTrestF, 'kCdcF' = kCdcF, 'kDalcF' = kDalcF, 'kPtcF' =  kPtcF,   'kDtcF' = kDtcF,
    'kPtcTu'=kPtcTu, 'kDalcTu' = kDalcTu, 'kDtcTu' = kDtcTu, 'kCdcTu' = kCdcTu, 
    'kLFLT'=kLFLT,  'kAFAT'=kAFAT, 
    'kRFRT'=kRFRT,
    'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF, 'kLuTLuAF'=kLuTLuAF, 'kSPFSPT' =kSPFSPT,
    'kSTFSTT' =kSTFSTT, 'kINFINT' =kINFINT, 'kHFHT' =kHFHT,
    'kBrFBrT' =kBrFBrT, 'kGoFGoT' =kGoFGoT,
    'kSKFSKT' =kSKFSKT, 'kBoFBoT'=kBoFBoT,
    
    'PeffK'=PeffK, 'PeffL'=PeffL, 
    'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'PeffLu' = PeffLu,
    'PeffSP'=PeffSP, 'PeffH'=PeffH, 'PeffBr'=PeffBr, 'PeffST' = PeffST,
    'PeffIN'=PeffIN, 'PeffGo'=PeffGo,
    'PeffSK' = PeffSK,  'PeffBo' = PeffBo
    
    
  ))
}  


################################################################################


setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")

MW <- 414.07 #g/mol
source("Goodness-of-fit-metrics.R")

# Read data
dzi_OR_Ftissues <- openxlsx::read.xlsx("Data/Dzierlenga_tissue_female_ORAL_2021.xlsx")

dzi_OR_Ftissues$Concentration_microM <- dzi_OR_Ftissues$Concentration_microM* MW/1000 #convert from uM to ug/g

setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/proximal_tubule/scenario22")

dataset <- list( "df8" = dzi_OR_Ftissues)


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
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000


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

# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))



# Put the experiments in a list
experiments <- list(experiment8 = experiment8)

colnames(preds_dzi_OR_Ftissues) <- c("Time","Liver","Kidney","Brain")

# Create a list containing the corresponding predictions
simulations <- list(predictions8 = preds_dzi_OR_Ftissues)


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
  
  print(final_plot)
}



