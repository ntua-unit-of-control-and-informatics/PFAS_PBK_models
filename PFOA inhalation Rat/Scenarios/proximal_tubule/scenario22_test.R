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
    PVINB <- ((0.0795+0.0651)+(0.0458+0.0375))/280 #mL/g BW weight
    VINB <- PVINB * BW #intestine  blood volume kg=L
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
    PVMB <- 0.04 #Brown et al. 1997
    VMB <- PVMB * PVM * BW #muscle blood volume kg=L
    PVMF <- 0.12 #pkSim
    VMF <- PVMF * PVM * BW #muscle interstitial fluid volume kg=L
    VMT <- VM - VMF #muscle tissue volume kg=L
    
    #Adipose
    PVA <- 7e-2 #Brown et al. 1997
    VA <- PVA * BW #adipose volume kg=L
    PVAB <- 0.02 #Brown et al. 1997
    VAB <- PVAB * PVA * BW #% adipose blood volume kg=L
    PVAF <- 0.174 #Ng, 2013
    VAF <- PVAF * PVA * BW #adipose interstitial fluid volume kg=L
    VAT <- VA - VAF #adipose tissue volume kg=L
    
    #Lung
    PVLu <- 0.48e-2 #Brown et al. 1997, p 418-19 mean of values for male rat
    VLu <- PVLu * BW
    PVLuB <- 9/100*PVLu
    VLuB <- PVLuB*BW #Brown et al. 1997, p 459 --> capillary blood occupied 9% of the lung volume
    PVLuF <- 0.263/280 #0.263 ml, Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VLuF <- PVLuF * BW #lung interstitial fluid volume
    PVLuAF <- 0.4/275 #0.4 mL Leslie et al, 1989 https://doi.org/10.1164/ajrccm/139.2.360 --> Watkins & Rannels 1979 https://doi.org/10.1152/jappl.1979.47.2.325  
    VLuAF <- PVLuAF * BW #lung alveolar lining fluid volume kg=LL
    VLuT <- VLu - VLuF - VLuAF #lung tissue volume kg=L
    
    #Spleen
    PVSP <- 0.2e-2  #Brown et al. 1997, p 416, Table 5
    VSP <- PVSP * BW
    PVSPB <- 0.22 #Brown et al. 1997, p 458, Table 30
    VSPB <- PVSPB * PVSP * BW #volume of the blood of spleen kg=L
    PVSPF <- 0.554/280 #ml/g BW-> Shah & Betts, 2012. https://doi.org/10.1007/s10928-011-9232-2
    VSPF <- PVSPF * BW #spleen interstitial fluid volume kg=L
    VSPT <- VSP - VSPF #spleen tissue volume kg=L
    
    #Heart
    PVH <- 0.33e-2  #Brown et al. 1997, p 416, Table 5
    VH <- PVH * BW
    PVHB <- 0.26 #Brown et al. 1997, p 458, Table 30
    VHB <- PVHB * PVH * BW #volume of the blood of heart kg=L
    PVHF <- 0.1 #pkSim
    VHF <- PVHF * PVH * BW #heart interstitial fluid volume kg=L
    VHT <- VH - VHF #heart tissue volume kg=L
    
    #Brain
    PVBr <- 0.57e-2  #Brown et al. 1997, p 416, Table 5
    VBr <- PVBr * BW
    PVBrB <- 0.03 #Brown et al. 1997, p 458, Table 30
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
    PVSKB <- 0.02 #Brown et al. 1997, p 458, Table 3
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
    
    ##Capillary surface area for each tissue (Ai) as percentage of body weight (m^2/kg),
    #values from pkSim "Endothelial Surface area", Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    
    #Peritubular capillary density in From Gazzard et al. (2024) [https://doi.org/10.1002/ar.25576] 
    d_peritubular <- 0.024 #um^2/um^3,
    VK_gazzard <- 0.00363#L
    Vcortex <- 2*1295* VK/VK_gazzard #mm^3, NEED TO RECHECK IF IT REFERS TO ONE OF BOTH KIDNEYS
    A_peritubular <- (d_peritubular*Vcortex*1e3/1e6) #m^2
    A_peritubular_PTC <- A_peritubular * VPTC/ (VPTC + VDALC)
    A_peritubular_DTC <- A_peritubular * VDALC/ (VPTC + VDALC)
    
    # We assume that surface areas scale to the power of 2/3, based on the way that total body SA scales [https://doi.org/10.1111/j.1476-5381.2009.00267.x]
    BW_PKSim <- 0.23 #kg
    SA_scaling_factor <- (BW/BW_PKSim)^(2/3)  
    
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
    
    #Enhancement factors due to villi and microvilli
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
    
    ###############################
    #-----------------------------#
    #   Reflection Coefficients   #
    #-----------------------------#
    ###############################
    # Pore diameters from Price & Gesquiere (2020), doi:https://doi.org/10.1126/sciadv.aax2642
    DpKi <- 200 #nm
    DpLi <- 280 #nm
    DpSt <- 80 #nm, assumption
    DpIn <- 80 #nm
    DpMu <- 80 #nm
    DpAd <- 80 #nm, assumption
    DpRe <- 80 #nm, assumption
    DpLu <- 27 #nm
    DpSp <- 5000 #nm
    DpHt <- 50 #nm
    DpBr <- 0.99 #nm
    DpGo <- 80 #nm, assumption
    DpSk <- 60 #nm
    DpBo <- 40000 #nm
    
    Dps <- c(DpKi, DpLi, DpSt, DpIn, DpMu, DpAd, DpRe, DpLu, DpSp, DpHt, DpBr, DpGo, DpSk, DpBo)
    s_r <- rep(NA, length(Dps))
    # (C-C) bond length is 0.154 nm ==> 7*0.154 = 1.078nm
    # For carboxyl group we assume 0.13nm, So the total size is around 1.2 nm
    np_size <- 1.2/2 #nm, PFOA equivalent radius
    
    for (i in 1:length(s_r)){
      a_r <- np_size/(Dps[i]/2)
      Phi = (1-a_r)^2
      F_r <- (((1-a_r^2)^(3/2))*Phi)/(1+0.2*(a_r^2)*(1-a_r^2)^16)
      G_r <- ((1- (2*a_r^2)/3 - 0.20217*a_r^5 )/ (1-0.75851*a_r^5)) - (0.0431*(1-(1-a_r^10)))
      s_r[i] <- 1-(1-(1-Phi)^2)*G_r+2*a_r^2*Phi*F_r
    }
    SKi <- s_r[1] 
    SLi <- s_r[2]
    SSt <- s_r[3]
    SIn <- s_r[4]
    SMu <- s_r[5]
    SAd <- s_r[6]
    SRe <- s_r[7]
    SLu <- s_r[8]
    SSp <- s_r[9]
    SHt <- s_r[10]
    SBr <- 1
    SGo <- s_r[12]
    SSk <- s_r[13]
    SBo <- s_r[14]
    
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
    
    #############################################
    #               Lymph flow rates            #
    #############################################
    #Paracellular flow as a fraction of organ blood flow, 
    #from Niederalt et al.(2017). https://doi.org/10.1007/s10928-017-9559-4
    fQparaKi <- 7.09E-4
    fQparaLi <- 1.99E-2
    fQparaSt <- 2.04E-3
    fQparaIn <- (1.95E-3+1.44E-2)/2
    fQparaMu <- 2.01E-3
    fQparaAd <- 7.54E-3 
    fQparaRe <- 2.0E-3 # Assumption based on 1/500 of flow (Dosgra et al. 2020, https://doi.org/10.1016/j.csbj.2020.02.014)
    fQparaLu <- 3.56E-5
    fQparaSp <- 1.99E-2
    fQparaHt <- 1.47E-3
    fQparaBr <- 7.27E-5
    fQparaGo <- 1.11E-2
    fQparaSk <- 3.52E-3
    fQparaBo <- 6.62E-4 
    
    #Estimation of paracellular flow rates:
    QparaKi <- fQparaKi*QBK
    QparaLi <- fQparaLi*QBL
    QparaSt <- fQparaSt*QBST
    QparaIn <- fQparaIn*QBIN
    QparaMu <- fQparaMu*QBM
    QparaAd <- fQparaAd*QBA
    QparaRe <- fQparaRe*QBR
    QparaLu <- fQparaLu*QBLu
    QparaSp <- fQparaSp*QBSP
    QparaHt <- fQparaHt*QBH
    QparaBr <- fQparaBr*QBBr
    QparaGo <- fQparaGo*QBGo
    QparaSk <- fQparaSk*QBSK
    QparaBo <- fQparaBo*QBBo
    
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
    QGE<- 1.41/BW^0.25 #gastric emptying time (1/(h*BW^0.25)); from Yang, 2013
    
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
    QTDL <- Q_scaling_factor * QTDL_ref #L/h, Thin descending limb (short +long)
    QTAL <- Q_scaling_factor * QTAL_ref #L/h, Thin+thick (short +long) ascending 
    QDT <- Q_scaling_factor * QDT_ref  #L/h, Distal tubule (short + long)
    QCD <- Q_scaling_factor * QCD_ref  #L/h, collecting duct (Cortical + medullary)
    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    #Surface areas Interstitial - Intracellular (m^2), from PKSim 
    BW_ref <- 0.23
    AcK_total= A_peritubular*15#437.16*BW/BW_ref
    AcK_PTC <- AcK_total*VPTC/VFil # surface area of decending/ascending limb cells (loop of Henle)
    AcK_DALC <- AcK_total*VDALC/VFil # surface area of decending/ascending limb cells (loop of Henle)
    AcK_DTC <- AcK_total*VDTC/VFil # surface area of decending/ascending limb cells (loop of Henle)
    AcK_CDC <- AcK_total*VCDC/VFil # surface area of collecting duct cells 
    AcKTrest <- AcK_total* VKTrest/VFil
    AcL=  AL*15*8#84.45*BW/BW_ref
    AcST= AST*15#1007.31*BW/BW_ref
    AcIN= AIN*15 #400.94+152.39) *BW/BW_ref  small+large intestine
    AcM = AM*15 #8.2*BW/BW_ref
    AcA= AA*15 #3.87*BW/BW_ref
    AcLu= ALu*15#0.05*BW/BW_ref
    AcSP= ASP*15#564.05*BW/BW_ref
    AcH= AH*15#5.60*BW/BW_ref
    AcBr= ABr*15#6.12e-4*BW/BW_ref
    AcGo= AGo*15#2.01*BW/BW_ref
    AcSK= ASK*15#0.11*BW/BW_ref
    AcBo= ABo*15#6.52*BW/BW_ref
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
    AcALF = ((78.8*2*5320*1e-6) + (125*2*123*1e-6))*BW/0.29  #Stone et al., 1992, BW_ref = 0.29, values for each lung , https://doi.org/10.1165/ajrcmb/6.2.235
    
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
      'VLuAF'=VLuAF, 'VLuT'=VLuT,
      'VSP'=VSP, 'VSPB'=VSPB, 'VSPF'=VSPF, 'VSPT'=VSPT,
      'VH'=VH, 'VHB'=VHB, 'VHF'=VHF, 'VHT'=VHT,
      'VBr'=VBr, 'VBrB'=VBrB, 'VBrF'=VBrF, 'VBrT'=VBrT,
      'VGo'=VGo, 'VGoB'=VGoB, 'VGoF'=VGoF, 'VGoT'=VGoT,
      'VIN'=VIN, 'VINB'=VINB, 'VINF'=VINF, 'VINT'=VINT,
      'VST'=VST, 'VSTB'=VSTB, 'VSTF'=VSTF, 'VSTT'=VSTT,
      'VSTL'=VSTL, 'VINL'=VINL,
      'VSK'=VSK,'VSKB'=VSKB, 'VSKF'=VSKF, 'VSKT'=VSKT,
      'VBo'=VBo,'VBoB'=VBoB, 'VBoF'=VBoF, 'VBoT'=VBoT,
      
      
      'A_peritubular_PTC' = A_peritubular_PTC, 
      'A_peritubular_DTC' = A_peritubular_DTC,
      'AL'=AL, 'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
      'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
      'AIN'=AIN, 'AGo'=AGo,
      'ASK'= ASK, 'ABo'=ABo,
      
      'AINL' = AINL, 'AcL' = AcL, 'AcM' = AcM, 'AcST' = AcST, 
      'AcIN' = AcIN, 'AcA' = AcA, 'AcLu' = AcLu, 'AcALF' = AcALF, 
      'AcSP' = AcSP, 'AcH' = AcH, 'AcBr' = AcBr, 'AcGo' = AcGo, 
      'AcSK' = AcSK, 'AcBo' = AcBo, 'AcR' = AcR, 'APT' = APT, 
      'ADAL' = ADAL, 'ADT' = ADT, 'ACD' = ACD, 'AcK_DALC' = AcK_DALC,  
      'AcK_CDC' = AcK_CDC, 'AcKTrest' = AcKTrest,
      'AcK_PTC' = AcK_PTC,   "AcK_DTC" = AcK_DTC, 
      
      "SKi" = SKi,"SLi" = SLi,"SSt" = SSt,"SIn" = SIn,
      "SMu" = SMu,"SAd" = SAd,"SRe" = SRe,"SLu" = SLu,
      "SSp" = SSp,"SHt" = SHt,"SBr" = SBr,"SGo" = SGo,
      "SSk" = SSk,"SBo" = SBo,
      
      'Qcardiac'=Qcardiac, 'QBK'=QBK, 
      'QBL'=QBL, 'QBLtot'=QBLtot,
      'QBM'=QBM, 'QBA'=QBA,
      'QBR'=QBR, 'QBLu'=QBLu, 'Qfeces'=Qfeces,  'feces_density'=feces_density,
      'Qbile'=Qbile, 'QGFR'=QGFR,'Qurine'=Qurine, 
      'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST,
      'QBIN'=QBIN, 'QGE'=QGE,
      'QBGo'=QBGo,
      'QBSK'=QBSK, 'QBBo'=QBBo, 'Hct' = Hct,
      
      "QparaKi" = QparaKi,"QparaLi" = QparaLi,"QparaSt" = QparaSt,"QparaIn" = QparaIn,
      "QparaMu" = QparaMu,"QparaAd" = QparaAd,"QparaRe" = QparaRe,"QparaLu" = QparaLu,
      "QparaSp" = QparaSp,"QparaHt" = QparaHt,"QparaBr" = QparaBr,"QparaGo" = QparaGo,
      "QparaSk" = QparaSk,"QparaBo" = QparaBo,
      
      "QPT" = QPT, "QTDL" = QTDL, "QTAL" = QTAL, "QDT" = QDT, "QCD" = QCD,
      
      'f_tubular' =  f_tubular,  'f_PTC_prot_to_tub_prot' = f_PTC_prot_to_tub_prot, 
      'f_DALC_prot_to_tub_prot' = f_DALC_prot_to_tub_prot, 
      'f_DTC_prot_to_tub_prot' = f_DTC_prot_to_tub_prot , 
      'f_CDC_prot_to_tub_prot' = f_CDC_prot_to_tub_prot,
      'kabST'=kabST, 
      
      
      
      "admin.time" = admin.time, "admin.dose" = admin.dose,
      "admin.type" = admin.type, "MW"=MW
      
    ))
    
  })
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



