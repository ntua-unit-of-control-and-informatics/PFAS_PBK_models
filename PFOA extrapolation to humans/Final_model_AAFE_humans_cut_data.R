# Check notes for description of scenario
#female increase in BW: 3.5 g/d
#Male increase in BW: 5.9 g.d

library(deSolve)

create_variable_params <- function(BW,sex, estimated_params, fixed_params){
  
  kabsUA <-0.1 #L/h/m^2
  kCLEua <- 0 #L/h/m^2
  
  
  RAFOatp_lu_ap <- 0#0.00188877 #4.992703e-01#4.060400e-01
  RAFOatp_lu_bas <- 0.42333523#RAFOatp_lu_ap
  k_desorption_fast <-  10.86622465 #1/h
  k_desorption_slow <-  0.01364133  #1/h
  f_fast <- plogis( 0.34862342 )# fraction of PFOA released fast from dust
  fdust_retained <-  plogis(0.09039221)
  
  kUAB <- kabsUA * fixed_params$Nasal_SA #absorption rate from upper airways to blood, in L/h
  CLEua <- kCLEua * fixed_params$Nasal_SA #clearance rate rate from upper airways to stomach, in L/h
  
  # if(sex == "M"){
  #   RAFOatp_k <-  estimated_params[1] #1.567845e+01
  # }else{
  #   RAFOatp_k <- estimated_params[2] #3.950214e-03 
  # }
  
  RAFOatp_k <- estimated_params[1] #1.567845e+01
  
  RAFOat3 <- estimated_params[2] #0#9.75e+03
  
  CL_int <- 1.961917e+00 #uL/min/million hepatocytes
  
  HEPGL <- 139  #million hepatocytes/gram of human liver (Sohlenius-Sternbeck, 2006, [doi: https://doi.org/10.1016/j.tiv.2006.06.003] )
  # Scaled hepatobiliary clearance (assume liver density ~1 kg/L so volume equals mass)
  CL_hepatobiliary <- CL_int*1e-6*HEPGL*(fixed_params$VLi_tot*1000) *60 #L/h
  
  RAFOatp_l <- 4.992703e-01
  RAFUrat <- RAFOatp_k
  RAFOat1 <- RAFOat3
  RAFOatp2_l <- RAFOatp_l
  
  RAFNtcp <- RAFOatp_l
  RAFOatp2_Int <- 6.158004e-06
  
  # The parameters below were part of tests involving hypothetical transporters in the kidneys.
  # In the final model these transporters are switched off.
  VmK_api <- 0
  VmK_baso <- 0
  KmK_baso <- 1e20
  KmK_api <-   5e4
  
  KLfabp <-  1.2e5#[L/mol]. From Sheng et al. (2018) [doi:10.1007/s00204-017-2055-1]
  Ka <- 5e5#2.391203e+04#[L/mol]. From Rue et al. (2024)
  
  CLfeces_unscaled <- estimated_params[3] #8.214905e-05 #in L/h/BW^(-0.25), scaling similar to Loccisano et al. (2012)
  CLfeces <- CLfeces_unscaled*BW^(-0.25) 
  
  f_alb_avail<-  1
  f_fabp_avail  <- 1
  
  koff_alb <-  1e-2*3600#[1/h] #assumed 
  koff_fabp <-  1e-2*3600#[1/h] #assumed
  
  reduction_factor <- 3.714130e-02 
  
 
  #===================================================================================================================================================================
  

  f_tubular <- 0.52 # https://doi.org/10.1016/j.kint.2018.11.028, Proximal tubule cells account for roughly 52% of the estimated 206 million tubule epithelial cells per kidney.
  f_PTC_prot_to_total_prot <- 0.65 #https://doi.org/10.1681/ASN.2019040415


  # For the section below, we assume that all measurement were made relative to rinsed organs that 
  # contained only capillary blood weight. Consequently, we multiply by the mass of the drained organ
  
  MW = 414.07 #g/mol, PFOA molecular weight
  
  #Kidney
  kidney_protein_total <- 17e-2*(1e6*fixed_params$VKi_tot) #mg https://www.icrp.org/publication.asp?id=ICRP%20Publication%2023
  PTC_protein <- f_PTC_prot_to_total_prot*kidney_protein_total #mg
  
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
  liver_protein_total <- 18e-2*(1e6*fixed_params$VLi_tot) #mg https://www.icrp.org/publication.asp?id=ICRP%20Publication%2023
  
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
  lung_protein_total <-  17.8e-2*(1e6*fixed_params$VLu_tot) #mg https://www.icrp.org/publication.asp?id=ICRP%20Publication%2023
  
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
  intestine_protein_total <- 13e-2*(1e6*fixed_params$VIn_tot) #mg https://www.icrp.org/publication.asp?id=ICRP%20Publication%2023
  #oatp2b1-intestine 
  VmIn_Oatp2_in_vitro= 456.63e-3 #nmol/mg protein/min  (Kimura et al., 2017) 
  #assuming that the mediated transport is performed only by this transporter
  VmIn_Oatp2_scaled = 60*VmIn_Oatp2_in_vitro*MW*intestine_protein_total/1000   #physiologically scaled to in vivo, ug/h
  VmIn_Oatp2 = VmIn_Oatp2_scaled*RAFOatp2_Int #in vivo value, in  ug/h, same RAF as in liver
  KmIn_Oatp2 = 8.3*MW #umol/L (Kimura et al., 2017) --> ug/L
  
  Mr_albumin <- 66500#g/mol
  Mr_fabp <- 12000 #g/mol, from Ockner and  Manning (1976) [doi:10.1172/JCI108510] 
  CalbB_init  <- f_alb_avail*4.25*10/Mr_albumin #g/dL --> mol/L, Calbumin in humans https://www.ncbi.nlm.nih.gov/books/NBK459198/

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
  
  ##From various studies, the interstitial fluid/plasma concentration ratio for albumin was 0.54 (Table 1 and Table 2). 
  # https://doi.org/10.1016/S0022-2275(20)38701-0
  #I replaced all the 0.5 values with 0.54
  
  IPR_K = 0.54
  IPR_L = 0.54
  IPR_ST = 0.54
  IPR_IN = 0.54 #(NEED TO BE CHANGED)
  IPR_M = 0.27 #mean interstitial ﬂuid-to-serum ratio, https://doi.org/10.1152/ajpendo.2000.278.2.E352
  IPR_A = 0.15 #mean interstitial ﬂuid-to-serum ratio, https://doi.org/10.1152/ajpendo.2000.278.2.E352
  IPR_Lu = 0.54
  IPR_Sp = 0.54
  IPR_H = 0.54
  IPR_SK = 0.62 #https://doi.org/10.1080/00365517409050824, subcutis
  IPR_Br = 0.54
  IPR_Go = 0.54 #assumption
  IPR_Bo = 0.54 #assumption
  IPR_R = (IPR_K+IPR_L+IPR_ST+IPR_IN+IPR_M+IPR_A+IPR_Lu+IPR_Sp+IPR_H+IPR_SK+IPR_Br+IPR_Go+IPR_Bo)/13 #average IPR of all the included organs (kg=L)
  
  if (sex == "M"){
    Hct = 0.47 #hematocrit for male human, https://www.ncbi.nlm.nih.gov/books/NBK259/
  }else if (sex == "F"){
    Hct= 0.42 #hematocrit for female human, https://www.ncbi.nlm.nih.gov/books/NBK259/
  }
  
  
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
  cytosolic_protein = 44.4 # mg/g https://doi.org/10.1124/dmd.123.001417
  
  L_FABP_L = cytosolic_protein*4.5/100*1000 #mg/L, L-FABP comprises 2–7% of human cytosolic protein, https://doi.org/10.1111/febs.12780
  CFabpLT_init = f_fabp_avail*(L_FABP_L*1e-3/14e3) #[mol/L]
  
  #LFABP concentration in kidney and liver tissue (mol/m^3) (difficult to find for humans!!!!!!!!!!!!!!!)
  CFabpKT_init <- f_fabp_avail*2.65*1e-6  #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
  
  
  kon_alb <- Ka * koff_alb #1/h/M
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
  kDalcTu <- ((2*Peff_monolayer/100) * fixed_params$ALH) *1000 #diffusion between descending/ascending cells and tubule filtrate
  kDtcTu <- ((2*Peff_monolayer/100) * fixed_params$ADT) *1000 #diffusion between distal tubule cells and tubule filtrate
  kCdcTu <- ((2*Peff_monolayer/100) * fixed_params$ACD) *1000 #diffusion between collecting duct cells and tubule filtrate
  
  #Diffusion rates in L/h between  tubule cells and interstitial space
  kDtcF <- ((2*Peff_monolayer/100) * fixed_params$AcK_DTC) *1000
  kPtcF <- ((2*Peff_monolayer/100) * fixed_params$AcK_PTC) *1000
  kDalcF <- ((2*Peff_monolayer) * fixed_params$AcK_LHC) *1000 #diffusion between proximal tubule cells and interstitial space
  kCdcF <- ((2*Peff_monolayer/100) * fixed_params$AcK_CDC) *1000 #diffusion between descending/ascending cells and interstitial space
  
  # Physiologic upper limit of pore size from Sarin et al. (2010) [10.1186/2040-2384-2-14]
  pore_diameters <- c(
    Ki = 9,
    Li = 105,     # use 135 as average for human liver
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
  f_kidney <- 0.35 #Haraldsson et al., 2008 [10.1152/physrev.00055.2006] 20-50% of the entire endothelial surface in human kidney
  f_liver <- 0.08 #Szafranska et al., 2021 [10.3389/fphys.2021.735573]2-20% of the human liver sinusoidal endothelial cells surface is covered by fenestrations
  
  #The spleen has discontinuous capillaries, which are a type of sinusoid,
  #not a fixed percentage of fenestrations in a capillary surface.
  #These capillaries are characterized by very large openings and gaps in the endothelium,
  #rather than discrete fenestrations like those in fenestrated capillaries,
  #which allow for a wide range of substances to pass through.
  #The exact "fraction of fenestrations" is not a standard metric, as the openings are part of the open,
  #discontinuous structure of the vessel wall itself.
  
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
    
    'CFabpKT_init'=CFabpKT_init,'CFabpLT_init'=CFabpLT_init, 
    
    'Ka'=Ka, 'KLfabp'=KLfabp,
    
    "koff_alb" = koff_alb, "koff_fabp" = koff_fabp,
    "kon_alb" = kon_alb, "kon_fabp" = kon_fabp,
    
    'f_tubular' =  f_tubular,  'f_PTC_prot_to_total_prot' = f_PTC_prot_to_total_prot,
    
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
    
    'kUAB'= kUAB, 'CLEua'=CLEua,
    'k_gut_in' = k_gut_in, 'k_gut_out' = k_gut_out,'kabST'=kabST,
    'kCdcF' = kCdcF, 'kDalcF' = kDalcF, 'kPtcF' = kPtcF, 'kDtcF' = kDtcF,
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
    
    "MW" = MW, "Hct" = Hct, "k_desorption_fast" = k_desorption_fast,  "k_desorption_slow" = k_desorption_slow,
    "f_fast" = f_fast, "fdust_retained" = fdust_retained
    
  ))
}  

create_fixed_params <- function(user.input){
  if (!"depfr_AF" %in% names(user.input)){
    user.input$depfr_AF <- 0
    
  }
  
  if (!"depfr_head" %in% names(user.input)){
    user.input$depfr_head <- 0
    
  }
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
    PVKi <- 0.44e-2 #Brown et al. 1997, Table 7
    VKi_tot <- PVKi * BW #kidney volume kg=L 
    PVKiB <- 0.23 #PK-Sim
    VKiB_cap <- PVKiB * PVKi * BW #kidney blood volume kg=L
    PVKiF <- 0.20 #PK-Sim
    VKiF <- PVKiF * PVKi * BW #kidney interstitial fluid volume kg=L
    
    #From Huang & Isoherranen, 2018, https://doi.org/10.1002/psp4.12321, Table 1
    #volumes of filtrate compartments for total of two kidneys,[L] 
    VPT_1 <- 0.0305 # Proximal tubule segment 1
    VPT_2 <- 0.0305 # Proximal tubule segment 2
    VPT_3 <- 0.0305 # Proximal tubule segment 3
    VPT_tot <- VPT_1 + VPT_2 + VPT_3 # Proximal tubule
    VLH_D <- 0.0027 # Loop of Henle, descending segment
    VLH_A <- 0.0027 # Loop of Henle, ascending segment
    VLH_tot <- VLH_D + VLH_A # Loop of Henle
    VDT_tot <- 0.0194 # Distal tubule 
    
    # The table multiplies every segment by the nephron count,
    # which is not fine for the collecting duct (CD), because CDs are a branching tree shared by many nephrons, #https://doi.org/10.1016/B978-0-12-391448-4.00011-3
    # calculations for collecting duct scaling the collecting duct by a branching factor
    RCD <- 0.1 #mm
    LCD <- 21 #mm
    tot_neph_per_CD <- 6 #An average of six nephrons drains into a collecting duct.
    N_CD <- 20*2 # 20 collecting ducts per kidney
    CD_neph <- N_CD * tot_neph_per_CD
    VCD_tot <- pi*(RCD^2)*LCD*CD_neph*1e-6 #L, Collecting duct volume
   
    VFil <-  VPT_tot + VLH_tot + VDT_tot + VCD_tot #L
    
    # We hypothesize that when kidney is weighted the renal tubule content remains inside 
    VKiT <- VKi_tot - VKiF - VKiB_cap #kidney tissue volume kg=L
    
    # Huang et al. (2018) [doi:10.1002/psp4.12321]
    # state that Table S1. The physiological volume of each tubular
    # segment was calculated as a cylinder and physiological
    # volume of the cellular compartment of each segment was assumed to be the same as the tubular volume
  
    VPTC <- VPT_tot 
    VLHC <- VLH_tot
    VDTC <- VDT_tot 
    VCDC <- VCD_tot 
    
    VBladder <- 0.350 # L https://doi.org/10.1111/j.1742-1241.2011.02763.x
    
 
    #Liver
    PVLi <- 2.57e-2 #Brown et al. 1997, Table 7
    VLi_tot <- PVLi * BW #liver volume kg=L
    PVLiB <- 0.17  #pkSim
    VLiB_cap <- PVLiB * PVLi * BW #liver blood volume kg=L
    PVLiF <- 0.16  #pkSim
    VLiF <- PVLiF * PVLi* BW #liver interstitial fluid volume kg=L
    PVBile <- 29*1e-3/65.26 #L/kg BW,  https://doi.org/10.18203/2320-6012.ijrms20241964
    VBile <- PVBile * BW 
    VLiT <- VLi_tot - VLiF - VLiB_cap - VBile #liver tissue volume kg=L
    
    
    #Intestine (small and large)
    PVIn <- (0.91+0.53)*1e-2 #Brown et al. 1997, Table 7
    VIn_tot <- PVIn * BW #intestine volume kg=L
    # In the following we add the plasma and blood cell volumes of the small and large intestine from Shah and betts, (2012)
    PVInB <- 0.02 #pkSim
    VInB_cap <- PVInB * PVIn * BW #intestine  blood volume kg=L
    PVInF <- 0.09 #pkSim
    VInF <- PVInF * PVIn * BW #intestine interstitial fluid volume kg=L
    VInT <- VIn_tot - VInF - VInB_cap #intestine tissue volume kg=L
    
    #Stomach
    PVSt <- 0.21e-2 #Brown et al. 1997, p 416, Table 7
    VSt_tot <- PVSt * BW #stomach volume kg=L
    PVStB <- 0.03 #from pkSim
    VStB_cap <- PVStB * PVSt * BW 
    PVStF <-  0.10 # from pkSim
    VStF <- PVStF * PVSt * BW 
    VStT <- VSt_tot - VStF - VStB_cap #stomach tissue volume kg=L
    
    #Stomach and intestine lumen
    # The values below are from Conell et al. (2008) [doi: 10.1211/jpp.60.1.0008] and represent measured average water content of fed and fasted state from Figure 4.
    PVStL <- 0.53/175 #mL/g BW
    VStL <- PVStL * BW 
    PVInL <- (5.4- 0.53)/175 # mL/g BW
    VInL <- PVInL * BW 
    
    #Stomach and intestine lumen (CHECK)
    
    # The small intestine is about 3.05 meters (10 feet) long in a living person (but about twice as long in a cadaver due to the loss of muscle tone)
    
    PVStL <- 0.768/64.1 #L/kg BW, max gastric volume/females, https://doi.org/10.1016/S0031-9384(01)00619-9
    VStL <- PVStL * BW #stomach lumen volume kg=L
    PVInL <- (0.894+0.792+0.678+0.598+0.442)/230 # mL/g BW, Funai et al., 2023 https://doi.org/10.1038/s41598-023-44742-y --> Figure 3C
    VInL <- PVInL * BW #intestine lumen volume kg=L
    
    #Muscle 
    PVMu <- 40e-2 #Brown et al. 1997, Table 7
    VMu_tot <- PVMu * BW #muscle volume kg=L
    PVMuB <- 0.03 #pkSim
    VMuB_cap <- PVMuB * PVMu * BW #muscle blood volume kg=L
    PVMuF <- 0.16 #pkSim
    VMuF <- PVMuF * PVMu * BW #muscle interstitial fluid volume kg=L
    VMuT <- VMu_tot - VMuF - VMuB_cap #muscle tissue volume kg=L
  
    #Adipose
    PVAd <- 21.42e-2 #Brown et al. 1997, Table 7
    VAd_tot <- PVAd * BW #adipose volume kg=L
    PVAdB <- 0.02 #pkSim
    VAdB_cap <- PVAdB * PVAd * BW #% adipose blood volume kg=L
    PVAdF <- 0.16 #pkSim
    VAdF <- PVAdF * PVAd * BW #adipose interstitial fluid volume kg=L
    VAdT <- VAd_tot - VAdF - VAdB_cap#adipose tissue volume kg=L
    
    # Lung
    #Upper airways (no body weight scaling)
    #PVUA <- 257e-3/288 #mm^3/g, for a 16-wk old male 288 g, Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    VUA <- 29.32e-3 # cm^3 --> L, https://doi.org/10.1016/j.jobcr.2025.01.009 
    PVLu <- 0.76e-2 #Brown et al. 1997, Table 7
    VLu_tot <- PVLu * BW
    PVLuB <- 0.58 #pkSim
    VLuB_cap <- PVLuB * PVLu * BW #Brown et al. 1997, p 459 --> capillary blood occupied 9% of the lung volume
    PVLuF <- 0.19 #pkSim
    VLuF <- PVLuF * PVLu * BW #lung interstitial fluid volume
    PVLuAF <- 36e-3/70 #70 kg is an assumption, the value is reported as an estimation, https://doi.org/10.3389/fphys.2012.00146
    VLuAF <- PVLuAF * BW #lung alveolar lining fluid volume kg=LL
    VLuT <- VLu_tot - VLuF - VLuAF - VLuB_cap #lung tissue volume kg=L
    
    #Spleen
    PVSp <- 0.26e-2  #Brown et al. 1997, p 416, Table 7
    VSp_tot <- PVSp * BW
    PVSpB <- 0.33 #pkSim
    VSpB_cap <- PVSpB * PVSp * BW #volume of the blood of spleen kg=L
    PVSpF <- 0.15 #pkSim
    VSpF <- PVSpF * PVSp * BW #spleen interstitial fluid volume kg=L
    VSpT <- VSp_tot - VSpF - VSpB_cap#spleen tissue volume kg=L
    
    #Heart
    PVHt <- 0.47e-2  #Brown et al. 1997, p 416, Table 7
    VHt_tot <- PVHt * BW
    PVHtB <- 0.14 #pkSim
    VHtB_cap <- PVHtB * PVHt * BW #volume of the blood of heart kg=L
    PVHtF <- 0.10 #pkSim
    VHtF <- PVHtF * PVHt * BW #heart interstitial fluid volume kg=L
    VHtT <- VHt_tot - VHtF - VHtB_cap #heart tissue volume kg=L
    
    #Brain
    PVBr <- 2e-2  #Brown et al. 1997, p 416, Table 7
    VBr_tot <- PVBr * BW
    PVBrB <- 0.04 #pkSim
    VBrB_cap <- PVBrB * PVBr * BW #volume of the blood of brain kg=L
    PVBrF <- 4e-3 #pkSim
    VBrF <- PVBrF * PVBr * BW 
    VBrT <- VBr_tot - VBrF - VBrB_cap #brain tissue volume kg=L
    
    #gonads
    PVGo <- 0.04/70 #pKsim, L/kg
    VGo_tot <- PVGo * BW
    PVGoB <- 0.06 #pKsim
    VGoB_cap <-PVGoB * PVGo * BW #volume of the blood of gonads kg=L
    PVGoF <- 0.07 #pKsim
    VGoF <- PVGoF * PVGo * BW #gonads interstitial fluid volume kg=L
    VGoT <- VGo_tot - VGoF - VGoB_cap#gonads tissue volume kg=L
    
    #Skin
    PVSk <- 3.71e-2 #Brown et al. 1997, p 416, Table 7
    VSk_tot <- PVSk * BW
    PVSkB <- 0.02 #pkSim
    VSkB_cap <- PVSkB * PVSk * BW #volume of the blood of skin kg=L
    PVSkF <- 0.3 #pkSim
    VSkF <- PVSkF * PVSk * BW #skin interstitial fluid volume kg=L
    VSkT <- VSk_tot - VSkF - VSkB_cap #skin tissue volume kg=L
    
    #Bones (ok)
    PVBo <- 11.82/70 #pkSim, L/kg
    VBo_tot <- PVBo * BW
    PVBoB <- 0.03 #pkSim
    VBoB_cap <-PVBoB * PVBo * BW #volume of the blood of bones kg=L
    PVBoF <- 0.1 #pkSim
    VBoF <- PVBoF * PVBo * BW #bones interstitial fluid volume kg=L
    VBoT <- VBo_tot - VBoF - VBoB_cap#bones tissue volume kg=L
    
    VBlood <- (0.06 * BW *1000+ 0.77)/1000    # Total blood. From Lee and Blaufox (1981) [PMID: 3965655]

      if (sex == "M"){
        PVB = 78e-3 #L/kg BW https://www.oxfordreference.com/display/10.1093/oi/authority.20110803095513769
        VBlood = PVB* BW  #L
      }else if (sex == "F"){
        PVB = 56e-3 #L/kg BW https://www.oxfordreference.com/display/10.1093/oi/authority.20110803095513769
        VBlood = PVB* BW  #L
      }
    
    PVRe <- 1 - PVB - PVKi - PVLi - PVMu - PVAd - PVSp - PVHt - PVBr - PVGo - PVSt - PVIn - PVSk - PVBo
    VRe_tot <- PVRe * BW #volume of the rest of the body kg=LL
    PVReB <-(PVKiB+PVLiB+PVLuB+PVMuB+PVAdB+PVSpB+PVHtB+PVBrB+PVGoB+PVStB+PVInB+PVSkB+PVBoB)/13 #average VF of all the included organs (kg=L)
    VReB_cap <- PVReB * PVRe * BW #volume of the blood of RoB kg=L
    PVReF <-(PVKiF+PVLiF+PVLuF+PVMuF+PVAdF+PVSpF+PVHtF+PVBrF+PVGoF+PVStF+PVInF+PVSkF+PVBoF)/13 #average VF of all the included organs (kg=L)
    VReF <- PVReF * PVRe * BW #RoB of the blood of rest of body kg=L
    VReT <- VRe_tot - VReF - VReB_cap #tissue volume of the rest of body kg=L
    

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
    BW_ref <- 70
    linear_scaling_factor <- BW/BW_ref
    nonlinear_scaling_factor <- (BW/BW_ref)^0.75
    ##Capillary surface area for each tissue (Ai)(m^2),
    
  
    # The following values from PK-Sim "Endothelial Surface area", Niederalt et al., 2018, [doi: 10.1007/s10928-017-9559-4]
    PAL <- 38.05 #m^2
    AL <- PAL *  linear_scaling_factor #liver surface area (m^2)
    PAST <- 0.51 #m^2
    AST <- PAST * linear_scaling_factor #stomach surface area (m^2)
    PAIN <- (1.32+0.75) #m^2 #large and small intestine
    AIN <- PAIN * linear_scaling_factor #intestine surface area (m^2)
    PAM <- 75.95 #m2
    AM <- PAM * linear_scaling_factor #muscle surface area (m^2)
    PAA <- 21.26 #m2
    AA <- PAA * linear_scaling_factor #adipose surface area (m^2)
    PAR <- 10 #m^2, assumption
    AR <- PAR * linear_scaling_factor #surface area of rest of body (m^2)
    PLu <- 66.52 #m^2 
    ALu <- PLu * linear_scaling_factor #lung surface area (m^2)
    PSP <- 6.47 #m^2
    ASP <- PSP * linear_scaling_factor #spleen surface area (m^2)
    PH <-  5.55 #m^2 
    AH <- PH * linear_scaling_factor #heart surface area (m^2)
    PBr <- 5.59 #m^2
    ABr <- PBr * linear_scaling_factor #brain surface area (m^2)
    PGo <- 0.21 #m^2
    AGo <- PGo * linear_scaling_factor #gonads surface area (m^2)
    PSK <- 16.09 #m^2
    ASK <- PSK * linear_scaling_factor #skin surface area (m^2)
    PBo <- 38.17 #m^2
    ABo <- PBo * linear_scaling_factor #bones surface area (m^2)
    
    
    # Gut surface areas from PkSim
    Duodenum <- 628.32 #cm^2
    Upper_jejunum <- 510.06 #cm^2
    Lower_jejunum <- 477.16 #cm^2
    Upper_ileum <- 635.91 #cm^2
    Lower_ileum <- 538.08 #cm^2
    Cecum <- 126.78 #cm^2
    colon_ascendens <- 354.52 #cm^2
    colon_transversum <- 547.72 #cm^2
    colon_descendens <-265.66 #cm^2
    colon_sigmoid <-328.74 #cm^2
    
    SA_PKSim <- Duodenum+Upper_jejunum+Lower_jejunum+Upper_ileum+Lower_ileum+Cecum+     
      colon_ascendens+colon_transversum+colon_descendens+colon_sigmoid
    
    #Enhancement factors
    EF_Duodenum <- 292.69
    EF_Upper_jejunum <- 447.99
    EF_Lower_jejunum <- 372.94
    EF_Upper_ileum <- 260.75
    EF_Lower_ileum <- 146.57
    EF_Cecum <- 1.80
    EF_colon_ascendens <- 2.5 
    EF_colon_transversum <- 2.5
    EF_colon_descendens <- 2.5
    EF_colon_sigmoid <- 3.56
    
    
    SA_PKSim_effective <- (Duodenum*EF_Duodenum+Upper_jejunum*EF_Upper_jejunum+Lower_jejunum*EF_Lower_jejunum+
                             Upper_ileum*EF_Upper_ileum+Lower_ileum*EF_Lower_ileum+Cecum*EF_Cecum+
                             colon_ascendens*EF_colon_ascendens+colon_transversum*EF_colon_transversum+
                             colon_descendens*EF_colon_descendens+colon_sigmoid*EF_colon_sigmoid)
    
    
    AINL <- SA_PKSim_effective * 1e-4 #m^2, no body weight scaling
    
    #Surface areas Interstitial - Intracellular (m^2), from PK-Sim 
    AcK_total= 22228.45*nonlinear_scaling_factor
    # for renal tubule, we assume that the total surface area is split proportionally to the volume of each segment.
    AcK_PTC <- AcK_total*VPTC/VFil # surface area of proximal tubule cells
    AcK_LHC <- AcK_total*VLHC/VFil # surface area of loop of Henle cells
    AcK_DTC <- AcK_total*VDTC/VFil # surface area of distal tubule cells
    AcK_CDC <- AcK_total*VCDC/VFil # surface area of collecting duct cells 
    AcL= 4930.84*nonlinear_scaling_factor
    AcST= 43538.81*nonlinear_scaling_factor
    AcIN= (19458.36+12682.70) *nonlinear_scaling_factor # small+large intestine
    AcM= 530.51*nonlinear_scaling_factor
    AcA= 804.69*nonlinear_scaling_factor
    AcLu= 10.89*nonlinear_scaling_factor
    AcSP= 44741.62*nonlinear_scaling_factor
    AcH= 606.98*nonlinear_scaling_factor
    AcBr= 0.1*nonlinear_scaling_factor
    AcGo= 16.08*nonlinear_scaling_factor
    AcSK= 3.35*nonlinear_scaling_factor
    AcBo= 926.21*nonlinear_scaling_factor
    # We don't have data for the surface area of IS-IC for the rest of the body, thus 
    # we naively assume an average:
    AcR= median(c(AcK_total,AcL,AcST,AcIN,AcM,AcA,AcLu,AcSP,AcH,AcBr,AcGo,AcSK,AcBo))
    
    n <- 5 #enlargement factor of the apical membrane of tubule cells
    
    #From Huang & Isoherranen, 2018, https://doi.org/10.1002/psp4.12321, Table 1
    # Surface areas of the different subcompartments of kidney filtrate, dm^2 --> m^2

    APT_1 <- 6107*1e-2 # Proximal tubule segment 1
    APT_2 <- 6107*1e-2  # Proximal tubule segment 2
    APT_3 <- 6107*1e-2  # Proximal tubule segment 3
    APT_tot <- APT_1 + APT_2 + APT_3 # Proximal tubule
    ALH_D <- 61*1e-2  # Loop of Henle, descending segment
    ALH_A <- 61*1e-2  # Loop of Henle, ascending segment
    ALH_tot <- ALH_D + ALH_A # Loop of Henle
    ADT_tot <- 156*1e-2  # Distal tubule 
    ACD_1 <- 6.7*1e-2 # Collecting duct segment 1
    ACD_2 <- 6.7*1e-2 # Collecting duct segment 2
    ACD_3 <- 6.7*1e-2 # Collecting duct segment 3
    ACD_4 <- 6.7*1e-2 # Collecting duct segment 4
    ACD_5 <- 6.7*1e-2 # Collecting duct segment 5
    ACD_tot <- ACD_1 + ACD_2 + ACD_3 + ACD_4 + ACD_5
    
    AFil <-  APT_tot + ALH_tot + ADT_tot + ACD_tot #m^2
    
    #Peritubular capillary density from Weijden et al. (2023) [doi: 10.1007/s40620-023-01734-5] 
    #Peritubular capillary density was assessed as number of PTCs per tubule (PTC/tubule)
    #and number of PTCs per surface area (PTC/50,000 μm2).
    
    d_peritubular <- 25.9 # PTCs per tubular surface area (PTC/50,000 μm2)
    PTC_diameter <- 14.5 #μm, assumption PTCs are simple cuboidal epithelial cells, epithelial cell diameters from 8-21 microns, 10.1016/s0014-4835(87)80052-0
    PTC_mean_capillary_length <- 75 #μm, assumption 
    A_PTC_mean <- pi*PTC_diameter*PTC_mean_capillary_length
    N_PTC_total <- d_peritubular * (AFil * 1e12 / 50000)
    A_peritubular <- N_PTC_total * A_PTC_mean * 1e-12   # m2
    

    #Alveolar cells surface area, m^2
    AcALF = 0.4* (BW/0.363)  #Stone et al., 1992 [doi:10.1165/ajrcmb/6.2.235]
    #nasal surface area, m^2
    Nasal_SA <-  18.5*1e-4 #Oller and Oberdörster (2010) [doi:https://doi.org/10.1016/j.yrtph.2010.02.006]
    
    #----------------------------------------------------------------------------------------------
    #-------------------------             Flow Rates           -----------------------------------
    #----------------------------------------------------------------------------------------------
    
    ########
    #Blood flow rates 
    ###
    
    #(QBi, in L/h) to different tissues expressed using the corresponding percentage of cardiac output
    
    Qcardiac <- 5200/1000 *(BW^0.75) *60 #[L/h]
    PQBK <- 19/100 #(10.1002/jps.22811--> 1064/5600*100)
    QBK <- PQBK * Qcardiac #L/h
    PQBL <- 6.5/100 #(10.1002/jps.22811--> 364/5600*100)
    QBL <- PQBL * Qcardiac #L/h Brown et al. 1997, p 438, Table 23
    PQBST <- 1/100 #[ICRP] 
    QBST <- PQBST * Qcardiac #L/h
    PQBIN <- 5/100 #[ICRP, mean value for small and large intestine] 
    QBIN <- PQBIN * Qcardiac #L/h
    PQBM <- 17/100 #(10.1002/jps.22811--> 952/5600*100)
    QBM <- PQBM * Qcardiac #L/h
    PQBA <- 5.2/100 #Brown et al. 1997, p 438, Table 23
    QBA <- PQBA * Qcardiac #L/h
    PQBLu <- 1 
    QBLu <- PQBLu * Qcardiac #L/h
    PQBSP <- 2/100 #[doi: 110.1002/prp2.305]  
    QBSP <- PQBSP * Qcardiac #L/h
    PQBH <- 4/100 #Brown et al. 1997, p 438, Table 23
    QBH <- PQBH * Qcardiac #L/h
    PQBBr <- 11.4/100 #Brown et al. 1997, p 438, Table 23
    QBBr <- PQBBr * Qcardiac #L/h 
    PQBGo <- 0.5/100 #assumption
    QBGo <- PQBGo * Qcardiac #L/h
    PQBSK <- 5.8/100 # Brown et al. (1997)
    QBSK <- PQBSK * Qcardiac #L/h
    PBBo <- 4.2/100 #Brown et al. 1997, p 438, Table 23
    QBBo <- PBBo * Qcardiac #L/h
    
    # Total blood outflow from liver
    QBLtot <- QBL+QBSP+QBIN+QBST
    
    PQBR = 1 - PQBK - PQBL - PQBST - PQBIN - PQBM - PQBA - PQBH - PQBSK - PQBSP - PQBGo - PQBBr - PBBo
    QBR <- PQBR * Qcardiac #L/h
    
    
    #########
    #     Other fluids flow rates    
    ###
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration rate (GFR), in L/h
    

    Qbile = 900*1e-3/24 #mL/day --> L/h, https://www.ncbi.nlm.nih.gov/books/NBK279386/
    Vbile = 55*1e-3 #L, https://www.ncbi.nlm.nih.gov/books/NBK279386/
 
    Qfeces <- 128*1e-3/24 #g/day --> L/h, https://doi.org/10.1080/10643389.2014.1000761
    feces_density <- 1.06 #g/cm^3 --> g/mL Brown et al., 1996 [https://doi.org/10.1016/0273-1223(96)00366-6]
    
    # if (sex == "M"){
    #   PQGFR <- 62.1  #L/h/kg kidney. From Corley et al., 2005 [doi: 10.1093/toxsci/kfi119]
    #   QGFR <- PQGFR * VKi_tot #L/h
    #   Qurine <- (60/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, [doi: 10.1002/nau.1006]
    # }else if(sex == "F"){
    #   PQGFR <- 41.04 #L/h/kg kidney. From Corley et al., 2005 [doi: 10.1093/toxsci/kfi119]
    #   QGFR <- PQGFR * VKi_tot#L/h
    #   Qurine <- (85/1000)*BW/24 #([ml/d/kg]/1000)*BW/24 --> L/h, from Schmidt et al., 2001, [doi: 10.1002/nau.1006] 
    # }
    
    QGFR <- 116.66*60/1000 #PK-Sim, mL/min --> L/h
    
    if (sex == "M"){
      PQurine <- 0.212  #L/h/kg kidney. From Corley et al., 2005 [doi: 10.1093/toxsci/kfi119]
      Qurine <- PQurine * VKi_tot #L/h
    }else if(sex == "F"){
      PQurine <- 0.152  #L/h/kg kidney. From Corley et al., 2005 [doi: 10.1093/toxsci/kfi119]
      Qurine <- PQurine * VKi_tot #L/h
    }
    
    # #From Abraham 2024 data
    # Qurine <- 2.21*1e-14 #ug/h --> L/h
    
    Qelim <- 28*1e-3/24 #ng/day --> ug/h
    
    QGE<- 3.5*BW^(-0.25) #gastric emptying time [1/(h/BW^(-0.25)]*BW^(-0.25)--> [1/h]. From Yang et al. (2014) [doi: 10.1371/journal.pone.0106101]
    
    #From Huang & Isoherranen, 2018, https://doi.org/10.1002/psp4.12321, Table 1
    #flows of filtrate compartments 
    # The tubular flow rate exiting each tubular subsegment is equal to the tubular flow rate entering the next tubular subsegment.
    # Units: [mL/min]*1e-3*60 ---> L/h 
    
    QPT_1 <- 120*1e-3*60 # Proximal tubule segment 1
    QPT_2 <- 94*1e-3*60  # Proximal tubule segment 2
    QPT_3 <- 68*1e-3*60  # Proximal tubule segment 3
    QPT_tot_ref <- QPT_1 + QPT_2 + QPT_3 # Proximal tubule
    QLH_D <- 43*1e-3*60  # Loop of Henle, descending segment
    QLH_Q <- 24*1e-3*60  # Loop of Henle, ascending segment
    QLH_tot_ref <- QLH_D + QLH_Q # Loop of Henle
    QDT_tot_ref <- 24*1e-3*60  # Distal tubule 
    QCD_1 <- 11*1e-3*60 # Collecting duct segment 1
    QCD_2 <- 9*1e-3*60 # Collecting duct segment 2
    QCD_3 <- 7*1e-3*60 # Collecting duct segment 3
    QCD_4 <- 5*1e-3*60 # Collecting duct segment 4
    QCD_5 <- 3*1e-3*60 # Collecting duct segment 5
    QCD_tot_ref <- QCD_1 + QCD_2 + QCD_3 + QCD_4 + QCD_5
    
    Q_scaling_factor = QGFR/QPT_tot_ref
    
    QPT_tot <- Q_scaling_factor * QPT_tot_ref
    QLH_tot <- Q_scaling_factor * QLH_tot_ref
    QDT_tot <- Q_scaling_factor * QDT_tot_ref
    QCD_tot <- Q_scaling_factor * QCD_tot_ref
    
    return(list(
      
      'VKiF'=VKiF, 'VKiT'=VKiT, 'VFil'=VFil,'VBladder' = VBladder, 
      'VPT_tot' = VPT_tot, 'VLH_tot' = VLH_tot, 'VDT_tot' = VDT_tot, 'VCD_tot' = VCD_tot,
      'VPTC' = VPTC, 'VLHC' = VLHC, 'VDTC' = VDTC, 'VCDC' = VCDC,
      
      'VKiB_cap'=VKiB_cap,'VLiB_cap'=VLiB_cap,'VLiF'=VLiF, 'VLiT'=VLiT, 'VBile'=VBile,
      'VInB_cap'=VInB_cap,'VInF'=VInF, 'VInT'=VInT,
      'VStB_cap'=VStB_cap, 'VStF'=VStF, 'VStT'=VStT,
      'VStL'=VStL, 'VInL'=VInL,
      'VMuB_cap'=VMuB_cap, 'VMuF'=VMuF, 'VMuT'=VMuT, 
      'VAdB_cap'=VAdB_cap,'VAdF'=VAdF, 'VAdT'=VAdT, 
      'VLuB_cap'=VLuB_cap,'VUA' = VUA, 'VLuF'=VLuF,'VLuAF'=VLuAF, 'VLuT'=VLuT, 
      'VSpB_cap'=VSpB_cap,'VSpF'=VSpF, 'VSpT'=VSpT,
      'VHtB_cap'=VHtB_cap,'VHtF'=VHtF, 'VHtT'=VHtT,
      'VBrB_cap'=VBrB_cap,'VBrF'=VBrF, 'VBrT'=VBrT,
      'VGoB_cap'=VGoB_cap,'VGoF'=VGoF, 'VGoT'=VGoT,
      'VSkB_cap'=VSkB_cap,'VSkF'=VSkF, 'VSkT'=VSkT,
      'VBoB_cap'=VBoB_cap, 'VBoF'=VBoF, 'VBoT'=VBoT,
      'VReB_cap'=VReB_cap,'VReF'=VReF, 'VReT'=VReT, 
      
      "VKi_tot" = VKi_tot, "VLi_tot" = VLi_tot,  "VIn_tot" = VIn_tot,
      "VSt_tot" = VSt_tot, "VMu_tot" = VMu_tot, 
      "VAd_tot" = VAd_tot, "VLu_tot" = VLu_tot, "VSp_tot" = VSp_tot,
      "VHt_tot" = VHt_tot, "VBr_tot" = VBr_tot, "VGo_tot" = VGo_tot,
      "VSk_tot" = VSk_tot, "VBo_tot" = VBo_tot, "VRe_tot" = VRe_tot,
      
      "VBlood" = VBlood, 'VVen' = VVen, 'VArt' = VArt,
      
      'A_peritubular' = A_peritubular, 
      'AL'=AL, 'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
      'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
      'AIN'=AIN, 'AGo'=AGo,
      'ASK'= ASK, 'ABo'=ABo,
      
      'AINL' = AINL, 'AcL' = AcL, 'AcM' = AcM, 'AcST' = AcST, 
      'AcIN' = AcIN, 'AcA' = AcA, 'AcLu' = AcLu, 'AcALF' = AcALF, 
      'AcSP' = AcSP, 'AcH' = AcH, 'AcBr' = AcBr, 'AcGo' = AcGo, 
      'AcSK' = AcSK, 'AcBo' = AcBo, 'AcR' = AcR,
      'AcK_PTC' = AcK_PTC, 'AcK_LHC' = AcK_LHC, 'AcK_DTC' = AcK_DTC, 'AcK_CDC' = AcK_CDC,
      'APT_tot' = APT_tot, 'ALH_tot' = ALH_tot, 'ADT_tot' = ADT_tot,
      'ACD_tot' = ACD_tot, 'AFil' = AFil,
      "Nasal_SA" = Nasal_SA,
      
      'Qcardiac'=Qcardiac, 'QBK'=QBK, 
      'QBL'=QBL, 'QBLtot'=QBLtot,
      'QBM'=QBM, 'QBA'=QBA,
      'QBR'=QBR, 'QBLu'=QBLu, 'Qfeces'=Qfeces,  'feces_density'=feces_density,
      'Qbile'=Qbile, 'Vbile'=Vbile, 'QGFR'=QGFR,'Qurine'=Qurine, 'Qelim'=Qelim,
      'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST,
      'QBIN'=QBIN, 'QGE'=QGE,
      'QBGo'=QBGo,
      'QBSK'=QBSK, 'QBBo'=QBBo,
      
      
      "QPT_tot" = QPT_tot, "QLH_tot" = QLH_tot, "QDT_tot" = QDT_tot, "QCD_tot" = QCD_tot,
      "admin.time" = admin.time, "admin.dose" = admin.dose,
      "admin.type" = admin.type, "MW"=MW, "sex" = sex, "BW" = BW,
      "depfr_head" = depfr_head, "depfr_AF" = depfr_AF
      
      
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
  #                        'Qfeces','Qbile', 'QGFR','Qurine','kabST',"QPT" , "QLH_tot", "QTAL" , "QDT", "QCD", 'CLfeces',
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
    
    CPTCf <- MPTCf/VPTC
    CPTCb <- MPTCb/VPTC
    CLHCf <- MLHCf/VLHC
    CLHCb <- MLHCb/VLHC
    CDTCf <- MDTCf/VDTC
    CDTCb <- MDTCb/VDTC
    CCDCf <- MCDCf/VCDC
    CCDCb <- MCDCb/VCDC
    
    MKTf <- MPTCf + MLHCf + MDTCf + MCDCf
    MKTb <- MPTCb + MLHCb + MDTCb + MCDCb
    MKT <- MKTf + MKTb
    Mfil <- MPT+MLH+MDT+MCD
    CKTb <- MKTb/VKiT
    CKTf <- MKTf/VKiT
    CKT <- MKT/VKiT # tissue concentration
    #renal filtrate
    CPT <- MPT/VPT_tot #proximal tubule
    CLH <- MLH/ VLH_tot # tubule
    CDT <- MDT/ VDT_tot #distal tubule
    CCD <- MCD/ VCD_tot #proximal tubule
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
    CUA = MUA/VUA # upper airways concentration
    
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
    #CLuAFdust <- MLuAFdust/VLuAF
    
    
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
    dCFabpPTCf <- koff_fabp*CPTCb/MW/1e6 - kon_fabp*CFabpPTCf*CPTCf/MW/1e6
    dCFabpLHCf <- koff_fabp*CLHCb/MW/1e6 - kon_fabp*CFabpLHCf*CLHCf/MW/1e6
    dCFabpDTCf <- koff_fabp*CDTCb/MW/1e6 - kon_fabp*CFabpDTCf*CDTCf/MW/1e6
    dCFabpCDCf <- koff_fabp*CCDCb/MW/1e6 - kon_fabp*CFabpCDCf*CCDCf/MW/1e6
    
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
    dMPTCb <- kon_fabp*CFabpPTCf*CPTCf*VPTC - koff_fabp*CPTCb*VPTC
    
    dMLHCb <- kon_fabp*CFabpLHCf*CLHCf*VLHC - koff_fabp*CLHCb*VLHC
    
    dMDTCb <- kon_fabp*CFabpDTCf*CDTCf*VDTC - koff_fabp*CDTCb*VDTC
    
    dMCDCb <- kon_fabp*CFabpCDCf*CCDCf*VCDC - koff_fabp*CCDCb*VCDC
    
    dMLTb <- kon_fabp*CFabpLTf*CLTf*VLiT - koff_fabp*CLTb*VLiT 
    
    #Alveolar lining fluid
    dMLuAFb <-  kon_alb*CalbLuAFf*CLuAFf*VLuAF -  koff_alb*CLuAFb*VLuAF 
    
    #====================================================================================================================
    
    #Arterial Blood
    dMArtf = QBLu*CLuBf - CArtf*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+
                                   QBST+QBIN+QBGo+QBSK+QBBo) - QGFR*CArtf +
      koff_alb*CArtb*VArt - kon_alb*CalbArtf*CArtf*VArt 
    
    #Venous Blood
    dMVenf = kUAB*CUA - CVenf*QBLu + QBK*CKBf + QBLtot*CLBf + QBM*CMBf + QBA*CABf + QBR*CRBf+
      QBH*CHBf + QBBr*CBrBf+ QBGo*CGoBf + QBSK*CSKBf + QBBo*CBoBf+
      koff_alb*CVenb*VVen - kon_alb*CalbVenf*CVenf*VVen 
    
    #Kidney
    #blood subcompartment
    dMKBf = QBK*CArtf - QBK*CKBf- (Ptrans_diff_K+PparaKi)*A_peritubular*(CKBf-CKFf) + 
      (VmK_baso*CPTCf/(KmK_baso+CPTCf))+
      koff_alb*CKBb*VKiB_cap - kon_alb*CalbKBf*CKBf*VKiB_cap
    
    #interstitial fluid subcompartment
    dMKFf =  (Ptrans_diff_K+PparaKi)*A_peritubular*(CKBf-CKFf)- 
      kPtcF*(CKFf-CPTCf) - kDalcF*(CKFf-CLHCf) -
      kDtcF*(CKFf-CDTCf) - kCdcF*(CKFf-CCDCf) + 
      koff_alb*CKFb*VKiF - kon_alb*CalbKFf*CKFf*VKiF-
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) - (VmK_Oat3*CKFf/(KmK_Oat3+CKFf))
    
    #proximal tubule  cells subcompartment
    dMPTCf = kPtcF*(CKFf-CPTCf) - kPtcTu*(CPTCf - CPT)  +
      (VmK_Oatp*CPT/(KmK_Oatp+CPT)) + (VmK_Urat*CPT/(KmK_Urat+CPT))+
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) + (VmK_Oat3*CKFf/(KmK_Oat3+CKFf)) - 
      (VmK_baso*CPTCf/(KmK_baso+CPTCf)) -(VmK_api*CPTCf/(KmK_api+CPTCf))-
      (kon_fabp*CFabpPTCf*CPTCf*VPTC - koff_fabp*CPTCb*VPTC) 
    
    #Tubule cells in Loop of Henle 
    dMLHCf =  kDalcF*(CKFf-CLHCf) - kDalcTu*(CLHCf - CLH)- 
      (kon_fabp*CFabpLHCf*CLHCf*VLHC -
         koff_fabp*CLHCb*VLHC)
    
    #Distal convoluted tubule cells 
    dMDTCf =   kDtcF*(CKFf-CDTCf)- kDtcTu*(CDTCf - CDT)-
      (kon_fabp*CFabpDTCf*CDTCf*VDTC -
         koff_fabp*CDTCb*VDTC)
    
    #Collecting duct cells 
    dMCDCf =  kCdcF*(CKFf-CCDCf) - kCdcTu*(CCDCf - CCD)- 
      (kon_fabp*CFabpCDCf*CCDCf*VCDC -
         koff_fabp*CCDCb*VCDC)
    
    #Proximal convoluted tubule
    dMPT =  QGFR*CArtf + kPtcTu*(CPTCf - CPT) - (VmK_Oatp*CPT/(KmK_Oatp+CPT)) - 
      (VmK_Urat*CPT/(KmK_Urat+CPT)) + (VmK_api*CPTCf/(KmK_api+CPTCf))- QLH_tot*CPT
    
    #Descending limb, Ascending limb (Loop of Henle )
    dMLH =  QLH_tot*CPT + kDalcTu*(CLHCf - CLH)  -  QDT_tot*CLH
    
    # Distal convoluted tubule 
    dMDT =   QDT_tot*CLH + kDtcTu*(CDTCf - CDT) -  QCD_tot*CDT
    
    #Collecting duct
    dMCD = QCD_tot*CDT + kCdcTu*(CCDCf - CCD) - Qurine*CCD
    
    # Bladder
    dMBladder = Qurine*CCD - Qelim*CBladder
    
    
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
    dMBile = CL_hepatobiliary*CLTf - CBile*Qbile
    
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
    dMSTL = - QGE*MSTL -kabST*CSTL+ CLEua*CUA
    
    
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
    dMINL = QGE*MSTL - (CLfeces*CINL) - k_gut_in*CINL + k_gut_out*CINTf + CBile*Qbile - 
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
    #Upper airways
    dMUA =  - kUAB * CUA - CLEua*CUA 
    
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
    dMLuAFf =  k_desorption_fast * MLuAFdust_fast + k_desorption_slow * MLuAFdust_slow+ kLuTLuAF*(CLuTf-CLuAFf) + 
      koff_alb*CLuAFb*VLuAF - kon_alb*CalbLuAFf*CLuAFf*VLuAF -
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf)) 
    dMLuAFdust_fast <- -k_desorption_fast * MLuAFdust_fast
    dMLuAFdust_slow <- -k_desorption_slow * MLuAFdust_slow
    
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
    dMurine <- Qelim*CBladder
    dVurine = Qurine
    dVfeces = Qfeces
    
    #Concentration calculation in each compartment 
    Cblood <- (MVen +MArt)/ (VVen+VArt)
    Cplasma <- Cblood/(1-Hct)
    
    background_plasma <- 1.48                 # µg/L 
    Cblood_background <- background_plasma * (1 - Hct)  # convert to whole-blood concentration
    
    Mblood <- MVen +MArt
    Mkidney <- MKB + MKF+ MKT + Mfil
    Mbrain <- MBrB + MBrF+ MBrT 
    Mliver <- MLB + MLF+ MLT + MBile
    Mdust <- MLuAFdust_slow + MLuAFdust_fast
    
    Ckidney <- (MKB + MKF+ MKT + Mfil)/VKi_tot
    Cliver <- (MLB + MLF+ MLT + MBile)/(VLi_tot)
    Cintestine <-  (MINB + MINF+ MINT+MINL)/VIn_tot
    Cstomach <-  (MSTB + MSTF+ MSTT + MSTL)/VSt_tot
    Cmuscle <-  (MMB + MMF+ MMT)/VMu_tot
    Cadipose <-  (MAB + MAF+ MAT)/VAd_tot
    Clungs <-  (MLuB + MLuF+ MLuT + MLuAF+Mdust)/VLu_tot
    Clungtissue <- (MLuB + MLuF+ MLuT+Mdust*fdust_retained)/(VLuB_cap+VLuF+VLuT)
    CUpperair <- MUA/VUA
    CalveolarLF <- (MLuAF+MLuAFdust_slow + MLuAFdust_fast)/VLuAF
    VBALF_Gustaffson <- 0.005 #L
    CBALF <-  (MLuAF+Mdust*(1-fdust_retained))/ VBALF_Gustaffson# This is for the Gustaffson et al. (2022) study which used 
    Cspleen <-  (MSPB + MSPF+ MSPT)/VSp_tot
    Cheart <-  (MHB + MHF+ MHT)/VHt_tot
    Cbrain <-  (MBrB + MBrF+ MBrT)/VBr_tot
    Cgonads <-  (MGoB + MGoF+ MGoT)/VGo_tot
    Cskin <-  (MSKB + MSKF+ MSKT)/VSk_tot
    Cbones <-  (MBoB + MBoF+ MBoT)/VBo_tot
    Crest <-  (MRB + MRF+ MRT)/VRe_tot
    Ccarcass <- (MMB+MMF+MMT+MAB+MAF+MAT+MRB+MRF+MRT+MBoB+MBoF+MBoT+MSKB+MSKF+MSKT)/(VMu_tot+VAd_tot+VRe_tot+VBo_tot+VSk_tot)
    
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
            'dCalbBoFf' = dCalbBoFf, 
            'dCFabpPTCf' = dCFabpPTCf, 'dCFabpLHCf' = dCFabpLHCf, 'dCFabpDTCf' = dCFabpDTCf,
            'dCFabpCDCf' = dCFabpCDCf, 
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
            'dMPTCb' = dMPTCb, 'dMLHCb' = dMLHCb, 'dMDTCb' = dMDTCb,
            'dMCDCb' = dMCDCb, 
            'dMLTb' = dMLTb, 'dMLuAFb'=dMLuAFb,
            
            
            'dMArtf'=dMArtf, 'dMVenf'=dMVenf, 'dMKBf'=dMKBf, 
            'dMKFf'=dMKFf, 
            'dMPTCf' = dMPTCf, 'dMLHCf' = dMLHCf,  
            'dMDTCf' = dMDTCf, 'dMCDCf' = dMCDCf,
            'dMPT' = dMPT,  'dMLH' = dMLH,
            'dMDT' =  dMDT, 'dMCD' = dMCD,
            
            'dMBladder' = dMBladder, 'dMLBf'=dMLBf, 
            'dMLFf'=dMLFf, 'dMLTf'=dMLTf, 'dMBile'=dMBile,
            
            'dMSTBf'=dMSTBf, 'dMSTFf'=dMSTFf, 'dMSTTf'=dMSTTf, 'dMSTL'=dMSTL,
            'dMINBf'=dMINBf, 'dMINFf'=dMINFf, 'dMINTf'=dMINTf,'dMINL'=dMINL,
            
            'dMMBf'=dMMBf, 'dMMFf'=dMMFf, 'dMMTf'=dMMTf,
            'dMABf'=dMABf, 'dMAFf'=dMAFf, 'dMATf'=dMATf, 
            'dMRBf'=dMRBf, 'dMRFf'=dMRFf,'dMRTf'=dMRTf,
            'dMUA' = dMUA, 'dMLuBf'=dMLuBf, 'dMLuFf'=dMLuFf,'dMLuTf'=dMLuTf,
            'dMLuAFf' = dMLuAFf, "dMLuAFdust_fast" = dMLuAFdust_fast, "dMLuAFdust_slow" = dMLuAFdust_slow,
            
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
    'CKTf'=CKTf, 'CKTb'=CKTb,  'CBladder' = CBladder,'CPTCf' = CPTCf,
    'CPTCb' = CPTCb, 'CLHCf' = CLHCf, 'CLHCb' = CLHCb,
    'CDTCf' = CDTCf, 'CDTCb' = CDTCb, 'CCDCf' = CCDCf,
    'CCDCb' = CCDCb, 'CPT' = CPT, 'CLH' = CLH, 'CDT' = CDT, 'CCD' = CCD,
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
    'Cblood_background'=Cblood_background,
    
    'Ckidney'=Ckidney, 'Cliver'=Cliver,
    'Cstomach'=Cstomach, 'Cintestine'=Cintestine, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose,
    'Clungs'=Clungs, 'Clungtissue'=Clungtissue, 'Crest'=Crest, 'Ccarcass'=Ccarcass,
    'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain, 
    'Cgonads'=Cgonads, 'Cskin'=Cskin, 'Cbones'=Cbones, 
    
    "CBile" = CBile, 'CalveolarLF' = CalveolarLF, "CUpperair" = CUpperair,
    "Cfeces" = Cfeces,  "Curine" = Curine, "CBALF" = CBALF
    
    
    )
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    CalbVenf<- CalbB_init; MVenf<- 0#0.13 * VVen;
    MVenb<- 0;
    CalbArtf<- CalbB_init; MArtf<- 0 #0.13 * VArt;
    MArtb<- 0; 
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
 
    CFabpPTCf <- CFabpKT_init; CFabpLHCf <- CFabpKT_init;
    CFabpDTCf <- CFabpKT_init; CFabpCDCf <- CFabpKT_init; 
    
    MPTCf <- 0;MLHCf <- 0;
    MDTCf <- 0; MCDCf <- 0;
    MPT <- 0;  MLH<- 0;
    MDT <- 0; MCD <- 0;
    
    CFabpLTf<- CFabpLT_init; CalbLuAFf<- CalbLuAF_init;
    
    MLTf<- 0; MLTb<- 0;  MPTCb <- 0; MLHCb <- 0; MDTCb <- 0
    MCDCb <- 0; MBile <-0; MSTTf <- 0;  MINTf <- 0; MLuTf <- 0; 
    MUA <- 0; MLuAFf <- 0; MLuAFb<- 0; MLuAFdust_fast <- 0; MLuAFdust_slow <- 0;
    MSPTf <- 0; MHTf <- 0;  MBrTf <- 0;
    MGoTf <- 0;  MSKTf <- 0; MBoTf <- 0; MMTf <- 0; MATf <- 0; MRTf <- 0;
    
    MBladder <- 0; Murine <-0;MSTL <-0;  MINL <-0;
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
             'CFabpPTCf' = CFabpPTCf, 'CFabpLHCf' = CFabpLHCf, 'CFabpDTCf' = CFabpDTCf,
             'CFabpCDCf' = CFabpCDCf, 
             'CFabpLTf' = CFabpLTf, 'CalbLuAFf' = CalbLuAFf,
             
             
             'MVenb' = MVenb, 'MArtb' = MArtb, 'MKBb' = MKBb, 
             'MLBb' = MLBb,'MSTBb' = MSTBb,'MINBb' = MINBb,'MMBb' = MMBb,
             'MABb' = MABb, 'MRBb' = MRBb,'MLuBb' = MLuBb, 'MSPBb' = MSPBb, 
             'MHBb' = MHBb,  'MBrBb' = MBrBb,  'MGoBb' = MGoBb, 
             'MSKBb' = MSKBb, 'MBoBb' = MBoBb,'MKFb' = MKFb, 'MLFb' = MLFb, 
             'MSTFb' = MSTFb,  'MINFb' = MINFb, 'MMFb' = MMFb, 
             'MAFb' = MAFb, 'MRFb' = MRFb, 'MLuFb' = MLuFb, 
             'MSPFb' = MSPFb,  'MHFb' = MHFb, 'MBrFb' = MBrFb, 
             'MGoFb' = MGoFb, 'MSKFb' = MSKFb, 'MBoFb' = MBoFb,
             'MPTCb' = MPTCb, 'MLHCb' = MLHCb, 'MDTCb' = MDTCb,
             'MCDCb' = MCDCb, 'MLTb' = MLTb, 'MLuAFb'=MLuAFb,
             
             
             'MArtf'=MArtf, 'MVenf'=MVenf, 'MKBf'=MKBf, 
             'MKFf'=MKFf,  'MPTCf' = MPTCf, 'MLHCf' = MLHCf,  
             'MDTCf' = MDTCf, 'MCDCf' = MCDCf,
             'MPT' = MPT,  'MLH' = MLH,
             'MDT' =  MDT, 'MCD' = MCD,
             
             'MBladder' = MBladder, 'MLBf'=MLBf, 
             'MLFf'=MLFf, 'MLTf'=MLTf, 'MBile'=MBile,
             
             'MSTBf'=MSTBf, 'MSTFf'=MSTFf, 'MSTTf'=MSTTf, 'MSTL'=MSTL,
             'MINBf'=MINBf, 'MINFf'=MINFf, 'MINTf'=MINTf,'MINL'=MINL,
             
             'MMBf'=MMBf, 'MMFf'=MMFf, 'MMTf'=MMTf,
             'MABf'=MABf, 'MAFf'=MAFf, 'MATf'=MATf, 
             'MRBf'=MRBf, 'MRFf'=MRFf,'MRTf'=MRTf,
             'MUA' = MUA, 'MLuBf'=MLuBf, 'MLuFf'=MLuFf,'MLuTf'=MLuTf,
             'MLuAFf' = MLuAFf, 'MLuAFdust_fast' = MLuAFdust_fast,
             'MLuAFdust_slow' = MLuAFdust_slow,
             
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
      }else if (admin.type == "inh") {
        
        events <- list(data = rbind(
          data.frame(
            var    = c("MLuAFdust_fast"),
            time   = admin.time,
            value  = admin.dose * depfr_AF * f_fast,
            method = c("add")
          ),
          data.frame(
            var    = c("MLuAFdust_slow"),
            time   = admin.time,
            value  = admin.dose * depfr_AF * (1 - f_fast),
            method = c("add")
          )
        ))
        
      } else if (admin.type == "nasal") {
        
        events <- list(data = rbind(
          data.frame(
            var    = c("MUA"),
            time   = admin.time,
            value  = admin.dose * depfr_head,
            method = c("add")
          ),
          data.frame(
            var    = c("MLuAFf"),
            time   = admin.time,
            value  = admin.dose * depfr_AF ,
            method = c("add")
          )
        ))
      }else if (admin.type == "dermal"){
        events <- list(data = rbind(data.frame(var = c("MSKTf"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
        
      }
    }
    return(events)
  })
}


create_all_fixed_params <- function(){
  params <- list()
  
  # Set up simulations for the arbritary values
  BW <- 82  # body weight (kg)
  admin.dose <- 3.96 #ug PFOA (Table 1)
  admin.time <- 0 # time when doses are administered, in hours
  admin.type <- "oral"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time,
                     "admin.type" = admin.type,
                     "sex" = sex)
  params <- create_fixed_params(user_input)
 
  
  return(params)
  
}

obj.func <- function(x, dataset, fixed_params){
  optimized_vector <- x
  print(optimized_vector)
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)
 
#=====================================================================================================   
#Abraham 2024 - oral exposure
#===================================================================================================== 
  
  BW <- 82  # body weight (kg)
  sex <- "M"
  variable_params <- create_variable_params(BW, sex, estimated_params, fixed_params)
  params <- c(fixed_params, variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  t_hours <- c(
    seq(0, 11, by = 0.25),         
    seq(11 + 6, 10800, by = 6)  
  )
  
  sample_time <- t_hours

  plasma_times_h <- sort(unique(dataset$df1$Time_h))
  urine_times_h  <- sort(unique(dataset$df2$Time_h))
  feces_times_h  <- sort(unique(dataset$df3$Time_h))
  
  #experimental_times <- c(plasma_times_h, urine_times_h, feces_times_h)
  experimental_times <- c(plasma_times_h, urine_times_h)
  
  simulation_time = sort(unique(c(sample_time, experimental_times)))

  solution <- data.frame(deSolve::ode(times = simulation_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-07, atol = 1e-07))

  #======================================df1=========================================================

  exp_data <- dataset$df1 # retrieve data of Abraham et al. 2024
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("Cplasma")

  preds_Abraham_2024_plasma <- list()
  compartment <- "Plasma" #unique(exp_data$Tissue)[i]
  exp_time <- exp_data[, 2]
  
  preds_Abraham_2024_plasma<- solution[solution$time %in% exp_time, "Cplasma"]
  

  obs_Abraham_2024_plasma <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
  

  score[1] <- AAFE(predictions = preds_Abraham_2024_plasma, observations = obs_Abraham_2024_plasma)

  #======================================df2=========================================================

  exp_data <- dataset$df2 # retrieve data of Abraham et al. 2024
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Murine")

  preds_Abraham_2024_urine <- list()
  compartment <- "Urine" #unique(exp_data$Tissue)[i]
  exp_time <- exp_data[, 2]

  preds_Abraham_2024_urine <- solution[solution$time %in% exp_time, "Murine"]

  obs_Abraham_2024_urine <- list(exp_data[exp_data$Tissue == "Urine", "mass"])
  obs_Abraham_2024_urine <- cumsum(unlist(obs_Abraham_2024_urine))
  

  score[2] <- AAFE(predictions = preds_Abraham_2024_urine, observations = obs_Abraham_2024_urine)
  
  #======================================df3=========================================================
  
  
  exp_data <- dataset$df3 # retrieve data of Abraham et al. 2024
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mfeces")

  preds_Abraham_2024_feces <- list()
  compartment <- "Feces" #unique(exp_data$Tissue)[i]
  exp_time <- exp_data[, 2]

  preds_Abraham_2024_feces <- solution[solution$time %in% exp_time, "Mfeces"]*1000

  obs_Abraham_2024_feces <- list(exp_data[exp_data$Tissue == "Feces", "mass"])
  obs_Abraham_2024_feces <- cumsum(unlist(obs_Abraham_2024_feces))

  score[3] <- AAFE(predictions = preds_Abraham_2024_feces, observations = obs_Abraham_2024_feces)


  
  
# #=====================================================================================================   
# #Abraham 2022 - 2024 - dermal exposure
# #=====================================================================================================
#   
#   BW <- 80  # body weight (kg)
#   sex <- "M"
#   variable_params <- create_variable_params(BW, sex, estimated_params, fixed_params[[2]])
#   params <- c(fixed_params[[2]], variable_params)
#   inits <- create.inits(params)
#   events <- create.events(params)
#   sample_time=seq(0,6900,15)
#   solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
#                                       y = inits, parms = params,
#                                       events = events,
#                                       method="lsodes",rtol = 1e-07, atol = 1e-07))
#   
#   #======================================df5=========================================================
#   
#   exp_data <- dataset$df5 # retrieve data of Abraham et al. 2022
#   colnames(exp_data)[c(2,3)] <- c("time", "concentration")
#   column_names <- c("Cplasma")
#   
#   preds_Abraham_2022_plasma_dermal <- list()
#   # loop over compartments with available data
#   for (i in 1:length(unique(exp_data$Tissue))) {
#     compartment <- unique(exp_data$Tissue)[i]
#     #Retrieve time points at which measurements are available for compartment i
#     exp_time <- exp_data[exp_data$Tissue == compartment, 2]
#     
#     preds_Abraham_2022_plasma_dermal[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
#   }
#   
#   preds_Abraham_2022_plasma_dermal[[1]] <- (preds_Abraham_2022_plasma_dermal[[1]] *1000) + 1400
#   
#   obs_Abraham_2022_plasma_dermal <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
#   
#   score[5] <- AAFE(predictions = preds_Abraham_2022_plasma_dermal, observations = obs_Abraham_2022_plasma_dermal)
#   
#   #======================================df6=========================================================
#   
#   BW <- 80  # body weight (kg)
#   sex <- "M"
#   variable_params <- create_variable_params(BW, sex, estimated_params, fixed_params[[3]])
#   params <- c(fixed_params[[3]], variable_params)
#   inits <- create.inits(params)
#   events <- create.events(params)
#   sample_time=seq(0,9806,1)
#   solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
#                                       y = inits, parms = params,
#                                       events = events,
#                                       method="lsodes",rtol = 1e-07, atol = 1e-07))
#   
#   
#   exp_data <- dataset$df6 # retrieve data of Abraham et al. 2024
#   colnames(exp_data)[c(2,3)] <- c("time", "concentration")
#   column_names <- c("Cplasma")
#   
#   preds_Abraham_2024_plasma_dermal <- list()
#   # loop over compartments with available data
#   for (i in 1:length(unique(exp_data$Tissue))) {
#     compartment <- unique(exp_data$Tissue)[i]
#     #Retrieve time points at which measurements are available for compartment i
#     exp_time <- exp_data[exp_data$Tissue == compartment, 2]
#     
#     preds_Abraham_2024_plasma_dermal[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
#   }
#   
#   preds_Abraham_2024_plasma_dermal[[1]] <- (preds_Abraham_2024_plasma_dermal[[1]]) + 0.13
#   
#   obs_Abraham_2024_plasma_dermal <- list(exp_data[exp_data$Tissue == "Plasma", "concentration"])
#   
#  score[6] <- AAFE(predictions = preds_Abraham_2024_plasma_dermal, observations = obs_Abraham_2024_plasma_dermal)
#   
  
 final_score <- mean(score, na.rm = TRUE)
 return(final_score)
 
}

  # ########################################################################################################################################
  
  
  setwd('/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat')
  
  MW <- 414.07 #g/mol
  source("Goodness-of-fit-metrics.R")
  
  setwd('/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/PFOA extrapolation to humans/Data')
  
  Abraham_2024_plasma <- openxlsx::read.xlsx("Abraham et al.2024_oral_plasma.xlsx") #ug/L
  Abraham_2024_urine <- openxlsx::read.xlsx("Abraham et al.2024_oral_urine_cut.xlsx") #ug
  Abraham_2024_feces <- openxlsx::read.xlsx("Abraham et al.2024_oral_feces_cut.xlsx") #ng
  
  Abraham_2024_plasma_6h <- openxlsx::read.xlsx("Abraham et al.2024_oral_plasma_6h.xlsx") #ug/L
  Abraham_2022_plasma_dermal<- openxlsx::read.xlsx("Abraham et al.2022_dermal_plasma.xlsx") #ng/L, substract background
  Abraham_2024_plasma_dermal<- openxlsx::read.xlsx("Abraham et al.2024_dermal_plasma.xlsx") #ug/L
  
  
  dataset <- list("df1" = Abraham_2024_plasma, "df2" = Abraham_2024_urine, "df3" = Abraham_2024_feces)
  
  setwd('/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/PFOA extrapolation to humans/Trials')
  
  #Initialise optimiser to NULL for better error handling later
  opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA"   "NLOPT_LN_SBPLX"
                "xtol_rel" = 1e-07,
                "ftol_rel" = 1e-07,
                "ftol_abs" = 0.0,
                "xtol_abs" = 0.0, 
                "maxeval" = 500, 
                "print_level" = 1)
  
  # Create initial conditions (zero initialisation)
  
  N_pars <- 2 # Number of parameters to be fitted
  #fit <- log(c(0.17, 10.9, 1.36e-2, 0.586, 0.3, 1, 1, 30, 0.16, 0.04, 3.1e-6, 4.3e4, 1.6e-5))
  #fit <- rep(log(c(1)), 6)
  
  fit <-  log(c(1e-2, 1, 1e-03))
  
  # lb = rep(log(c(1e-12)), 2)
  # ub = rep(log(c(1e5)), 2)
 
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
  save.image("Final_model_AAFE_humans.RData")
  
  # ########################################################################################################################################
  
  BW <- 82  # body weight (kg)
  sex <- "M"
  variable_params <- create_variable_params(BW, sex, estimated_params, fixed_params)
  params <- c(fixed_params, variable_params)
  inits <- create.inits(params)
  events <- create.events(params)
  t_hours <- c(
    seq(0, 11, by = 0.25),         
    seq(11 + 6, 10800, by = 6)  
  )
  
  sample_time <- t_hours
  
  plasma_times_h <- sort(unique(dataset$df1$Time_h))
  urine_times_h  <- sort(unique(dataset$df2$Time_h))
  feces_times_h  <- sort(unique(dataset$df3$Time_h))
  
  experimental_times <- c(plasma_times_h, urine_times_h, feces_times_h)
  
  simulation_time = sort(unique(c(sample_time, experimental_times)))
  
  solution <- data.frame(deSolve::ode(times = simulation_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-07, atol = 1e-07))
  
  preds_Abraham_2024_plasma <- solution[solution$time %in% plasma_times_h, c("time", "Cplasma")]
  preds_Abraham_2024_urine <- solution[solution$time %in% urine_times_h, c("time", "Murine")]
  preds_Abraham_2024_feces <- solution[solution$time %in% feces_times_h, c("time", "Mfeces")]
  
  preds_Abraham_2024_feces$Mfeces <- preds_Abraham_2024_feces$Mfeces * 1000
  

  # ######################################################################################
  #Plot the predictions against the observations
  library(ggplot2) 
  
  # Function that creates a plot given a compartment name and the respective predictions and observations
  create.plots <- function(predictions, observations, compartment){  
    #Colours of observations and predictions
    cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
    x_rng <- range(observations$Time, na.rm = TRUE) 
    
    ggplot(data = predictions)+
      geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                            color = '"predictions"'),  size=1.5,alpha = 0.7) +
      geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                               color='"Observations"'), size=4)+
      coord_cartesian(xlim = x_rng) +
      labs(title = rlang::expr(!!compartment), 
           y = expression("PFOA concentration (ug/L)" ),
           x = "Time (h)")+
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
  
  
  # Convert data from long to wide format using reshape
  experiment1 <- reshape(Abraham_2024_plasma[c("Tissue" ,"Time_h",
                                               "Concentration_ug/L")],
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
  colnames(experiment1) <- c("Time",unique(Abraham_2024_plasma$Tissue))
  
  experiment2 <- reshape(Abraham_2024_urine[c("Tissue" ,"Time_h",
                                               "Mass_(ug)")],
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
  colnames(experiment2) <- c("Time",unique(Abraham_2024_urine$Tissue))
  experiment2$Urine <- cumsum(experiment2$Urine)
  
  experiment3 <- reshape(Abraham_2024_feces[c("Tissue" ,"Time_h",
                                              "Mass_(ng)")],
                         idvar = "Time_h", timevar = "Tissue", direction = "wide")
  colnames(experiment3) <- c("Time",unique(Abraham_2024_feces$Tissue))
  experiment3$Feces <- cumsum(experiment3$Feces)
  
  # Put the experiments in a list
  experiments <- list(experiment1 = experiment1,
                      experiment2 = experiment2,
                      experiment3 = experiment3)

  
  # Rename predictions so that they share the same name as the names of the experimental data 
  colnames(preds_Abraham_2024_plasma) <- c( "Time", "Plasma")
  colnames(preds_Abraham_2024_urine) <- c( "Time", "Urine")
  colnames(preds_Abraham_2024_feces) <- c( "Time", "Feces")
  
  # Create a list containing the corresponding predictions
  simulations <- list(predictions1 = preds_Abraham_2024_plasma,
                      predictions2 = preds_Abraham_2024_urine,
                      predictions3 = preds_Abraham_2024_feces)
  
  
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
  