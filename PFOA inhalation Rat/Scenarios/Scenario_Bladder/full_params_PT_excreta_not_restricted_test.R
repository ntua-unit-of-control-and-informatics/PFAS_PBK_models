setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Scenario_Bladder/Training/AAFE/full_params_PT_excreta_not_restricted_test")
load("full_params_PT_excreta_not_restricted.RData")


create.params <- function(user.input){
  with(as.list(user.input),{
    # BW in kg
    # Cheng and Ng 2017 Table S1
    # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
    #assuming the density of tissue is 1 g/mL.
    # Estimated parameters
    if(sex == "M"){
      RAFOatp_k <- estimated_params[1]
      RAFOat1 <- estimated_params[3]
      RAFbile_transp <- estimated_params[5]
    }else{
      RAFOatp_k <- 1#estimated_params[2]
      RAFOat1 <-1#estimated_params[4]
      RAFbile_transp <- 0#estimated_params[6]
    }
    RAFOatp_l <- estimated_params[7]
    RAFUrat <- RAFOatp_k
    RAFOat3 <- RAFOat1
    RAFOatp2_l <- RAFOatp_l
    RAFOatp_lu_ap <- estimated_params[8]
    RAFOatp_lu_bas <- RAFOatp_lu_ap
    RAFNtcp <- RAFOatp_l
    RAFOatp2_Int <- estimated_params[9]
    
    RAF_papp <- estimated_params[10]
    
    f_fabp_avail <- estimated_params[11]
    f_alb_avail <- estimated_params[12]
    
    koff_alb <-100
    koff_fabp <-  koff_alb
    koff_a2u <- koff_alb
    
    VmK_api <- 0
    VmK_baso <- 0
    KmK_baso <- 1e20
    KmK_api <-   1e20
    KLfabp <- (1.2e5+4e4+1.9e4)  #[L/mol]*1e-3 , value from Cheng et al. (2017)
    Ka <- estimated_params[13] # 5.8e05 from Rue et al. (2024)#mol/L
    
    
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
    PVKF <- 0.2 #pkSim
    VKF <- PVKF * PVK * BW #kidney interstitial fluid volume kg=L
    VKT <- VK - VKF #kidney tissue volume kg=L
    
    #volumes of filtrate compartments, https://doi.org/10.1152/ajprenal.00219.2018, [μm^3]*1e-15 --> L
    
    VPT <- pi*((22.9/2)^2*9886*26980+(23.1/2)^2*11077*11020)*2*1e-15 # Proximal tubule (short- and long-looped)
    VTDL <- pi*((15/2)^2*1500*26980+(15/2)^2*6200*11020)*2*1e-15 #Thin descending limb (short- and long-looped)
    VThinAL <-  pi*((15/2)^2*4700*11020)*2*1e-15 #Thin ascending limb (long-looped)
    VThickAL <- pi*((25.4/2)^2*1450*38000+(29/2)^2*2100*38000)*2*1e-15 #Thick ascending limb (Cortical and Medullary)
    VDT <- pi*((39/2)*1452*26980+(43/2)^2*1650*11020)*2*1e-15 #Distal tubule (superficial+deep)
    Vduct <- pi*((24/2)^2*2900*6000+(24/2)^2*2100*6000)*2*1e-15 #collecting duct (Cortical+Outer)
    VFil_rest <- VTDL+VThinAL+VThickAL+VDT+Vduct
    #VFil <- 0.25/1000 #renal filtrate volume in L,  from Cheng et al., 2017 (from Arthur, 1986; Bonvalet, 1981)
    VBladder <- 0.001
    
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
    
    PAK <- 231.7e-4/0.23 #m^2/kg
    AK <- PAK * BW #kidney surface area (m^2)
    PAKG <- 68.90e-4 #Kirkman and Stowell, 1942, https://doi.org/10.1002/ar.1090820310
    AKG <- PAKG * VK * 1e3 #the surface area of glomerular capillary (m^2)
    PAL <- 1136e-4/0.23 #m^2/kg
    AL <- PAL * BW #liver surface area (m^2)
    
    PAST <- 33.77e-4/0.23 #m^2/kg
    AST <- PAST * BW #stomach surface area (m^2)
    PASTL<- 33.77e-4/0.23 #m^2/kg
    ASTL<- PASTL * VSTL #stomach lumen surface area (m^2)
    
    PAIN <- 74.82e-4/0.23 #m^2/kg
    AIN <- PAIN * BW #intestine surface area (m^2)
    
    #Calculations of rat intestinal lumen surface area based on Kothari et al. (2020),https://doi.org/10.1002/btm2.10146
    # Lengths
    L_duodenum <- 9.6# cm
    L_jejunum <- 26 #cm
    L_ileum <- 34.4# cm
    # Inner diameters
    d_duodenum <- 2.21 #cm
    d_jejunum <- 2.56 #cm
    d_ileum <- 3.36 #cm
    
    SA <- 2*pi*(d_duodenum/2)*L_duodenum + 2*pi*(d_jejunum/2)*L_jejunum + 2*pi*(d_ileum/2)*L_ileum #cm^2
    n <- 5 #enlargement factor of apical membrane of proximal tubule or enterocytes
    PAINL <- n * SA * 1e-4/0.195 #m^2/kg, scaled to reference body weight
    AINL <- PAINL*BW #m^2
    
    PAM <- 3042e-4/0.23 #m2/kg
    AM <- PAM * BW #muscle surface area (m^2)
    PAA <- 95.93e-4/0.23 #m2/kg
    AA <- PAA * BW #adipose surface area (m^2)
    
    PAR <- 100e-4#m2/g tissue, assumption
    AR <- PAR * VR *1e03#surface area of rest of body (m^2)
    
    PLu <- 600.5e-4/0.23 #m2/kg 
    ALu <- PLu * BW #lung surface area (m^2)
    PSP <- 162.3e-4/0.23 #m2/kg
    ASP <- PSP * BW #spleen surface area (m^2)
    PH <-  201.1e-4/0.23 #m2/kg 
    AH <- PH * BW #heart surface area (m^2)
    PBr <- 60.34e-4/0.23 #m2/kg
    ABr <- PBr * BW #brain surface area (m^2)
    PGo <- 335.8e-4/0.23#m2/kg
    AGo <- PGo * BW #gonads surface area (m^2)
    PSK <- 729.1e-4/0.23#m2/kg
    ASK <- PSK * BW #skin surface area (m^2)
    PBo <- 621.4e-4/0.23#m2/kg
    ABo <- PBo * BW #skin surface area (m^2)
    
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
    
    #Estimation of lymph flow rates:
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
    Qfeces <- (8.18/0.21)*BW #g/kg BW, based on Cui et al.(2010)
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
    QGE<- 0.54/BW^0.25 #gastric emptying time (1/(h*BW^0.25)); from Yang, 2013
    
    #flows of filtrate compartments, https://doi.org/10.1152/ajprenal.00219.2018, nL/min ---> L/h 
    
    QPT <- (40.9*26980+39.7*11020)*2*1e-9*60 #proximal tubule flow (short +long)
    
    Q_scaling_factor = QGFR/QPT
    #an den ta xreiastoume pouthena ta svinoume auta
    
    QTDL_reference <- (26.4*26980+26.3*11020)*2*1e-9*60 #Thin descending limb (short +long)
    QTDL <- QTDL_reference*Q_scaling_factor
    # QTAL <- 26.3*11020*2*1e-9*60 #Thin descending limb (short +long)
    # QThinAL <- 8*11020*2*1e-9*60 #Thin ascending limb (long-looped)
    # QThickAL <- (8.1*38000+8*38000)*2*1e-9*60 #Thick ascending limb (Cortical and Medullary)
    # QDT <- 6*(26980+11020)*2*1e-9*60 #Distal tubule (short +long)
    # Qduct <- Qurine
    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    #Surface areas Interstitial - Intracellular PKSim (m^2)
    BW_ref <- 0.23
    AcK= 437.16*BW/BW_ref
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
    AcR= (AcK+AcL+AcST+AcIN+AcM+AcA+AcLu+AcSP+AcH+AcBr+AcGo+AcSK+AcBo)/13
    # Surface areas of the different subcompartments of kidney filtrate
    APT <-  2*pi*((22.9/2)*9886*26980+(23.1/2)*11077*11020)*2*n*1e-12 # Proximal tubule (short- and long-looped)
    ATDL <- 2*pi*((15/2)*1500*26980+(15/2)*6200*11020)*2*n*1e-12 #Thin descending limb (short- and long-looped)
    AThinAL <-  2*pi*((15/2)*4700*11020)*2*n*1e-12 #Thin ascending limb (long-looped)
    AThickAL <- 2*pi*((25.4/2)*1450*38000+(29/2)*2100*38000)*2*n*1e-12 #Thick ascending limb (Cortical and Medullary)
    ADT <- 2*pi*((39/2)*1452*26980+(43/2)*1650*11020)*2*n*1e-12 #Distal tubule (superficial+deep)
    Aduct <- 2*pi*((24/2)*2900*6000+(24/2)*2100*6000)*2*n*1e-12 #collecting duct (Cortical+Outer)
    AFil_rest <- ATDL+AThinAL+AThickAL+ADT+Aduct
    
    #Alveolar cells surface area (Type I and II), m^2
    AcALF = ((78.8*2*5320*1e-6) + (125*2*123*1e-6))*BW/0.29  #Stone et al., 1992, BW_ref = 0.29, values for each lung , https://doi.org/10.1165/ajrcmb/6.2.235
    
    #canalicular surface area, m^2
    rat_hep_surf_area = 22.95 * 1e2 # 22.95*1e6 cm2 --> m2,  https://doi.org/10.1074/jbc.271.12.6702
    AcLBilec = 0.01 * 22.95 * 1e2 # m2 , canalicular membrane 1% of the surface area of the hepatocyte,https://www.ncbi.nlm.nih.gov/books/NBK470209/
    
    # Following the calculations  of Lin et al. (2023) for Caco-2 cells
    ClINFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
    muscle_protein <- 158.45 #mg/g muscle (protein data from Cheek et al.,1971 
    #and muscle mass from Caster et al.,1956)
    intestine_protein <- muscle_protein
    intestine_protein_total <- intestine_protein*(1000*VINT)
    muscle_protein <- 158.45 #mg/g muscle (protein data from Cheek et al.,1971 
    #and muscle mass from Caster et al.,1956)
    Awell = 9 #cm^2 (for a 35 mm culture dish)
    Swell = 1.12 #cm^2
    well_protein = 0.346 #mg protein
    P = (well_protein * Awell)/Swell #mg protein/well
    Papp = RAF_papp*(ClINFT_unscaled*60*1e-03*P)/(Awell *2) #cm/h
    P_passive = ( (Papp/100) * AINL)*1000 #L/h
    
    
    #passive diffusion rates
    
    kKFKT = ((Papp/100) * AcK)*1000 #m^3/h * 1000 --> L/h
    kLFLT = ((Papp/100) * AcL)*1000 #m^3/h * 1000 --> L/h
    kLTLbile = ((Papp/100) * AcLBilec)*1000 #m^3/h * 1000 --> L/h
    kMFMT = ((Papp/100) * AcM)*1000 #m^3/h * 1000 --> L/h
    kSTFSTT = ((Papp/100) * AcST)*1000 #m^3/h * 1000 --> L/h 
    kINFINT = ((Papp/100) * AcIN)*1000 #m^3/h * 1000 --> L/h 
    kAFAT = ((Papp/100) * AcA)*1000 #m^3/h * 1000 --> L/h 
    kLuTLuF = ((Papp/100) * AcLu)*1000 #m^3/h * 1000 --> L/h
    kLuTLuAF = ((Papp/100) * AcALF)*1000 #m^3/h * 1000 --> L/h
    kSPFSPT = ((Papp/100) * AcSP)*1000 #m^3/h * 1000 --> L/h 
    kHFHT = ((Papp/100) * AcH)*1000 #m^3/h * 1000 --> L/h 
    kBrFBrT = ((Papp/100) * AcBr)*1000 #m^3/h * 1000 --> L/h 
    kGoFGoT = ((Papp/100) * AcGo)*1000 #m^3/h * 1000 --> L/h 
    kSKFSKT = ((Papp/100) * AcSK)*1000 #m^3/h * 1000 --> L/h
    kBoFBoT = ((Papp/100) * AcBo)*1000 #m^3/h * 1000 --> L/h
    kRFRT = ((Papp/100) * AcR)*1000 #m^3/h*1000 --> L/h 
    kPTKT <- ((Papp/100) * APT) *1000 
    kFilKT <- ((Papp/100) * AFil_rest)*1000 #m^3/h *1000 ---> L/h
    
    
    #For all CMTs
    MW = 414.07 #g/mol, PFOA molecular weight
    Acell = 4000 #um^2/cell
    #Kidney
    kidney_protein_per_rat <- 1000*(0.218+0.225+0.212)/3#mg of protein per rat  (Addis 1936)
    rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
    rat_kidney_weight_addis <- rat_weight_addis*0.0073 # kidney fraction to BW, Brown (1997)
    kidney_protein_per_gram <- kidney_protein_per_rat/rat_kidney_weight_addis #mg of protein/g kidney
    
    kidney_cells = 1.47e07 #cells/g https://doi.org/10.1038/s41598-024-53270-2
    kidney_cells_total <- kidney_cells* (1000*VKT)
    kidney_protein_total <- kidney_protein_per_gram* (1000*VKT) #mg
    
    #Oatp kidney
    VmK_Oatp_in_vitro <- 9.3 #nmol/mg protein/min (Weaver et al. 2010)
    VmK_Oatp_scaled <- 60*VmK_Oatp_in_vitro*MW*kidney_protein_total/1000  #physiologically scaled to in vivo, ug/h
    VmK_Oatp <- VmK_Oatp_scaled*RAFOatp_k #in vivo value, in  ug/h
    KmK_Oatp=  126.4*MW# [umol/L] * g/mol  --> ug/L, from Weaver et al. (2010)
    
    #oat1 kidney
    VmK_Oat1_in_vitro= 2.6 #nmol/mg protein/min (Weaver et al. 2010)
    VmK_Oat1_scaled = 60*VmK_Oat1_in_vitro*MW*kidney_protein_total/1000 #physiologically scaled to in vivo, ug/h
    VmK_Oat1= VmK_Oat1_scaled*RAFOat1 #in vivo value, in   ug/h
    KmK_Oat1= 43.2 * MW #umol/L (Weaver et al. 2010) --> ug/L
    
    #oat3 kidney
    VmK_Oat3_in_vitro= 3.8 #nmol/mg protein/min  (Weaver et al. 2010)
    VmK_Oat3_scaled = 60*VmK_Oat3_in_vitro*MW*kidney_protein_total/1000 #physiologically scaled to in vivo, ug/h
    VmK_Oat3 = VmK_Oat3_scaled*RAFOat3 #in vivo value, in   ug/h
    KmK_Oat3= 65.7 * MW #umol/L (Weaver et al. 2010) --> ug/L
    
    #Urat1 kidney
    VmK_Urat_in_vitro= 1520e-3 #nmol/mg protein/min  (Lin et al. 2023)
    VmK_Urat_scaled = 60*VmK_Urat_in_vitro*MW*kidney_protein_total/1000 #physiologically scaled to in vivo, ug/h
    VmK_Urat = VmK_Urat_scaled*RAFUrat #in vivo value, in   ug/h
    KmK_Urat = 820.04 * MW #umol/L (Lin et al. 2023) --> ug/L
    
    #Liver
    liver_protein_per_rat <- 1000*(1.52+1.53+1.52)/3#mg of protein per rat  (Addis 1936)
    rat_weight_addis <- 200 #g, average rat weight in Addis, 1936
    rat_liver_weight_addis <- rat_weight_addis*0.0366 # liver fraction to BW, Brown (1997)
    liver_protein_per_gram <- liver_protein_per_rat/rat_liver_weight_addis #mg or protein/g liver
    liver_cells = 117*10^6 #hepatocytes per g of liver (Sohlenius-Sternbeck et al. 2006) (2e09 cells: https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=110895)
    
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
    
    #bile_transp_liver
    Vmbile_transp_in_vitro= mean(c(VmL_Oatp_in_vitro,VmL_Oatp2_in_vitro, VmL_Ntcp_in_vitro )) #nmol/mg protein/min  (mean value from other liver transporters)
    Vmbile_transp_scaled = 60*Vmbile_transp_in_vitro*MW*liver_protein_per_gram*(VL*1000)/1000   #physiologically scaled to in vivo, ug/h
    Vmbile_transp = Vmbile_transp_scaled*RAFbile_transp #in vivo value, in  ug/h
    Kmbile_transp = mean(c(KmL_Ntcp,KmL_Oatp2, KmL_Oatp )) #same as kidney
    
    
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
    
    #Stomach
    # For identifiability reasons we assume that absorption takes place only through the intestines
    kabST <- (kabs_st* ASTL)*1000 #L/h
    
    #Effective permeability (Peff, in mm/h) for blood (B), liver(L), kidney(K),
    #stomach(ST),intestine (IN), adipose(A), muscle(M), spleen (SP), heart (H), 
    #brain (Br), gonads (Go), rest of body(R)
    
    PeffK <- Papp*10 #mm/h
    PeffL <- Papp*10 #mm/h
    PeffST <- Papp*10 #mm/h
    PeffIN <- Papp*10 #mm/h
    PeffA <- Papp*10 #mm/h
    PeffM <- Papp*10 #mm/h
    PeffR <- Papp*10 #mm/h
    PeffLu <- Papp*10 #mm/h
    PeffSP <- Papp*10 #mm/h
    PeffH <- Papp*10 #mm/h
    PeffBr <- Papp*10 #mm/h
    PeffGo <- Papp*10 #mm/h
    PeffSK <- Papp*10 #mm/h
    PeffBo <- Papp*10 #mm/h
    
    #Albumin concentration in blood and interstitial fluid compartments(mol/m^3 = 1e-6* nmol/g)
    
    CalbB_init <- f_alb_avail*486*1e-06 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    
    
    CalbKF_init <- f_alb_avail*243*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbLF_init <- f_alb_avail*243*1e-6 # #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbSTF_init <- f_alb_avail*146*1e-6 # [umol/L]*1e-6 -->(mol/L), same as Gut (assumption)
    CalbINF_init <- f_alb_avail*146*1e-6 #[umol/L]*1e-6 -->(mol/L), same as Gut (assumption)
    CalbMF_init <- f_alb_avail*146*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbAF_init <- f_alb_avail*73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbRF_init <- f_alb_avail*73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbBoF_init <-f_alb_avail* 73*1e-6 #[umol/L]*1e-6 -->(mol/L), from Cheng et al. (2017)
    CalbLuF_init <- f_alb_avail*CalbINF_init #assumption 
    CalbLuAF_init <- f_alb_avail*10/100 * CalbB_init #based on Woods et al. 2015 statement https://doi.org/10.1016/j.jconrel.2015.05.269
    
    MW_albumin <- 66500#g/mol
    CalbSPF_init <- f_alb_avail*243e-6 #[umol/L]*1e-6 -->(mol/L), same as liver (assumption)
    CalbGoF_init <- f_alb_avail*41/MW_albumin #mg/mL-->  (mol/L) from https://doi.org/10.1210/endo-116-5-1983 --> 41 mg/mL, MW=65 kg/mol
    CalbHF_init <- f_alb_avail* 65/MW_albumin##mg/mL--> (mol/L) https://doi.org/10.1007/s12291-010-0042-x --> 6.5 g/100 g tissue, MW=65 kg/mol 
    CalbBrF_init <- f_alb_avail*8e-2/MW_albumin  ##mg/mL--> (mol/L) https://doi.org/10.1016/0014-4886(90)90158-O --> 0.08 g/L, MW=65 kg/mol 
    CalbSKF_init <- f_alb_avail*21/MW_albumin ##mg/mL-->  (mol/L) https://doi.org/10.1111/j.1748-1716.1973.tb05464.x -->Table 2: 2.1 g/100 mL
    
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
    
    
    CalbKB_init <- CalbKF_init*(1/IPR_K) 
    CalbLB_init <- CalbLF_init*(1/IPR_L) 
    CalbSTB_init <- CalbLF_init*(1/IPR_ST)
    CalbINB_init <- CalbINF_init*(1/IPR_IN)
    CalbMB_init <- CalbMF_init*(1/IPR_M)
    CalbAB_init <- CalbAF_init*(1/IPR_A)
    CalbRB_init <- CalbRF_init*(1/IPR_R)
    CalbBoB_init <- CalbBoF_init*(1/IPR_Bo)
    CalbLuB_init <- CalbLuF_init*(1/IPR_Lu)
    CalbSPB_init <- CalbSPF_init*(1/IPR_Sp)
    CalbGoB_init <- CalbGoF_init*(1/IPR_Go)
    CalbHB_init <- CalbHF_init*(1/IPR_H)
    CalbBrB_init <- CalbBrF_init*(1/IPR_Br)
    CalbSKB_init <- CalbSKF_init*(1/IPR_SK)
    
    #Alpha2mu-globulin concentration in kidney tissue (mol/L)
    if (sex == "M"){
      a2u_globulin_k = 8.77*kidney_protein_total*1e-3/VKT #mg/L, 8.77 mg/g kidney protein from https://doi.org/10.1016/0300-483X(86)90197-6 
      Ca2uKT_init <- f_alb_avail*(a2u_globulin_k*1e-3/15.5e3) #[mol/L]
      
      #Ca2uKT_init <- 321.51*1e-3 #[umol/L]*1e-3 -->(mol/m3), from Cheng et al. (2017)
      
    }else if(sex == "F"){
      Ca2uKT_init <- 0 #mol/L
    }
    
    #LFABP concentration in kidney and liver tissue (mol/m^3)
    L_FABP_L = 28.2e-3*liver_protein_per_rat/VLT #mg/L, 28.2 ug/mg cytosolic protein from https://doi.org/10.1016/S0021-9258(18)34463-6
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
    
    
    return(list('VB'=VB, 'Vplasma'=Vplasma, 'VK'=VK, 'VKB'=VKB, 
                'VKF'=VKF, 'VKT'=VKT, 'VFil_rest'=VFil_rest, 'VPT' = VPT, 'VBladder' = VBladder,
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
                
                "VPT" = VPT, "VTDL"= VTDL, "VThinAL"=VThinAL, "VThickAL"=VThickAL,
                "VDT" = VDT, "Vduct"=Vduct,
                
                'AK'=AK, 'AKG'=AKG, 'AL'=AL, 
                'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
                'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
                'AIN'=AIN, 'AGo'=AGo,
                'ASK'= ASK, 'ABo'=ABo,
                
                "SKi" = SKi,"SLi" = SLi,"SSt" = SSt,"SIn" = SIn,
                "SMu" = SMu,"SAd" = SAd,"SRe" = SRe,"SLu" = SLu,
                "SSp" = SSp,"SHt" = SHt,"SBr" = SBr,"SGo" = SGo,
                "SSk" = SSk,"SBo" = SBo,
                
                
                'PeffK'=PeffK, 'PeffL'=PeffL, 
                'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'PeffLu' = PeffLu,
                'PeffSP'=PeffSP, 'PeffH'=PeffH, 'PeffBr'=PeffBr, 'PeffST' = PeffST,
                'PeffIN'=PeffIN, 'PeffGo'=PeffGo,
                'PeffSK' = PeffSK,  'PeffBo' = PeffBo,  
                
                'Qcardiac'=Qcardiac, 'QBK'=QBK, 
                'QBL'=QBL, 'QBLtot'=QBLtot,
                'QBM'=QBM, 'QBA'=QBA,
                'QBR'=QBR, 'QBLu'=QBLu, 'Qfeces'=Qfeces, 'feces_density'=feces_density,
                'Qbile'=Qbile, 'QGFR'=QGFR,'Qurine'=Qurine, 
                'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST,
                'QBIN'=QBIN, 'QGE'=QGE,
                'QBGo'=QBGo,
                'QBSK'=QBSK, 'QBBo'=QBBo, 'Hct' = Hct,
                
                "QparaKi" = QparaKi,"QparaLi" = QparaLi,"QparaSt" = QparaSt,"QparaIn" = QparaIn,
                "QparaMu" = QparaMu,"QparaAd" = QparaAd,"QparaRe" = QparaRe,"QparaLu" = QparaLu,
                "QparaSp" = QparaSp,"QparaHt" = QparaHt,"QparaBr" = QparaBr,"QparaGo" = QparaGo,
                "QparaSk" = QparaSk,"QparaBo" = QparaBo,
                
                "QPT" = QPT, 'QTDL' = QTDL,
                # "QTDL"= QTDL, "QTAL"= QTAL, "QThinAL"=QThinAL, "QThickAL"=QThickAL,
                # "QDT" = QDT, "Qduct"=Qduct,
                
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
                
                'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'= KmL_Oatp, 'VmL_Ntcp'= VmL_Ntcp,
                'VmL_Oatp2'=VmL_Oatp2, 'KmL_Oatp2'= KmL_Oatp2, 
                'VmIn_Oatp2'=VmIn_Oatp2, 'KmIn_Oatp2'= KmIn_Oatp2,
                'KmL_Ntcp'= KmL_Ntcp,'VmK_Oatp'= VmK_Oatp, 
                'Vmbile_transp' = Vmbile_transp,    'Kmbile_transp' = Kmbile_transp,
                'KmK_Oatp'=KmK_Oatp, 'VmLu_Oatp_ap'= VmLu_Oatp_ap,
                'KmLu_Oatp_ap'=KmLu_Oatp_ap, 'VmLu_Oatp_bas'= VmLu_Oatp_bas,
                'KmLu_Oatp_bas'=KmLu_Oatp_bas,
                
                
                'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 
                'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
                'VmK_Urat'=VmK_Urat, 'KmK_Urat'=KmK_Urat, 
                
                'KmK_baso' = KmK_baso, 'KmK_api' = KmK_api,
                'VmK_baso' = VmK_baso,'VmK_api' = VmK_api,
                
                'Papp' = Papp, 'P_passive' = P_passive,
                'kKFKT'=kKFKT, 'kFilKT'=kFilKT, 'kPTKT' = kPTKT,  
                'kLFLT'=kLFLT, 'kLTLbile'=kLTLbile,  'kAFAT'=kAFAT, 
                'kRFRT'=kRFRT,
                'kabST'=kabST, 
                'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF, 'kLuTLuAF'=kLuTLuAF, 'kSPFSPT' =kSPFSPT,
                'kSTFSTT' =kSTFSTT, 'kINFINT' =kINFINT, 'kHFHT' =kHFHT,
                'kBrFBrT' =kBrFBrT, 'kGoFGoT' =kGoFGoT,
                'kSKFSKT' =kSKFSKT, 'kBoFBoT'=kBoFBoT,
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW
    ))
    
  })
}  


setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")

MW <- 414.07 #g/mol

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
Kemp_OR_Ffeces <- openxlsx::read.xlsx("Data/PFOA_Feces_female_oral_25_mg_per_kg-Loc.xlsx")
Kemp_OR_Mfeces <- openxlsx::read.xlsx("Data/PFOA_Feces_male_oral_25_mg_per_kg-Loc.xlsx")
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

dataset <- list("df1" = kudo_high_dose, "df2" = kudo_low_dose, "df3" = kim_IV_Mtissues, "df4" = kim_OR_Mtissues,
                "df5" = kim_IV_Ftissues, "df6" = kim_OR_Ftissues, "df7" = dzi_OR_Mtissues, "df8" = dzi_OR_Ftissues,
                "df9" = kim_OR_Mblood, "df10" = kim_IV_Mblood, "df11" = Lup_OR_Ftissues, "df12" = Kemp_OR_Ffeces,
                "df13" = Kemp_OR_Mfeces, "df14" = Kemp_OR_Furine_low, "df15" = Kemp_OR_Furine_med, "df16" = Kemp_OR_Furine_high,
                "df17" = Kemp_OR_Murine_low, "df18" = Kemp_OR_Murine_med, "df19" = Kemp_OR_Murine_high, 
                "df20" = dzi_IV_Mserum, "df21" = dzi_OR_Mserum_low, "df22" = dzi_OR_Mserum_medium,
                "df23" = dzi_OR_Mserum_high, "df24" = dzi_IV_Fserum, "df25" = dzi_OR_Fserum_low, "df26" = dzi_OR_Fserum_medium,
                "df27" = dzi_OR_Fserum_high, "df28" = kim_OR_Fblood, "df29" = kim_IV_Fblood, "df30" = gus_OR_Mblood,
   "df31" = gus_OR_Mtissues)


setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Scenario_Bladder/Training/AAFE/full_params_PT_excreta_not_restricted_test")

# Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_IV_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

# Set up simulations for the 6th case, i.e. kim (2016) ORAL female tissues
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)



sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]


# Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
BW <- 0.2  # body weight (kg) not reported
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)



sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]




##########################
#-------------------------
# Kemper 2003 (Loccisano)
#-------------------------
##########################

# Set up simulations for the 12th case, i.e.Kemper 2003 (Loccisano) ORAL female feces

sex <- "F"
BW <- 0.2 #kg
sample_time <- seq(0,672,1)
admin.type <-"oral"
admin.dose <- 25 * BW*1000 #ug
admin.time <- 0

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


# sample_time: a vector of time points to solve the ODEs
sample_time=c(0, 0.25, 0.5, 1, 1.5, 2, seq(4,100,2), seq(104,672,4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_Kemp_OR_Ffeces <-  solution[, c("time", "Mfeces")]


##########################
#-------------------------
# Kemper 2003 (Worley)
#-------------------------
##########################

# Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL female urine LOW

sex <- "F"
BW <- 0.2 #kg
sample_time <- c(0, 0.25, 0.5, 1, 1.5, 2, seq(4,168,1), 169)
admin.type <-"oral"
admin.time <- 0

#Female, oral 1mg/kg dose
admin.dose <- 1 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))


events <- create.events(parameters)
inits <- create.inits (parameters)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))


preds_Kemp_OR_Furine_low <-  solution[, c("time", "Murine")]


# Set up simulations for the 15th case, i.e.Kemper 2003 (Worley) ORAL female urine MEDIUM

admin.dose <- 5 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))


preds_Kemp_OR_Furine_med <-  solution[, c("time", "Murine")]


# Set up simulations for the 16th case, i.e.Kemper 2003 (Worley) ORAL female urine HIGH

admin.dose <- 25 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))


events <- create.events(parameters)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))


preds_Kemp_OR_Furine_high <-  solution[, c("time", "Murine")]



# Set up simulations for the 24th case, i.e. Dzierlenga 2021, IV female serum 
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)


inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, seq(4, 192, 2))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_IV_Fserum <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 25th case, i.e. Dzierlenga 2021, ORAL female serum low
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)


inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0, 96, 1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 26th case, i.e. Dzierlenga 2021, ORAL female serum medium
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)


inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.25, 0.5, 1, 1.5,2, seq(4, 192, 2))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_medium <-  solution[, c("time", "Cplasma")]



# Set up simulations for the 27th case, i.e. Dzierlenga 2021, ORAL female serum high
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 320 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)


inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0, 96, 1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_high <-  solution[, c("time", "Cplasma")]



#################################################################################
#--------------------------------------------------------------------------------
#                                Kim 2016 female
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 28th case, i.e. Kim (2016) ORAL male blood
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)


inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_OR_Fblood <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 29th case, i.e. Kim (2016) IV male blood
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)


inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]






#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets

preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] <- preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] /1000
preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] <- preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000
preds_Kemp_OR_Ffeces[,2:dim(preds_Kemp_OR_Ffeces)[2]] <- preds_Kemp_OR_Ffeces[,2:dim(preds_Kemp_OR_Ffeces)[2]] /1000
preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] <- preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] /1000
preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] <- preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] /1000
preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] <- preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] /1000
preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] <- preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] /1000
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] <- preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] /1000
preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] <- preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] /1000
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000


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


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))
# Convert Kemper ORAL female feces from long to wide format using reshape
experiment12 <- reshape(Kemp_OR_Ffeces[c("Tissue" ,"Time_h", 
                                         "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment12) <- c("Time",unique(Kemp_OR_Ffeces$Tissue))
experiment12$Feces = (experiment12$Feces/100)*0.2*25



# Convert Kemper ORAL female urine low from long to wide format using reshape
experiment14 <- reshape(Kemp_OR_Furine_low [c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment14) <- c("Time",unique(Kemp_OR_Furine_low$Tissue))
experiment14$Urine = (experiment14$Urine/100)*0.2*1


# Convert Kemper ORAL female urine med from long to wide format using reshape
experiment15 <- reshape(Kemp_OR_Furine_med[c("Tissue" ,"Time_h", 
                                             "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment15) <- c("Time",unique(Kemp_OR_Furine_med$Tissue))
experiment15$Urine = (experiment15$Urine/100)*0.2*5

# Convert Kemper ORAL female urine high from long to wide format using reshape
experiment16 <- reshape(Kemp_OR_Furine_high[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment16) <- c("Time",unique(Kemp_OR_Furine_high$Tissue))
experiment16$Urine = (experiment16$Urine/100)*0.2*25

#Convert Dzierlenga 2021, IV female serum from long to wide format using reshape
experiment24 <- reshape(dzi_IV_Fserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment24) <- c("Time",unique(dzi_IV_Fserum$Tissue))

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment25 <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment25) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL female serum medium from long to wide format using reshape
experiment26 <- reshape(dzi_OR_Fserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment26) <- c("Time",unique(dzi_OR_Fserum_medium$Tissue))

#Convert Dzierlenga 2021, ORAL female serum high from long to wide format using reshape
experiment27 <- reshape(dzi_OR_Fserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment27) <- c("Time",unique(dzi_OR_Fserum_high$Tissue))

#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment28 <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment28) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment29<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment29) <- c("Time",unique(kim_IV_Fblood$Tissue))


# Put the experiments in a list
experiments <- list(
                    experiment5 = experiment5, experiment6 = experiment6,  experiment8 = experiment8,
                     experiment12 = experiment12,
                   experiment14 = experiment14, experiment15 = experiment15, experiment16 = experiment16,
                    experiment24 = experiment24,
                    experiment25 = experiment25, experiment26 = experiment26, experiment27 = experiment27,
                    experiment28 = experiment28, experiment29 = experiment29)



colnames(preds_kim_IV_Ftissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")
colnames(preds_kim_OR_Ftissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")

colnames(preds_dzi_OR_Ftissues) <- c("Time","Liver","Kidney","Brain")

colnames(preds_Kemp_OR_Ffeces) <- c ("Time", "Feces")

colnames(preds_Kemp_OR_Furine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_high) <- c ("Time", "Urine")

colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")



# Create a list containing the corresponding predictions
simulations <- list( predictions5 = preds_kim_IV_Ftissues, predictions6 = preds_kim_OR_Ftissues,
                   predictions8 = preds_dzi_OR_Ftissues,  predictions12 = preds_Kemp_OR_Ffeces,
       predictions14 = preds_Kemp_OR_Furine_low, predictions15 = preds_Kemp_OR_Furine_med,
                    predictions16 =preds_Kemp_OR_Furine_high, 
                    predictions24 =preds_dzi_IV_Fserum, predictions24 =preds_dzi_OR_Fserum_low, 
       predictions26 =preds_dzi_OR_Fserum_medium,
                    predictions27 =preds_dzi_OR_Fserum_high, predictions28 = preds_kim_OR_Fblood, 
                    predictions29 = preds_kim_IV_Fblood)


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

solution[,69:80]

