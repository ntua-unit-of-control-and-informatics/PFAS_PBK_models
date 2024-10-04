#female increase in BW: 3.5 g/d
#Male increase in BW: 5.9 g.d
options(max.print=1000000)
library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    # BW in kg
    # Cheng and Ng 2017 Table S1
    # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3),
    #assuming the density of tissue is 1 g/mL.
    # Estimated parameters
    if (sex == "M"){
      RAFOatp_k <- 2.683341e+00
      RAFOat1 <- 3.007561e+04
      RAFOatp_l <- 1.852061e+02
      
      
    }else if(sex == "F"){
      RAFOatp_k <- 2.776635e-19
      RAFOat1 <- 2.207686e+05
      RAFOatp_l <- 2.202331e-06
      
      
    }
    RAFOat3 <- RAFOat1
    RAFUrat <- RAFOatp_k
    RAFOatp2_l <- RAFOatp_l
    RAFOatp_lu_ap <- 8.526492e+00
    RAFOatp_lu_bas <- RAFOatp_lu_ap
    RAFNtcp <- RAFOatp_l
    RAFOatp2_Int <- 6.922875e+01
    f_fabp_avail <- 1
    f_alb_avail <- 1
    
   
    RAF_papp <- 1
    
    koff_alb <- 0.9
    koff_fabp <-  koff_alb
    koff_a2u <- koff_alb
    
    VmK_api <- 0
    VmK_baso <- 0
    KmK_baso <- 1e20
    KmK_api <-   1e20
    KLfabp <- (1.2e5+4e4+1.9e4)  #[L/mol]*1e-3 , value from Cheng et al. (2017)
    Ka <- 5.323317e+04 # 5.8e05 #mol/L
    
    kabsUA <- estimated_params[1]
    kCLEal <- estimated_params[2]
    kCLEua <- estimated_params[3]
    
    
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
    VFil <- 0.25/1000 #renal filtrate volume in L,  from Cheng et al., 2017 (from Arthur, 1986; Bonvalet, 1981)
    
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
    PVINF <- (0.867+0.5)/280 #mL/g BW weight
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
    
    #Upper airways
    PVUA <- 257e-3/288 #mm^3/g, for a 16-wk old male 288 g, Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    VUA <- PVUA * BW  #total volume of nasal cavity 
    
    #respiratory area of nasal cavity
    RA_area = 623 #mm^2 Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    RA_capillary_density = 362 #capillaries/mm^2
    capillary_area = RA_area/(RA_area*RA_capillary_density) #mm^2, area of one capillary
    
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
    PVR <- 1 - PVB - PVK - PVL - PVM - PVA - PVSP - PVH - PVBr - PVGo - PVST - PVIN - PVSK - PVBo - PVUA
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
    
    #Organ surface area
    
    PAUA <- ((1343.5+76.8)*1e-06)/0.288 #m2/kg, total surface area 1343.5 mm^2 for a 16-wk old male 288 g, Gross et al., 1982, https://pubmed.ncbi.nlm.nih.gov/7130058/
    #plus nasopharynx area 76.8 mm^2 calculated by Ménache et al., 1997, https://doi.org/10.1080/00984109708984003
    AUA <- PAUA*BW
    PALF <- 27.2e-4/0.3 #27.2 cm2, BW=0.3 in the study--> 27.2e-4 m2/kg, Mercer et al., 1994 https://doi.org/10.1152/jappl.1987.62.4.1480
    ALF <- PALF*BW
    
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
    DpUa <- 27 #nm ???, assumption
    
    Dps <- c(DpKi, DpLi, DpSt, DpIn, DpMu, DpAd, DpRe, DpLu, DpSp, DpHt, DpBr, DpGo, DpSk, DpBo, DpUa)
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
    SUa <- s_r[15]
    
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
    #fQparaUa <- 2.0E-3 # Assumption based on 1/500 of flow (Dosgra et al. 2020, https://doi.org/10.1016/j.csbj.2020.02.014)
    
    
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
    
    PQbile <- 90/1000/24 # [mL/d/kg BW]/1000/24 --> L/h/kg BW
    Qbile <- PQbile * BW #L/h
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
    
    #Alveolar cells surface area (Type I and II), m^2
    AcALF = ((78.8*2*5320*1e-6) + (125*2*123*1e-6))*BW/0.29  #Stone et al., 1992, BW_ref = 0.29, values for each lung , https://doi.org/10.1165/ajrcmb/6.2.235
    
    #canalicular surface area, m^2
    rat_hep_surf_area = 22.95 * 1e2 # 22.95*1e6 cm2 --> m2,  https://doi.org/10.1074/jbc.271.12.6702
    AcLBilec = 0.01 * 22.95 * 1e2 # m2 , canalicular membrane 1% of the surface area of the hepatocyte,https://www.ncbi.nlm.nih.gov/books/NBK470209/
    
    
    # Following the calculations of Lin et al. (2023) for Caco-2 cells
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
    kFKT <- ((Papp/100) * AK) * n*1000 #m^3/h *1000 ---> L/h
    
    
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
    
    #Lung
    lung_protein_per_gram <- 134 # 134 mg/mL tissue --> 134 mg/g tissue, Figure 2, https://doi.org/10.1007/s00580-021-03242-z 
    kUAB <- kabsUA * RA_area * 1000#absorption rate from upper airways to blood, RA_area in m^2
    CLEal <- kCLEal * ALF* 1000 #clearance rate from alveolar lining fluid to stomach, ALF in m^2
    CLEua <- kCLEua * AUA * 1000#clearance rate rate from upper airways to stomach, AUA in m^2
    #kLuTLuAF <- kLuAF * ALF #transport rate from alveolar lining fluid to lung tissue, ALF in m^2
    
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
    PeffUA <- Papp*10 #mm/h
    
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
                'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,
                'VL'=VL, 'VLB'=VLB, 'VLF'=VLF, 'VLT'=VLT, 'VLbile'=VLbile,
                'VM'=VM, 'VMB'=VMB, 'VMF'=VMF, 'VMT'=VMT, 'VA'=VA, 'VAB'=VAB, 
                'VAF'=VAF, 'VAT'=VAT, 'VUA'=VUA, 'VR'=VR, 'VRB'=VRB, 
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
                
                'AK'=AK, 'AKG'=AKG, 'AL'=AL, 
                'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
                'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
                'AIN'=AIN, 'AGo'=AGo,
                'ASK'= ASK, 'ABo'=ABo,
                
                "SKi" = SKi,"SLi" = SLi,"SSt" = SSt,"SIn" = SIn,
                "SMu" = SMu,"SAd" = SAd,"SRe" = SRe,"SLu" = SLu,
                "SSp" = SSp,"SHt" = SHt,"SBr" = SBr,"SGo" = SGo,
                "SSk" = SSk,"SBo" = SBo, "SUa" =SUa,
                
                
                'PeffK'=PeffK, 'PeffL'=PeffL, 
                'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'PeffLu' = PeffLu,
                'PeffSP'=PeffSP, 'PeffH'=PeffH, 'PeffBr'=PeffBr, 'PeffST' = PeffST,
                'PeffIN'=PeffIN, 'PeffGo'=PeffGo,
                'PeffSK' = PeffSK,  'PeffBo' = PeffBo, 'PeffUA'=PeffUA, 
                
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
                'KmK_Oatp'=KmK_Oatp, 'VmLu_Oatp_ap'= VmLu_Oatp_ap,
                'KmLu_Oatp_ap'=KmLu_Oatp_ap, 'VmLu_Oatp_bas'= VmLu_Oatp_bas,
                'KmLu_Oatp_bas'=KmLu_Oatp_bas,
                
                
                'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 
                'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
                'VmK_Urat'=VmK_Urat, 'KmK_Urat'=KmK_Urat, 
                
                'KmK_baso' = KmK_baso, 'KmK_api' = KmK_api,
                'VmK_baso' = VmK_baso,'VmK_api' = VmK_api,
                
                'Papp' = Papp, 'P_passive' = P_passive,
                'kKFKT'=kKFKT, 'kFKT'=kFKT,  
                'kLFLT'=kLFLT, 'kLTLbile'=kLTLbile, 'kAFAT'=kAFAT, 
                'kRFRT'=kRFRT,
                'kabST'=kabST, 
                'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF, 'kLuTLuAF'=kLuTLuAF, 'kSPFSPT' =kSPFSPT,
                'kSTFSTT' =kSTFSTT, 'kINFINT' =kINFINT, 'kHFHT' =kHFHT,
                'kBrFBrT' =kBrFBrT, 'kGoFGoT' =kGoFGoT,
                'kSKFSKT' =kSKFSKT, 'kBoFBoT'=kBoFBoT,
                'kUAB'= kUAB, 'CLEal'=CLEal, 'CLEua'=CLEua, 
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW, "capillary_area"=capillary_area,
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
    MKT <- MKTf + MKTb
    CKB <- MKB/VKB # blood concentration
    CKBf <- MKBf/VKB
    CKBb <- MKBb/VKB
    CKF <- MKF/VKF  #interstitial fluid concentration
    CKFf <- MKFf/VKF
    CKFb <- MKFb/VKF
    CKT <- MKT/VKT # tissue concentration
    CKTf <- MKTf/VKT
    CKTb <- MKTb/VKT
    CFil <- MFil/VFil# filtrate concentration
    
    
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
    CLbile <- MLbile/VLbile #Bile  canaliculi
    
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
    CSTT <- MSTT/VSTT # tissue concentration
    
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
    CINT <- MINT/VINT # tissue concentration
    
    
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
    CMT <-  MMT/VMT # tissue concentration
    
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
    CAT <- MAT/VAT # tissue concentration
    
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
    CRT <- MRT/VRT # tissue concentration
    
    #Upper airways
    CUA = MUA/VUA # upper airways concentration
    
    #Lung
    
    MLuB <- MLuBf + MLuBb
    MLuF <- MLuFf + MLuFb
    MLuT <- MLuTf
    MLuAF <- MLuAFf + MLuAFb
    CLuB <- MLuB/VLuB # blood concentration
    CLuBf <- MLuBf/VLuB
    CLuBb <- MLuBb/VLuB
    CLuF <- MLuF/VLuF  #Interstitial fluid concentration
    CLuFf <- MLuFf/VLuF
    CLuFb <- MLuFb/VLuF
    CLuT <- MLuT/VLuT #tissue concentration
    CLuAF <- MLuAF/VLuAF #alveolar lining fluid concentration
    CLuAFf <- MLuAFf/VLuAF
    CLuAFb <- MLuAFb/VLuAF
    
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
    CSPT <-  MSPT/VSPT # tissue concentration
    
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
    CHT <- MHT/VHT # tissue concentration
    
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
    CBrT <-  MBrT/VBrT # tissue concentration
    
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
    CGoT <-  MGoT/VGoT # tissue concentration
    
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
    CSKT <- MSKT/VSKT # tissue concentration
    
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
    CBoT <- MBoT/VBoT # tissue concentration
    
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
    dCa2uKTf <- koff_a2u*CKTb/MW/1e6 - kon_a2u*Ca2uKTf*CKTf/MW/1e6
    dCFabpKTf <- koff_fabp*CKTb/MW/1e6 - kon_fabp*CFabpKTf*CKTf/MW/1e6
    
    
    #Calculation of free and bound PFOA in liver tissue
    dCFabpLTf <- koff_fabp*CLTb/MW/1e6 - kon_fabp*CFabpLTf*CLTf/MW/1e6
    
    #Calculation of free and bound PFOA in alveolar lining fluid
    dCalbLuAFf = koff_alb*CLuAFb/MW/1e6 - kon_alb*CalbLuAFf*CLuAFf/MW/1e6
    
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
    dMKTb <- kon_a2u*Ca2uKTf*CKTf*VKT + kon_fabp*CFabpKTf*CKTf*VKT -
      koff_fabp*CKTb*VKT - koff_a2u*CKTb*VKT
    dMLTb <-  kon_fabp*CFabpLTf*CLTf*VLT - koff_fabp*CLTb*VLT
    
    #Alveolar lining fluid
    dMLuAFb <-  kon_alb*CalbLuAFf*CLuAFf*VLuAF -  koff_alb*CLuAFb*VLuAF
    
    #====================================================================================================================
    
    #Arterial Blood
    dMArtf = QBLu*CLuBf - CArtf*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+
                                   QBST+QBIN+QBGo+QBSK+QBBo) - QGFR*CArtf +
      koff_alb*CArtb*VArt - kon_alb*CalbArtf*CArtf*VArt
    
    #Venous Blood
    dMVenf = kUAB * CUA - CVenf*QBLu + QBK*CKBf + QBLtot*CLBf + QBM*CMBf + QBA*CABf + QBR*CRBf+
      QBH*CHBf + QBBr*CBrBf+ QBGo*CGoBf + QBSK*CSKBf + QBBo*CBoBf +
      koff_alb*CVenb*VVen - kon_alb*CalbVenf*CVenf*VVen
    
    #Kidney
    #blood subcompartment
    dMKBf = QBK*CArtf - QBK*CKBf - PeffK*AK*(CKBf-CKFf) - QparaKi*(1-SKi)*CKBf +
      koff_alb*CKBb*VKB - kon_alb*CalbKBf*CKBf*VKB
    #interstitial fluid subcompartment
    dMKFf = QparaKi*(1-SKi)*CKBf+ PeffK*AK*(CKBf-CKFf) - kKFKT*(CKFf-CKTf) -
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) - (VmK_Oat3*CKFf/(KmK_Oat3+CKFf))  +
      (VmK_baso*CKTf/(KmK_baso+CKTf)) +  koff_alb*CKFb*VKF - kon_alb*CalbKFf*CKFf*VKF
    #Kidney proximal tubule cells subcompartment
    dMKTf = kKFKT*(CKFf-CKTf) - kFKT*(CKTf - CFil) + (VmK_Oatp*CFil/(KmK_Oatp+CFil)) +
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) + (VmK_Oat3*CKFf/(KmK_Oat3+CKFf)) + 
      (VmK_Urat*CFil/(KmK_Urat+CFil))  - (VmK_baso*CKTf/(KmK_baso+CKTf)) -
      (VmK_api*CKTf/(KmK_api+CKTf))+
      koff_fabp*CKTb*VKT + koff_a2u*CKTb*VKT -kon_fabp*CFabpKTf*CKTf*VKT - kon_a2u*Ca2uKTf*CKTf*VKT
    
    dMFil =  QGFR*CArtf + kFKT*(CKTf - CFil) - (VmK_Oatp*CFil/(KmK_Oatp+CFil)) - 
      (VmK_Urat*CFil/(KmK_Urat+CFil)) + (VmK_api*CKTf/(KmK_api+CKTf))- (Qurine*CFil)
    
    #Liver
    #blood subcompartment
    dMLBf = QBL*CArtf + QBSP*CSPBf + QBIN*CINBf + QBST*CSTBf - 
      QBLtot*CLBf - PeffL*AL*(CLBf-CLFf) - QparaLi*(1-SLi)*CLBf +
      koff_alb*CLBb*VLB -kon_alb*CalbLBf*CLBf*VLB
    #interstitial fluid subcompartment 
    dMLFf =  QparaLi*(1-SLi)*CLBf + PeffL*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) - 
      (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) - (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf)) -
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_alb*CLFb*VLF -kon_alb*CalbLFf*CLFf*VLF
    #Liver tissue subcompartment
    dMLTf = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) + (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf))+
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_fabp*CLTb*VLT-
      kon_fabp*CFabpLTf*CLTf*VLT - kLTLbile*(CLTf-CLbile)
    
    #Bile  canaliculi subcompartment
    
    dMLbile = kLTLbile*(CLTf-CLbile) - CLbile*Qbile
    
    
    #Stomach
    #blood subcompartment
    dMSTBf = QBST*CArtf - QBST*CSTBf - PeffST*AST*(CSTBf-CSTFf) -  QparaSt*(1-SSt)*CSTBf +
      koff_alb*CSTBb*VSTB-kon_alb*CalbSTBf*CSTBf*VSTB
    #interstitial fluid subcompartment 
    dMSTFf = QparaSt*(1-SSt)*CSTBf + PeffST*AST*(CSTBf-CSTFf) - kSTFSTT*(CSTFf-CSTT) +
      koff_alb*CSTFb*VSTF - kon_alb*CalbSTFf*CSTFf*VSTF
    #Stomach tissue subcompartment
    dMSTTf = kSTFSTT*(CSTFf-CSTT) + kabST*CSTL
    #Stomach lumen
    dMSTL = - QGE*CSTL -kabST*CSTL + CLEal*CLuAFf + CLEua*CUA
    
    
    #Intestine
    #blood subcompartment
    dMINBf = QBIN*CArtf - QBIN*CINBf - PeffIN*AIN*(CINBf-CINFf) - QparaIn*(1-SIn)*CINBf +
      koff_alb*CINBb*VINB - kon_alb*CalbINBf*CINBf*VINB
    #interstitial fluid subcompartment 
    dMINFf = QparaIn*(1-SIn)*CINBf + PeffIN*AIN*(CINBf-CINFf) - kINFINT*(CINFf-CINT) +
      koff_alb*CINFb*VINF - kon_alb*CalbINFf*CINFf*VINF
    #Intestine tissue subcompartment
    dMINTf = kINFINT*(CINFf-CINT) + P_passive*CINL + (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))
    #Intestine lumen
    dMINL = QGE*CSTL - (Qfeces*CINL) - P_passive*CINL +  CLbile*Qbile - 
      (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))
    
    
    #Muscle
    #blood subcompartment
    dMMBf = QBM*CArtf - QBM*CMBf - PeffM*AM*(CMBf-CMFf) - QparaMu*(1-SMu)*CMBf +
      koff_alb*CMBb*VMB - kon_alb*CalbMBf*CMBf*VMB
    #interstitial fluid subcompartment 
    dMMFf = QparaMu*(1-SMu)*CMBf + PeffM*AM*(CMBf-CMFf) - kMFMT*(CMFf- CMT) +
      koff_alb*CMFb*VMF - kon_alb*CalbMFf*CMFf*VMF
    #Muscle tissue subcompartment 
    dMMTf = kMFMT*(CMFf- CMT)
    
    
    #Adipose
    #blood subcompartment
    dMABf = QBA*CArtf - QBA*CABf - PeffA*AA*(CABf-CAFf) - QparaAd*(1-SAd)*CABf +
      koff_alb*CABb*VAB - kon_alb*CalbABf*CABf*VAB
    #interstitial fluid subcompartment 
    dMAFf = QparaAd*(1-SAd)*CABf + PeffA*AA*(CABf-CAFf) - kAFAT*(CAFf-CAT) +
      koff_alb*CAFb*VAF - kon_alb*CalbAFf*CAFf*VAF
    #Adipose tissue subcompartment 
    dMATf =  kAFAT*(CAFf-CAT) 
    
    
    #Rest of body
    #blood subcompartment
    dMRBf = QBR*CArtf - QBR*CRBf - PeffR*AR*(CRBf-CRFf) - QparaRe*(1-SRe)*CRBf +
      koff_alb*CRBb*VRB - kon_alb*CalbRBf*CRBf*VRB
    #interstitial fluid subcompartment 
    dMRFf = QparaRe*(1-SRe)*CRBf + PeffR*AR*(CRBf-CRFf) - kRFRT*(CRFf -CRT) +
      koff_alb*CRFb*VRF - kon_alb*CalbRFf*CRFf*VRF
    #Rest of body tissue subcompartment 
    dMRTf = kRFRT*(CRFf -CRT) 
    
    #Upper airways
    dMUA = - CLEua*CUA  - kUAB * CUA 
    
    #Lung 
    #blood subcompartment
    dMLuBf = CVenf*QBLu - QBLu*CLuBf - PeffLu*ALu*(CLuBf-CLuFf) - QparaLu*(1-SLu)*CLuBf +
      koff_alb*CLuBb*VLuB - kon_alb*CalbLuBf*CLuBf*VLuB
    #interstitial fluid subcompartment
    dMLuFf = QparaLu*(1-SLu)*CLuBf + PeffLu*ALu*(CLuBf-CLuFf) + kLuTLuF*(CLuT-CLuFf) + 
      koff_alb*CLuFb*VLuF - kon_alb*CalbLuFf*CLuFf*VLuF - (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf))
    #Lung tissue
    dMLuTf =  - kLuTLuF*(CLuT-CLuFf) -  kLuTLuAF*(CLuT-CLuAFf) + (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) +
                (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf))
    #Alveolar lining fluid
    dMLuAFf =  kLuTLuAF*(CLuT-CLuAFf) + koff_alb*CLuAFb*VLuAF - kon_alb*CalbLuAFf*CLuAFf*VLuAF -
               (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf)) - CLEal*CLuAFf
    
    
    #Spleen
    #blood subcompartment
    dMSPBf = QBSP*CArtf - QBSP*CSPBf - PeffSP*ASP*(CSPBf-CSPFf) - QparaSp*(1-SSp)*CSPBf + 
      koff_alb*CSPBb*VSPB - kon_alb*CalbSPBf*CSPBf*VSPB 
    #interstitial fluid subcompartment 
    dMSPFf = QparaSp*(1-SSp)*CSPBf + PeffSP*ASP*(CSPBf-CSPFf) - kSPFSPT*(CSPFf -CSPT) +
      koff_alb*CSPFb*VSPF - kon_alb*CalbSPFf*CSPFf*VSPF
    #Spleen tissue subcompartment 
    dMSPTf = kSPFSPT*(CSPFf -CSPT) 
    
    
    #Heart
    #blood subcompartment
    dMHBf = QBH*CArtf - QBH*CHBf - PeffH*AH*(CHBf-CHFf) - QparaHt*(1-SHt)*CHBf + 
      koff_alb*CHBb*VHB - kon_alb*CalbHBf*CHBf*VHB
    #interstitial fluid subcompartment 
    dMHFf = QparaHt*(1-SHt)*CHBf + PeffH*AH*(CHBf-CHFf) - kHFHT*(CHFf -CHT) + 
      koff_alb*CHFb*VHF - kon_alb*CalbHFf*CHFf*VHF
    #Heart tissue subcompartment 
    dMHTf = kHFHT*(CHFf -CHT) 
    
    
    #Brain
    #blood subcompartment
    dMBrBf = QBBr*CArtf - QBBr*CBrBf - PeffBr*ABr*(CBrBf-CBrFf) - QparaBr*(1-SBr)*CBrBf + 
      koff_alb*CBrBb*VBrB - kon_alb*CalbBrBf*CBrBf*VBrB
    #interstitial fluid subcompartment 
    dMBrFf = QparaBr*(1-SBr)*CBrBf + PeffBr*ABr*(CBrBf-CBrFf) - kBrFBrT*(CBrFf -CBrT) +
      koff_alb*CBrFb*VBrF - kon_alb*CalBrFf*CBrFf*VBrF
    #Brain tissue subcompartment 
    dMBrTf = kBrFBrT*(CBrFf -CBrT) 
    
    
    #Gonads
    #blood subcompartment
    dMGoBf = QBGo*CArtf - QBGo*CGoBf - PeffGo*AGo*(CGoBf-CGoFf) - QparaGo*(1-SGo)*CGoBf +
      koff_alb*CGoBb*VGoB - kon_alb*CalbGoBf*CGoBf*VGoB
    #interstitial fluid subcompartment 
    dMGoFf = QparaGo*(1-SGo)*CGoBf + PeffGo*AGo*(CGoBf-CGoFf) - kGoFGoT*(CGoFf -CGoT) +
      koff_alb*CGoFb*VGoF - kon_alb*CalbGoFf*CGoFf*VGoF
    #gonads tissue subcompartment 
    dMGoTf = kGoFGoT*(CGoFf -CGoT) 
    
    
    #Skin
    #blood subcompartment
    dMSKBf = QBSK*CArtf - QBSK*CSKBf - PeffSK*ASK*(CSKBf-CSKFf) - QparaSk*(1-SSk)*CSKBf +
      koff_alb*CSKBb*VSKB - kon_alb*CalbSKBf*CSKBf*VSKB
    #interstitial fluid subcompartment
    dMSKFf = QparaSk*(1-SSk)*CSKBf + PeffSK*ASK*(CSKBf-CSKFf) - kSKFSKT*(CSKFf -CSKT) +
      koff_alb*CSKFb*VSKF - kon_alb*CalbSKFf*CSKFf*VSKF
    #Skin tissue subcompartment
    dMSKTf = kSKFSKT*(CSKFf -CSKT)
    
    
    #Bones
    #blood subcompartment
    dMBoBf = QBBo*CArtf - QBBo*CBoBf - PeffBo*ABo*(CBoBf-CBoFf) - QparaBo*(1-SBo)*CBoBf +
      koff_alb*CBoBb*VBoB - kon_alb*CalbBoBf*CBoBf*VBoB
    #interstitial fluid subcompartment
    dMBoFf = QparaBo*(1-SBo)*CBoBf + PeffBo*ABo*(CBoBf-CBoFf) - kBoFBoT*(CBoFf -CBoT) +
      koff_alb*CBoFb*VBoF -  kon_alb*CalbBoFf*CBoFf*VBoF
    #Bones tissue subcompartment
    dMBoTf = kBoFBoT*(CBoFf -CBoT)
    
    #Excreta#
    dMfeces <- Qfeces*CINL
    dMurine <- Qurine*CFil
    dVurine = Qurine
    dVfeces = Qfeces
    
    #Concentration calculation in each compartment 
    
    Cblood <- (MVen +MArt)/ (VVen+VArt)
    Mblood <- MVen +MArt
    Cplasma <- Cblood/(1-Hct)
    
    Ckidney <- (MKB + MKF+ MKT)/(VKB+VKF+VKT)
    Mkidney <- MKB + MKF+ MKT
    
    Cliver <- (MLB + MLF+ MLT + MLbile )/(VLB+VLF+VLT+VLbile)
    Mliver <- MLB + MLF+ MLT + MLbile
    
    Cstomach <-  (MSTB + MSTF+ MSTT + MSTL)/(VSTB+VSTF+VSTT+VSTL)
    Cintestine <-  (MINB + MINF+ MINT+MINL)/(VINB+VINF+VINT+VINL)
    Cmuscle <-  (MMB + MMF+ MMT)/(VMB+VMF+VMT)
    Cadipose <-  (MAB + MAF+ MAT)/(VAB+VAF+VAT)
    
    CUpperair <- MUA/VUA
    CalveolarLF <- (MLuAFf+MLuAFb)/VLuAF
    Clungs <-  (MLuB + MLuF+ MLuT + MLuAF)/(VLuB+VLuF+VLuT+VLuAF)
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
    
    CalveolarLF <- (MLuAFf+MLuAFb)/VLuAF
    
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
            'dCalbBoFf' = dCalbBoFf, 'dCa2uKTf' = dCa2uKTf, 'dCFabpKTf' = dCFabpKTf,
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
            'dMKTb' = dMKTb, 'dMLTb' = dMLTb, 'dMLuAFb'=dMLuAFb,
            
            
            'dMArtf'=dMArtf, 'dMVenf'=dMVenf, 'dMKBf'=dMKBf, 
            'dMKFf'=dMKFf, 'dMKTf'=dMKTf,
            'dMFil'=dMFil,  'dMLBf'=dMLBf, 
            'dMLFf'=dMLFf, 'dMLTf'=dMLTf, 'dMLbile'=dMLbile,
            
            'dMSTBf'=dMSTBf, 'dMSTFf'=dMSTFf, 'dMSTTf'=dMSTTf, 'dMSTL'=dMSTL,
            'dMINBf'=dMINBf, 'dMINFf'=dMINFf, 'dMINTf'=dMINTf,'dMINL'=dMINL,
            
            'dMMBf'=dMMBf, 'dMMFf'=dMMFf, 'dMMTf'=dMMTf,
            'dMABf'=dMABf, 'dMAFf'=dMAFf, 'dMATf'=dMATf, 
            'dMRBf'=dMRBf, 'dMRFf'=dMRFf,'dMRTf'=dMRTf,
            'dMLuBf'=dMLuBf, 'dMLuFf'=dMLuFf,'dMLuTf'=dMLuTf,
            'dMLuAFf' = dMLuAFf,'dMUA'=dMUA,
            
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
    'CKTf'=CKTf, 'CKTb'=CKTb, 'CFil'=CFil, 'CLB'=CLB, 'CLBf'=CLBf, 'CLBb'=CLBb, 
    'CLF'=CLF, 'CLFf'=CLFf, 'CLFb'=CLFb, 'CLT'=CLT, 'CLTf'=CLTf, 'CLTb'=CLTb, 
    'CSTB'=CSTB, 'CSTBf'=CSTBf, 'CSTBb'=CSTBb, 'CSTF'=CSTF, 'CSTFf'=CSTFf, 'CSTFb'=CSTFb,
    'CSTT'=CSTT, 'CINB'=CINB, 'CINBf'=CINBf, 'CINBb'=CINBb, 'CINF'=CINF, 'CINFf'=CINFf,
    'CINFb'=CINFb, 'CINT'=CINT, 'CSTL'=CSTL, 'CINL'=CINL, 'CMB'=CMB, 'CMBf'=CMBf,
    'CMBb'=CMBb, 'CMF'=CMF, 'CMFf'=CMFf, 'CMFb'=CMFb, 'CMT'=CMT, 'CAB'=CAB, 'CABf'=CABf,
    'CABb'=CABb, 'CAF'=CAF, 'CAFf'=CAFf, 'CAFb'=CAFb, 'CAT'=CAT, 'CRB'=CRB, 'CRBf'=CRBf,
    'CRBb'=CRBb, 'CRF'=CRF, 'CRFf'=CRFf, 'CRFb'=CRFb, 'CRT'=CRT, 'CUA'=CUA, 'CLuB'=CLuB,
    'CLuBf'=CLuBf, 'CLuBb'=CLuBb, 'CLuF'=CLuF, 'CLuFf'=CLuFf, 'CLuFb'=CLuFb, 'CLuT'=CLuT,
    'CLuAF'=CLuAF, 'CLuAFf'=CLuAFf, 'CLuAFb'=CLuAFb, 'CSPB'=CSPB, 'CSPBf'=CSPBf, 
    'CSPBb'=CSPBb, 'CSPF'=CSPF, 'CSPFf'=CSPFf, 'CSPFb'=CSPFb, 'CSPT'=CSPT, 'CHB'=CHB,
    'CHBf'=CHBf, 'CHBb'=CHBb, 'CHF'=CHF, 'CHFf'=CHFf, 'CHFb'=CHFb, 'CHT'=CHT,
    'CBrB'=CBrB, 'CBrBf'=CBrBf, 'CBrBb'=CBrBb, 'CBrF'=CBrF, 'CBrFf'=CBrFf, 'CBrFb'=CBrFb,
    'CBrT'=CBrT, 'CGoB'=CGoB, 'CGoBf'=CGoBf, 'CGoBb'=CGoBb, 'CGoF'=CGoF, 'CGoFf'=CGoFf, 
    'CGoFb'=CGoFb, 'CGoT'=CGoT, 'CSKB'=CSKB, 'CSKBf'=CSKBf, 'CSKBb'=CSKBb, 'CSKF'=CSKF,
    'CSKFf'=CSKFf, 'CSKFb'=CSKFb, 'CSKT'=CSKT, 'CBoB'=CBoB, 'CBoBf'=CBoBf, 'CBoBb'=CBoBb,
    'CBoF'=CBoF, 'CBoFf'=CBoFf, 'CBoFb'=CBoFb, 'CBoT'=CBoT,
    
    'Cblood'=Cblood, 'Mblood'=Mblood, 'Cplasma'=Cplasma, 
    'Ckidney'=Ckidney, 'Mkidney'=Mkidney, 'Cliver'=Cliver, 'Mliver'=Mliver, 
    'Cstomach'=Cstomach, 'Cintestine'=Cintestine, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose,
    'CUpperair'=CUpperair, 'CalveolarLF'=CalveolarLF, "Clungtissue" = Clungtissue,
    'Clungs'=Clungs, 'Crest'=Crest, 'Ccarcass'=Ccarcass, 'Cfeces'=Cfeces,
    'Curine'=Curine, 'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain, 
    'Mbrain'=Mbrain, 'Cgonads'=Cgonads, 'Cskin'=Cskin, 'Cbones'=Cbones 
    
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
    CalBrFf<- CalbBrF_init; MBrFf<- 0; MBrFb<- 0;
    CalbGoFf<- CalbGoF_init; MGoFf<- 0; MGoFb<- 0; 
    CalbSKFf<- CalbSKF_init; MSKFf<- 0; MSKFb<- 0;
    CalbBoFf<- CalbBoF_init; MBoFf<- 0; MBoFb<- 0;
    
    Ca2uKTf<- Ca2uKT_init; CFabpKTf<- CFabpKT_init; MKTf<- 0; MKTb<- 0; 
    CFabpLTf<- CFabpLT_init; CalbLuAFf<- CalbLuAF_init;
    MLTf<- 0; MLTb<- 0; MLbile <-0; MSTTf <- 0; MINTf <- 0; MMTf <- 0; MATf <- 0; MRTf <- 0; MLuTf <- 0;
    MLuAFf <- 0; MUA <- 0;  MLuAFb<- 0; MSPTf <- 0; MHTf <- 0; MBrTf <- 0;
    MGoTf <- 0; MSKTf <- 0; MBoTf <- 0;  
    
    MFil <-0; Murine <-0;MSTL <-0;  MINL <-0;
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
             'CalbBoFf' = CalbBoFf, 'Ca2uKTf' = Ca2uKTf, 'CFabpKTf' = CFabpKTf,
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
             'MKTb' = MKTb, 'MLTb' = MLTb, 'MLuAFb'=MLuAFb,
             
             
             'MArtf'=MArtf, 'MVenf'=MVenf, 'MKBf'=MKBf, 
             'MKFf'=MKFf, 'MKTf'=MKTf,
             'MFil'=MFil,  'MLBf'=MLBf, 
             'MLFf'=MLFf, 'MLTf'=MLTf, 'MLbile'=MLbile,
             
             'MSTBf'=MSTBf, 'MSTFf'=MSTFf, 'MSTTf'=MSTTf, 'MSTL'=MSTL,
             'MINBf'=MINBf, 'MINFf'=MINFf, 'MINTf'=MINTf,'MINL'=MINL,
             
             'MMBf'=MMBf, 'MMFf'=MMFf, 'MMTf'=MMTf,
             'MABf'=MABf, 'MAFf'=MAFf, 'MATf'=MATf, 
             'MRBf'=MRBf, 'MRFf'=MRFf,'MRTf'=MRTf,
             'MLuBf'=MLuBf, 'MLuFf'=MLuFf,'MLuTf'=MLuTf,'MLuAFf' = MLuAFf,
             'MUA'=MUA,
             
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
        events <- list(data = rbind(data.frame(var = c("MLuAFf"),  time = admin.time, 
                                               value = admin.dose*depfr_AF, method = c("add")) ))
      }else if (admin.type == "nasal"){
        events <- list(data = rbind(data.frame(var = c("MUA"),  time = admin.time, 
                                               value = c(admin.dose*depfr_head), method = c("add")),
                                    data.frame(var = c("MLuAFf"),  time = admin.time, 
                                               value = c(admin.dose*depfr_AF), method = c("add")) ))
       
        
      }
    }
    return(events)
  })
}


#############################
#--------------------------
# Objective function
#-------------------------
obj.func <- function(x, dataset){
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)
  depfr_head <- 0
  
  ##########################
  #-------------------------
  # Gustafsson Inhalation male blood
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. Gustafsson (2022) Inhalation male blood
  BW <- 0.5125  #kg, from Gustafsson et al., 2022
  sex <- "M"
  inhalation_params=estimate_BFi_TVi(sex, BW)
  BFi = inhalation_params["BFi"]# breaths/h
  TVi = inhalation_params["TVi"]# L per breath
  duration <- 0.375 #hours, 22.5 min
  admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
  #admin.dose_mg_per_L <- 1000*335 * 0.3/1000 # [ug PFOA/g dust] * [g dust/L]  administered dose in mg/m^3 0.25-0.35 mg/L
  #depfr_AF <- (0.6107+0.0543)
  k = duration*24 #partition of administration packages
  admin.dose <- rep((admin.dose_per_g*BW*1000)/k, length.out = k) #ug PFOA, for 22.5 min inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "inh"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = 0, "depfr_AF" = 1 )
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,48,0.2)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Gustafsson (2022)
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
  
  score[1] <- AAFE(predictions = preds_gus_INH_Mblood, observations = obs_gus_INH_Mblood)
  
  
  
  ##########################
  #-------------------------
  # Gustafsson Inhalation male tissues
  #-------------------------
  ##########################
 
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Gustafsson (2022) Inhalation male tissues
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  column_names <- c("CalveolarLF","Cliver","Clungtissue", "Ckidney")
  
  preds_gus_INH_Mtissues <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, "time"]
    
    preds_gus_INH_Mtissues [[i]] <- solution[solution$time %in% exp_time, column_names[i]]/1000
  }
  
  
  obs_gus_INH_Mtissues <- list(exp_data[exp_data$Tissue == "ALF", "concentration"],
                               exp_data[exp_data$Tissue == "Liver", "concentration"],
                               exp_data[exp_data$Tissue == "Lung", "concentration"], 
                               exp_data[exp_data$Tissue == "Kidney", "concentration"]) 
  
  
  
  score[2] <- AAFE(predictions = preds_gus_INH_Mtissues, observations = obs_gus_INH_Mtissues)
  
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single low
  #-------------------------
  ##########################
  # Set up simulations for the 3rd case, i.e. Hinderliter Inhalation male single low
  BW <- 0.225  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
  sex <- "M"
  inhalation_params=estimate_BFn_TVn(sex, BW)
  BFn = inhalation_params["BFn"]# 1/h
  TVn = inhalation_params["TVn"]# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 1.2 # administered dose in mg/m^3
  depfr_head <- 0.2864
  depfr_AF <- (0.1440+0.0254)
  k = 6*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of Hinderliter Inhalation male single low
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
  
  
  score[3] <- AAFE(predictions = preds_hind_INH_Mblood_low, observations = obs_hind_INH_Mblood_low)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single medium
  #-------------------------
  ##########################
  # Set up simulations for the 4th case, i.e. Hinderliter Inhalation male single medium
  BW <- 0.225  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
  sex <- "M"
  inhalation_params=estimate_BFn_TVn(sex, BW)
  BFn = inhalation_params["BFn"]# 1/h
  TVn = inhalation_params["TVn"]# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 9.8 # administered dose in mg/m^3
  depfr_head <- 0.3057
  depfr_AF <- (0.1195+0.0243)
  k = 6*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of Hinderliter Inhalation male single medium
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
  
  
  score[4] <- AAFE(predictions = preds_hind_INH_Mblood_medium, observations = obs_hind_INH_Mblood_medium)
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation male single high
  #-------------------------
  ##########################
  # Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single high
  BW <- 0.225  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
  sex <- "M"
  inhalation_params=estimate_BFn_TVn(sex, BW)
  BFn = inhalation_params["BFn"]# 1/h
  TVn = inhalation_params["TVn"]# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 27 # administered dose in mg/m^3
  depfr_head <- 0.3573
  depfr_AF <- (0.1618+0.0241)
  k = 6*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of Hinderliter Inhalation male single high
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
  
  
  score[5] <- AAFE(predictions = preds_hind_INH_Mblood_high, observations = obs_hind_INH_Mblood_high)
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single low
  #-------------------------
  ##########################
  # Set up simulations for the 6th case, i.e. Hinderliter Inhalation female single low
  BW <- 0.21  #kg, not reported in the study - 180-240 g average BW of female CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
  sex <- "F"
  inhalation_params=estimate_BFn_TVn(sex, BW)
  BFn = inhalation_params["BFn"]# 1/h
  TVn = inhalation_params["TVn"]# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 1.2 # administered dose in mg/m^3
  depfr_head <- 0.2822
  depfr_AF <- (0.1148+0.0177)
  k = 6*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,9,0.04)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of Hinderliter Inhalation male single low
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
  
  
  score[6] <- AAFE(predictions = preds_hind_INH_Fblood_low, observations = obs_hind_INH_Fblood_low)
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single medium
  #-------------------------
  ##########################
  # Set up simulations for the 7th case, i.e. Hinderliter Inhalation female single low
  BW <- 0.21  #kg, not reported in the study - 180-240 g average BW of female CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
  sex <- "F"
  inhalation_params=estimate_BFn_TVn(sex, BW)
  BFn = inhalation_params["BFn"]# 1/h
  TVn = inhalation_params["TVn"]# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 9.8 # administered dose in mg/m^3
  depfr_head <- 0.3101
  depfr_AF <- (0.0939+0.0165)
  k = 6*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df7=========================================================
  
  exp_data <- dataset$df7 # retrieve data of Hinderliter Inhalation male single medium
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
  
  
  score[7] <- AAFE(predictions = preds_hind_INH_Fblood_medium, observations = obs_hind_INH_Fblood_medium)
  
  
  
  ##########################
  #-------------------------
  # Hinderliter Inhalation female single high
  #-------------------------
  ##########################
  # Set up simulations for the 8th case, i.e. Hinderliter Inhalation female single high
  BW <- 0.21  #kg, not reported in the study - 180-240 g average BW of female CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
  sex <- "F"
  inhalation_params=estimate_BFn_TVn(sex, BW)
  BFn = inhalation_params["BFn"]# 1/h
  TVn = inhalation_params["TVn"]# L
  duration <- 6 #hours
  admin.dose_mg_per_m3 <- 27 # administered dose in mg/m^3
  depfr_head <- 0.3372
  depfr_AF <- (0.1327+0.0177)
  k = 6*duration
  admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
  admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
  admin.type <- "nasal"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )
  
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time= seq(0,30,0.1)
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df8=========================================================
  
  exp_data <- dataset$df8 # retrieve data of Hinderliter Inhalation male single high
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
  
  
  score[8] <- AAFE(predictions = preds_hind_INH_Fblood_high, observations = obs_hind_INH_Fblood_high)
  
  
  ########################################################################################
  # Estimate final score
  
  final_score <- mean(score, na.rm = TRUE)
  return(final_score)
  
}


################################################################################

#setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")
#setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")

MW <- 414.07 #g/mol
source("Goodness-of-fit-metrics.R")

# Read data
gus_INH_Mblood <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Plasma Male rats_Inhalation.xlsx")
gus_INH_Mtissues <- openxlsx::read.xlsx("Inhalation_data/Gustafsson 2022_PFOA_Tissues Male rats_Inhalation.xlsx")
hind_INH_Mblood_low <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_Low_dose.xlsx")
hind_INH_Mblood_medium <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_Medium_dose.xlsx")
hind_INH_Mblood_high <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_male_plasma_single_High_dose.xlsx")
hind_INH_Fblood_low <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_Low_dose.xlsx")
hind_INH_Fblood_medium <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_Medium_dose.xlsx")
hind_INH_Fblood_high <- openxlsx::read.xlsx("Inhalation_data/Hinderliter_2006_female_plasma_single_High_dose.xlsx")

#setwd("C:/Users/user/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Scenario_2/Training/AAFE/NoStomachAbs")
setwd("C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/New_Inhalation_models/Training/AAFE/Inhalation_koff_09_lung_same_RAF_bile")

dataset <- list("df1" = gus_INH_Mblood, "df2" = gus_INH_Mtissues, 
                "df3" = hind_INH_Mblood_low, "df4" = hind_INH_Mblood_medium,
                "df5" = hind_INH_Mblood_high,"df6" = hind_INH_Fblood_low,
                "df7" = hind_INH_Fblood_medium, "df8" = hind_INH_Fblood_high)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX",#"NLOPT_LN_NEWUOA","NLOPT_LN_SBPLX"
              "xtol_rel" = 1e-07,
              "ftol_rel" = 0.0,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = 500, 
              "print_level" = 1)

# Create initial conditions (zero initialisation)
#Parameter names:
# Male RAFOatp_k, Male RAFOat1, Male RAFOat3, Male RAFOatp_l,Male RAFNtcp
# Female RAFOatp_k, Female RAFOat1, Female RAFOat3, Female RAFOatp_l,female RAFNtcp
#  CF_Peff


N_pars <- 3 # Number of parameters to be fitted
fit <-  c(log(1e2), log(1),log(1))

lb	= c(log(1e-5), log(1e-5),log(1e-5))
ub = c(log(1e5), log(1e5),log(1e5))


# Run the optimization algorithmm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             lb	= lb,
                             ub = ub,
                             opts = opts,
                             dataset = dataset)

#estimated_params <- exp(optimizer$solution)
estimated_params <- exp(optimizer$solution)



#################################################################################
#--------------------------------------------------------------------------------
#                                Gustafsson Inhalation male blood
#-------------------------------------------------------------------------------
#################################################################################

depfr_head <- 0
# Set up simulations for the 1st case, i.e. Gustafsson Inhalation male blood
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
inhalation_params=estimate_BFi_TVi(sex, BW)
BFi = inhalation_params["BFi"]# 1/h
TVi = inhalation_params["TVi"]# L
duration <- 0.375 #hours, 22.5 min
admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
admin.dose_mg_per_m3 <- 300 # administered dose in mg/m^3 0.25-0.35 mg/L
depfr_AF <- (0.6107+0.0543)
k = duration*24 #partition of administration packages
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFi*TVi)/k, length.out = k) #ug PFOA, for 22.5 min inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "inh"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,48,0.2)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_gus_INH_Mblood <-  solution[, c("time", "Cplasma")]



#################################################################################
#--------------------------------------------------------------------------------
#                                Gustafsson Inhalation male tissues
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 2nd case, i.e. Gustafsson Inhalation male tissues
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
inhalation_params=estimate_BFi_TVi(sex, BW)
BFi = inhalation_params["BFi"]# 1/h
TVi = inhalation_params["TVi"]# L
duration <- 0.375 #hours, 22.5 min
admin.dose_per_g <- 0.164 # administered dose in mg PFOA/kg BW 
admin.dose_mg_per_m3 <- 300 # administered dose in mg/m^3 0.25-0.35 mg/L
depfr_AF <- (0.6107+0.0543)
k = duration*24 #partition of administration packages
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFi*TVi)/k, length.out = k) #ug PFOA, for 22.5 min inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "inh"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,48,0.2)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_gus_INH_Mtissues <-  solution[, c("time", "CalveolarLF","Cliver", "Clungtissues", "Ckidney")]



#################################################################################
#--------------------------------------------------------------------------------
#                        Hinderliter Inhalation male single low
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 3d case, i.e. Hinderliter Inhalation male single low
BW <- 0.225  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "M"
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 1.2 # administered dose in mg/m^3
depfr_head <- 0.2864
depfr_AF <- (0.1440+0.0254)
k = 6*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )



params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_hind_INH_Mblood_low <-  solution[, c("time", "Cplasma")]


#################################################################################
#--------------------------------------------------------------------------------
#                        Hinderliter Inhalation male single medium
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 4th case, i.e. Hinderliter Inhalation male single medium
BW <- 0.225  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "M" 
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 10 # administered dose in mg/m^3
depfr_head <- 0.3057
depfr_AF <- (0.1195+0.0243)
k = 6*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )



params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_hind_INH_Mblood_medium <-  solution[, c("time", "Cplasma")]


#################################################################################
#--------------------------------------------------------------------------------
#                        Hinderliter Inhalation male single high
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 5th case, i.e. Hinderliter Inhalation male single high
BW <- 0.225  #kg, not reported in the study - 200-250 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "M"
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 25 # administered dose in mg/m^3
depfr_head <- 0.3573
depfr_AF <- (0.1618+0.0241)
k = 6*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )



params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_hind_INH_Mblood_high <-  solution[, c("time", "Cplasma")]



#################################################################################
#--------------------------------------------------------------------------------
#                        Hinderliter Inhalation female single low
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 6th case, i.e. Hinderliter Inhalation female single low
BW <- 0.21  #kg, not reported in the study - 180-240 g average BW of female CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "F" 
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 1.2 # administered dose in mg/m^3
depfr_head <- 0.2822
depfr_AF <- (0.1148+0.0177)
k = 6*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )



params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,9,0.04)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_hind_INH_Fblood_low <-  solution[, c("time", "Cplasma")]



#################################################################################
#--------------------------------------------------------------------------------
#                        Hinderliter Inhalation female single medium
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 7th case, i.e. Hinderliter Inhalation female single medium
BW <- 0.21  #kg, not reported in the study - 180-240 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "F" 
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 10 # administered dose in mg/m^3
depfr_head <- 0.3101
depfr_AF <- (0.0939+0.0165)
k = 6*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )



params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_hind_INH_Fblood_medium <-  solution[, c("time", "Cplasma")]


#################################################################################
#--------------------------------------------------------------------------------
#                        Hinderliter Inhalation female single high
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 8th case, i.e. Hinderliter Inhalation female single high
BW <- 0.21  #kg, not reported in the study - 180-240 g average BW of male CD® IGS (SD) rats at 6 to 8 weekshttps://animalab.eu/cd-sprague-dawley-igs-rat-crl-cd-sd
sex <- "F" 
inhalation_params=estimate_BFn_TVn(sex, BW)
BFn = inhalation_params["BFn"]# 1/h
TVn = inhalation_params["TVn"]# L
duration <- 6 #hours
admin.dose_mg_per_m3 <- 25 # administered dose in mg/m^3
depfr_head <- 0.3372
depfr_AF <- (0.1327+0.0177)
k = 6*duration
admin.dose <- rep((admin.dose_mg_per_m3*duration*BFn*TVn)/k,length.out = k) #ug PFOA, for 6h inhalation
admin.time <- seq(0,duration ,length.out = k) #time when doses are administered, in hours
admin.type <- "nasal"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex, "depfr_head" = depfr_head, "depfr_AF" = depfr_AF )



params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)


sample_time= seq(0,30,0.1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-3, atol = 1e-3))

preds_hind_INH_Fblood_high <-  solution[, c("time", "Cplasma")]



#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets

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



# Convert Gustafsson Inhalation male blood from long to wide format using reshape
experiment_inh_1 <- reshape(gus_INH_Mblood[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_1) <- c("Time",unique(gus_INH_Mblood$Tissue))


# Convert Gustafsson Inhalation male tissues from long to wide format using reshape
experiment_inh_2 <- reshape(gus_INH_Mtissues[c("Tissue" ,"Time_hours", 
                                               "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_2) <- c("Time",unique(gus_INH_Mtissues$Tissue))


# Convert Hinderliter Inhalation male single low from long to wide format using reshape
experiment_inh_3 <- reshape(hind_INH_Mblood_low[c("Tissue" ,"Time_hours", 
                                                  "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_3) <- c("Time",unique(hind_INH_Mblood_low$Tissue))


# Convert Hinderliter Inhalation male single medium from long to wide format using reshape
experiment_inh_4 <- reshape(hind_INH_Mblood_medium[c("Tissue" ,"Time_hours", 
                                                     "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_4) <- c("Time",unique(hind_INH_Mblood_medium$Tissue))


# Convert Hinderliter Inhalation male single high from long to wide format using reshape
experiment_inh_5 <- reshape(hind_INH_Mblood_high[c("Tissue" ,"Time_hours", 
                                                   "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_5) <- c("Time",unique(hind_INH_Mblood_high$Tissue))


# Convert Hinderliter Inhalation female single low from long to wide format using reshape
experiment_inh_6 <- reshape(hind_INH_Fblood_low[c("Tissue" ,"Time_hours", 
                                                  "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_6) <- c("Time",unique(hind_INH_Fblood_low$Tissue))


# Convert Hinderliter Inhalation female single medium from long to wide format using reshape
experiment_inh_7 <- reshape(hind_INH_Fblood_medium[c("Tissue" ,"Time_hours", 
                                                     "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_7) <- c("Time",unique(hind_INH_Fblood_medium$Tissue))


# Convert Hinderliter Inhalation female single high from long to wide format using reshape
experiment_inh_8 <- reshape(hind_INH_Fblood_high[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment_inh_8) <- c("Time",unique(hind_INH_Fblood_high$Tissue))


# Put the experiments in a list
experiments_inh <- list(experiment_inh_1 = experiment_inh_1, experiment_inh_2 = experiment_inh_2,
                        experiment_inh_3 = experiment_inh_3, experiment_inh_4 = experiment_inh_4,
                        experiment_inh_5=experiment_inh_5, experiment_inh_6=experiment_inh_6,
                        experiment_inh_7=experiment_inh_7, experiment_inh_8=experiment_inh_8)


# Rename predictions so that they share the same name as the names of the experimental dataframe


colnames(preds_gus_INH_Mblood) <- c ("Time", "Plasma")
colnames(preds_gus_INH_Mtissues) <- c ("Time", "ALF", "Liver", "Lung", "Kidney")
colnames(preds_hind_INH_Mblood_low) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Mblood_medium) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Mblood_high) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Fblood_low) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Fblood_medium) <- c ("Time", "Plasma")
colnames(preds_hind_INH_Fblood_high) <- c ("Time", "Plasma")


# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_gus_INH_Mblood, predictions2 = preds_gus_INH_Mtissues,
                    predictions3 = preds_hind_INH_Mblood_low, predictions4 = preds_hind_INH_Mblood_medium,
                    predictions5 = preds_hind_INH_Mblood_high, predictions6 = preds_hind_INH_Fblood_low,
                    predictions7 = preds_hind_INH_Fblood_medium, predictions8 = preds_hind_INH_Fblood_high)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments_inh)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments_inh[[i]]
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
  ggsave(paste0("experiment_inh", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}
save.image("Inhalation_koff_09_lung_same_RAF_bile.RData")





