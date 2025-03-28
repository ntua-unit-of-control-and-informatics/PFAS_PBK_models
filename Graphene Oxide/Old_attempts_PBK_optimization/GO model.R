#Graphene oxide model

library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{
    
    #Volumes of organs as percent of BW
    #blood, kidneys, liver, stomach, small intestine, large intestine, lungs, 
    #spleen, heart, brain, RoB  
    
    #Blood
    PVB <- 1.7/0.02e-3 #Davies et al. 1993, 1.7 for BW = 0.02 kg
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 1/0.02e-3 #Davies et al. 1993, 1.0 for BW = 0.02 kg
    Vplasma <- PVplasma * BW #plasma volume kg=L
    VVen <- BW*0.5/20 	#volume of venous plasma (L); from doi:10.1007/bf02353860
    VArt <- BW*0.22/20	#volume of arterial plasma (L); from doi:10.1007/bf02353860
    
    #Kidney
    PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
    VKi <- PVKi * BW #kidney volume kg=L
    PVKiB <- 0.24 #Brown et al. 1997, Table 30
    VKiB <- PVKiB * PVKi * BW #kidney blood volume kg=L
    PVKiF <- 0.2 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VKiF <- PVKF * PVKi * BW #kidney interstitial fluid volume kg=L
    VKiT <- VKi - VKiF #kidney tissue volume kg=L
    
      
    #Liver
    PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
    VLi <- PVLi * BW #kidney volume kg=L
    PVLiB <- 0.31  #Brown et al. 1997. Table 30
    VLiB <- PVLiB * PVL * BW #liver blood volume kg=L
    PVLiF <- 0.2 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VLiF <- PVLiF * PVLi * BW #liver interstitial fluid volume kg=L
    VLiT <- VLi - VLiF #liver tissue volume kg=L
    
    #Stomach
    PVSt <- 6e-3 #Brown et al. 1997, Table 4
    VSt <- PVSt * BW #Stomach volume kg=L
    PVStB <- 0.31  #Brown et al. 1997, Table 30
    VStB <- PVStB * PVSt * BW #Stomach blood volume kg=L
    PVStF <- 0.1 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VStF <- PVStF * PVSt * BW #Stomach interstitial fluid volume kg=L
    VStT <- VSt - VStF #Stomach tissue volume kg=L
    
    #Small intestine
    PVSmIn <- 2.53e-2 #Brown et al. 1997, Table 4
    VSmIn <- PVSmIn * BW #Small intestine volume kg=L
    PVSmInB <- 0.02 #pKSim, - 
    VSmInB <- PVSmInB * BW #Small intestine volume kg=L
    PVSmInF <- 0.09 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VSmInF <- PVSmInF * PVSmIn * BW #Small intestine interstitial fluid volume kg=L
    VSmInT <- VSmIn - VSmInF - VSmInB #Small intestine tissue volume kg=L
    
    #Large intestine
    PVLgIn <- 1.09e-2 #Brown et al. 1997, Table 4
    VLgIn <- PVLgIn * BW #Large intestine volume kg=L
    PVLgInB <- 0.02  #pKSim, -
    VLgInB <- PVLgInB * BW #kidney volume kg=L
    PVLgInF <- 0.09 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VLgInF <- PVLgInF * PVLgIn * BW #Large intestine interstitial fluid volume kg=L
    VLgInT <- VLgIn - VLgInF #Large intestine tissue volume kg=L
    
    #Lung
    PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
    VLn <- PVLn * BW #kidney volume kg=L
    PVLnB <- 0.5 #Brown et al. 1997, Table 30
    VLnB <- PVLnB * BW #Lung volume kg=L
    PVLnF <- 0.19 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VLnF <- PVLnF * PVLn * BW #Lung interstitial fluid volume kg=L
    VLnT <- VLn - VLnF #Lung tissue volume kg=L
    
    #Spleen
    PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
    VSpl <- PVSpl * BW #Spleen volume kg=L
    PVSplB <- 0.17 #Brown et al. 1997, Table 30
    VSplB <- PVSplB * BW #Spleen volume kg=L
    PVSplF <- 0.15 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VSplF <- PVSplF * PVSpl * BW #Spleen interstitial fluid volume kg=L
    VSplT <- VSpl - VSplF #Spleen tissue volume kg=L
    
    #Heart
    PVH <- 5e-3 #Brown et al. 1997, Table 4
    VH <- PVH * BW #Heart volume kg=L
    PVHB <- 0.26  #pKSim
    VHB <- PVHB * BW #Heart volume kg=L
    PVHF <- 0.1 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VHF <- PVHF * PVH * BW #Heart interstitial fluid volume kg=L
    VHT <- VH - VHF#Heart tissue volume kg=L
    
    
    #Brain
    PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
    VBr <- PVBr * BW #Brain volume kg=L
    PVBrB <- 0.04 #Brown et al. 1997, Table 30
    VBrB <- PVBrB * BW #Brain volume kg=L
    PVBrF <- 4e-3 #pKSim, Kawai et al., 1994   https://doi.org/10.1007/BF02353860
    VBrF <- PVBrF * PVBr * BW #Brain interstitial fluid volume kg=L
    VBrT <- VBr - VBrF - VBrB #Brain tissue volume kg=L
    
    
    #RoB
    PVR <- 1 - PVB - PVKi - PVLi - PVSt - PVSmIn - PVLgIn - PVLn - PVSpl - PVH - PVBr
    VR <- PVR * BW #volume of the rest of the body kg=L
    PVRB <- PVB - PVKiB - PVLiB - PVStB - PVSmInB - PVLgInB - PVLnB - PVSplB - PVHB - PVBrB
    VRB <- PVRB * PVR * BW #volume of the blood of the rest of body kg=L
  
    
    
    PVRF <- (0.1*0.00158)+((0.14*BW)/0.001)+((0.07*BW)/0.00025)+
            ((0.12*BW)/0.00013)+((0.12*BW)/0.01)+((0.3*BW)/0.0029)  #weighted average VF of bones,fat,gonads, pancreas,muscle,skin
    VRF <- PVRF * PVR * BW #interstitial fluid volume of the rest of body kg=L
    VRT <- VR - VRF #tissue volume of the rest of body kg=L
    
    
    #Capillary surface area for each tissue (Ai) as percentage of body weight (m^2/kg), values from pkSim "Endothelial Surface area"
    
    PAKi <- 33.92e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AKi <- PAKi * BW #the surface area of kidney (m^2)
    PAKiGl <-NA
    AKiGl <- PAKiGl * BW #the surface area of glomerular capillary (m^2)
    PALi <- 142.0e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALi <- PALi * BW #liver surface area (m^2)
    
    PASt <- 3.34e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASt <- PASt * VSt * 1e3 #stomach surface area (m^2)
    PASTL<- 3.34e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASTL<- PASTL * VSmInL #stomach lumen surface area (m^2)
    
    PASmIn <- 9.62e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASmIn <- PASmIn * VLgIn * 1e3 #small intestine surface area (m^2)
    
    PALgIn <- 5.88e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALgIn <- PALgIn * VLgIn * 1e3 #small intestine surface area (m^2)
    

    # #Calculations of mouse small intestinal lumen surface area based on Casteleyn et al. (2010), https://doi.org/10.1258/la.2009.009112
    # # Lengths
    # L_duodenum <- 7# cm
    # L_jejunum <- 32.5 #cm
    # L_ileum <- 4.3# cm
    # # Inner diameters
    # d_duodenum <- 2.9 #cm
    # d_jejunum <- 3 #cm
    # d_ileum <- 2.5 #cm
    # 
    # SA_Sm <- 2*pi*(d_duodenum/2)*L_duodenum + 2*pi*(d_jejunum/2)*L_jejunum + 2*pi*(d_ileum/2)*L_ileum #cm^2
    # n <- 5 #enlargement factor of apical membrane of proximal tubule or enterocytes
    # PASmInL <- n * SA_Sm * 1e-4/0.037 #m^2/kg
    # ASmInL <- PASmInL*BW #m^2
    # 
    # #Calculations of mouse large intestinal lumen surface area based on Casteleyn et al. (2010), https://doi.org/10.1258/la.2009.009112
    # # Length
    # L_colon <- 8.2# cm
    # # Inner diameter
    # d_colon <- 2.9 #cm
    # 
    # SA_Ln <- 2*pi*(d_colon/2)*L_colon  #cm^2
    # n <- 5 #enlargement factor of apical membrane of proximal tubule or enterocytes
    # PALgInL <- n * SA_Ln * 1e-4/0.037 #m^2/kg
    # ALgInL <- PALgInL*BW #m^2
    # 
    
    PLn <- 59.47e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALn <- PLn* BW #lung surface area (m^2)
    PSpl <- 26.79e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASpl <- PSpl* BW #spleen surface area (m^2), same as muscle #assumption
    PH <-  23.65e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AH <- PH* VH #heart surface area (m^2)
    PBr <- 5.98e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ABr <- PBr* VBr *1e03 #brain surface area (m^2)
    PAR <- 10e-4#m2/g tissue, assumption ???
    AR <- PAR * VR *1e3#surface area of rest of body (m^2)
    AR <- (AKi+ALi+ASt+ASmIn+ALgIn+ALn+ASpl+AH+ABr)/9
    
    
    ###############################
    #-----------------------------#
    #   Reflection Coefficients   #
    #-----------------------------#
    ###############################
    # Pore diameters from Price & Gesquiere (2020), doi:https://doi.org/10.1126/sciadv.aax2642
    DpKi <- 200 #nm
    DpLi <- 280 #nm
    DpSt <- 80 #nm, assumption
    DpSmIn <- 80 #nm
    DpLgIn <- 80 #nm
    DpR <- 80 #nm, assumption
    DpLn <- 27 #nm
    DpSpl <- 5000 #nm
    DpH <- 50 #nm
    DpBr <- 0.99 #nm
   
    
    Dps <- c(DpKi, DpLi, DpSt, DpSmIn, DpLgIn, DpR, DpLn, DpSpl, DpH, DpBr)
    s_r <- rep(NA, length(Dps))
    # (C-C) bond length is 0.154 nm ==> 7*0.154 = 1.078nm
    # For carboxyl group we assume 0.13nm, So the total size is around 1.2 nm
    np_size_small <- 3/2 #nm, small GO equivalent radius
    np_size_large <- 300/2 #nm, large GO equivalent radius
    
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
      
    
    #(QBi, in L/h) to different tissues (i=L, K, G, A, M, R)
    #as a percentage of cardiac output (Qcardiac L/h), which itself is a function
    #of body weight (BW)
    
    Qcardiac <- 20.4e-3*60*(BW^0.75)#mL/min->*60/1000-> L/h 
    PQBKi <- 9.1/100 #Brown et al. 1997, p 438, Table 23
    QBKi <- PQBKi * Qcardiac #L/h
    PQBLi <- 16.1/100 
    QBLi <- PQBLi * Qcardiac #L/h Brown et al. 1997, p 438, Table 23
    PQBSt <- NA
    QBSt <- PQBSt * Qcardiac #L/h
    PQBSmIn <- 1e-3
    QBSmIn <- PQBSmIn * Qcardiac #L/h
    PQBLgIn <- NA
    QBLgIn <- PQBLgIn * Qcardiac #L/h
    PQBLn <- 1 
    QBLn <- PQBLn * Qcardiac #L/h
    PQBSpl <- NA 
    QBSpl <- PQBSpl * Qcardiac #L/h
    PQBH <- 6.6/100 #Brown et al. 1997, p 438, Table 23
    QBH <- PQBH * Qcardiac #L/h
    PQBBr <- 3.3/100 #Brown et al. 1997, p 438, Table 23
    QBBr <- PQBBr * Qcardiac #L/h 
    
    
    # Total blood outflow from liver
    QBLitot <- QBLi+QBSpl+QBSmIn+QBLgIn+QBSt
    
    PQBR = 1 - PQBKi - PQBLi - PQBSt - PQBSmIn - PQBLgIn - PQBLn- PQBSpl - PQBH - PQBBr
    QBR <- PQBR * Qcardiac #L/h
    
    
    #############################################
    #               Lymph flow rates            #
    #############################################
    #Paracellular flow as a fraction of organ blood flow, L/h
    #from Niederalt et al.(2017). https://doi.org/10.1007/s10928-017-9559-4
    fQparaKi <- 7.09E-4
    fQparaLi <- 1.99E-2
    fQparaSt <- 2.04E-3
    fQparaSmIn <- 1.95E-3
    fQparaLgIn <- 1.44E-2
    fQparaRe <- 2.0E-3 # Assumption based on 1/500 of flow (Dosgra et al. 2020, https://doi.org/10.1016/j.csbj.2020.02.014 )
    fQparaLn <- 3.56E-5
    fQparaSpl <- 1.99E-2
    fQparaH <- 1.47E-3
    fQparaBr <- 7.27E-5
    
    
    #Estimation of lymph flow rates:
    QparaKi <- fQparaKi*QBKi
    QparaLi <- fQparaLi*QBLi
    QparaSt <- fQparaSt*QBSt
    QparaSmIn <- fQparaSmIn*QBSmIn
    QparaLgIn <- fQparaLgIn*QBLgIn
    QparaRe <- fQparaRe*QBR
    QparaLn <- fQparaLn*QBLn
    QparaSpl <- fQparaSpl*QBSpl
    QparaH <- fQparaH*QBH
    QparaBr <- fQparaBr*QBBr
    
    
    #Surface areas Interstitial - Intracellular PKSim (m^2)
    BW_ref <- 0.02 #kg
    AcKi= 99.71*BW/BW_ref
    AcLi= 17.11*BW/BW_ref
    AcSt= 171.37*BW/BW_ref
    AcSmIn= 82.94*BW/BW_ref
    AcLgIn= 57.32*BW/BW_ref
    AcLn= 9.11e-3*BW/BW_ref
    AcSpl= 140.76*BW/BW_ref
    AcH= 1.08*BW/BW_ref
    AcBr= 1.04e-4*BW/BW_ref
    AcR= (AcKi+AcLi+AcSt+AcSmIn+AcLgIn+AcLn+AcSpl+AcH+AcBr)/9

    
    # Following the calculations  of Lin et al. (2023) for Caco-2 cells
    ClINFT_unscaled= 18.1 #uL/min/mg protein, Kimura et al. 2017
    Awell = 9 #cm^2 (for a 35 mm culture dish)
    Swell = 1.12 #cm^2
    well_protein = 0.346 #mg protein
    P = (well_protein * Awell)/Swell #mg protein/well
    Papp = RAF_papp*(ClINFT_unscaled*60*1e-03*P)/(Awell *2) #cm/h
    P_passive = ( (Papp/100) * AINL)*1000 #L/h
    
    
    #passive diffusion rates
    
    kKiFKiT = ((Papp/100) * AcKi)*1000 #m^3/h * 1000 --> L/h
    kLiFLiT = ((Papp/100) * AcLi)*1000 #m^3/h * 1000 --> L/h
    kStFStT = ((Papp/100) * AcSt)*1000 #m^3/h * 1000 --> L/h 
    kSmInFSmInT = ((Papp/100) * AcSmIn)*1000 #m^3/h * 1000 --> L/h 
    kLgInFLgInT = ((Papp/100) * AcLgIn)*1000 #m^3/h * 1000 --> L/h 
    kLnTLnF = ((Papp/100) * AcLn)*1000 #m^3/h * 1000 --> L/h
    #kLnTLnAF = ((Papp/100) * AcALF)*1000 #m^3/h * 1000 --> L/h
    kSplFSplT = ((Papp/100) * AcSpl)*1000 #m^3/h * 1000 --> L/h 
    kHFHT = ((Papp/100) * AcH)*1000 #m^3/h * 1000 --> L/h 
    kBrFBrT = ((Papp/100) * AcBr)*1000 #m^3/h * 1000 --> L/h 
    kRFRT = ((Papp/100) * AcR)*1000 #m^3/h*1000 --> L/h 
    
    
    #Stomach
    # For identifiability reasons we assume that absorption takes place only through the intestines
    kabST <- (kabs_st* ASTL)*1000 #L/h
    
    #Effective permeability (Peff, in mm/h) for blood (B), liver(L), kidney(K),
    #stomach(ST),intestine (IN), adipose(A), muscle(M), spleen (SP), heart (H), 
    #brain (Br), gonads (Go), rest of body(R)
    
    PeffKi <- Papp*10 #mm/h
    PeffLi <- Papp*10 #mm/h
    PeffSt <- Papp*10 #mm/h
    PeffSmIn <- Papp*10 #mm/h
    PeffLgIn <- Papp*10 #mm/h
    PeffR <- Papp*10 #mm/h
    PeffLn <- Papp*10 #mm/h
    PeffSpl <- Papp*10 #mm/h
    PeffH <- Papp*10 #mm/h
    PeffBr <- Papp*10 #mm/h
    
    
    
    return(list('VB'=VB,'VKi'=VKi,'VKiB'=VKiB,'VKiF'=VKiF,'VKiT'=VKiT,'VFil'=VFil,
                'VLi'=VLi,'VLiB'=VLiB,'VLiF'=VLiF,'VLiT'=VLiT, 
                'VSt'=VSt,'VStB'=VStB,'VStF'=VStF,'VStT'=VStT, 
                'VSmIn'=VSmIn,'VSmInB'=VSmInB,'VSmInF'=VSmInF,'VSmInT'=VSmInT, 
                'VLgIn'=VLgIn,'VLgInB'=VLgInB,'VLgInF'=VLgInF,'VLgInT'=VLgInT,
                'VLn'=VLn,'VLnB'=VLnB,'VLnF'=VLnF,'VLnT'=VLnT,
                'VSpl'=VSpl,'VSplB'=VSplB,'VSplF'=VSplF,'VSplT'=VSplT, 
                'VH'=VH,'VHB'=VHB,'VHF'=VHF,'VHT'=VHT,
                'VBr'=VBr,'VBrB'=VBrB,'VBrF'=VBrF,'VBrT'=VBrT,
                
                'AKi'=AKi,'ALi'=ALi,'ASt'=ASt,'ASmIn'=ASmIn,'ALgIn'=ALgIn,
                'ASmInL'=ASmInL,'ALgInL'=ALgInL,'ALn'=ALn,'ASpl'=ASpl,'AH'=AH,
                'ABr'=ABr,'AR'=AR,
                
                'Qcardiac'=Qcardiac, 'QBKi'=QBKi, 
                'QBLi'=QBLi, 'QBLitot'=QBLitot,
                'QBR'=QBR, 'QBLn'=QBLn, 'QBSpl'=QBSpl, 'QBH'=QBH,
                'QBBr'=QBBr, 'QBSt'=QBSt,'QSmIn'=QSmIn, 'QLgIn'=QLgIn,
                
                "QparaKi" = QparaKi,"QparaLi" = QparaLi,"QparaSt" = QparaSt,"QparaSmIn" = QparaSmIn,
                "QparaLgIn" = QparaLgIn, "QparaR" = QparaR,"QparaLn" = QparaLn,
                "QparaSpl" = QparaSpl,"QparaH" = QparaH,"QparaBr" = QparaBr,
                
                
                'Papp' = Papp, 'P_passive' = P_passive,
                'kKiFKiT'=kKiFKiT, 'kLiFLiT'=kLiFLiT, 
                'kRFRT'=kRFRT,'kabSt'=kabSt, 
                'kLnTLnF' =kLnTLnF,
                'kLnTLnAF'=kLuTLuAF,
                'kSplFSplT' =kSplFSplT, 'kStFStT' =kStFStT, 
                'kSmInFSmInT' =kSmInFSmInT, 'kLgInFLgInT'=kLgInFLgInT,
                'kHFHT' =kHFHT, 'kBrFBrT' =kBrFBrT, 'kRFRT'=kRFRT,
                
                'PeffKi'=PeffKi, 'PeffLi'=PeffLi, 
                'PeffR'=PeffR, 'PeffLn' = PeffLn,'PeffSt' = PeffSt,
                'PeffSmIn' = PeffSmIn, 'PeffLgmIn' = PeffLgmIn,
                'PeffSpl'=PeffSpl, 'PeffH'=PeffH, 'PeffBr'=PeffBr, 
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "MW"=MW
                
                
    
  ))

  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
    # Blood concentration
    CBVen <- Mven/Vven
    CBart <- Mart/Vart
  
    # Kidney 
    CKiB = MBKi/VKiB # blood concentration
    CKiF = MKiF/VKiF #interstitial fluid concentration
    CKiT = MKiT/VKiT # tissue concentration
    
    #Liver
    CLiB = MBiL/VLiB # blood concentration
    CLiF = MiLF/VLiF #interstitial fluid concentration
    CLiT = MLiT/VLiT # tissue concentration
    
    #Stomach
    CStB = MBSt/VStB # blood concentration
    CStF = MStF/VStF  #interstitial fluid concentration
    CStT = MStT/VStT # tissue concentration
    
    #Small Intestine
    CSmInB = MBSmIn/VSmInB # blood concentration
    CSmInF = MSmInF/VSmInF  #Interstitial fluid concentration
    CSmInT = MSmInT/VSmInT # tissue concentration
    
    #Large Intestine
    CLgInB = MBLgIn/VLgInB # blood concentration
    CLgInF = MLgInF/VLgInF  #Interstitial fluid concentration
    CLgInT = MLgInT/VLgInT # tissue concentration

    #Lungs
    CLnB = MBLn/VLnB # blood concentration
    CLnF = MLnF/VLnF  #interstitial fluid concentration
    CLnT = MLnT/VLnT # tissue concentration
    #CLnAF = MLnAF/VLnAF #alveolar lining fluid concentration
    
    #Spleen
    CSplB = MBSpl/VSplB # blood concentration
    CSplF = MSplF/VSplF  #interstitial fluid concentration
    CSplT = MSplT/VSplT # tissue concentration
    
    #Heart
    CHB = MBH/VHB # blood concentration
    CHF = MHF/VHF  #interstitial fluid concentration
    CHT = MHT/VHT # tissue concentration
    
    #Brain
    CBrB = MBBr/VBrB # blood concentration
    CBrF = MBrF/VBrF  #interstitial fluid concentration
    CBrT = MBrT/VBrT # tissue concentration
    
    #Rest-of-the-body
    CRB = MBR/VRB # blood concentration
    CRF = MRF/VRF  #interstitial fluid concentration
    CRT = MRT/VRT # tissue concentration
    
    
    #Arterial Blood
    dMBart = QBLn*CLnB - CBart*(QBKi+QBLi+QBR+QBSpl+QBH+QBBr+QBSmIn+QBSt+QBLgIn) - QGFR*CBfart
    
    #Venous Blood
    dMBven = - QBLn*CBven + QBKi*CKiB + QBLitot*CLiB + QBR*CRB + QBH*CHB + QBBr*CBrB
           
    
    #Kidney
    #blood subcompartment
    dMBKi = QBKi*CBart - QBKi*CKiB - PeffKi*AKi*(CKiB-CKiF) - QparaKi*(1-SKi)*CKiB     
    #interstitial fluid subcompartment
    dMKiF = QparaKi*(1-SKi)*CKiB + PeffKi*AKi*(CKiB-CKiF) - kKiFT*(CKiF-CKiT)
    #Kidney proximal tubule cells subcompartment
    dMKiT = kKiFT*(CKiF-CKiT) - kFKiT*(CKiT - CFil)
    dMFil =  QGFR*CBfart + kFKiT*(CKiT - CFil) - (Qurine*CFil)
    dMurine = Qurine*CFil
    
    
    #Liver
    #blood subcompartment
    dMBLi = QBLi*CBart + QBSpl*CSplB + QBSmIn*CSmInB + QBLgIn*CLgInB + QBSt*CStB - 
      QBLitot*CLBf - PeffLi*ALi*(CLiB-CLiF) - QparaLi*(1-SLi)*CLiB
    #interstitial fluid subcompartment 
    dMLiF =  QparaLi*(1-SLi)*CLiB + PeffLi*ALi*(CLiB-CLiF) - kLiFT*(CLiF-CLiT) 
    #Liver tissue subcompartment
    dMLiT = kLiFT*(CLiF-CLiT) -  P_liver_bile*Qbile*CLiT
    
    # Feces
    dMfeces = Qfeces*CLgInL 
    
    
    #Stomach
    #blood subcompartment
    dMBSt = QBSt*CBart - QBSt*CStB - PeffSt*ASt*(CStB-CStF) -  QparaSt*(1-SSt)*CStB
    #interstitial fluid subcompartment 
    dMSTF = QparaSt*(1-SSt)*CStB + PeffSt*ASt*(CStB-CStF) - kStFT*(CStF-CSTt)
    #Stomach tissue subcompartment
    dMSTT = kStFT*(CStF-CSTt) + kabSt*CStL
    #Stomach lumen
    dMSTL = - QGE*CStL -kabSt*CStL 
    
    
    #Small Intestine
    #blood subcompartment
    dMBSmIn = QBSmIn*CBart - QBSmIn*CSmInB - PeffSmIn*ASmIn*(CSmInB-CSmInF) - QparaSmIn*(1-SSmIn)*CSmInB
    #interstitial fluid subcompartment 
    dMSmInF = QparaSmIn*(1-SSmIn)*CSmInB + PeffSmIn*ASmIn*(CSmInB-CSmInF) - kSmInFT*(CSmInF-CSmInT) 
    #Intestine tissue subcompartment
    dMSmInT = kSmInFT*(CSmInF-CSmInT) + P_passive*CSmInL 
    #Intestine lumen
    dMSmInL = QGE*CStL - P_passive*CSmInL + P_liver_bile*Qbile*CLiT 
      
    
    #Large Intestine
    #blood subcompartment
    dMBLgIn = QBLgIn*CBart - QBLgIn*CLgInB - PeffLgIn*ALgIn*(CLgInB-CLgInF) - QparaLgIn*(1-SLgIn)*CLgInB
    #interstitial fluid subcompartment 
    dMLgInF = QparaLgIn*(1-SLgIn)*CLgInB + PeffLgIn*ALgIn*(CLgInB-CLgInF) - kLgInFT*(CLgInF-CLgInT) 
    #Intestine tissue subcompartment
    dMLgInT = kLgInFT*(CLgInF-CLgInT) 
    #Intestine lumen
    dMLgInL = - (Qfeces*CLgInL)
    
    
    #Lung 
    #blood subcompartment
    dMBLn = CBven*QBLn - QBLn*CLnB - PeffLn*ALn*(CLnB-CLnF) - QparaLn*(1-SLn)*CLnB
    #interstitial fluid subcompartment
    dMLnF = QparaLn*(1-SLn)*CLnB+ PeffLn*ALn*(CLnB-CLnF) + kLnTF*(CLnT-CLnF)
    #Lung tissue
    dMLnT =  - kLnTF*(CLnT-CLnF)  
    #Alveolar lining fluid
    dMLnAF = 0 
    
    
    #Spleen
    #blood subcompartment
    dMBSpl = QBSpl*CBart - QBSpl*CSplB - PeffSpl*ASpl*(CSplB-CSplF) - QparaSpl*(1-SSpl)*CSplB
    #interstitial fluid subcompartment 
    dMSplF = QparaSpl*(1-SSpl)*CSplB + PeffSpl*ASpl*(CSplB-CSplF) - kSplFT*(CSplF -CSplT) 
    #Spleen tissue subcompartment 
    dMSplT = kSplFT*(CSplF -CSplT) 
    
    
    #Heart
    #blood subcompartment
    dMBH = QBH*CBart - QBH*CHB - PeffH*AH*(CHB-CHF) - QparaHt*(1-SHt)*CHB
    #interstitial fluid subcompartment 
    dMHF = QparaHt*(1-SHt)*CHB + PeffH*AH*(CHB-CHF) - kHFT*(CHF -CHT) 
    #Heart tissue subcompartment 
    dMHT = kHFT*(CHF -CHT) 
    
    
    #Brain
    #blood subcompartment
    dMBBr = QBBr*CBart - QBBr*CBrB - PeffBr*ABr*(CBrB-CBrF) - QparaBr*(1-SBr)*CBrB
    #interstitial fluid subcompartment 
    dMBrF = QparaBr*(1-SBr)*CBrB + PeffBr*ABr*(CBrB-CBrF) - kBrFT*(CBrF -CBrT) 
    #Brain tissue subcompartment 
    dMBrT = kBrFT*(CBrF -CBrT)
    
    
    #Rest of body
    #blood subcompartment
    dMBR = QBR*CBart - QBR*CRB - PeffR*AR*(CRB-CRF) - QparaRe*(1-SRe)*CRB
    #interstitial fluid subcompartment 
    dMRF = QparaRe*(1-SRe)*CRB + PeffR*AR*(CRB-CRF) - kRFT*(CRF -CRT) 
    #Rest of body tissue subcompartment 
    dMRT = kRFT*(CRF -CRT)
    
    #Vurine, Vfeces
    dVurine = Qurine
    dVfeces = Qfeces
    
    #Concentration calculation in each compartment 
    Cven <- CBven
    Cart <- CBart
    Cblood <- (MBven + MBart)/ (Vven+Vart)
    Ckidney <- (MBKi + MKiF+ MKiT)/(VKiB+VKiF+VKiT)
    Cliver <- (MBLi + MLiF+ MLiT)/(VLiB+VLiF+VLiT)
    Cstomach <-  (MBST + MSTF+ MSTT)/(VSTB+VSTF+VSTT)
    Cintestine <-  (MBIN + MINF+ MINT+MINL)/(VINB+VINF+VINT+VINL)
    Cmuscle <-  (MBM + MMF+ MMT)/(VMB+VMF+VMT)
    Cadipose <-  (MBA + MAF+ MAT)/(VAB+VAF+VAT)
    Clungs <-  (MBLu + MLuF+ MLuT)/(VLuB+VLuF+VLuT)
    Crest <-  (MBR + MRF+ MRT)/(VRB+VRF+VRT)
    Ccarcass <- (MBM+MMF+MMT+MBA+MAF+MAT+MBR+MRF+MRT+MBBo+MBoF+MBoT+MBSK+MSKF+MSKT)/(VM+VA+VR+VBo+VSK)
    Cfeces <- Mfeces/(Vfeces*feces_density)
    Curine <- Murine/Vurine
    Cspleen <-  (MBSP + MSPF+ MSPT)/(VSPB+VSPF+VSPT)
    Cheart <-  (MBH + MHF+ MHT)/(VHB+VHF+VHT)
    
    Cbrain <-  (MBBr + MBrF+ MBrT)/(VBrB+VBrF+VBrT)
    Mbrain <- MBBr + MBrF+ MBrT
    
    Cgonads <-  (MBGo + MGoF+ MGoT)/(VGoB+VGoF+VGoT)
    Cskin <-  (MBSK + MSKF+ MSKT)/(VSKB+VSKF+VSKT)

    list(c(
    
  )
    )
  })
}