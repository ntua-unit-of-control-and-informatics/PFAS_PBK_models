library(deSolve)
getwd()


create.params <- function(user.input){
  with(as.list(user.input),{
    # Cheng and Ng 2017 Table S1
    # Volume of tissue i as percentage of body weight (PVi, unitless) and % volume (Vi, m^3), assuming the density of tissue is 1 g/mL.
    
    #units conversion from Cheng 2017R, time-> h, PFOA mass->ng, tissues mass-> g
    
    #======Table S1=======#    
    
    PVB <- 54e-3 #13.5 mL/244 g=0.055 mL/g~55e-3 mL/g (kg=L)
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 31.2e-3 
    Vplasma <- PVplasma * BW #plasma volume kg=L
    PVK <- 7.3e-3 
    VK <- PVK * BW #kidney volume kg=L
    PVKB <- 0.16 
    VKB <- PVKB * PVK * BW #kidney blood volume kg=L
    PVKF <- 0.13 
    VKF <- PVKF * PVK * BW #kidney interstitial fluid volume kg=L
    VKT <- VK - VKF #kidney tissue volume kg=L
    VFil <- 0.25 #renal filtrate volume kg=L
    PVL <- 3.66e-2 
    VL <- PVL * BW #liver volume kg=L
    PVLB <- 0.21 
    VLB <- PVLB * PVL * BW #liver blood volume kg=L
    PVLF <- 0.049 
    VLF <- PVLF * PVL* BW #liver interstitial fluid volume kg=L
    VLT <- VL - VLF #liver tissue volume kg=L
    PVbile <- 0.004 
    Vbile <- PVbile * PVL * BW #bile volume kg=L
    PVG <- 2.69e-2 
    VG <- PVG * BW #gut volume kg=L
    PVGB <- 0.034
    VGB <- PVGB * PVG * BW #gut blood volume kg=L
    PVGF <- 0.28
    VGF <- PVGF * PVG * BW #gut interstitial fluid volume kg=L
    VGT <- VG - VGF #gut tissue volume kg=L
    PVGL <- 4.5e-2 
    VGL <- PVGL * BW #gut lumen volume kg=L
    PVM <- 40.43e-2 
    VM <- PVM * BW #muscle volume kg=L
    PVMB <- 0.04 
    VMB <- PVMB * PVM * BW #muscle blood volume kg=L
    PVMF <- 0.054
    VMF <- PVMF * PVM * BW #muscle interstitial fluid volume kg=L
    VMT <- VM - VMF #muscle tissue volume kg=L
    PVA <- 7e-2 
    VA <- PVA * BW #adipose volume kg=L
    PVAB <- 0.02 
    VAB <- PVAB * PVA * BW #% adipose blood volume kg=L
    PVAF <- 0.174
    VAF <- PVAF * PVA * BW #adipose interstitial fluid volume kg=L
    VAT <- VA - VAF #adipose tissue volume kg=L
    PVR <- 1 - PVB - PVK - PVL - PVG - PVM - PVA
    VR <- PVR * BW #volume of the rest of the body kg=LL
    PVRB <- 0.036 
    VRB <- PVRB * PVR * BW #volume of the blood of the rest of body kg=L
    PVRF <- 0.18  
    VRF <- PVRF * PVR * BW #interstitial fluid volume of the rest of body kg=L
    VRT <- VR - VRF #tissue volume of the rest of body kg=L
    
    #Capillary surface area for each tissue (Ai) as percentage of body weight
    #or weight of corresponding tissue (PAi, m^2/g) and surface area (m^2)
    
    PAK <- 350e-4 
    AK <- PAK * VK * 1e-3 #kidney surface area (m^2)
    PAKG <- 68.90e-4 
    AKG <- PAKG * VK * 1e-3 #the surface area of glomerular capillary (m^2)
    PAL <- 250e-4 
    AL <- PAL * VL * 1e-3 #liver surface area (m^2)
    PAG <- 100e-4 
    AG <- PAG * VG * 1e-3 #gut surface area (m^2)
    PAGL <- 4.14 #"The surface area of gut lumen would be 4.14 m2/kg" !!!!???
    AGL <- PAGL * BW #gut lumen surface area (m^2)
    PAM <- 70e-4 
    AM <- PAM * VM * 1e-3 #muscle surface area (m^2)
    PAA <- 70e-4
    AA <- PAA * VA * 1e-3 #adipose surface area (m^2)
    PAR <- 100e-4
    AR <- PAR * VR * 1e-3 #surface area of rest of body (m^2)
    
    #Effective permeability (Peff, in m/h) for blood (B), liver(L), kidney(K),
    #gut(G),adipose(A), muscle(M), rest of body(R)
    
    PeffB <- 4.98e-8/3600
    PeffK <- 4.38e-8/3600
    PeffL <- 5.15e-8/3600
    PeffG <- 2.65e-8/3600
    PeffA <- 2.65e-8/3600
    PeffM <- 2.65e-8/3600
    PeffR <- 2.65e-8/3600
    
    
    #Blood flow rates (QBi, in mL/h) to different tissues (i=L, K, G, A, M, R)
    #as a percentage of cardiac output (Qcardiac L/h), which itself is a function
    #of body weight (BW)
    
    Qcardiac <- 0.235 * (BW^0.75) *60 #L/min->*60-> L/h
    PQBK <- 14.1/100 
    QBK <- PQBK * Qcardiac #L/h
    PQBG <- 15.1/100 
    QBG <- PQBG * Qcardiac #L/h
    PQBL <- 2.4/100 
    QBL <- (PQBL+PQBG) * Qcardiac #L/h
    PQBM <- 27.8/100 
    QBM <- PQBM * Qcardiac #L/h
    PQBA <- 7/100 
    QBA <- PQBA * Qcardiac #L/h
    PQBR = 1 - PQBK - PQBG - PQBL - PQBM - PQBA
    QBR = PQBR * Qcardiac #L/h
    
    #Flow rate of fluids including feces, bile, urine and glomerular filtration rate (GFR), in L/h
    
    Qfeces <- 5.63/24 #L per day -> L/h
    PQbile <- 90/24 #L/h/kg BW
    Qbile <- PQbile * BW #L/h
    PQurine <- 200/24 #L/h/kg BW
    Qurine <- PQurine * BW #L/h
    PQGFR <- 10.74*60 #L/h/kg BW
    QGFR <- PQGFR * BW #L/h
    
    
    #======Table S2=======#
    #Albumin concentration in blood and interstitial fluid compartments(mol/m^3 = 1e-6* nmol/g)

    CalbB <- 281e-9*7.8 #n=7.8 binding sites (nmol/g)
    CalbKF <- 243e-9*7.8 #n=7.8 binding sites (nmol/g)
    CalbLF <- 243e-9*7.8 #n=7.8 binding sites (nmol/g)
    CalbGF <- 146e-9*7.8 #n=7.8 binding sites (nmol/g)
    CalbMF <- 146e-9*7.8 #n=7.8 binding sites (nmol/g)
    CalbAF <- 73e-9*7.8 #n=7.8 binding sites (nmol/g)
    CalbRF <- 73e-9*7.8 #n=7.8 binding sites (nmol/g)
    
    
    #Alpha2mu-globulin concentration in kidney tissue (nmol/g)

    Ca2uKT <- 110

    #LFABP concentration in kidney and liver tissue (mol/m^3)

    CLfabpKT <- 2.65 * 3 #n=3 binding sites (nmol/g)
    CLfabpLT <- 133e-3 * 3 #n=3 binding sites (nmol/g)
  

    #======Table S2=======#
    #Equilibrium association constant (m^3/mol= 10^-3*M-1) for albumin(Ka), LFABP(KL_fabp),
    #and alpha2mu-globulin(Ka2u). See SI section S2-2 for details

    Ka <-  24.18    # 3.1*7.8 m3/mol (M-1) multiplying by number of binding sites (Cheng et al. 2021)
    KLfabp <- 135  # 45.0*3 geo_mean of 3 binding affinity, may try normal mean (Cheng et al. 2021)
    Ka2u <- 0.5

    
    #Overall mass transfer coefficients between subcompartments and passive
    #diffusion rate constants. See SI section S3-1 for details
    
    #passive diffusion rates
    
    ClLFT= 3.71e+04 #μL/min/mg protein, Han et al. 2008
    ClKFT= 17.5 #μL/min/mg protein, Yang et al. 2010
    ClGFT= 18.1 #μL/min/mg protein, Kimura et al. 2017
    ClMFT= 17.8 #mean of ClKFT and ClGFT
    ClAFT= 17.8 #mean of ClKFT and ClGFT
    ClRFT= 17.8 #mean of ClKFT and ClGFT
    
      
    #Kidney
    kidney_protein <- 3.02e-06	#amount of protein in proximal tubule cells (mg protein/proximal tubule cell) (Addis 1936)
    kidney_cells <- 2.14e07
    RAFOatp_k <-
    RAFOat1 <-
    RAFOat3 <- 
    
    kBKF <- ((1/QBK) + 1/(PeffB * AK))^(-1)
    kKFKT <- ClKFT*(kidney_protein*kidney_cells)*VKT/60
    n <- 5 #enlargement factor of apical membrane of proximal tubule
    kFKT <- PeffK * AK * n
    kBF <- PeffB * AKG
   
    VmK_Oatp= 9.3*kidney_protein*kidney_cells*RAFOatp_k #nmol/min (Weaver et al. 2010)
    KmK_Oatp= 5.2e-2 #g/mL (Weaver et al. 2010)
    #VmK_Osta= 9.3*2.0e-6*75960 #nmol/min (Weaver et al. 2010 - same as Oatp1)
    #KmK_Osta= 52.3 #mg/mL (Weaver et al. 2010 - mean of Oat1,3,Oatp)
    VmK_Oat1= 2.6*kidney_protein*kidney_cells*RAFOat1 #nmol/min (Weaver et al. 2010)
    KmK_Oat1= 1.8e-2 #g/mL (Weaver et al. 2010)
    VmK_Oat3= 3.8*kidney_protein*kidney_cells*RAFOat3 #nmol/min (Weaver et al. 2010)
    KmK_Oat3= 2.7e-2 #g/mL (Weaver et al. 2010)

      
    #Liver
    liver_protein <- 2.89e-06 #mount of protein in hepatocytes (mg protein/hepatocyte)(Addis 1936)
    liver_cells <- 5.27e08   
    RAFOatp_l <-
    RAFNtcp <-
    
    kBLF <- ((1/QBL) + 1/(PeffB * AL))^(-1)
    kLFLT <- ClLFT*(liver_protein*liver_cells)*VLT/60
    kbileLT <- PeffL * AL
    
    VmL_Oatp= 9.3*4.0e-5*liver_protein*liver_cells*RAFOatp_l #nmol/min, Weaver et al. 2010 - same as kidney
    KmL_Oatp= 5.2e-2 #g/mL, Weaver et al. 2010 - same as kidney
    VmL_Ntcp= 3*4.0e-5*liver_protein*liver_cells*RAFNtcp #nmol/min, Ruggiero et al. 2021
    KmL_Ntcp= 8.3e-3 #g/mL, Ruggiero et al. 2021
    
   
    #Gut
    kBGF <- ((1/QBG) + 1/(PeffB * AG))^(-1)
    gut_protein <- #5034
    kGFGT <- ClGFT*gut_protein*VGT/60
    kGLGT <- PeffG * AGL
    
    #Muscle
    kBMF <- ((1/QBM) + 1/(PeffB * AM))^(-1)
    muscle_protein <- 24368 #mg (protein data from Cheek et al.,1971 and muscle mass from Caster et al.,1956)
    kMFMT <- ClMFT*muscle_protein*VMT/60
    
    #Adipose
    kBAF <- ((1/QBA) + 1/(PeffB * AA))^(-1)
    adipose_protein <- #10507
    kAFAT <- ClAFT*adipose_protein*VAT/60
    
    #Rest of body
    kBRF <- ((1/QBR) + 1/(PeffB *AR))^(-1)
    RoB_protein <- #18456
    kRFRT <- ClRFT*RoB_protein*VRT/60
    
    
    
    return(list('admin.dose'=admin.dose,'PVB'=PVB, 'VB'=VB, 'PVplsma'=PVplasma, 'Vplasma'=Vplasma, 'PVK'=PVK, 'VK'=VK, 'PVKB'=PVKB, 'VKB'=VKB, 'PVKF'=PVKF, 'VKF'=VKF, 'VKT'=VKT, 'VFil'=VFil,
                'PVL'=PVL, 'VL'=VL, 'PVLB'=PVLB, 'VLB'=VLB, 'PVLF'=PVLF, 'VLF'=VLF, 'VLT'=VLT, 'PVbile'=PVbile, 'Vbile'=Vbile, 'PVG'=PVG, 'VG'=VG, 'PVGB'=PVGB, 'VGB'=VGB,
                'PVGF'=PVGF,'VGF'=VGF, 'VGT'=VGT, 'PVGL'=PVGL, 'VGL'=VGL, 'PVM'=PVM, 'VM'=VM, 'PVMB'=PVMB, 'VMB'=VMB, 'PVMF'=PVMF, 'VMF'=VMF, 'VMT'=VMT,
                'PVA'=PVA, 'VA'=VA, 'PVAB'=PVAB, 'VAB'=VAB, 'PVAF'=PVAF, 'VAF'=VAF, 'VAT'=VAT, 'PVR'=PVR, 'VR'=VR, 'PVRB'=PVRB, 'VRB'=VRB, 'PVRF'=PVRF, 'VRF'=VRF, 'VRT'=VRT,
                
                'AKG'=AKG, 'PAL'=PAL, 'AL'=AL, 'PAG'=PAG, 'AG'=AG, 'PAGL'=PAGL, 'AGL'=AGL, 'PAM'=PAM, 'AM'=AM, 'PAA'=PAA, 'AA'=AA, 'PAR'=PAR, 'AR'=AR,
                
                'PeffB'=PeffB, 'PeffK'=PeffK, 'PeffL'=PeffL, 'PeffG'=PeffG, 'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR,
                
                'Qcardiac'=Qcardiac, 'PQBK'=PQBK, 'QBK'=QBK, 'PQBG'=PQBG, 'QBG'=QBG, 'PQBL'=PQBL, 'QBL'=QBL, 'PQBM'=PQBM, 'QBM'=QBM, 'PQBA'=PQBA, 'QBA'=QBA, 'PQBR'=PQBR, 'QBR'=QBR,
                'Qfeces'=Qfeces, 'PQbile'=PQbile, 'Qbile'=Qbile, 'PQurine'=PQurine, 'PQGFR'=PQGFR, 'QGFR'=QGFR,'Qurine'=Qurine, 'QGFR'=QGFR,
                
                'CalbB'=CalbB, 'CalbKF'=CalbKF, 'CalbLF'=CalbLF, 'CalbGF'=CalbGF, 'CalbMF'=CalbMF, 'CalbAF'=CalbAF, 'CalbRF'=CalbRF,'Ca2uKT'=Ca2uKT,
                'CLfabpKT'=CLfabpKT,'CLfabpLT'=CLfabpLT,
                
                'Ka'=Ka, 'Ka2u'=Ka2u, 'KLfabp'=KLfabp,
                
                # 'kidney_protein'=kidney_protein, 'liver_protein'=liver_protein, 'gut_protein'=gut_protein, 'muscle_protein'=muscle_protein, 'adipose_protein'=adipose_protein, 'RoB_protein'=RoB_protein,
                # 'kidney_cells'-kidney_cells, 'liver_cells'=liver_cells,
                
                'ClLFT'=ClLFT, 'ClKFT'=ClKFT, 'ClGFT'=ClGFT, 'ClMFT'=ClMFT, 'ClAFT'=ClAFT, 'ClRFT'=ClRFT,
                
                'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'= KmL_Oatp, 'VmL_Ntcp'= VmL_Ntcp, 'KmL_Ntcp'= KmL_Ntcp,'VmK_Oatp'= VmK_Oatp, 'VmK_Oatp'=VmK_Oatp,'KmK_Oatp'=KmK_Oatp,
                #'VmK_Osta'=VmK_Osta,'KmK_Osta'=KmK_Osta, 
                'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
                  
                
                'kBKF'=kBKF, 'kKFKT'=kKFKT, 'kFKT'=kFKT, 'kBF'=kBF, 'kBLF'=kBLF,'kLFLT'=kLFLT, 'kbileLT'=kbileLT, 'kBAF'=kBAF, 'kAFAT'=kAFAT, 
                'kBRF'=kBRF, 'kRFRT'=kRFRT, 'kBGF'=kBGF, 'kGFGT'=kGFGT,'kGLGT'=kGLGT, 'kBMF'=kBMF, 'kMFMT'=kMFMT
                
                
                
                
    ))
    
    
  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
    
    #====================PFOA mass balance at each tissue or fluid compartment==============================     
      

   
    #Blood subcompartment
    
    dCB = kBKF*(CKFf-CBf) + kBLF*(CLFf-CBf) + kBGF*(CGFf-CBf) + kBMF*(CMFf-CBf) + kBAF*(CAFf-CBf) + kBRF*(CRFf-CBf) + kBF*(CFil-CBf) 
    
    #Kidney
    
    #interstitial fluid subcompartment
    #dCKF = kBKF*(CBf-CKFf) + kKFKT*(CKTf-CKFf) + (VmK_Osta*CKTf/KmK_Osta+CKTf) - (VmK_Oat1*CKFf/KmK_Oat1+CKFf) - (VmK_Oat3*CKFf/KmK_Oat3+CKFf)
    dCKF = kBKF*(CBf-CKFf) + kKFKT*(CKTf-CKFf) - (VmK_Oat1*CKFf/KmK_Oat1+CKFf) - (VmK_Oat3*CKFf/KmK_Oat3+CKFf)
    #Kidney proximal tubule cells subcompartment
    #dCKT = kKFKT*(CKFf-CKTf) + kFKT*(CFil-CKTf) + (VmK_Oatp*CFil/KmK_Oatp+CFil) - (VmK_Osta*CKTf/KmK_Osta+CKTf) + (VmK_Oat1*CKFf/KmK_Oat1+CKFf) + (VmK_Oat3*CKFf/KmK_Oat3+CKFf)
    dCKT = kKFKT*(CKFf-CKTf) + kFKT*(CFil-CKTf) + (VmK_Oatp*CFil/KmK_Oatp+CFil) + (VmK_Oat1*CKFf/KmK_Oat1+CKFf) + (VmK_Oat3*CKFf/KmK_Oat3+CKFf)
    dCFil = kBF*(CBf-CFil) + kFKT*(CKTf-CFil) - (VmK_Oatp*CFil/KmK_Oatp+CFil) - (Qurine/VFil)*CFil
    dCurine = (Qurine/VFil)*CFil
    
    #Liver
    
    #interstitial fluid subcompartment 
    dCLF = kBLF*(CBf-CLFf) + kLFLT*(CLTf-CLFf) - (VmL_Oatp*CLFf/KmL_Oatp+CLFf) - (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf)
    #Liver tissue subcompartment
    dCLT = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/KmL_Oatp+CLFf) + (VmL_Ntcp*CLFf/KmL_Ntcp+CLFf) + kbileLT*(Cbile-CLTf)
    dCbile = kbileLT*(CLTf-Cbile) - (Qbile/Vbile)*Cbile
    
    
    #Gut
    
    #interstitial fluid subcompartment 
    dCGF = kBGF*(CBf-CGFf) + kGFGT*(CGT-CGFf) 
    #Gut tissue subcompartment
    dCGT = kGFGT*(CGFf-CGT) + kGLGT*(CGL-CGT)
    dCGL = kGLGT*(CGT-CGL) + (Qbile/Vbile)*Cbile - (Qfeces/VGL)*CGL
    dCfeces = (Qfeces/VGL)*CGL
    
    #Muscle
    
    #interstitial fluid subcompartment 
    dCMF = kBMF*(CBf-CMFf) + kMFMT*(CMT-CMFf) 
    #Muscle tissue subcompartment 
    dCMT = kMFMT*(CMFf-CMT)
    
    #Adipose
    
    #interstitial fluid subcompartment 
    dCAF = kBAF*(CBf-CAFf) + kAFAT*(CAT-CAFf) 
    #Adipose tissue subcompartment 
    dCAT = kAFAT*(CAFf-CAT)
    
    #Rest of body
    
    #interstitial fluid subcompartment 
    dCRF = kBRF*(CBf-CRFf) + kRFRT*(CRT-CRFf) 
    #Rest of body tissue subcompartment 
    dCRT = kRFRT*(CRFf-CRT)
    
    
    #Lung Tissue subcompartment
    
    #dLNf = bALFLN*MALFf - bLNALF*MLNf + bBLN*MBf - bLNB*MLNf
    
    
    #Cfree calculation using the expression of free fraction ff
    
    CBf = CB * 1.0 / (1.0 + CalbB * Ka)
    CKFf = CKF * 1.0 / (1.0 + CalbKF * Ka)
    CLFf = CLF * 1.0 / (1.0 + CalbLF * Ka)
    CGFf = CGF * 1.0 / (1.0 + CalbGF * Ka)
    CMFf = CMF * 1.0 / (1.0 + CalbMF * Ka)
    CAFf = CAF * 1.0 / (1.0 + CalbAF * Ka)
    CRFf = CRF * 1.0 / (1.0 + CalbRF * Ka)
    CKTf = CKT * 1.0 / (1.0 + Ca2uKT * Ka2u + CLfabpKT * KLfabp)
    CLTf = CLT * 1.0 / (1.0 + CLfabpLT * KLfabp)
    
    
    
    Cblood <- CB * VB /Vplasma
    Ckidney <- ((CB*VKB) + (CKF*VKF) + (CKT*VKT) + (CFil*VFil))/(VB+VKF+VKT+VFil)
    Cliver <- ((CB*VLB) + (CLF*VLF) + (CLT*VLT))/(VB+VLF+VLT)
    Cgut <- ((CB*VGB) + (CGF*VGF) + (CGT*VGT))/(VB+VGF+VGT)
    Cmuscle <- ((CB*VMB) + (CMF*VMF) + (CMT*VMT))/(VB+VMF+VMT)
    Cadipose <- ((CB*VAB) + (CAF*VAF) + (CAT*VAT))/(VB+VAF+VAT)
    Crest <- ((CB*VRB) + (CRF*VRF) + (CRT*VRT))/(VB+VRF+VRT)
    Cfeces <- CGL
    Cbile <- Cbile
    Curine <- CFil
    
    return(list(c('dCB'=dCB,'dCKF'=dCKF,'dCLF'=dCLF,'dCGF'=dCGF,'dCMF'=dCMF,'dCAF'=dCAF, 'dCRF'=dCRF,  
                  'dCAT'=dCAT, 'dCMT'=dCMT, 'dCRT'=dCRT, 'dCKT'=dCKT, 'dCFil'=dCFil, 'dCurine'=dCurine,
                  'dCLT'=dCLT, 'dCbile'=dCbile, 'dCGT'=dCGT, 'dCGL'=dCGL, 'dCfeces'=dCfeces,
                  'CBf'=CBf, 'CKFf'=CKFf, 'CLFf'=CLFf, 'CGFf'=CGFf, 'CMFf'=CMFf, 'CAFf'=CAFf, 'CRFf'=CRFf, 'CKTf'=CKTf, 'CLTf'=CLTf,
                  'CLT'=CLT, 'CGT'=CGT,'CMT'=CMT, 'CAT'=CAT, 'CRT'=CRT, 'CGL'=CGL, 'Cbile'=Cbile, 'CFil'=CFil, 'Cfeces'=Cfeces),
                  'Cblood'=Cblood, 'Ckidney'=Ckidney, 'Cliver'=Cliver, 'Cgut'=Cgut, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose, 'Crest'=Crest,
                  'Cfeces'=Cfeces, 'Cbile'=Cbile, 'Curine'=Curine
               
                
                ))
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    CB <- admin.dose/VB; CKF <-0;  CLF <-0; CGF <-0; CMF <-0; CAF <-0; CRF <-0; CKT <-0; CLT <-0
    CKT <-0; CFil <-0; CLT <-0; CLT <-0; CGT <-0; CMT <-0
    CAT <-0; CRT <-0; CGL <-0; Cbile <-0; CFil <-0; Cfeces <-0; Curine <-0
    CBf <-0; CKFf <-0; CLFf <-0; CGFf <-0; CMFf <-0; CAFf <-0; CRFf <-0; CKTf <-0; CLTf <-0;
    CLT <-0; CGT <-0; CMT <-0; CAT <-0; CRT <-0; CGL <-0; Cbile <-0; CFil <-0; Cfeces <-0 
    
    return(c('CB'=CB,'CKF'=CKF,'CLF'=CLF,'CGF'=CGF,'CMF'=CMF,'CAF'=CAF, 'CRF'=CRF,  
             'CAT'=CAT, 'CMT'=CMT, 'CRT'=CRT, 'CKT'=CKT, 'CFil'=CFil, 'Curine'=Curine,
             'CLT'=CLT, 'Cbile'=Cbile, 'CGT'=CGT, 'CGL'=CGL, 'Cfeces'=Cfeces,
             'CBf'=CBf, 'CKFf'=CKFf, 'CLFf'=CLFf, 'CGFf'=CGFf, 'CMFf'=CMFf, 'CAFf'=CAFf, 'CRFf'=CRFf, 'CKTf'=CKTf, 'CLTf'=CLTf,
             'CLT'=CLT, 'CGT'=CGT,'CMT'=CMT, 'CAT'=CAT, 'CRT'=CRT, 'CGL'=CGL, 'Cbile'=Cbile, 'CFil'=CFil, 'Cfeces'=Cfeces
              
                
             
    ))
    
    
  })
}


################################################################################

BW <- 0.244  # body weight (kg)
admin.dose_per_g <- 1e-6 # administered dose in kg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e09 #ng PFOA
admin.time <- 0 # time when doses are administered, in days
user_input <- list('BW'=BW,
                   "admin.dose"=admin.dose)

params <- create.params(user_input)
inits <- create.inits(params)


sample_time=seq(0,20*24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
rowSums(solution[,c(2:18)])
######################################################################################


# Plot with ggplot2
library (ggplot2)
compartment <- c('Cblood')
color_codes <- scales::hue_pal()(length(compartment))

plot <- ggplot(data = solution)+
  geom_line( aes(x = time, y = Cblood, color='Cblood'), size=1.3)+
  scale_y_log10()+       
  scale_x_continuous()+
 
  labs(title = 'Predicted values',
       y = 'Concentration (ng/g)', x = "Time (hours)")+
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25),
        legend.text=element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
  
  scale_color_manual("Compartment", values=color_codes)+
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))
print(plot)
