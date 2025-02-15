

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
                         'VLuAF','VSP', 'VSPB', 'VSPF','VH', 'VHB', 'VHF','VBr', 'VBrB', 'VBrF','VGo', 'VGoB', 'VGoF',
                         'VIN', 'VINB', 'VINF', 'VST', 'VSTB', 'VSTF','VSTL', 'VINL','VSK','VSKB', 'VSKF',
                         'VBo','VBoB', 'VBoF','VLT', 'VKTrest','VINT', 'VSTT','VMT', 'VAT', 
                         'VLuT', 'VSPT','VHT', 'VBrT','VGoT', 'VSKT','VBoT', 'VRT', 'VUA',
                         'VKT', 'A_peritubular_PTC', 'A_peritubular_DTC','AL', 'AM', 'AA', 'AR', 'ALu','ASP', 'AH', 'ABr', 'AST',
                         'AIN', 'AGo','ASK', 'ABo','AINL', 'AcL' , 'AcM' , 'AcST' ,'AcIN', 'AcA' , 'AcLu' , 'AcALF' , 'AUA', 'ALF',
                         'AcSP', 'AcH' , 'AcBr' , 'AcGo','AcSK', 'AcBo' , 'AcR' , 'APT' , 'ADAL' , 'ADT', 'ACD' , 'AcK_DALC','AcK_CDC' , 'AcKTrest',
                         'Qfeces','Qbile', 'QGFR','Qurine','kabST',"QPT" , "QTDL", "QTAL" , "QDT", "QCD", 'CLfeces', "CL_hepatobiliary", 
                         'VmL_Oatp', 'VmL_Ntcp','VmL_Oatp2',  'VmIn_Oatp2', 'VmK_Oatp','VmLu_Oatp_ap', 'VmLu_Oatp_bas',
                         'VmK_Oat1', 'VmK_Oat3','VmK_Urat','kKTrestF', 'kCdcF' , 'kDalcF' , 'kPtcF' , 'kDtcF' ,
                         'kPtcTu', 'kDalcTu' , 'kDtcTu' , 'kCdcTu' , 'kLFLT',  'kAFAT', 'kRFRT','kMFMT', 'kLuTLuF', 'kLuTLuAF', 'kSPFSPT' ,
                         'kSTFSTT' , 'kINFINT' , 'kHFHT' ,'kBrFBrT' , 'kGoFGoT' ,'kSKFSKT' , 'kBoFBoT', "k_gut_in", 'k_gut_out')
  Q_scaled_params <- c('QBK', 'QBL', 'QBLtot','QBM', 'QBA',
                       'QBR', 'QBLu','QBSP', 'QBH', 'QBBr', 'QBST','QBIN', 'QGE','QBGo','QBSK', 'QBBo', "QparaKi" ,"QparaLi" ,"QparaSt" ,"QparaIn" ,
                       "QparaMu" ,"QparaAd" ,"QparaRe" ,"QparaLu" ,"QparaSp","QparaHt" ,"QparaBr" ,"QparaGo","QparaSk" ,"QparaBo")
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
    CLuTf <- MLuTf/VLuT #tissue concentration
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
    dMLuAFb <-  kon_alb*CalbLuAFf*CLuAFf*VLuAF -  koff_alb*CLuAFb*VLuAF 
    
    #====================================================================================================================
    
    #Arterial Blood
    dMArtf = QBLu*CLuBf - CArtf*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+
                                   QBST+QBIN+QBGo+QBSK+QBBo) - QGFR*CArtf +
      koff_alb*CArtb*VArt - kon_alb*CalbArtf*CArtf*VArt 
    
    #Venous Blood
    dMVenf =  kUAB*CUA  - CVenf*QBLu + QBK*CKBf + QBLtot*CLBf + QBM*CMBf + QBA*CABf + QBR*CRBf+
      QBH*CHBf + QBBr*CBrBf+ QBGo*CGoBf + QBSK*CSKBf + QBBo*CBoBf+
      koff_alb*CVenb*VVen - kon_alb*CalbVenf*CVenf*VVen 
    #Kidney
    #blood subcompartment
    dMKBf = QBK*CArtf - QBK*CKBf   - PeffK*A_peritubular_PTC*(CKBf-CPTCf) -
      PeffK*A_peritubular_DTC*(CKBf-CDTCf) - QparaKi*(1-SKi)*CKBf -
      (VmK_Oat1*CKBf/(KmK_Oat1+CKBf)) - (VmK_Oat3*CKBf/(KmK_Oat3+CKBf))+ (VmK_baso*CPTCf/(KmK_baso+CPTCf))+
      koff_alb*CKBb*VKB - kon_alb*CalbKBf*CKBf*VKB
    
    #interstitial fluid subcompartment
    dMKFf =  - kPtcF*(CKFf-CPTCf) - kDalcF*(CKFf-CDALCf) -
      kDtcF*(CKFf-CDTCf) - kCdcF*(CKFf-CCDCf)   -  kKTrestF*(CKFf-CKTrestf) + 
      koff_alb*CKFb*VKF - kon_alb*CalbKFf*CKFf*VKF
    
    #proximal tubule  cells subcompartment
    dMPTCf =  QparaKi*(1-SKi)*CKBf + PeffK*A_peritubular_PTC*(CKBf-CPTCf) + 
      kPtcF*(CKFf-CPTCf) - kPtcTu*(CPTCf - CPT)  +
      (VmK_Oatp*CPT/(KmK_Oatp+CPT)) + (VmK_Urat*CPT/(KmK_Urat+CPT))+
      (VmK_Oat1*CKBf/(KmK_Oat1+CKBf)) + (VmK_Oat3*CKBf/(KmK_Oat3+CKBf)) - 
      (VmK_baso*CPTCf/(KmK_baso+CPTCf)) -(VmK_api*CPTCf/(KmK_api+CPTCf))-
      (kon_a2u*Ca2uPTCf*CPTCf*VPTC + kon_fabp*CFabpPTCf*CPTCf*VPTC -
         koff_fabp*CPTCb*VPTC - koff_a2u*CPTCb*VPTC) 
    
    #Tubule cells in Loop of Henle 
    dMDALCf =  kDalcF*(CKFf-CDALCf) - kDalcTu*(CDALCf - CDAL)- 
      (kon_a2u*Ca2uDALCf*CDALCf*VDALC + kon_fabp*CFabpDALCf*CDALCf*VDALC -
         koff_fabp*CDALCb*VDALC - koff_a2u*CDALCb*VDALC )
    
    #Distal convoluted tubule cells 
    dMDTCf =   PeffK*A_peritubular_DTC*(CKBf-CDTCf) +kDtcF*(CKFf-CDTCf)- kDtcTu*(CDTCf - CDT)-
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
    dMPT = kPtcTu*(CPTCf - CPT) - (VmK_Oatp*CPT/(KmK_Oatp+CPT)) - 
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
      QBLtot*CLBf - PeffL*AL*(CLBf-CLFf) - QparaLi*(1-SLi)*CLBf +
      koff_alb*CLBb*VLB -kon_alb*CalbLBf*CLBf*VLB 
    #interstitial fluid subcompartment 
    dMLFf =  QparaLi*(1-SLi)*CLBf + PeffL*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) - 
      (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) - (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf)) -
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_alb*CLFb*VLF -kon_alb*CalbLFf*CLFf*VLF 
    #Liver tissue subcompartment
    dMLTf = kLFLT*(CLFf-CLTf) + (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) + (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf))+
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) + koff_fabp*CLTb*VLT-
      kon_fabp*CFabpLTf*CLTf*VLT - CL_hepatobiliary*CLTf
    #Bile  canaliculi subcompartment
    dMLbile = CL_hepatobiliary*CLTf - CLbile*Qbile
    
    
    
    #Stomach
    #blood subcompartment
    dMSTBf = QBST*CArtf - QBST*CSTBf - PeffST*AST*(CSTBf-CSTFf) -  QparaSt*(1-SSt)*CSTBf +
      koff_alb*CSTBb*VSTB-kon_alb*CalbSTBf*CSTBf*VSTB 
    #interstitial fluid subcompartment 
    dMSTFf = QparaSt*(1-SSt)*CSTBf + PeffST*AST*(CSTBf-CSTFf) - kSTFSTT*(CSTFf-CSTTf) +
      koff_alb*CSTFb*VSTF - kon_alb*CalbSTFf*CSTFf*VSTF 
    #Stomach tissue subcompartment
    dMSTTf = kSTFSTT*(CSTFf-CSTTf) + kabST*CSTL 
    #Stomach lumen
    dMSTL = - QGE*CSTL -kabST*CSTL + CLEal*CLuAFf + CLEua*CUA    
    
    #Intestine
    #blood subcompartment
    dMINBf = QBIN*CArtf - QBIN*CINBf - PeffIN*AIN*(CINBf-CINFf) - QparaIn*(1-SIn)*CINBf +
      koff_alb*CINBb*VINB - kon_alb*CalbINBf*CINBf*VINB
    #interstitial fluid subcompartment 
    dMINFf = QparaIn*(1-SIn)*CINBf + PeffIN*AIN*(CINBf-CINFf) - kINFINT*(CINFf-CINTf) +
      koff_alb*CINFb*VINF - kon_alb*CalbINFf*CINFf*VINF
    #Intestine tissue subcompartment
    dMINTf = kINFINT*(CINFf-CINTf) + k_gut_in*CINL - k_gut_out*CINTf + (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL)) 
    #Intestine lumen
    dMINL = QGE*CSTL - (CLfeces*CINL) - k_gut_in*CINL + k_gut_out*CINTf + CLbile*Qbile - 
      (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))
    
    
    #Muscle
    #blood subcompartment
    dMMBf = QBM*CArtf - QBM*CMBf - PeffM*AM*(CMBf-CMFf) - QparaMu*(1-SMu)*CMBf +
      koff_alb*CMBb*VMB - kon_alb*CalbMBf*CMBf*VMB 
    #interstitial fluid subcompartment 
    dMMFf = QparaMu*(1-SMu)*CMBf + PeffM*AM*(CMBf-CMFf) - kMFMT*(CMFf- CMTf) +
      koff_alb*CMFb*VMF - kon_alb*CalbMFf*CMFf*VMF
    #Muscle tissue subcompartment 
    dMMTf = kMFMT*(CMFf- CMTf) 
    
    
    #Adipose
    #blood subcompartment
    dMABf = QBA*CArtf - QBA*CABf - PeffA*AA*(CABf-CAFf) - QparaAd*(1-SAd)*CABf +
      koff_alb*CABb*VAB - kon_alb*CalbABf*CABf*VAB
    #interstitial fluid subcompartment 
    dMAFf = QparaAd*(1-SAd)*CABf + PeffA*AA*(CABf-CAFf) - kAFAT*(CAFf-CATf) +
      koff_alb*CAFb*VAF - kon_alb*CalbAFf*CAFf*VAF 
    #Adipose tissue subcompartment 
    dMATf =  kAFAT*(CAFf-CATf)
    
    
    #Rest of body
    #blood subcompartment
    dMRBf = QBR*CArtf - QBR*CRBf - PeffR*AR*(CRBf-CRFf) - QparaRe*(1-SRe)*CRBf +
      koff_alb*CRBb*VRB - kon_alb*CalbRBf*CRBf*VRB
    #interstitial fluid subcompartment 
    dMRFf = QparaRe*(1-SRe)*CRBf + PeffR*AR*(CRBf-CRFf) - kRFRT*(CRFf -CRTf) +
      koff_alb*CRFb*VRF - kon_alb*CalbRFf*CRFf*VRF 
    #Rest of body tissue subcompartment 
    dMRTf = kRFRT*(CRFf -CRTf) 
    
    #Upper airways
    dMUA = - CLEua*CUA  - kUAB * CUA 
    
    #Lung 
    #blood subcompartment
    dMLuBf = CVenf*QBLu - QBLu*CLuBf - PeffLu*ALu*(CLuBf-CLuFf) - QparaLu*(1-SLu)*CLuBf +
      koff_alb*CLuBb*VLuB - kon_alb*CalbLuBf*CLuBf*VLuB 
    #interstitial fluid subcompartment
    dMLuFf = QparaLu*(1-SLu)*CLuBf + PeffLu*ALu*(CLuBf-CLuFf) + kLuTLuF*(CLuTf-CLuFf) + 
      koff_alb*CLuFb*VLuF - kon_alb*CalbLuFf*CLuFf*VLuF - (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) 
    #Lung tissue
    dMLuTf =  - kLuTLuF*(CLuTf-CLuFf) -  kLuTLuAF*(CLuTf-CLuAFf) + (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) +
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf)) 
    #Alveolar lining fluid
    dMLuAFf =  kLuTLuAF*(CLuTf-CLuAFf) + koff_alb*CLuAFb*VLuAF - kon_alb*CalbLuAFf*CLuAFf*VLuAF -
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf)) - CLEal*CLuAFf
    
    
    #Spleen
    #blood subcompartment
    dMSPBf = QBSP*CArtf - QBSP*CSPBf - PeffSP*ASP*(CSPBf-CSPFf) - QparaSp*(1-SSp)*CSPBf + 
      koff_alb*CSPBb*VSPB - kon_alb*CalbSPBf*CSPBf*VSPB
    #interstitial fluid subcompartment 
    dMSPFf = QparaSp*(1-SSp)*CSPBf + PeffSP*ASP*(CSPBf-CSPFf) - kSPFSPT*(CSPFf -CSPTf) +
      koff_alb*CSPFb*VSPF - kon_alb*CalbSPFf*CSPFf*VSPF 
    #Spleen tissue subcompartment 
    dMSPTf = kSPFSPT*(CSPFf -CSPTf)
    
    
    #Heart
    #blood subcompartment
    dMHBf = QBH*CArtf - QBH*CHBf - PeffH*AH*(CHBf-CHFf) - QparaHt*(1-SHt)*CHBf + 
      koff_alb*CHBb*VHB - kon_alb*CalbHBf*CHBf*VHB 
    #interstitial fluid subcompartment 
    dMHFf = QparaHt*(1-SHt)*CHBf + PeffH*AH*(CHBf-CHFf) - kHFHT*(CHFf -CHTf) + 
      koff_alb*CHFb*VHF - kon_alb*CalbHFf*CHFf*VHF 
    #Heart tissue subcompartment 
    dMHTf = kHFHT*(CHFf -CHTf) 
    
    
    #Brain
    #blood subcompartment
    dMBrBf = QBBr*CArtf - QBBr*CBrBf - PeffBr*ABr*(CBrBf-CBrFf) - QparaBr*(1-SBr)*CBrBf + 
      koff_alb*CBrBb*VBrB - kon_alb*CalbBrBf*CBrBf*VBrB 
    #interstitial fluid subcompartment 
    dMBrFf = QparaBr*(1-SBr)*CBrBf + PeffBr*ABr*(CBrBf-CBrFf) - kBrFBrT*(CBrFf -CBrTf) +
      koff_alb*CBrFb*VBrF - kon_alb*CalBrFf*CBrFf*VBrF 
    #Brain tissue subcompartment 
    dMBrTf = kBrFBrT*(CBrFf -CBrTf) 
    
    
    #Gonads
    #blood subcompartment
    dMGoBf = QBGo*CArtf - QBGo*CGoBf - PeffGo*AGo*(CGoBf-CGoFf) - QparaGo*(1-SGo)*CGoBf +
      koff_alb*CGoBb*VGoB - kon_alb*CalbGoBf*CGoBf*VGoB 
    #interstitial fluid subcompartment 
    dMGoFf = QparaGo*(1-SGo)*CGoBf + PeffGo*AGo*(CGoBf-CGoFf) - kGoFGoT*(CGoFf -CGoTf) +
      koff_alb*CGoFb*VGoF - kon_alb*CalbGoFf*CGoFf*VGoF 
    #gonads tissue subcompartment 
    dMGoTf = kGoFGoT*(CGoFf -CGoTf) 
    
    
    #Skin
    #blood subcompartment
    dMSKBf = QBSK*CArtf - QBSK*CSKBf - PeffSK*ASK*(CSKBf-CSKFf) - QparaSk*(1-SSk)*CSKBf +
      koff_alb*CSKBb*VSKB - kon_alb*CalbSKBf*CSKBf*VSKB 
    #interstitial fluid subcompartment
    dMSKFf = QparaSk*(1-SSk)*CSKBf + PeffSK*ASK*(CSKBf-CSKFf) - kSKFSKT*(CSKFf -CSKTf) +
      koff_alb*CSKFb*VSKF - kon_alb*CalbSKFf*CSKFf*VSKF
    #Skin tissue subcompartment
    dMSKTf = kSKFSKT*(CSKFf -CSKTf) 
    
    
    #Bones
    #blood subcompartment
    dMBoBf = QBBo*CArtf - QBBo*CBoBf - PeffBo*ABo*(CBoBf-CBoFf) - QparaBo*(1-SBo)*CBoBf +
      koff_alb*CBoBb*VBoB - kon_alb*CalbBoBf*CBoBf*VBoB 
    #interstitial fluid subcompartment
    dMBoFf = QparaBo*(1-SBo)*CBoBf + PeffBo*ABo*(CBoBf-CBoFf) - kBoFBoT*(CBoFf -CBoTf) +
      koff_alb*CBoFb*VBoF -  kon_alb*CalbBoFf*CBoFf*VBoF 
    #Bones tissue subcompartment
    dMBoTf = kBoFBoT*(CBoFf -CBoTf) 
    
    #Excreta#
    dMfeces <- CLfeces*CINL
    dMurine <- Qurine*CBladder +  QGFR*CArtf
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
            'dMRBf'=dMRBf, 'dMRFf'=dMRFf,'dMRTf'=dMRTf, 'dMUA'=dMUA,
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
    'CRBb'=CRBb, 'CRF'=CRF, 'CRFf'=CRFf, 'CRFb'=CRFb, 'CRTf'=CRTf, 'CUA'=CUA, 'CLuB'=CLuB,
    'CLuBf'=CLuBf, 'CLuBb'=CLuBb, 'CLuF'=CLuF, 'CLuFf'=CLuFf, 'CLuFb'=CLuFb, 'CLuTf'=CLuTf,
    'CLuAF'=CLuAF, 'CLuAFf'=CLuAFf, 'CLuAFb'=CLuAFb, 'CSPB'=CSPB, 'CSPBf'=CSPBf, 
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
    'CUpperair'=CUpperair, 'CalveolarLF'=CalveolarLF, "Clungtissue" = Clungtissue, 'Clungs'=Clungs,
    'Crest'=Crest, 'Ccarcass'=Ccarcass, 'Cfeces'=Cfeces,
    'Curine'=Curine, 'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain, 
    'Mbrain'=Mbrain, 'Cgonads'=Cgonads, 'Cskin'=Cskin, 'Cbones'=Cbones, 'CalveolarLF' = CalveolarLF,
    "CLbile" = CLbile
    
    )
    
  })
}