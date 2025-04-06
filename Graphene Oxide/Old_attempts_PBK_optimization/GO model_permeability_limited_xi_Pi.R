#Graphene oxide model

library(deSolve)
setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide")


create.params <- function(user.input){
  with(as.list(user.input),{
    
    
    #Volumes of organs as percent of BW
    #blood, kidneys, liver, stomach, Small intestine, large intestine, lungs, 
    #spleen, heart, brain, RoB  
    
    #Blood
    PVB <- 1.7e-3/0.02 #Davies et al. 1993, 1.7 for BW = 0.02 kg
    VB <- PVB * BW #blood volume kg=L
    PVplasma <- 1e-3 /0.02 #Davies et al. 1993, 1.0 for BW = 0.02 kg
    Vplasma <- PVplasma * BW #plasma volume kg=L
    VBven <- BW*0.5/20 	#volume of venous plasma (L); from doi:10.1007/bf02353860
    VBart <- BW*0.22/20	#volume of arterial plasma (L); from doi:10.1007/bf02353860
    
    #Kidney
    PVKi <- 1.67e-2 #Brown et al. 1997, Table 4
    VKi <- PVKi * BW #kidney volume kg=L
    PVKiB <- 0.24 #Brown et al. 1997, Table 30
    VKiB <- PVKiB * PVKi * BW #kidney blood volume kg=L
    
    
    #Liver
    PVLi <- 5.49e-2 #Brown et al. 1997, Table 4
    VLi <- PVLi * BW #liver volume kg=L
    PVLiB <- 0.31  #Brown et al. 1997. Table 30
    VLiB <- PVLiB * PVLi * BW #liver blood volume kg=L
    V_macro_Li = 27.5/100*VLi # https://doi.org/10.3892/etm.2019.7450
    
    #Stomach
    PVSt <- 6e-3 #Brown et al. 1997, Table 4
    VSt <- PVSt * BW #Stomach volume kg=L
    PVStB <- 0.31  #Brown et al. 1997, Table 30
    VStB <- PVStB * PVSt * BW #Stomach blood volume kg=L
    
    #Small intestine
    PVSIn <- 2.53e-2 #Brown et al. 1997, Table 4
    VSIn <- PVSIn * BW #Small intestine volume kg=L
    PVSInB <- 0.02 #pKSim, - 
    VSInB <- PVSInB * BW #Small intestine volume kg=L
   
    #Large intestine
    PVLIn <- 1.09e-2 #Brown et al. 1997, Table 4
    VLIn <- PVLIn * BW #Large intestine volume kg=L
    PVLInB <- 0.02  #pKSim, -
    VLInB <- PVLInB * BW #Large intestine volume kg=L
    
    #length
    L_cecum <- 3.5 #cm
    L_colon_ascendens <- 2.05 #cm
    L_colon_transversum <- 2.05 #cm
    L_colon_descendens <-1.37 #cm
    L_colon_sigmoid <-1.37 #cm
    L_rectum <-1.37 #cm
    
    #radius
    r_cecum <- (0.15+0.27)/2 #cm
    r_colon_ascendens <- 0.15 #cm
    r_colon_transversum <- 0.15 #cm
    r_colon_descendens <-0.15 #cm
    r_colon_sigmoid <-0.15 #cm
    r_rectum <-0.15 #cm
    
    VLInlumen <-(pi*(r_cecum^2*L_cecum+r_colon_ascendens^2*L_colon_ascendens+
                    r_colon_transversum^2*L_colon_transversum+r_colon_descendens^2*L_colon_descendens+
                    r_colon_sigmoid^2*L_colon_sigmoid+r_rectum^2*L_rectum))/1000
    
    #Lung
    PVLn <- 7.3e-3 #Brown et al. 1997, Table 4
    VLn <- PVLn * BW #Lung volume kg=L
    PVLnB <- 0.5 #Brown et al. 1997, Table 30
    VLnB <- PVLnB * BW #Lung volume kg=L
    
    #Spleen
    PVSpl <- 3.5e-3 #Brown et al. 1997, Table 4
    VSpl <- PVSpl * BW #Spleen volume kg=L
    PVSplB <- 0.17 #Brown et al. 1997, Table 30
    VSplB <- PVSplB * BW #Spleen volume kg=L
    V_macro_Spl = 6.94/100*VSpl #https://doi.org/10.1038/s41374-018-0137-1
   
    #Heart
    PVH <- 5e-3 #Brown et al. 1997, Table 4
    VH <- PVH * BW #Heart volume kg=L
    PVHB <- 0.26  #pKSim
    VHB <- PVHB * BW #Heart volume kg=L
   
    #Brain
    PVBr <- 1.65e-2 #Brown et al. 1997, Table 4
    VBr <- PVBr * BW #Brain volume kg=L
    PVBrB <- 0.04 #Brown et al. 1997, Table 30
    VBrB <- PVBrB * BW #Brain volume kg=L
    
    #RoB
    PVRe <- 1 - PVB - PVKi - PVLi - PVSt - PVSIn - PVLIn - PVLn - PVSpl - PVH - PVBr
    VRe <- PVRe * BW #volume of the rest of the body kg=L
    PVReB <- (PVKiB+PVLiB+PVStB+PVSInB+PVLInB+PVLnB+PVSplB+PVHB+PVBrB)/9 #average PViB of all the included organs (kg=L)
    VReB <- PVReB * PVRe * BW #volume of the blood of the rest of body kg=L
  

    #Capillary surface area for each tissue (Ai) as percentage of body weight (m^2/kg), values from pkSim "Endothelial Surface area"
    
    PAKi <- 33.92e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AKi <- PAKi * BW #the surface area of kidney (m^2)
    
    PALi <- 142.0e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALi <- PALi * BW #liver surface area (m^2)
    
    PASt <- 3.34e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASt <- PASt * VSt * 1e3 #stomach surface area (m^2)
    
    PASIn <- 9.62e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASIn <- PASIn * VLIn * 1e3 #Small intestine surface area (m^2)
    
    PALIn <- 5.88e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALIn <- PALIn * VLIn * 1e3 #Small intestine surface area (m^2)
    
    PLn <- 59.47e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ALn <- PLn* BW #lung surface area (m^2)
    
    PSpl <- 26.79e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ASpl <- PSpl* BW #spleen surface area (m^2), same as muscle #assumption
    
    PH <-  23.65e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    AH <- PH* VH #heart surface area (m^2)
    
    PBr <- 5.98e-4/0.02 #m^2/kg, pKSim, Niederalt et al., 2018, https://doi.org/10.1007/s10928-017-9559-4
    ABr <- PBr* VBr *1e03 #brain surface area (m^2)
    
    ARe <- (AKi+ALi+ASt+ASIn+ALIn+ALn+ASpl+AH+ABr)/9 #assumption ???
    
    np_size_small <- 148/2 #nm,  #3/2
    np_size_large <- 556/2 #nm, #300/2
    
    ###############################
    #-----------------------------#
    #   Reflection Coefficients   #
    #-----------------------------#
    ###############################
    # Pore diameters from Price & Gesquiere (2020), doi:https://doi.org/10.1126/sciadv.aax2642
    DpKi <- 200 #nm
    DpLi <- 280 #nm
    DpSt <- 80 #nm, assumption
    DpSIn <- 80 #nm
    DpLIn <- 80 #nm
    DpRe <- 80 #nm, assumption
    DpLn <- 27 #nm
    DpSpl <- 5000 #nm
    DpH <- 50 #nm
    DpBr <- 0.99 #nm
   
    
    Dps <- c(DpKi, DpLi, DpSt, DpSIn, DpLIn, DpLn, DpRe, DpSpl, DpH, DpBr)
    s_r <- rep(NA, length(Dps))
    for (i in 1:length(s_r)){
      a_r <- np_size /(Dps[i]/2)
      Phi = (1-a_r)^2
      
      if (1 - a_r^2 >= 0) {
        F_r <- (((1 - a_r^2)^(3/2)) * Phi) / (1 + 0.2 * (a_r^2) * (1 - a_r^2)^16)
      } else {
        F_r <- 0
      }
      
      G_r <- ((1- (2*a_r^2)/3 - 0.20217*a_r^5 )/ (1-0.75851*a_r^5)) - (0.0431*(1-(1-a_r^10)))
      s_r[i] <- 1-(1-(1-Phi)^2)*G_r+2*a_r^2*Phi*F_r
    }
    SKi <- s_r[1] 
    SLi <- s_r[2]
    SSt <- s_r[3]
    SSIn <- s_r[4]
    SLIn <- s_r[5]
    SLn <- s_r[6]
    SRe <- s_r[7]
    SSpl <- s_r[8]
    SHt <- s_r[9]
    SBr <- 1
   
      
    
    #(QBi, in L/min) to different tissues (i=L, K, G, A, M, R)
    #Hall et al., 2012, https://doi.org/10.1002/jps.22811
    
    BW_ref <- 0.03 #kg
    
    QBgut <- 1.5/1000 * BW/BW_ref # mL/min --> L/min
    QBKi <- 1.3/1000 * BW/BW_ref # mL/min --> L/min
    QBLi <- 1.8/1000 * BW/BW_ref # mL/min --> L/min
    QBSt <- QBgut
    QBSIn <- QBgut
    QBLIn <- QBgut
    QBLn <- 1
    QBSpl <- 0.09/1000 * BW/BW_ref # mL/min --> L/min
    QBH <- 0.28/1000 * BW/BW_ref # mL/min --> L/min
    QBBr <- 0.26/1000 * BW/BW_ref # mL/min --> L/min
    
    Cardiac_output <- 11.4 #mL/min
    Qorgans <- 70.18/100 * Cardiac_output
    Qorgans_scaled <- Qorgans/1000*BW/BW_ref # mL/min --> L/min
    QBRe <- Qorgans_scaled - (QBKi+QBLi+QBSpl+QBH+QBBr+QBgut)
    
    
    QGE<- 1.233/60*BW^0.25 #gastric emptying time (1/(min*BW^0.25)); https://doi.org/10.1002/mrm.10207
    feces_density<- 0.92 #g/cm³ #need to check https://doi.org/10.1038/s41598-022-22626-x
    # Total blood outflow from liver
    QBLitot <- QBLi+QBSpl+QBSIn+QBLIn+QBSt
    
    Qtotal <- QBKi+QBLi+QBRe+QBSpl+QBH+QBBr+QBSIn+QBSt+QBLIn
    Qurine = (5.4+3.7+2.7)/3/60/0.02 #mL/h/kg, unknown BW, https://doi.org/10.1007/BF02035147
    Qfeces = 6.68/20.6 #mg/g --> g/kg BW, https://doi.org/10.3390/toxins10050204
      
    #############################################
    #               Lymph flow rates            #
    #############################################
    #Paracellular flow as a fraction of organ blood flow, L/min
    #from Niederalt et al.(2017). https://doi.org/10.1007/s10928-017-9559-4
    fQparaKi <- 7.09E-4/60
    fQparaLi <- 1.99E-2/60
    fQparaSt <- 2.04E-3/60
    fQparaSIn <- 1.95E-3/60
    fQparaLIn <- 1.44E-2/60
    fQparaRe <- 2.0E-3/60 # Assumption based on 1/500 of flow (Dosgra et al. 2020, https://doi.org/10.1016/j.csbj.2020.02.014 )
    fQparaLn <- 3.56E-5/60
    fQparaSpl <- 1.99E-2/60
    fQparaH <- 1.47E-3/60
    fQparaBr <- 7.27E-5/60
    
    
    #Estimation of lymph flow rates:
    QparaKi <- fQparaKi*QBKi
    QparaLi <- fQparaLi*QBLi
    QparaSt <- fQparaSt*QBSt
    QparaSIn <- fQparaSIn*QBSIn
    QparaLIn <- fQparaLIn*QBLIn
    QparaRe <- fQparaRe*QBRe
    QparaLn <- fQparaLn*QBLn
    QparaSpl <- fQparaSpl*QBSpl
    QparaH <- fQparaH*QBH
    QparaBr <- fQparaBr*QBBr
    
    
    #Surface areas Interstitial - Intracellular PKSim (m^2)
    BW_ref <- 0.02 #kg
    AcKi= 99.71*BW/BW_ref
    AcLi= 17.11*BW/BW_ref
    AcSt= 171.37*BW/BW_ref
    AcSIn= 82.94*BW/BW_ref
    AcLIn= 57.32*BW/BW_ref
    AcLn= 9.11e-3*BW/BW_ref
    AcSpl= 140.76*BW/BW_ref
    AcH= 1.08*BW/BW_ref
    AcBr= 1.04e-4*BW/BW_ref
    AcRe= (AcKi+AcLi+AcSt+AcSIn+AcLIn+AcLn+AcSpl+AcH+AcBr)/9

    
    # #Liver cells
    # liver_cells = 35*1e6 #cells/per g of liver, https://doi.org/10.1016/j.tiv.2006.06.003
    # macrophages = 12.5/100*liver_cells #10-15% of total liver cells, https://doi.org/10.1136/bmjgast-2016-000079
  
    
    return(list('VB'=VB,'Vplasma'=Vplasma,'VBven'=VBven,'VBart'=VBart,
                'VKi'=VKi,'VKiB'=VKiB,
                'VLi'=VLi,'VLiB'=VLiB,'VSt'=VSt,'VStB'=VStB, 
                'VSIn'=VSIn,'VSInB'=VSInB,'VLIn'=VLIn,'VLInB'=VLInB,'VLInlumen'=VLInlumen,
                'VLn'=VLn,'VLnB'=VLnB,'VSpl'=VSpl,'VSplB'=VSplB,'VH'=VH,
                'VHB'=VHB,'VBr'=VBr,'VReB'=VReB,'VRe'=VRe,'VBrB'=VBrB,
                'V_macro_Li'=V_macro_Li, 'V_macro_Spl'=V_macro_Spl,
                
                'AKi'=AKi,'ALi'=ALi,'ASt'=ASt,'ASIn'=ASIn,'ALIn'=ALIn,
                'ALn'=ALn,'ASpl'=ASpl,'AH'=AH,'ABr'=ABr,'ARe'=ARe,
                
                "SKi" = SKi,"SLi" = SLi,"SSt" = SSt,"SSIn" = SSIn,
                "SLIn" = SLIn,"SRe" = SRe,"SLn" = SLn,
                "SSpl" = SSpl,"SHt" = SHt,"SBr" = SBr,
                
                'QBKi'=QBKi, 
                'QBLi'=QBLi, 'QBLitot'=QBLitot,
                'QBRe'=QBRe, 'QBLn'=QBLn, 'QBSpl'=QBSpl, 'QBH'=QBH,
                'QBBr'=QBBr, 'QBSt'=QBSt,'QBSIn'=QBSIn, 'QBLIn'=QBLIn,
                'QGE'=QGE,'feces_density'=feces_density,
                
                "QparaKi" = QparaKi,"QparaLi" = QparaLi,"QparaSt" = QparaSt,"QparaSIn" = QparaSIn,
                "QparaLIn" = QparaLIn, "QparaRe" = QparaRe,"QparaLn" = QparaLn,
                "QparaSpl" = QparaSpl,"QparaH" = QparaH,"QparaBr" = QparaBr,
                
                
                "admin.time" = admin.time, "admin.dose" = admin.dose,
                "admin.type" = admin.type, "np_size"=np_size,
                "np_size_small"=np_size_small, "np_size_large"=np_size_large,
                "Qtotal"=Qtotal,"Qurine"=Qurine,"Qfeces"=Qfeces, 
                "sex"=sex, "estimated_params"=estimated_params
          
                
    
  ))

  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
#coef_i is the ratio of permeability coefficient (xi) between capillary blood and organ
#and partition coefficient of nanoparticles between tissue and blood (Pi), coef_i=xi/Pi
# based on the work of Li et al., 2013, https://doi.org/10.3109/17435390.2013.863406
  
    Pi_liver <- estimated_params[1]
    Pi_spleen <- estimated_params[2]
    Pi_kidney <- estimated_params[3]
    Pi_heart <- estimated_params[4]
    Pi_lung <- estimated_params[5]
    Pi_rob <- estimated_params[6]
    Pi_stomach <- estimated_params[7]
    Pi_smallIn <- estimated_params[8]
    Pi_largeIn <- estimated_params[9]
    Pi_brain <- estimated_params[10]
    
    xi_liver <- estimated_params[11]
    xi_spleen <- xi_liver
    xi_kidney <- xi_liver
    xi_heart <- xi_liver
    xi_lung <- xi_liver
    xi_rob <- xi_liver
    xi_stomach <- xi_liver
    xi_smallIn <- xi_liver
    xi_largeIn <- xi_liver
    xi_brain <- xi_liver
    CLurine <- estimated_params[12]
    CLfeces <- estimated_params[13]
    
      
    # Blood concentration
    CBven <- MBven/VBven
    CBart <- MBart/VBart
    
    # Kidney 
    CKiB = MBKi/VKiB # blood concentration
    CKi = MKi/VKi # tissue concentration
    
    #Liver
    CLiB = MBLi/VLiB # blood concentration
    CLi = MLi/VLi # tissue concentration
    
    
    #Stomach
    CStB = MBSt/VStB # blood concentration
    CSt = MSt/VSt # tissue concentration
    
    #Small Intestine
    CSInB = MBSIn/VSInB # blood concentration
    CSIn = MSIn/VSIn # tissue concentration
    
    #Large Intestine
    CLInB = MBLIn/VLInB # blood concentration
    CLIn = MLIn/VLIn # tissue concentration
    
    #Lungs
    CLnB = MBLn/VLnB # blood concentration
    CLn = MLn/VLn # tissue concentration
    
    #Spleen
    CSplB = MBSpl/VSplB # blood concentration
    CSpl = MSpl/VSpl # tissue concentration
    
    #Heart
    CHB = MBH/VHB # blood concentration
    CH = MH/VH # tissue concentration
    
    #Brain
    CBrB = MBBr/VBrB # blood concentration
    CBr = MBr/VBr # tissue concentration
    
    #Rest-of-the-body
    CReB = MBRe/VReB # blood concentration
    CRe = MRe/VRe # tissue concentration
    
    
    #Arterial Blood
    dMBart = QBLn*CLnB - CBart*(QBKi+QBLi+QBRe+QBSpl+QBH+QBBr+QBSIn+QBSt+QBLIn)
    
    
    #Venous Blood
    dMBven = - QBLn*CBven + QBKi*CKiB + (QBLi+QBSpl+QBSIn+QBLIn+QBSt)*CLiB +
      QBRe*CReB + QBH*CHB + QBBr*CBrB
    
    
    
    #Kidney
    #blood subcompartment
    dMBKi = QBKi*CBart - QBKi*CKiB - QparaKi*(1-SKi)*CKiB - xi_kidney*QBKi*CKiB + (xi_kidney/Pi_kidney)*QBKi*CKi  - CLurine*CKiB
    #tissue subcompartment
    dMKi = QparaKi*(1-SKi)*CKiB + xi_kidney*QBKi*CKiB - (xi_kidney/Pi_kidney)*QBKi*CKi
    
    
    #Liver
    #blood subcompartment
    dMBLi = QBLi*CBart - (QBLi+QBSpl+QBSIn+QBLIn+QBSt)*CLiB + QBSpl*CSplB  + QBSIn*CSInB  + 
      QBLIn*CLInB + QBSt*CStB  - QparaLi*(1-SLi)*CLiB - 
      xi_liver*(QBLi+QBSpl+QBSIn+QBLIn+QBSt)*CLiB + (xi_liver/Pi_liver)*(QBLi+QBSpl+QBSIn+QBLIn+QBSt)*CLi -
      CLfeces*CLiB
    #tissue subcompartment
    dMLi =  QparaLi*(1-SLi)*CLiB + xi_liver*(QBLi+QBSpl+QBSIn+QBLIn+QBSt)*CLiB - 
            (xi_liver/Pi_liver)*(QBLi+QBSpl+QBSIn+QBLIn+QBSt)*CLi
    
    
    #Stomach
    #blood subcompartment
    dMBSt = QBSt*CBart - QBSt*CStB - QparaSt*(1-SSt)*CStB - xi_stomach*QBSt*CStB + (xi_stomach/Pi_stomach)*QBSt*CSt
    #tissue subcompartment 
    dMSt = QparaSt*(1-SSt)*CStB + xi_stomach*QBSt*CStB - (xi_stomach/Pi_stomach)*QBSt*CSt
    #lumen
    dMStlumen = - QGE*MStlumen
    
    #Small Intestine
    #blood subcompartment
    dMBSIn = QBSIn*CBart - QBSIn*CSInB - QparaSIn*(1-SSIn)*CSInB - xi_smallIn*QBSIn*CSInB + (xi_smallIn/Pi_smallIn)*QBSIn*CSIn 
    #tissue subcompartment 
    dMSIn = QparaSIn*(1-SSIn)*CSInB + QGE*MStlumen + xi_smallIn*QBSIn*CSInB - (xi_smallIn/Pi_smallIn)*QBSIn*CSIn 
    
    
    #Large Intestine
    #blood subcompartment
    dMBLIn = QBLIn*CBart - QBLIn*CLInB - QparaLIn*(1-SLIn)*CLInB - xi_largeIn*QBLIn*CLInB + (xi_largeIn/Pi_largeIn)*QBLIn*CLIn   
    #tissue subcompartment 
    dMLIn = QparaLIn*(1-SLIn)*CLInB + xi_largeIn*QBLIn*CLInB - (xi_largeIn/Pi_largeIn)*QBLIn*CLIn  
    
    
    #Lung 
    #blood subcompartment
    dMBLn = QBLn*CBven- QBLn*CLnB - QparaLn*(1-SLn)*CLnB  - xi_lung*QBLn*CLnB + (xi_lung/Pi_lung)*QBLn*CLn 
    #tissue subcompartment
    dMLn = QparaLn*(1-SLn)*CLnB + xi_lung*QBLn*CLnB - (xi_lung/Pi_lung)*QBLn*CLn 
    
    
    #Spleen
    #blood subcompartment
    dMBSpl = QBSpl*CBart - QBSpl*CSplB - QparaSpl*(1-SSpl)*CSplB - xi_spleen*QBSpl*CSplB + (xi_spleen/Pi_spleen)*QBSpl*CSpl 
    #tissue subcompartment 
    dMSpl = QparaSpl*(1-SSpl)*CSplB + xi_spleen*QBSpl*CSplB - (xi_spleen/Pi_spleen)*QBSpl*CSpl 
    
    
    #Heart
    #blood subcompartment
    dMBH = QBH*CBart - QBH*CHB - QparaH*(1-SHt)*CHB - xi_heart*QBH*CHB + (xi_heart/Pi_heart)*QBH*CH 
    #tissue subcompartment 
    dMH = QparaH*(1-SHt)*CHB +  xi_heart*QBH*CHB - (xi_heart/Pi_heart)*QBH*CH 
    
    
    #Brain
    #blood subcompartment
    dMBBr = QBBr*CBart - QBBr*CBrB - QparaBr*(1-SBr)*CBrB - xi_brain*QBBr*CBrB + (xi_brain/Pi_brain)*QBBr*CBr 
    #Tissue subcompartment 
    dMBr = QparaBr*(1-SBr)*CBrB + xi_brain*QBBr*CBrB - (xi_brain/Pi_brain)*QBBr*CBr 
    
    
    #Rest of body
    #blood subcompartment
    dMBRe = QBRe*CBart - QBRe*CReB - QparaRe*(1-SRe)*CReB - xi_rob*QBRe*CReB + (xi_rob/Pi_rob)*QBRe*CRe 
    #interstitial fluid subcompartment 
    dMRe = QparaRe*(1-SRe)*CReB + xi_rob*QBRe*CReB - (xi_rob/Pi_rob)*QBRe*CRe 
    
    
    # Urine
    dMurine = CLurine*CKiB
    # Feces
    dMfeces = CLfeces*CLiB
    
    
    dVurine = Qurine
    dVfeces = Qfeces
   
    
    #Concentration calculation in each compartment 
    Cven <- CBven
    Cart <- CBart
    Cblood <- (MBven + MBart)/ (VBven + VBart)
    Ckidneys <- (MBKi + MKi)/(VKiB + VKi)
    Cliver <- (MBLi + MLi)/(VLiB + VLi)
    Cstomach <-  (MBSt + MSt + MStlumen)/(VStB + VSt)
    Csmall_intestine <-  (MBSIn + MSIn)/(VSInB+VSIn)
    Clarge_intestine <-  (MBLIn + MLIn)/(VLInB+VLIn)
    Clungs <-  (MBLn + MLn)/(VLnB+VLn)
    Crest <-  (MBRe + MRe)/(VReB+VRe)
    Cfeces <- Mfeces/(Vfeces*feces_density)
    Curine <- Murine/Vurine
    Cspleen <-  (MBSpl + MSpl)/(VSplB+VSpl)
    Cheart <-  (MBH + MH)/(VHB+VH)
    Cbrain <-  (MBBr + MBr)/(VBrB+VBr)
    
    
    Mven <- MBven
    Mart <- MBart
    Mblood <- MBven + MBart
    Mkidneys <- MBKi + MKi
    Mliver <- MBLi + MLi 
    Mstomach <- MBSt + MSt + MStlumen
    Msmall_intestine <- MBSIn + MSIn
    Mlarge_intestine <- MBLIn + MLIn
    Mlungs <- MBLn + MLn
    Mrest <- MBRe + MRe
    Mfeces <- Mfeces
    Murine <- Murine
    Mspleen <-  MBSpl + MSpl 
    Mheart <-  MBH + MH
    Mbrain <-  MBBr + MBr
    
    
    list(c( "dMBart"=dMBart, "dMBven"=dMBven, "dMBKi"=dMBKi, "dMKi"=dMKi,
            "dMBLi"=dMBLi, "dMLi"=dMLi, "dMBSt"=dMBSt,
            "dMSt"=dMSt, "dMStlumen"=dMStlumen,
            "dMBSIn"=dMBSIn, "dMSIn"=dMSIn, "dMBLIn"=dMBLIn, "dMLIn"=dMLIn, 
            "dMBLn"=dMBLn, "dMLn"=dMLn, "dMBSpl"=dMBSpl,
            "dMSpl"=dMSpl, "dMBH"=dMBH, "dMH"=dMH, "dMBBr"=dMBBr, "dMBr"=dMBr,
            "dMBRe"=dMBRe, "dMRe"=dMRe, "dMurine"=dMurine, 
            "dMfeces"=dMfeces, "dVurine"=dVurine, "dVfeces"=dVfeces), 
         
           'Cblood'=Cblood, 'Ckidneys'=Ckidneys, 'Cliver'=Cliver, 'Cstomach'=Cstomach,
           'Csmall_intestine'=Csmall_intestine, 'Clarge_intestine'=Clarge_intestine,
           'Clungs'=Clungs, 'Crest'=Crest, 'Cfeces'=Cfeces, 'Curine'=Curine, 
           'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain,
            
           'Mblood'=Mblood, 'Mkidneys'=Mkidneys, 'Mliver'=Mliver, 'Mstomach'=Mstomach,
           'Msmall_intestine'=Msmall_intestine, 'Mlarge_intestine'=Mlarge_intestine,
           'Mlungs'=Mlungs, 'Mrest'=Mrest, 
           'Mspleen'=Mspleen, 'Mheart'=Mheart, 'Mbrain'=Mbrain,
         
           'CBven'=CBven, 'CBart'=CBart, 'CKiB'=CKiB, 'CKi'=CKi, 'CLiB'=CLiB,
           'CLi'=CLi, 'CStB'=CStB, 'CSt'=CSt, 'CSInB'=CSInB,
           'CSIn'=CSIn, 'CLInB'=CLInB,
           'CLIn'=CLIn, 'CLnB'=CLnB, 'CLn'=CLn, 'CSplB'=CSplB,
           'CSpl'=CSpl, 'CHB'=CHB, 'CH'=CH, 'CBrB'=CBrB, 'CBr'=CBr, 'CReB'=CReB, 'CRe'=CRe)
        
            
  
  })
}


create.inits <- function(parameters){
  with(as.list(parameters),{
    
    MBart<-0; MBven<-0; MBKi<-0; MKi<-0; MBLi<-0; MLi<-0; MBSt<-0; MSt<-0;
    MStlumen<-0; MBSIn<-0; MSIn<-0; MBLIn<-0; MLIn<-0; 
    MBLn<-0; MLn<-0; MBSpl<-0; MSpl<-0; MBH<-0; MH<-0; MBBr<-0; MBr<-0;
    MBRe<-0; MRe<-0; Murine<-0; Mfeces<-0; Vurine <-0; Vfeces <-0
    
    
    return(c("MBart"=MBart, "MBven"=MBven, "MBKi"=MBKi, "MKi"=MKi,
             "MBLi"=MBLi, "MLi"=MLi, "MBSt"=MBSt,
             "MSt"=MSt, "MStlumen"=MStlumen,
             "MBSIn"=MBSIn, "MSIn"=MSIn, "MBLIn"=MBLIn, "MLIn"=MLIn, 
             "MBLn"=MBLn, "MLn"=MLn,
             "MBSpl"=MBSpl,"MSpl"=MSpl,"MBH"=MBH, "MH"=MH,
             "MBBr"=MBBr, "MBr"=MBr,
             "MBRe"=MBRe, "MRe"=MRe, "Murine"=Murine,
             "Mfeces"=Mfeces, "Vurine"=Vurine, "Vfeces"=Vfeces
           
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
        events <- list(data = rbind(data.frame(var = c("MBven"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }else if (admin.type == "oral"){
        events <- list(data = rbind(data.frame(var = c("MStlumen"),  time = admin.time, 
                                               value = admin.dose, method = c("add")) ))
      }
    }
    return(events)
  })
}

obj.func <- function(x, dataset){
  N_data <- length(dataset)
  score <- rep(NA, N_data)
  
  # x: a vector with the values of the optimized parameters (it is not the x
  # from the odes!!!)
  estimated_params <- exp(x)
  
  ##########################
  #-------------------------
  # Liu et al., 2012
  #-------------------------
  ##########################
  # Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df1=========================================================
  
  exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  preds_Liu_1_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                  "Mliver","Mspleen", "Mstomach",
                                                                                                  "Mkidneys", "Mlungs", "Mbrain",
                                                                                                  "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_1_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100)
  
  score[1] <- AAFE(predictions = preds_Liu_1_small_tissues, observations = obs_Liu_1_small_tissues)
  
  
  # Set up simulations for the 2nd case, Liu (2012) 1 mg/kg small_tissues_diff_time_points
  
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df2=========================================================
  
  exp_data <- dataset$df2 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mheart", "Mliver","Mspleen", "Mstomach","Mkidneys", "Mlungs", "Mbrain",
                    "Msmall_intestine", "Mlarge_intestine")
  preds_Liu_1_small_diftp_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_small_diftp_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  obs_Liu_1_small_diftp_tissues <- list( exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  
  score[2] <- AAFE(predictions = preds_Liu_1_small_diftp_tissues, observations = obs_Liu_1_small_diftp_tissues)
  
  
  
  # Set up simulations for the 3rd case, i.e. Liu (2012) 1 mg/kg large_tissues_diff_time_points
  
  
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_large <- 556/2
  np_size <- np_size_large #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df3=========================================================
  
  exp_data <- dataset$df3 # retrieve data of Liu et al. 2012 tissues large p.s., 1 mg/kg diff_time_points
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mheart", "Mliver","Mspleen", "Mstomach","Mkidneys", "Mlungs", "Mbrain",
                    "Msmall_intestine", "Mlarge_intestine")
  preds_Liu_1_large_diftp_tissues <- list()
  
  
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_large_diftp_tissues[[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  obs_Liu_1_large_diftp_tissues <- list( exp_data[exp_data$Tissue == "Heart", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Liver", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Spleen", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Stomach", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Kidneys", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Lungs", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Brain", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Small_intestine", "mass"]*admin.dose/100,
                                         exp_data[exp_data$Tissue == "Large_intestine", "mass"]*admin.dose/100) 
  
  score[3] <- AAFE(predictions = preds_Liu_1_large_diftp_tissues, observations = obs_Liu_1_large_diftp_tissues)
  
  
  # Set up simulations for the 4th case, i.e. Liu (2012) 2 mg/kg small_tissues
  
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 2 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  
  #======================================df4=========================================================
  
  exp_data <- dataset$df4 # retrieve data of Liu et al. 2012 tissues Small p.s., 2 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_Liu_2_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                  "Mliver","Mspleen", "Mstomach",
                                                                                                  "Mkidneys", "Mlungs", "Mbrain",
                                                                                                  "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_2_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Liver", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Spleen", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Stomach", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Kidneys", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Lungs", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Brain", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Small_intestine", "concentration"]*admin.dose/100,
                                  exp_data[exp_data$Tissue == "Large_intestine", "concentration"]*admin.dose/100)
  
  score[4] <- AAFE(predictions = preds_Liu_2_small_tissues, observations = obs_Liu_2_small_tissues)
  
  
  # # Set up simulations for the 5th case, Liu (2012) 10 mg/kg small_tissues
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M" 
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,10,0.05) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #======================================df5=========================================================
  
  exp_data <- dataset$df5 # retrieve data of Liu et al. 2012 tissues Small p.s., 10 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "concentration")
  preds_Liu_10_small_tissues <- as.data.frame(solution[solution$time %in% unique(exp_data$time), c("Mheart",
                                                                                                   "Mliver","Mspleen", "Mstomach",
                                                                                                   "Mkidneys", "Mlungs", "Mbrain",
                                                                                                   "Msmall_intestine", "Mlarge_intestine")])
  
  
  
  obs_Liu_10_small_tissues <- list(exp_data[exp_data$Tissue == "Heart", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Liver", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Spleen", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Stomach", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Kidneys", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Lungs", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Brain", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Small_intestine", "concentration"]*admin.dose/100,
                                   exp_data[exp_data$Tissue == "Large_intestine", "concentration"]*admin.dose/100)
  
  score[5] <- AAFE(predictions = preds_Liu_10_small_tissues, observations = obs_Liu_10_small_tissues)
  
  
  # Set up simulations for the 6th case, i.e. Liu (2012) 1 mg/kg small_blood diff_time_points
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_small <- 148/2
  np_size <- np_size_small #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M"
  
  
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df6=========================================================
  
  exp_data <- dataset$df6 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  
  preds_Liu_1_small_diftp_blood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_small_diftp_blood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  
  obs_Liu_1_small_diftp_blood <- list(exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100)
  
  score[6] <- AAFE(predictions = preds_Liu_1_small_diftp_blood, observations = obs_Liu_1_small_diftp_blood)
  
  
  # Set up simulations for the 7th case, i.e. Liu (2012) 1 mg/kg large_blood_diff_time_points
  
  
  BW <- 0.04  # body weight (kg)
  admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
  admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
  np_size_large <- 556/2
  np_size <- np_size_large #nm, Small GO equivalent radius
  admin.time <- 0 # time when doses are administered, in mins
  admin.type <- "iv"
  sex <- "M"
  user_input <- list('BW'=BW,
                     "admin.dose"= admin.dose,
                     "np_size"=np_size,
                     "admin.time" = admin.time, 
                     "admin.type" = admin.type,
                     "estimated_params" = estimated_params,
                     "sex" = sex)
  
  params <- create.params(user_input)
  inits <- create.inits(params)
  events <- create.events(params)
  
  # sample_time: a vector of time points to solve the ODEs
  sample_time=seq(0,180,1) #min
  
  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                      y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  # We need to keep only the predictions for the relevant compartments for the time points 
  # at which we have available data. 
  
  #======================================df7=========================================================
  
  exp_data <- dataset$df7 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
  colnames(exp_data)[c(2,3)] <- c("time", "mass")
  column_names <- c("Mblood")
  
  preds_Liu_1_large_diftp_blood <- list()
  # loop over compartments with available data
  for (i in 1:length(unique(exp_data$Tissue))) {
    compartment <- unique(exp_data$Tissue)[i]
    #Retrieve time points at which measurements are available for compartment i
    exp_time <- exp_data[exp_data$Tissue == compartment, 2]
    
    preds_Liu_1_large_diftp_blood [[i]] <- solution[solution$time %in% exp_time, column_names[i]]
  }
  
  
  
  obs_Liu_1_large_diftp_blood <- list(exp_data[exp_data$Tissue == "Blood", "mass"]*admin.dose/100)
  
  score[7] <- AAFE(predictions = preds_Liu_1_large_diftp_blood, observations = obs_Liu_1_large_diftp_blood)
  
  # Estimate final score
  if (sum(is.na(score))>0){
    final_score <- 100
    
  }else{
    final_score <- mean(score)
    
  }
  return(final_score)
  
}

################################################################################


setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide")

MW <- 124.91 #g/mol
source("Goodness-of-fit-metrics.R")
# Read data
Liu_1_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues.xlsx")
Liu_1_small_diftp_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_tissues_dif_times.xlsx")
Liu_1_large_diftp_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_large_1_tissues_dif_times.xlsx")
Liu_2_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_2_tissues.xlsx")
Liu_10_small_tissues <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_10_tissues.xlsx")
Liu_1_small_diftp_blood <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_small_1_blood_dif_times.xlsx")
Liu_1_large_diftp_blood <- openxlsx::read.xlsx("Data/Liu_2012_GO_male_large_1_blood_dif_times.xlsx")


setwd("/Users/eviepapakyriakopoulou/Documents/GitHub/PFAS_PBK_models/Graphene Oxide/Training/AAFE/GO_model")


dataset <- list("df1" = Liu_1_small_tissues,"df2" = Liu_1_small_diftp_tissues,
                "df3" = Liu_1_large_diftp_tissues,"df4" = Liu_2_small_tissues,
                "df5" = Liu_10_small_tissues,"df6" = Liu_1_small_diftp_blood,
                "df7" = Liu_1_large_diftp_blood)


#Initialise optimiser to NULL for better error handling later
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA"
              "xtol_rel" = 1e-03,
              "ftol_rel" = 0.0,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0, 
              "maxeval" = 2000, 
              "print_level" = 1)

# Create initial conditions (zero initialisation)
#Parameter names:

N_pars <- 13 # Number of parameters to be fitted
fit <-  c(rep(log(1), 13))

lb = c(rep(log(1e-20), 13))
ub = c(rep(log(1e20), 13))


# Run the optimization algorithm to estimate the parameter values
optimizer <- nloptr::nloptr( x0= fit,
                             eval_f = obj.func,
                             # lb	= lb,
                             # ub = ub,
                             opts = opts,
                             dataset = dataset)


estimated_params <- exp(optimizer$solution)
save.image("GO_model.RData")



# Set up simulations for the 1st case, i.e. Liu (2012) 1 mg/kg small_tissues
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df1=========================================================

exp_data <- dataset$df1 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_tissues <- solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                          "Mkidneys", "Mlungs", "Mbrain",
                                          "Msmall_intestine", "Mlarge_intestine")]



# Set up simulations for the 2nd case, Liu (2012) 1 mg/kg small_tissues_diff_time_points
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))
#======================================df2=========================================================

exp_data <- dataset$df2 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_diftp_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                                 "Mkidneys", "Mlungs", "Mbrain",
                                                 "Msmall_intestine", "Mlarge_intestine")]



# Set up simulations for the 3rd case, i.e. Liu (2012) 1 mg/kg large_tissues_diff_time_points

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_large <- 556/2
np_size <- np_size_large #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df3=========================================================

exp_data <- dataset$df3 # retrieve data of Liu et al. 2012 tissues large p.s., 1 mg/kg diff_time_points
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_large_diftp_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                                 "Mkidneys", "Mlungs", "Mbrain",
                                                 "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 4th case, i.e. Liu (2012) 2 mg/kg small_tissues

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 2 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))


#======================================df4=========================================================

exp_data <- dataset$df4 # retrieve data of Liu et al. 2012 tissues Small p.s., 2 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
preds_Liu_2_small_tissues <- solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                          "Mkidneys", "Mlungs", "Mbrain",
                                          "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 5th case, Liu (2012) 10 mg/kg small_tissues
BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 10 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M" 

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,10,0.05) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

#======================================df5=========================================================

exp_data <- dataset$df5 # retrieve data of Liu et al. 2012 tissues Small p.s., 10 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "concentration")
preds_Liu_10_small_tissues <-  solution[, c("time","Mheart","Mliver","Mspleen", "Mstomach",
                                            "Mkidneys", "Mlungs", "Mbrain",
                                            "Msmall_intestine", "Mlarge_intestine")]


# Set up simulations for the 6th case, i.e. Liu (2012) 1 mg/kg small_blood diff_time_points

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_small <- 148/2
np_size <- np_size_small #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #h

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df6=========================================================

exp_data <- dataset$df6 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_small_diftp_blood <- solution[, c("time", "Mblood")]


# Set up simulations for the 7th case, i.e. Liu (2012) 1 mg/kg large_blood_diff_time_points

BW <- 0.04  # body weight (kg)
admin.dose_per_kg <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_kg*BW*1e03 #ug PFOA
np_size_large <- 556/2
np_size <- np_size_large #nm, Small GO equivalent radius
admin.time <- 0 # time when doses are administered, in mins
admin.type <- "iv"
sex <- "M"
user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "np_size"=np_size,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,180,1) #min

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

# We need to keep only the predictions for the relevant compartments for the time points 
# at which we have available data. 

#======================================df7=========================================================

exp_data <- dataset$df7 # retrieve data of Liu et al. 2012 tissues Small p.s., 1 mg/kg
colnames(exp_data)[c(2,3)] <- c("time", "mass")
preds_Liu_1_large_diftp_blood <- solution[, c("time", "Mblood")]



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
         y = expression("GO mass (" * mu* "g)" ),
         x = "Time (min)")+
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


# Convert Liu 2012, male small_1_tissues from long to wide format using reshape
experiment1 <- reshape(Liu_1_small_tissues[c("Tissue" ,"Time_min", 
                                             "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment1) <- c("Time",Liu_1_small_tissues$Tissue )

# Convert Liu 2012, male small_1_tissues_dif_times from long to wide format using reshape
experiment2 <- reshape(Liu_1_small_diftp_tissues[c("Tissue" ,"Time_min", 
                                                   "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment2) <- c("Time",unique(Liu_1_small_diftp_tissues$Tissue))

# Convert Liu 2012, male large_1_tissues_dif_times from long to wide format using reshape
experiment3 <- reshape(Liu_1_large_diftp_tissues[c("Tissue" ,"Time_min", 
                                                   "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment3) <- c("Time",unique(Liu_1_large_diftp_tissues$Tissue))

# Convert Liu 2012, male small_2_tissues from long to wide format using reshape
experiment4 <- reshape(Liu_2_small_tissues[c("Tissue" ,"Time_min", 
                                             "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment4) <- c("Time",Liu_2_small_tissues$Tissue)

# Convert Liu 2012, male small_10_tissues from long to wide format using reshape
experiment5 <- reshape(Liu_10_small_tissues[c("Tissue" ,"Time_min", 
                                              "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",Liu_10_small_tissues$Tissue)

# Convert Liu 2012, male small_1_blood_dif_times from long to wide format using reshape
experiment6 <- reshape(Liu_1_small_diftp_blood[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",unique(Liu_1_small_diftp_blood$Tissue))

# Convert Liu 2012, male large_1_blood_dif_times from long to wide format using reshape
experiment7 <- reshape(Liu_1_large_diftp_blood[c("Tissue" ,"Time_min", 
                                                 "%ID")], 
                       idvar = "Time_min", timevar = "Tissue", direction = "wide")
colnames(experiment7) <- c("Time",unique(Liu_1_large_diftp_blood$Tissue))



# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment3,
                    experiment4 = experiment4, experiment5 = experiment5, experiment6 = experiment6,
                    experiment7 = experiment7)


# Rename predictions so that they share the same name as the names of the experimental data dataframe
colnames(preds_Liu_1_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                         "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine")
colnames(preds_Liu_1_small_diftp_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                               "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_1_large_diftp_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                               "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_2_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                         "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_10_small_tissues) <- c("Time", "Heart", "Liver","Spleen", "Stomach",
                                          "Kidneys", "Lungs", "Brain","Small_intestine", "Large_intestine") 
colnames(preds_Liu_1_small_diftp_blood) <- c("Time", "Blood")
colnames(preds_Liu_1_large_diftp_blood) <- c("Time", "Blood") 


# Create a list containing the corresponding predictions
simulations <- list(predictions1 = preds_Liu_1_small_tissues, predictions2 = preds_Liu_1_small_diftp_tissues,
                    predictions3 = preds_Liu_1_large_diftp_tissues, predictions4 = preds_Liu_2_small_tissues,
                    predictions5 = preds_Liu_10_small_tissues, predictions6 = preds_Liu_1_small_diftp_blood,
                    predictions7 = preds_Liu_1_large_diftp_blood)


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


