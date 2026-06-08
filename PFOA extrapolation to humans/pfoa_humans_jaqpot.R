library(deSolve)

# =============================================================================
# 1. Parameters
# =============================================================================

create.params <- function(user.input) {
  if (!"depfr_AF"   %in% names(user.input)) user.input$depfr_AF   <- 0
  if (!"depfr_head" %in% names(user.input)) user.input$depfr_head <- 0

  with(as.list(user.input), {

    # -------------------------------------------------------------------------
    # FIXED PARAMS  (anatomy, physiology, volumes, flows, surface areas)
    # -------------------------------------------------------------------------

    #Kidney
    PVKi        <- 0.44e-2
    VKi_tot     <- PVKi * BW
    PVKiB       <- 0.23
    VKiB_cap    <- PVKiB * PVKi * BW
    PVKiF       <- 0.20
    VKiF        <- PVKiF * PVKi * BW
    VPT_1       <- 0.0305; VPT_2 <- 0.0305; VPT_3 <- 0.0305
    VPT_tot     <- VPT_1 + VPT_2 + VPT_3
    VLH_D       <- 0.0027; VLH_A <- 0.0027
    VLH_tot     <- VLH_D + VLH_A
    VDT_tot     <- 0.0194
    RCD         <- 0.1; LCD <- 21; tot_neph_per_CD <- 6; N_CD <- 20*2
    CD_neph     <- N_CD * tot_neph_per_CD
    VCD_tot     <- pi*(RCD^2)*LCD*CD_neph*1e-6
    VFil        <- VPT_tot + VLH_tot + VDT_tot + VCD_tot
    VKiT        <- VKi_tot - VKiF - VKiB_cap
    VPTC        <- VPT_tot; VLHC <- VLH_tot; VDTC <- VDT_tot; VCDC <- VCD_tot
    VBladder    <- 0.350

    #Liver
    PVLi        <- 2.57e-2
    VLi_tot     <- PVLi * BW
    PVLiB       <- 0.17
    VLiB_cap    <- PVLiB * PVLi * BW
    PVLiF       <- 0.16
    VLiF        <- PVLiF * PVLi * BW
    PVBile      <- 29*1e-3/65.26
    VBile       <- PVBile * BW
    VLiT        <- VLi_tot - VLiF - VLiB_cap - VBile

    #Intestine
    PVIn        <- (0.91+0.53)*1e-2
    VIn_tot     <- PVIn * BW
    PVInB       <- 0.02
    VInB_cap    <- PVInB * PVIn * BW
    PVInF       <- 0.09
    VInF        <- PVInF * PVIn * BW
    VInT        <- VIn_tot - VInF - VInB_cap

    #Stomach
    PVSt        <- 0.21e-2
    VSt_tot     <- PVSt * BW
    PVStB       <- 0.03
    VStB_cap    <- PVStB * PVSt * BW
    PVStF       <- 0.10
    VStF        <- PVStF * PVSt * BW
    VStT        <- VSt_tot - VStF - VStB_cap
    PVStL       <- 0.768/64.1
    VStL        <- PVStL * BW
    PVInL       <- (0.894+0.792+0.678+0.598+0.442)/230
    VInL        <- PVInL * BW

    #Muscle
    PVMu        <- 40e-2
    VMu_tot     <- PVMu * BW
    PVMuB       <- 0.03
    VMuB_cap    <- PVMuB * PVMu * BW
    PVMuF       <- 0.16
    VMuF        <- PVMuF * PVMu * BW
    VMuT        <- VMu_tot - VMuF - VMuB_cap

    #Adipose
    PVAd        <- 21.42e-2
    VAd_tot     <- PVAd * BW
    PVAdB       <- 0.02
    VAdB_cap    <- PVAdB * PVAd * BW
    PVAdF       <- 0.16
    VAdF        <- PVAdF * PVAd * BW
    VAdT        <- VAd_tot - VAdF - VAdB_cap

    #Lung
    VUA         <- 29.32e-3
    PVLu        <- 0.76e-2
    VLu_tot     <- PVLu * BW
    PVLuB       <- 0.58
    VLuB_cap    <- PVLuB * PVLu * BW
    PVLuF       <- 0.19
    VLuF        <- PVLuF * PVLu * BW
    PVLuAF      <- 36e-3/70
    VLuAF       <- PVLuAF * BW
    VLuT        <- VLu_tot - VLuF - VLuAF - VLuB_cap

    #Spleen
    PVSp        <- 0.26e-2
    VSp_tot     <- PVSp * BW
    PVSpB       <- 0.33
    VSpB_cap    <- PVSpB * PVSp * BW
    PVSpF       <- 0.15
    VSpF        <- PVSpF * PVSp * BW
    VSpT        <- VSp_tot - VSpF - VSpB_cap

    #Heart
    PVHt        <- 0.47e-2
    VHt_tot     <- PVHt * BW
    PVHtB       <- 0.14
    VHtB_cap    <- PVHtB * PVHt * BW
    PVHtF       <- 0.10
    VHtF        <- PVHtF * PVHt * BW
    VHtT        <- VHt_tot - VHtF - VHtB_cap

    #Brain
    PVBr        <- 2e-2
    VBr_tot     <- PVBr * BW
    PVBrB       <- 0.04
    VBrB_cap    <- PVBrB * PVBr * BW
    PVBrF       <- 4e-3
    VBrF        <- PVBrF * PVBr * BW
    VBrT        <- VBr_tot - VBrF - VBrB_cap

    #Gonads
    PVGo        <- 0.04/70
    VGo_tot     <- PVGo * BW
    PVGoB       <- 0.06
    VGoB_cap    <- PVGoB * PVGo * BW
    PVGoF       <- 0.07
    VGoF        <- PVGoF * PVGo * BW
    VGoT        <- VGo_tot - VGoF - VGoB_cap

    #Skin
    PVSk        <- 3.71e-2
    VSk_tot     <- PVSk * BW
    PVSkB       <- 0.02
    VSkB_cap    <- PVSkB * PVSk * BW
    PVSkF       <- 0.3
    VSkF        <- PVSkF * PVSk * BW
    VSkT        <- VSk_tot - VSkF - VSkB_cap

    #Bones
    PVBo        <- 11.82/70
    VBo_tot     <- PVBo * BW
    PVBoB       <- 0.03
    VBoB_cap    <- PVBoB * PVBo * BW
    PVBoF       <- 0.1
    VBoF        <- PVBoF * PVBo * BW
    VBoT        <- VBo_tot - VBoF - VBoB_cap

    #Blood
    if (sex == "M") {
      PVB    <- 78e-3
    } else {
      PVB    <- 56e-3
    }
    VBlood      <- PVB * BW
    PVRe        <- 1 - PVB - PVKi - PVLi - PVMu - PVAd - PVSp - PVHt - PVBr - PVGo - PVSt - PVIn - PVSk - PVBo
    VRe_tot     <- PVRe * BW
    PVReB       <- (PVKiB+PVLiB+PVLuB+PVMuB+PVAdB+PVSpB+PVHtB+PVBrB+PVGoB+PVStB+PVInB+PVSkB+PVBoB)/13
    VReB_cap    <- PVReB * PVRe * BW
    PVReF       <- (PVKiF+PVLiF+PVLuF+PVMuF+PVAdF+PVSpF+PVHtF+PVBrF+PVGoF+PVStF+PVInF+PVSkF+PVBoF)/13
    VReF        <- PVReF * PVRe * BW
    VReT        <- VRe_tot - VReF - VReB_cap
    Vcap_tot    <- VKiB_cap + VLiB_cap + VInB_cap + VStB_cap + VMuB_cap + VAdB_cap + VLuB_cap +
                   VSpB_cap + VHtB_cap + VBrB_cap + VGoB_cap + VSkB_cap + VBoB_cap + VReB_cap
    VB_central  <- VBlood - Vcap_tot
    f_ven       <- 3/4; f_art <- 1/4
    VVen        <- f_ven * VB_central
    VArt        <- f_art * VB_central

    #Surface areas (capillary)
    BW_ref                  <- 70
    linear_scaling_factor   <- BW / BW_ref
    nonlinear_scaling_factor <- (BW / BW_ref)^0.75

    AL   <- 38.05  * linear_scaling_factor
    AST  <- 0.51   * linear_scaling_factor
    AIN  <- (1.32+0.75) * linear_scaling_factor
    AM   <- 75.95  * linear_scaling_factor
    AA   <- 21.26  * linear_scaling_factor
    AR   <- 10     * linear_scaling_factor
    ALu  <- 66.52  * linear_scaling_factor
    ASP  <- 6.47   * linear_scaling_factor
    AH   <- 5.55   * linear_scaling_factor
    ABr  <- 5.59   * linear_scaling_factor
    AGo  <- 0.21   * linear_scaling_factor
    ASK  <- 16.09  * linear_scaling_factor
    ABo  <- 38.17  * linear_scaling_factor

    #Gut effective surface area
    Duodenum          <- 628.32; Upper_jejunum <- 510.06; Lower_jejunum <- 477.16
    Upper_ileum       <- 635.91; Lower_ileum   <- 538.08; Cecum          <- 126.78
    colon_ascendens   <- 354.52; colon_transversum <- 547.72
    colon_descendens  <- 265.66; colon_sigmoid <- 328.74
    EF_Duodenum       <- 292.69; EF_Upper_jejunum <- 447.99; EF_Lower_jejunum <- 372.94
    EF_Upper_ileum    <- 260.75; EF_Lower_ileum   <- 146.57; EF_Cecum          <- 1.80
    EF_colon_ascendens <- 2.5;   EF_colon_transversum <- 2.5
    EF_colon_descendens <- 2.5;  EF_colon_sigmoid <- 3.56
    AINL <- (Duodenum*EF_Duodenum + Upper_jejunum*EF_Upper_jejunum +
               Lower_jejunum*EF_Lower_jejunum + Upper_ileum*EF_Upper_ileum +
               Lower_ileum*EF_Lower_ileum + Cecum*EF_Cecum +
               colon_ascendens*EF_colon_ascendens + colon_transversum*EF_colon_transversum +
               colon_descendens*EF_colon_descendens + colon_sigmoid*EF_colon_sigmoid) * 1e-4

    #Surface areas interstitial-intracellular
    AcK_total <- 22228.45 * nonlinear_scaling_factor
    AcK_PTC   <- AcK_total * VPTC / VFil
    AcK_LHC   <- AcK_total * VLHC / VFil
    AcK_DTC   <- AcK_total * VDTC / VFil
    AcK_CDC   <- AcK_total * VCDC / VFil
    AcL       <- 4930.84   * nonlinear_scaling_factor
    AcST      <- 43538.81  * nonlinear_scaling_factor
    AcIN      <- (19458.36+12682.70) * nonlinear_scaling_factor
    AcM       <- 530.51    * nonlinear_scaling_factor
    AcA       <- 804.69    * nonlinear_scaling_factor
    AcLu      <- 10.89     * nonlinear_scaling_factor
    AcSP      <- 44741.62  * nonlinear_scaling_factor
    AcH       <- 606.98    * nonlinear_scaling_factor
    AcBr      <- 0.1       * nonlinear_scaling_factor
    AcGo      <- 16.08     * nonlinear_scaling_factor
    AcSK      <- 3.35      * nonlinear_scaling_factor
    AcBo      <- 926.21    * nonlinear_scaling_factor
    AcR       <- median(c(AcK_total, AcL, AcST, AcIN, AcM, AcA, AcLu, AcSP,
                           AcH, AcBr, AcGo, AcSK, AcBo))
    AcALF     <- 0.4 * (BW / 0.363)
    Nasal_SA  <- 18.5e-4

    #Tubule filtrate surface areas
    APT_tot   <- (6107 + 6107 + 6107) * 1e-2
    ALH_tot   <- (61 + 61) * 1e-2
    ADT_tot   <- 156 * 1e-2
    ACD_tot   <- (6.7*5) * 1e-2
    AFil      <- APT_tot + ALH_tot + ADT_tot + ACD_tot

    #Peritubular capillary surface area
    d_peritubular          <- 25.9
    PTC_diameter           <- 14.5
    PTC_mean_capillary_length <- 75
    A_PTC_mean             <- pi * PTC_diameter * PTC_mean_capillary_length
    N_PTC_total            <- d_peritubular * (AFil * 1e12 / 50000)
    A_peritubular          <- N_PTC_total * A_PTC_mean * 1e-12

    #Blood flows
    Qcardiac <- 5200/1000 * (BW^0.75) * 60
    QBK      <- 0.19  * Qcardiac
    QBL      <- 0.065 * Qcardiac
    QBST     <- 0.01  * Qcardiac
    QBIN     <- 0.05  * Qcardiac
    QBM      <- 0.17  * Qcardiac
    QBA      <- 0.052 * Qcardiac
    QBLu     <- 1     * Qcardiac
    QBSP     <- 0.02  * Qcardiac
    QBH      <- 0.04  * Qcardiac
    QBBr     <- 0.114 * Qcardiac
    QBGo     <- 0.005 * Qcardiac
    QBSK     <- 0.058 * Qcardiac
    QBBo     <- 0.042 * Qcardiac
    QBLtot   <- QBL + QBSP + QBIN + QBST
    PQBR     <- 1 - 0.19 - 0.065 - 0.01 - 0.05 - 0.17 - 0.052 - 0.04 - 0.058 - 0.02 - 0.005 - 0.114 - 0.042
    QBR      <- PQBR * Qcardiac

    #Other flows
    Qbile          <- 900e-3 / 24
    Vbile          <- 55e-3
    Qfeces         <- 128e-3 / 24
    feces_density  <- 1.06
    QGFR           <- 116.66 * 60 / 1000
    if (sex == "M") {
      Qurine <- 0.212 * VKi_tot
    } else {
      Qurine <- 0.152 * VKi_tot
    }
    Qelim <- 28e-3 / 24
    QGE   <- 3.5 * BW^(-0.25)

    #Tubular filtrate flows (scaled from Huang & Isoherranen 2018)
    QPT_tot_ref <- (120 + 94 + 68) * 1e-3 * 60
    QLH_tot_ref <- (43 + 24)       * 1e-3 * 60
    QDT_tot_ref <- 24              * 1e-3 * 60
    QCD_tot_ref <- (11+9+7+5+3)    * 1e-3 * 60
    Q_scaling   <- QGFR / QPT_tot_ref
    QPT_tot     <- Q_scaling * QPT_tot_ref
    QLH_tot     <- Q_scaling * QLH_tot_ref
    QDT_tot     <- Q_scaling * QDT_tot_ref
    QCD_tot     <- Q_scaling * QCD_tot_ref

    # -------------------------------------------------------------------------
    # ESTIMATED PARAMS  (fitted; set 1e-05 as placeholder)
    # -------------------------------------------------------------------------
    # Estimated parameters are retreived from "PFOA extrapolation to humans/Trials/Trial9_urine and feces cut data/Final_model_AAFE_humans.RData"
    estimated_params <- c(0.008039367, 0.631940503, 0.001779902)

    # -------------------------------------------------------------------------
    # VARIABLE PARAMS  (transport kinetics)
    # -------------------------------------------------------------------------
    kabsUA            <- 0.1
    kCLEua            <- 0
    RAFOatp_lu_ap     <- 0
    RAFOatp_lu_bas    <- 0.42333523
    k_desorption_fast <- 10.86622465
    k_desorption_slow <- 0.01364133
    f_fast            <- plogis(0.34862342)
    fdust_retained    <- plogis(0.09039221)

    kUAB  <- kabsUA * Nasal_SA
    CLEua <- kCLEua * Nasal_SA

    RAFOatp_k <- estimated_params[1]
    RAFOat3   <- estimated_params[2]
    RAFOatp_l <- 4.992703e-01
    RAFUrat   <- RAFOatp_k
    RAFOat1   <- RAFOat3
    RAFOatp2_l   <- RAFOatp_l
    RAFNtcp      <- RAFOatp_l
    RAFOatp2_Int <- 6.158004e-06

    VmK_api   <- 0; VmK_baso <- 0; KmK_baso <- 1e20; KmK_api <- 5e4

    KLfabp <- 1.2e5
    Ka     <- 5e5

    CLfeces_unscaled <- estimated_params[3]
    CLfeces          <- CLfeces_unscaled * BW^(-0.25)

    f_alb_avail  <- 1; f_fabp_avail <- 1
    koff_alb     <- 1e-2 * 3600
    koff_fabp    <- 1e-2 * 3600

    reduction_factor <- 3.714130e-02

    f_tubular            <- 0.52
    f_PTC_prot_to_total_prot <- 0.65

    MW <- 414.07

    # Kidney transporters
    kidney_protein_total <- 17e-2 * (1e6 * VKi_tot)
    PTC_protein          <- f_PTC_prot_to_total_prot * kidney_protein_total

    VmK_Oatp_in_vitro <- 9.3
    VmK_Oatp_scaled   <- 60 * VmK_Oatp_in_vitro * MW * PTC_protein / 1000
    VmK_Oatp          <- VmK_Oatp_scaled * RAFOatp_k
    KmK_Oatp          <- 126.4 * MW

    VmK_Oat1_in_vitro <- 2.6
    VmK_Oat1_scaled   <- 60 * VmK_Oat1_in_vitro * MW * PTC_protein / 1000
    VmK_Oat1          <- VmK_Oat1_scaled * RAFOat1
    KmK_Oat1          <- 43.2 * MW

    VmK_Oat3_in_vitro <- 3.8
    VmK_Oat3_scaled   <- 60 * VmK_Oat3_in_vitro * MW * PTC_protein / 1000
    VmK_Oat3          <- VmK_Oat3_scaled * RAFOat3
    KmK_Oat3          <- 65.7 * MW

    VmK_Urat_in_vitro <- 1520e-3
    VmK_Urat_scaled   <- 60 * VmK_Urat_in_vitro * MW * PTC_protein / 1000
    VmK_Urat          <- VmK_Urat_scaled * RAFUrat
    KmK_Urat          <- 820.04 * MW

    # Liver transporters
    liver_protein_total  <- 18e-2 * (1e6 * VLi_tot)

    VmL_Oatp_in_vitro <- 9.3
    VmL_Oatp_scaled   <- 60 * VmL_Oatp_in_vitro * MW * liver_protein_total / 1000
    VmL_Oatp          <- VmL_Oatp_scaled * RAFOatp_l
    KmL_Oatp          <- KmK_Oatp

    VmL_Oatp2_in_vitro <- 1493e-3
    VmL_Oatp2_scaled   <- 60 * VmL_Oatp2_in_vitro * MW * liver_protein_total / 1000
    VmL_Oatp2          <- VmL_Oatp2_scaled * RAFOatp2_l
    KmL_Oatp2          <- 148.68 * MW

    VmL_Ntcp_in_vitro <- 3
    VmL_Ntcp_scaled   <- 60 * VmL_Ntcp_in_vitro * MW * liver_protein_total / 1000
    VmL_Ntcp          <- VmL_Ntcp_scaled * RAFNtcp
    KmL_Ntcp          <- 20 * MW

    CL_int           <- 1.961917e+00
    HEPGL            <- 139
    CL_hepatobiliary <- CL_int * 1e-6 * HEPGL * (VLi_tot * 1000) * 60

    # Lung transporters
    lung_protein_total     <- 17.8e-2 * (1e6 * VLu_tot)
    VmLu_Oatp_ap_in_vitro <- 9.3
    VmLu_Oatp_ap_scaled   <- 60 * VmLu_Oatp_ap_in_vitro * MW * lung_protein_total / 1000
    VmLu_Oatp_ap          <- VmLu_Oatp_ap_scaled * RAFOatp_lu_ap
    KmLu_Oatp_ap          <- KmK_Oatp

    VmLu_Oatp_bas_in_vitro <- 9.3
    VmLu_Oatp_bas_scaled   <- 60 * VmLu_Oatp_bas_in_vitro * MW * lung_protein_total / 1000
    VmLu_Oatp_bas          <- VmLu_Oatp_bas_scaled * RAFOatp_lu_bas
    KmLu_Oatp_bas          <- KmK_Oatp

    # Intestine transporters
    intestine_protein_total <- 13e-2 * (1e6 * VIn_tot)
    VmIn_Oatp2_in_vitro    <- 456.63e-3
    VmIn_Oatp2_scaled      <- 60 * VmIn_Oatp2_in_vitro * MW * intestine_protein_total / 1000
    VmIn_Oatp2             <- VmIn_Oatp2_scaled * RAFOatp2_Int
    KmIn_Oatp2             <- 8.3 * MW

    # Albumin
    Mr_albumin   <- 66500; Mr_fabp <- 12000
    CalbB_init   <- f_alb_avail * 4.25 * 10 / Mr_albumin

    IPR_K  <- 0.54; IPR_L  <- 0.54; IPR_ST <- 0.54; IPR_IN <- 0.54
    IPR_M  <- 0.27; IPR_A  <- 0.15; IPR_Lu <- 0.54; IPR_Sp <- 0.54
    IPR_H  <- 0.54; IPR_SK <- 0.62; IPR_Br <- 0.54; IPR_Go <- 0.54
    IPR_Bo <- 0.54
    IPR_R  <- (IPR_K+IPR_L+IPR_ST+IPR_IN+IPR_M+IPR_A+IPR_Lu+IPR_Sp+
               IPR_H+IPR_SK+IPR_Br+IPR_Go+IPR_Bo)/13

    if (sex == "M") { Hct <- 0.47 } else { Hct <- 0.42 }

    CalbKF_init  <- (CalbB_init/(1-Hct)) * IPR_K
    CalbLF_init  <- (CalbB_init/(1-Hct)) * IPR_L
    CalbSTF_init <- (CalbB_init/(1-Hct)) * IPR_ST
    CalbINF_init <- (CalbB_init/(1-Hct)) * IPR_IN
    CalbMF_init  <- (CalbB_init/(1-Hct)) * IPR_M
    CalbAF_init  <- (CalbB_init/(1-Hct)) * IPR_A
    CalbRF_init  <- (CalbB_init/(1-Hct)) * IPR_R
    CalbBoF_init <- (CalbB_init/(1-Hct)) * IPR_Bo
    CalbLuF_init <- (CalbB_init/(1-Hct)) * IPR_Lu
    CalbSPF_init <- (CalbB_init/(1-Hct)) * IPR_Sp
    CalbGoF_init <- (CalbB_init/(1-Hct)) * IPR_Go
    CalbHF_init  <- (CalbB_init/(1-Hct)) * IPR_H
    CalbBrF_init <- (CalbB_init/(1-Hct)) * IPR_Br
    CalbSKF_init <- (CalbB_init/(1-Hct)) * IPR_SK
    CalbLuAF_init <- (10/100) * (CalbB_init/(1-Hct))

    cytosolic_protein <- 44.4
    L_FABP_L          <- cytosolic_protein * 4.5/100 * 1000
    CFabpLT_init      <- f_fabp_avail * (L_FABP_L * 1e-3 / 14e3)
    CFabpKT_init      <- f_fabp_avail * 2.65e-6

    kon_alb  <- Ka     * koff_alb
    kon_fabp <- KLfabp * koff_fabp

    # Gut passive diffusion
    ClINFT_unscaled <- 18.1; Awell <- 9; Swell <- 1.12; well_protein <- 0.346
    protein_per_well <- (well_protein * Awell) / Swell
    RAF_papp         <- 1
    Peff_monolayer   <- 1.46e-6 * 3600

    k_gut_in  <- ((2*Peff_monolayer/100) * AINL) * 1000
    k_gut_out <- k_gut_in
    kabST     <- 0

    # Tissue-IF diffusion rates
    kLFLT   <- ((2*Peff_monolayer/100) * AcL)  * 1000
    kMFMT   <- ((2*Peff_monolayer/100) * AcM)  * 1000
    kSTFSTT <- ((2*Peff_monolayer/100) * AcST) * 1000
    kINFINT <- ((2*Peff_monolayer/100) * AcIN) * 1000
    kAFAT   <- ((2*Peff_monolayer/100) * AcA)  * 1000
    kLuTLuF <- ((2*Peff_monolayer/100) * AcLu) * 1000
    kLuTLuAF <- ((2*Peff_monolayer/100) * AcALF) * 1000
    kSPFSPT <- ((2*Peff_monolayer/100) * AcSP) * 1000
    kHFHT   <- ((2*Peff_monolayer/100) * AcH)  * 1000
    kBrFBrT <- ((2*Peff_monolayer/100) * AcBr) * 1000
    kGoFGoT <- ((2*Peff_monolayer/100) * AcGo) * 1000
    kSKFSKT <- ((2*Peff_monolayer/100) * AcSK) * 1000
    kBoFBoT <- ((2*Peff_monolayer/100) * AcBo) * 1000
    kRFRT   <- ((2*Peff_monolayer/100) * AcR)  * 1000

    # Tubule cell - filtrate diffusion
    kPtcTu  <- ((2*Peff_monolayer/100) * APT_tot) * 1000
    kDalcTu <- ((2*Peff_monolayer/100) * ALH_tot) * 1000
    kDtcTu  <- ((2*Peff_monolayer/100) * ADT_tot) * 1000
    kCdcTu  <- ((2*Peff_monolayer/100) * ACD_tot) * 1000

    # Tubule cell - interstitial diffusion
    kPtcF   <- ((2*Peff_monolayer/100) * AcK_PTC) * 1000
    kDalcF  <- ((2*Peff_monolayer)     * AcK_LHC) * 1000
    kDtcF   <- ((2*Peff_monolayer/100) * AcK_DTC) * 1000
    kCdcF   <- ((2*Peff_monolayer/100) * AcK_CDC) * 1000

    # Capillary permeability (Renkin + fenestration model)
    dif      <- 5.46e-6
    kboltzman <- 1.38e-23
    Temp     <- 37 + 273
    dyn_visc <- 6.9e-4
    R_H      <- kboltzman * Temp / (6*pi*dif*1e-4*dyn_visc) * 1e9

    pore_diameters <- c(Ki=9, Li=105, St=9, In=9, Mu=5, Ad=5, Re=5,
                        Lu=5, Sp=5000, Ht=5, Br=1, Go=9, Sk=9, Bo=5)
    lambda           <- R_H / (pore_diameters/2)
    lambda["Br"]     <- 1
    renkin_reduction <- (1-lambda)^2 * (1 - 2.104*lambda + 2.09*lambda^3 - 0.95*lambda^5)
    wall_width       <- 0.5e-6 * 100
    basement_membrane <- 0.5e-6 * 100
    Pgap             <- reduction_factor * dif * renkin_reduction / (wall_width + basement_membrane)
    Pgap_sinusoidal  <- reduction_factor * dif * renkin_reduction / wall_width

    f_kidney        <- 0.35; f_liver  <- 0.08; f_spleen    <- 0.08
    f_intestine     <- 0.095; f_non_fenestrated <- 0.0048

    Ptrans_diff_K  <- Peff_monolayer*10 * (1-f_kidney)
    Ptrans_diff_L  <- Peff_monolayer*10 * (1-f_liver)
    Ptrans_diff_ST <- Peff_monolayer*10 * (1-f_intestine)
    Ptrans_diff_IN <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_A  <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_M  <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_R  <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_Lu <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_SP <- Peff_monolayer*10 * (1-f_spleen)
    Ptrans_diff_H  <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_Br <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_Go <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_SK <- Peff_monolayer*10 * (1-f_non_fenestrated)
    Ptrans_diff_Bo <- Peff_monolayer*10 * (1-f_non_fenestrated)

    PparaKi <- Pgap["Ki"]*3600*10*f_kidney
    PparaLi <- Pgap_sinusoidal["Li"]*3600*10*f_liver
    PparaSt <- Pgap["St"]*3600*10*f_non_fenestrated
    PparaIn <- Pgap["In"]*3600*10*f_intestine
    PparaMu <- Pgap["Mu"]*3600*10*f_non_fenestrated
    PparaAd <- Pgap["Ad"]*3600*10*f_non_fenestrated
    PparaRe <- Pgap["Re"]*3600*10*f_non_fenestrated
    PparaLu <- Pgap["Lu"]*3600*10*f_non_fenestrated
    PparaSp <- Pgap_sinusoidal["Sp"]*3600*10*f_spleen
    PparaHt <- Pgap["Ht"]*3600*10*f_non_fenestrated
    PparaBr <- 0
    PparaGo <- Pgap["Go"]*3600*10*f_non_fenestrated
    PparaSk <- Pgap["Sk"]*3600*10*f_non_fenestrated
    PparaBo <- Pgap["Bo"]*3600*10*f_non_fenestrated

    return(list(
      'BW'=BW, 'sex'=sex, 'MW'=MW, 'Hct'=Hct,
      'admin.dose'=admin.dose, 'admin.time'=admin.time,
      'admin.type'=admin.type, 'depfr_AF'=depfr_AF, 'depfr_head'=depfr_head,

      # Volumes
      'VKiB_cap'=VKiB_cap, 'VKiF'=VKiF, 'VKiT'=VKiT, 'VKi_tot'=VKi_tot,
      'VFil'=VFil, 'VBladder'=VBladder,
      'VPT_tot'=VPT_tot, 'VLH_tot'=VLH_tot, 'VDT_tot'=VDT_tot, 'VCD_tot'=VCD_tot,
      'VPTC'=VPTC, 'VLHC'=VLHC, 'VDTC'=VDTC, 'VCDC'=VCDC,
      'VLiB_cap'=VLiB_cap, 'VLiF'=VLiF, 'VLiT'=VLiT, 'VLi_tot'=VLi_tot, 'VBile'=VBile,
      'VInB_cap'=VInB_cap, 'VInF'=VInF, 'VInT'=VInT, 'VIn_tot'=VIn_tot,
      'VStB_cap'=VStB_cap, 'VStF'=VStF, 'VStT'=VStT, 'VSt_tot'=VSt_tot,
      'VStL'=VStL, 'VInL'=VInL,
      'VMuB_cap'=VMuB_cap, 'VMuF'=VMuF, 'VMuT'=VMuT, 'VMu_tot'=VMu_tot,
      'VAdB_cap'=VAdB_cap, 'VAdF'=VAdF, 'VAdT'=VAdT, 'VAd_tot'=VAd_tot,
      'VLuB_cap'=VLuB_cap, 'VUA'=VUA, 'VLuF'=VLuF, 'VLuAF'=VLuAF, 'VLuT'=VLuT, 'VLu_tot'=VLu_tot,
      'VSpB_cap'=VSpB_cap, 'VSpF'=VSpF, 'VSpT'=VSpT, 'VSp_tot'=VSp_tot,
      'VHtB_cap'=VHtB_cap, 'VHtF'=VHtF, 'VHtT'=VHtT, 'VHt_tot'=VHt_tot,
      'VBrB_cap'=VBrB_cap, 'VBrF'=VBrF, 'VBrT'=VBrT, 'VBr_tot'=VBr_tot,
      'VGoB_cap'=VGoB_cap, 'VGoF'=VGoF, 'VGoT'=VGoT, 'VGo_tot'=VGo_tot,
      'VSkB_cap'=VSkB_cap, 'VSkF'=VSkF, 'VSkT'=VSkT, 'VSk_tot'=VSk_tot,
      'VBoB_cap'=VBoB_cap, 'VBoF'=VBoF, 'VBoT'=VBoT, 'VBo_tot'=VBo_tot,
      'VReB_cap'=VReB_cap, 'VReF'=VReF, 'VReT'=VReT, 'VRe_tot'=VRe_tot,
      'VBlood'=VBlood, 'VVen'=VVen, 'VArt'=VArt,

      # Surface areas
      'A_peritubular'=A_peritubular,
      'AL'=AL, 'AST'=AST, 'AIN'=AIN, 'AM'=AM, 'AA'=AA, 'AR'=AR,
      'ALu'=ALu, 'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AGo'=AGo, 'ASK'=ASK, 'ABo'=ABo,
      'AINL'=AINL,
      'AcL'=AcL, 'AcST'=AcST, 'AcIN'=AcIN, 'AcM'=AcM, 'AcA'=AcA, 'AcLu'=AcLu,
      'AcALF'=AcALF, 'AcSP'=AcSP, 'AcH'=AcH, 'AcBr'=AcBr, 'AcGo'=AcGo,
      'AcSK'=AcSK, 'AcBo'=AcBo, 'AcR'=AcR,
      'AcK_PTC'=AcK_PTC, 'AcK_LHC'=AcK_LHC, 'AcK_DTC'=AcK_DTC, 'AcK_CDC'=AcK_CDC,
      'APT_tot'=APT_tot, 'ALH_tot'=ALH_tot, 'ADT_tot'=ADT_tot, 'ACD_tot'=ACD_tot,

      # Flows
      'Qcardiac'=Qcardiac, 'QBK'=QBK, 'QBL'=QBL, 'QBLtot'=QBLtot,
      'QBM'=QBM, 'QBA'=QBA, 'QBR'=QBR, 'QBLu'=QBLu,
      'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST, 'QBIN'=QBIN,
      'QBGo'=QBGo, 'QBSK'=QBSK, 'QBBo'=QBBo,
      'Qbile'=Qbile, 'Vbile'=Vbile, 'Qfeces'=Qfeces, 'feces_density'=feces_density,
      'QGFR'=QGFR, 'Qurine'=Qurine, 'Qelim'=Qelim, 'QGE'=QGE,
      'QPT_tot'=QPT_tot, 'QLH_tot'=QLH_tot, 'QDT_tot'=QDT_tot, 'QCD_tot'=QCD_tot,

      # Transporters
      'VmK_Oatp'=VmK_Oatp, 'KmK_Oatp'=KmK_Oatp,
      'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1,
      'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3,
      'VmK_Urat'=VmK_Urat, 'KmK_Urat'=KmK_Urat,
      'VmK_baso'=VmK_baso, 'KmK_baso'=KmK_baso,
      'VmK_api'=VmK_api,   'KmK_api'=KmK_api,
      'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'=KmL_Oatp,
      'VmL_Oatp2'=VmL_Oatp2, 'KmL_Oatp2'=KmL_Oatp2,
      'VmL_Ntcp'=VmL_Ntcp,   'KmL_Ntcp'=KmL_Ntcp,
      'VmLu_Oatp_ap'=VmLu_Oatp_ap, 'KmLu_Oatp_ap'=KmLu_Oatp_ap,
      'VmLu_Oatp_bas'=VmLu_Oatp_bas, 'KmLu_Oatp_bas'=KmLu_Oatp_bas,
      'VmIn_Oatp2'=VmIn_Oatp2, 'KmIn_Oatp2'=KmIn_Oatp2,
      'CL_hepatobiliary'=CL_hepatobiliary, 'CLfeces'=CLfeces,

      # Albumin / FABP
      'CalbB_init'=CalbB_init,
      'CalbKF_init'=CalbKF_init, 'CalbLF_init'=CalbLF_init,
      'CalbSTF_init'=CalbSTF_init, 'CalbINF_init'=CalbINF_init,
      'CalbMF_init'=CalbMF_init, 'CalbAF_init'=CalbAF_init,
      'CalbRF_init'=CalbRF_init, 'CalbBoF_init'=CalbBoF_init,
      'CalbLuF_init'=CalbLuF_init, 'CalbSPF_init'=CalbSPF_init,
      'CalbGoF_init'=CalbGoF_init, 'CalbHF_init'=CalbHF_init,
      'CalbBrF_init'=CalbBrF_init, 'CalbSKF_init'=CalbSKF_init,
      'CalbLuAF_init'=CalbLuAF_init,
      'CFabpKT_init'=CFabpKT_init, 'CFabpLT_init'=CFabpLT_init,
      'Ka'=Ka, 'KLfabp'=KLfabp,
      'koff_alb'=koff_alb, 'koff_fabp'=koff_fabp,
      'kon_alb'=kon_alb, 'kon_fabp'=kon_fabp,

      # Diffusion rates
      'kUAB'=kUAB, 'CLEua'=CLEua,
      'k_gut_in'=k_gut_in, 'k_gut_out'=k_gut_out, 'kabST'=kabST,
      'kLFLT'=kLFLT, 'kMFMT'=kMFMT, 'kSTFSTT'=kSTFSTT, 'kINFINT'=kINFINT,
      'kAFAT'=kAFAT, 'kLuTLuF'=kLuTLuF, 'kLuTLuAF'=kLuTLuAF, 'kSPFSPT'=kSPFSPT,
      'kHFHT'=kHFHT, 'kBrFBrT'=kBrFBrT, 'kGoFGoT'=kGoFGoT, 'kSKFSKT'=kSKFSKT,
      'kBoFBoT'=kBoFBoT, 'kRFRT'=kRFRT,
      'kPtcTu'=kPtcTu, 'kDalcTu'=kDalcTu, 'kDtcTu'=kDtcTu, 'kCdcTu'=kCdcTu,
      'kPtcF'=kPtcF, 'kDalcF'=kDalcF, 'kDtcF'=kDtcF, 'kCdcF'=kCdcF,

      # Capillary permeability
      'Ptrans_diff_K'=Ptrans_diff_K, 'Ptrans_diff_L'=Ptrans_diff_L,
      'Ptrans_diff_ST'=Ptrans_diff_ST, 'Ptrans_diff_IN'=Ptrans_diff_IN,
      'Ptrans_diff_A'=Ptrans_diff_A, 'Ptrans_diff_M'=Ptrans_diff_M,
      'Ptrans_diff_R'=Ptrans_diff_R, 'Ptrans_diff_Lu'=Ptrans_diff_Lu,
      'Ptrans_diff_SP'=Ptrans_diff_SP, 'Ptrans_diff_H'=Ptrans_diff_H,
      'Ptrans_diff_Br'=Ptrans_diff_Br, 'Ptrans_diff_Go'=Ptrans_diff_Go,
      'Ptrans_diff_SK'=Ptrans_diff_SK, 'Ptrans_diff_Bo'=Ptrans_diff_Bo,
      'PparaKi'=PparaKi, 'PparaLi'=PparaLi, 'PparaSt'=PparaSt, 'PparaIn'=PparaIn,
      'PparaMu'=PparaMu, 'PparaAd'=PparaAd, 'PparaRe'=PparaRe, 'PparaLu'=PparaLu,
      'PparaSp'=PparaSp, 'PparaHt'=PparaHt, 'PparaBr'=PparaBr, 'PparaGo'=PparaGo,
      'PparaSk'=PparaSk, 'PparaBo'=PparaBo,

      # Desorption (dust/inhalation)
      'k_desorption_fast'=k_desorption_fast, 'k_desorption_slow'=k_desorption_slow,
      'f_fast'=f_fast, 'fdust_retained'=fdust_retained
    ))
  })
}

# =============================================================================
# 2. Initial conditions
# =============================================================================

create.inits <- function(parameters) {
  with(as.list(parameters), {

    CalbVenf <- CalbB_init;   MVenf  <- 0; MVenb  <- 0
    CalbArtf <- CalbB_init;   MArtf  <- 0; MArtb  <- 0
    CalbKBf  <- CalbB_init;   MKBf   <- 0; MKBb   <- 0
    CalbLBf  <- CalbB_init;   MLBf   <- 0; MLBb   <- 0
    CalbSTBf <- CalbB_init;   MSTBf  <- 0; MSTBb  <- 0
    CalbINBf <- CalbB_init;   MINBf  <- 0; MINBb  <- 0
    CalbMBf  <- CalbB_init;   MMBf   <- 0; MMBb   <- 0
    CalbABf  <- CalbB_init;   MABf   <- 0; MABb   <- 0
    CalbRBf  <- CalbB_init;   MRBf   <- 0; MRBb   <- 0
    CalbLuBf <- CalbB_init;   MLuBf  <- 0; MLuBb  <- 0
    CalbSPBf <- CalbB_init;   MSPBf  <- 0; MSPBb  <- 0
    CalbHBf  <- CalbB_init;   MHBf   <- 0; MHBb   <- 0
    CalbBrBf <- CalbB_init;   MBrBf  <- 0; MBrBb  <- 0
    CalbGoBf <- CalbB_init;   MGoBf  <- 0; MGoBb  <- 0
    CalbSKBf <- CalbB_init;   MSKBf  <- 0; MSKBb  <- 0
    CalbBoBf <- CalbB_init;   MBoBf  <- 0; MBoBb  <- 0

    CalbKFf  <- CalbKF_init;  MKFf   <- 0; MKFb   <- 0
    CalbLFf  <- CalbLF_init;  MLFf   <- 0; MLFb   <- 0
    CalbSTFf <- CalbSTF_init; MSTFf  <- 0; MSTFb  <- 0
    CalbINFf <- CalbINF_init; MINFf  <- 0; MINFb  <- 0
    CalbMFf  <- CalbMF_init;  MMFf   <- 0; MMFb   <- 0
    CalbAFf  <- CalbAF_init;  MAFf   <- 0; MAFb   <- 0
    CalbRFf  <- CalbRF_init;  MRFf   <- 0; MRFb   <- 0
    CalbLuFf <- CalbLuF_init; MLuFf  <- 0; MLuFb  <- 0
    CalbSPFf <- CalbSPF_init; MSPFf  <- 0; MSPFb  <- 0
    CalbHFf  <- CalbHF_init;  MHFf   <- 0; MHFb   <- 0
    CalbBrFf <- CalbBrF_init; MBrFf  <- 0; MBrFb  <- 0
    CalbGoFf <- CalbGoF_init; MGoFf  <- 0; MGoFb  <- 0
    CalbSKFf <- CalbSKF_init; MSKFf  <- 0; MSKFb  <- 0
    CalbBoFf <- CalbBoF_init; MBoFf  <- 0; MBoFb  <- 0

    CFabpPTCf <- CFabpKT_init; CFabpLHCf <- CFabpKT_init
    CFabpDTCf <- CFabpKT_init; CFabpCDCf <- CFabpKT_init

    MPTCf <- 0; MLHCf <- 0; MDTCf <- 0; MCDCf <- 0
    MPT   <- 0; MLH   <- 0; MDT   <- 0; MCD   <- 0

    CFabpLTf  <- CFabpLT_init; CalbLuAFf <- CalbLuAF_init

    MLTf  <- 0; MLTb  <- 0; MPTCb <- 0; MLHCb <- 0; MDTCb <- 0
    MCDCb <- 0; MBile <- 0; MSTTf <- 0; MINTf <- 0; MLuTf <- 0
    MUA   <- 0; MLuAFf <- 0; MLuAFb <- 0
    MLuAFdust_fast <- 0; MLuAFdust_slow <- 0
    MSPTf <- 0; MHTf  <- 0; MBrTf <- 0
    MGoTf <- 0; MSKTf <- 0; MBoTf <- 0; MMTf <- 0; MATf <- 0; MRTf <- 0

    MBladder <- 0; Murine <- 0; MSTL <- 0; MINL <- 0
    Mfeces   <- 0; Vurine <- 0; Vfeces <- 0

    return(c(
      'CalbVenf'=CalbVenf, 'CalbArtf'=CalbArtf,
      'CalbKBf'=CalbKBf, 'CalbLBf'=CalbLBf, 'CalbSTBf'=CalbSTBf,
      'CalbINBf'=CalbINBf, 'CalbMBf'=CalbMBf, 'CalbABf'=CalbABf,
      'CalbRBf'=CalbRBf, 'CalbLuBf'=CalbLuBf, 'CalbSPBf'=CalbSPBf,
      'CalbHBf'=CalbHBf, 'CalbBrBf'=CalbBrBf, 'CalbGoBf'=CalbGoBf,
      'CalbSKBf'=CalbSKBf, 'CalbBoBf'=CalbBoBf,
      'CalbKFf'=CalbKFf, 'CalbLFf'=CalbLFf, 'CalbSTFf'=CalbSTFf,
      'CalbINFf'=CalbINFf, 'CalbMFf'=CalbMFf, 'CalbAFf'=CalbAFf,
      'CalbRFf'=CalbRFf, 'CalbLuFf'=CalbLuFf, 'CalbSPFf'=CalbSPFf,
      'CalbHFf'=CalbHFf, 'CalbBrFf'=CalbBrFf, 'CalbGoFf'=CalbGoFf,
      'CalbSKFf'=CalbSKFf, 'CalbBoFf'=CalbBoFf,
      'CFabpPTCf'=CFabpPTCf, 'CFabpLHCf'=CFabpLHCf,
      'CFabpDTCf'=CFabpDTCf, 'CFabpCDCf'=CFabpCDCf,
      'CFabpLTf'=CFabpLTf, 'CalbLuAFf'=CalbLuAFf,
      'MVenb'=MVenb, 'MArtb'=MArtb, 'MKBb'=MKBb,
      'MLBb'=MLBb, 'MSTBb'=MSTBb, 'MINBb'=MINBb, 'MMBb'=MMBb,
      'MABb'=MABb, 'MRBb'=MRBb, 'MLuBb'=MLuBb, 'MSPBb'=MSPBb,
      'MHBb'=MHBb, 'MBrBb'=MBrBb, 'MGoBb'=MGoBb,
      'MSKBb'=MSKBb, 'MBoBb'=MBoBb,
      'MKFb'=MKFb, 'MLFb'=MLFb, 'MSTFb'=MSTFb, 'MINFb'=MINFb,
      'MMFb'=MMFb, 'MAFb'=MAFb, 'MRFb'=MRFb, 'MLuFb'=MLuFb,
      'MSPFb'=MSPFb, 'MHFb'=MHFb, 'MBrFb'=MBrFb,
      'MGoFb'=MGoFb, 'MSKFb'=MSKFb, 'MBoFb'=MBoFb,
      'MPTCb'=MPTCb, 'MLHCb'=MLHCb, 'MDTCb'=MDTCb,
      'MCDCb'=MCDCb, 'MLTb'=MLTb, 'MLuAFb'=MLuAFb,
      'MArtf'=MArtf, 'MVenf'=MVenf, 'MKBf'=MKBf,
      'MKFf'=MKFf, 'MPTCf'=MPTCf, 'MLHCf'=MLHCf,
      'MDTCf'=MDTCf, 'MCDCf'=MCDCf,
      'MPT'=MPT, 'MLH'=MLH, 'MDT'=MDT, 'MCD'=MCD,
      'MBladder'=MBladder, 'MLBf'=MLBf,
      'MLFf'=MLFf, 'MLTf'=MLTf, 'MBile'=MBile,
      'MSTBf'=MSTBf, 'MSTFf'=MSTFf, 'MSTTf'=MSTTf, 'MSTL'=MSTL,
      'MINBf'=MINBf, 'MINFf'=MINFf, 'MINTf'=MINTf, 'MINL'=MINL,
      'MMBf'=MMBf, 'MMFf'=MMFf, 'MMTf'=MMTf,
      'MABf'=MABf, 'MAFf'=MAFf, 'MATf'=MATf,
      'MRBf'=MRBf, 'MRFf'=MRFf, 'MRTf'=MRTf,
      'MUA'=MUA, 'MLuBf'=MLuBf, 'MLuFf'=MLuFf, 'MLuTf'=MLuTf,
      'MLuAFf'=MLuAFf, 'MLuAFdust_fast'=MLuAFdust_fast, 'MLuAFdust_slow'=MLuAFdust_slow,
      'MSPBf'=MSPBf, 'MSPFf'=MSPFf, 'MSPTf'=MSPTf,
      'MHBf'=MHBf, 'MHFf'=MHFf, 'MHTf'=MHTf,
      'MBrBf'=MBrBf, 'MBrFf'=MBrFf, 'MBrTf'=MBrTf,
      'MGoBf'=MGoBf, 'MGoFf'=MGoFf, 'MGoTf'=MGoTf,
      'MSKBf'=MSKBf, 'MSKFf'=MSKFf, 'MSKTf'=MSKTf,
      'MBoBf'=MBoBf, 'MBoFf'=MBoFf, 'MBoTf'=MBoTf,
      'Mfeces'=Mfeces, 'Murine'=Murine, 'Vfeces'=Vfeces, 'Vurine'=Vurine
    ))
  })
}

# =============================================================================
# 3. Dosing events
# =============================================================================

create.events <- function(parameters) {
  with(as.list(parameters), {
    ldose  <- length(admin.dose)
    ltimes <- length(admin.time)
    if (ltimes != ldose) stop("admin.time and admin.dose must have the same length")

    if (admin.type == "iv") {
      events <- list(data = data.frame(var='MVenf', time=admin.time,
                                       value=admin.dose, method='add'))
    } else if (admin.type == "oral") {
      events <- list(data = data.frame(var='MSTL', time=admin.time,
                                       value=admin.dose, method='add'))
    } else if (admin.type == "inh") {
      events <- list(data = rbind(
        data.frame(var='MLuAFdust_fast', time=admin.time,
                   value=admin.dose * depfr_AF * f_fast, method='add'),
        data.frame(var='MLuAFdust_slow', time=admin.time,
                   value=admin.dose * depfr_AF * (1-f_fast), method='add')
      ))
    } else if (admin.type == "nasal") {
      events <- list(data = rbind(
        data.frame(var='MUA',    time=admin.time,
                   value=admin.dose * depfr_head, method='add'),
        data.frame(var='MLuAFf', time=admin.time,
                   value=admin.dose * depfr_AF,   method='add')
      ))
    } else if (admin.type == "dermal") {
      events <- list(data = data.frame(var='MSKTf', time=admin.time,
                                       value=admin.dose, method='add'))
    }
    return(events)
  })
}

# =============================================================================
# 4. Custom function (placeholder)
# =============================================================================

custom.func <- function() {
  return()
}

# =============================================================================
# 5. ODE system
# =============================================================================

ode.func <- function(time, inits, params, custom.func) {
  with(as.list(c(inits, params)), {

    MVen  <- MVenf + MVenb
    MArt  <- MArtf + MArtb
    CVen  <- MVen  / VVen;   CVenb <- MVenb / VVen;   CVenf <- MVenf / VVen
    CArt  <- MArt  / VArt;   CArtf <- MArtf / VArt;   CArtb <- MArtb / VArt

    MKB   <- MKBf  + MKBb;   MKF  <- MKFf  + MKFb
    CKB   <- MKB   / VKiB_cap; CKBf <- MKBf / VKiB_cap; CKBb <- MKBb / VKiB_cap
    CKF   <- MKF   / VKiF;     CKFf <- MKFf / VKiF;     CKFb <- MKFb / VKiF
    CPTCf <- MPTCf / VPTC;  CPTCb <- MPTCb / VPTC
    CLHCf <- MLHCf / VLHC;  CLHCb <- MLHCb / VLHC
    CDTCf <- MDTCf / VDTC;  CDTCb <- MDTCb / VDTC
    CCDCf <- MCDCf / VCDC;  CCDCb <- MCDCb / VCDC
    MKTf  <- MPTCf + MLHCf + MDTCf + MCDCf
    MKTb  <- MPTCb + MLHCb + MDTCb + MCDCb
    MKT   <- MKTf  + MKTb
    Mfil  <- MPT   + MLH   + MDT   + MCD
    CKTb  <- MKTb  / VKiT; CKTf <- MKTf / VKiT; CKT <- MKT / VKiT
    CPT   <- MPT   / VPT_tot; CLH <- MLH / VLH_tot
    CDT   <- MDT   / VDT_tot; CCD <- MCD / VCD_tot
    CBladder <- MBladder / VBladder

    MLB   <- MLBf  + MLBb;   MLF <- MLFf + MLFb;   MLT <- MLTf + MLTb
    CLB   <- MLB   / VLiB_cap; CLBf <- MLBf / VLiB_cap; CLBb <- MLBb / VLiB_cap
    CLF   <- MLF   / VLiF;     CLFf <- MLFf / VLiF;     CLFb <- MLFb / VLiF
    CLT   <- MLT   / VLiT;     CLTf <- MLTf / VLiT;     CLTb <- MLTb / VLiT
    CBile <- MBile / VBile

    MSTB  <- MSTBf + MSTBb;  MSTF <- MSTFf + MSTFb;  MSTT <- MSTTf
    CSTB  <- MSTB  / VStB_cap; CSTBf <- MSTBf / VStB_cap; CSTBb <- MSTBb / VStB_cap
    CSTF  <- MSTF  / VStF;     CSTFf <- MSTFf / VStF;     CSTFb <- MSTFb / VStF
    CSTTf <- MSTTf / VStT
    CSTL  <- MSTL  / VStL;     CINL  <- MINL  / VInL

    MINB  <- MINBf + MINBb;  MINF <- MINFf + MINFb;  MINT <- MINTf
    CINB  <- MINB  / VInB_cap; CINBf <- MINBf / VInB_cap; CINBb <- MINBb / VInB_cap
    CINF  <- MINF  / VInF;     CINFf <- MINFf / VInF;     CINFb <- MINFb / VInF
    CINTf <- MINTf / VInT

    MMB   <- MMBf  + MMBb;   MMF  <- MMFf  + MMFb;   MMT  <- MMTf
    CMB   <- MMB   / VMuB_cap; CMBf <- MMBf / VMuB_cap; CMBb <- MMBb / VMuB_cap
    CMF   <- MMF   / VMuF;     CMFf <- MMFf / VMuF;     CMFb <- MMFb / VMuF
    CMTf  <- MMTf  / VMuT

    MAB   <- MABf  + MABb;   MAF  <- MAFf  + MAFb;   MAT  <- MATf
    CAB   <- MAB   / VAdB_cap; CABf <- MABf / VAdB_cap; CABb <- MABb / VAdB_cap
    CAF   <- MAF   / VAdF;     CAFf <- MAFf / VAdF;     CAFb <- MAFb / VAdF
    CATf  <- MATf  / VAdT

    MRB   <- MRBf  + MRBb;   MRF  <- MRFf  + MRFb;   MRT  <- MRTf
    CRB   <- MRB   / VReB_cap; CRBf <- MRBf / VReB_cap; CRBb <- MRBb / VReB_cap
    CRF   <- MRF   / VReF;     CRFf <- MRFf / VReF;     CRFb <- MRFb / VReF
    CRTf  <- MRTf  / VReT

    CUA   <- MUA   / VUA
    MLuB  <- MLuBf + MLuBb;  MLuF <- MLuFf + MLuFb;  MLuT <- MLuTf
    MLuAF <- MLuAFf + MLuAFb
    CLuB  <- MLuB  / VLuB_cap; CLuBf <- MLuBf / VLuB_cap; CLuBb <- MLuBb / VLuB_cap
    CLuF  <- MLuF  / VLuF;     CLuFf <- MLuFf / VLuF;     CLuFb <- MLuFb / VLuF
    CLuTf <- MLuTf / VLuT
    CLuAF <- MLuAF / VLuAF; CLuAFf <- MLuAFf / VLuAF; CLuAFb <- MLuAFb / VLuAF

    MSPB  <- MSPBf + MSPBb;  MSPF <- MSPFf + MSPFb;  MSPT <- MSPTf
    CSPB  <- MSPB  / VSpB_cap; CSPBf <- MSPBf / VSpB_cap; CSPBb <- MSPBb / VSpB_cap
    CSPF  <- MSPF  / VSpF;     CSPFf <- MSPFf / VSpF;     CSPFb <- MSPFb / VSpF
    CSPTf <- MSPTf / VSpT

    MHB   <- MHBf  + MHBb;   MHF  <- MHFf  + MHFb;   MHT  <- MHTf
    CHB   <- MHB   / VHtB_cap; CHBf <- MHBf / VHtB_cap; CHBb <- MHBb / VHtB_cap
    CHF   <- MHF   / VHtF;     CHFf <- MHFf / VHtF;     CHFb <- MHFb / VHtF
    CHTf  <- MHTf  / VHtT

    MBrB  <- MBrBf + MBrBb;  MBrF <- MBrFf + MBrFb;  MBrT <- MBrTf
    CBrB  <- MBrB  / VBrB_cap; CBrBf <- MBrBf / VBrB_cap; CBrBb <- MBrBb / VBrB_cap
    CBrF  <- MBrF  / VBrF;     CBrFf <- MBrFf / VBrF;     CBrFb <- MBrFb / VBrF
    CBrTf <- MBrTf / VBrT

    MGoB  <- MGoBf + MGoBb;  MGoF <- MGoFf + MGoFb;  MGoT <- MGoTf
    CGoB  <- MGoB  / VGoB_cap; CGoBf <- MGoBf / VGoB_cap; CGoBb <- MGoBb / VGoB_cap
    CGoF  <- MGoF  / VGoF;     CGoFf <- MGoFf / VGoF;     CGoFb <- MGoFb / VGoF
    CGoTf <- MGoTf / VGoT

    MSKB  <- MSKBf + MSKBb;  MSKF <- MSKFf + MSKFb;  MSKT <- MSKTf
    CSKB  <- MSKB  / VSkB_cap; CSKBf <- MSKBf / VSkB_cap; CSKBb <- MSKBb / VSkB_cap
    CSKF  <- MSKF  / VSkF;     CSKFf <- MSKFf / VSkF;     CSKFb <- MSKFb / VSkF
    CSKTf <- MSKTf / VSkT

    MBoB  <- MBoBf + MBoBb;  MBoF <- MBoFf + MBoFb;  MBoT <- MBoTf
    CBoB  <- MBoB  / VBoB_cap; CBoBf <- MBoBf / VBoB_cap; CBoBb <- MBoBb / VBoB_cap
    CBoF  <- MBoF  / VBoF;     CBoFf <- MBoFf / VBoF;     CBoFb <- MBoFb / VBoF
    CBoTf <- MBoTf / VBoT

    # Albumin binding ODEs (blood)
    dCalbVenf  <- koff_alb*CVenb/MW/1e6 - kon_alb*CalbVenf*CVenf/MW/1e6
    dCalbArtf  <- koff_alb*CArtb/MW/1e6 - kon_alb*CalbArtf*CArtf/MW/1e6
    dCalbKBf   <- koff_alb*CKBb/MW/1e6  - kon_alb*CalbKBf*CKBf/MW/1e6
    dCalbLBf   <- koff_alb*CLBb/MW/1e6  - kon_alb*CalbLBf*CLBf/MW/1e6
    dCalbSTBf  <- koff_alb*CSTBb/MW/1e6 - kon_alb*CalbSTBf*CSTBf/MW/1e6
    dCalbINBf  <- koff_alb*CINBb/MW/1e6 - kon_alb*CalbINBf*CINBf/MW/1e6
    dCalbMBf   <- koff_alb*CMBb/MW/1e6  - kon_alb*CalbMBf*CMBf/MW/1e6
    dCalbABf   <- koff_alb*CABb/MW/1e6  - kon_alb*CalbABf*CABf/MW/1e6
    dCalbRBf   <- koff_alb*CRBb/MW/1e6  - kon_alb*CalbRBf*CRBf/MW/1e6
    dCalbLuBf  <- koff_alb*CLuBb/MW/1e6 - kon_alb*CalbLuBf*CLuBf/MW/1e6
    dCalbSPBf  <- koff_alb*CSPBb/MW/1e6 - kon_alb*CalbSPBf*CSPBf/MW/1e6
    dCalbHBf   <- koff_alb*CHBb/MW/1e6  - kon_alb*CalbHBf*CHBf/MW/1e6
    dCalbBrBf  <- koff_alb*CBrBb/MW/1e6 - kon_alb*CalbBrBf*CBrBf/MW/1e6
    dCalbGoBf  <- koff_alb*CGoBb/MW/1e6 - kon_alb*CalbGoBf*CGoBf/MW/1e6
    dCalbSKBf  <- koff_alb*CSKBb/MW/1e6 - kon_alb*CalbSKBf*CSKBf/MW/1e6
    dCalbBoBf  <- koff_alb*CBoBb/MW/1e6 - kon_alb*CalbBoBf*CBoBf/MW/1e6

    # Albumin binding ODEs (interstitial fluid)
    dCalbKFf   <- koff_alb*CKFb/MW/1e6  - kon_alb*CalbKFf*CKFf/MW/1e6
    dCalbLFf   <- koff_alb*CLFb/MW/1e6  - kon_alb*CalbLFf*CLFf/MW/1e6
    dCalbSTFf  <- koff_alb*CSTFb/MW/1e6 - kon_alb*CalbSTFf*CSTFf/MW/1e6
    dCalbINFf  <- koff_alb*CINFb/MW/1e6 - kon_alb*CalbINFf*CINFf/MW/1e6
    dCalbMFf   <- koff_alb*CMFb/MW/1e6  - kon_alb*CalbMFf*CMFf/MW/1e6
    dCalbAFf   <- koff_alb*CAFb/MW/1e6  - kon_alb*CalbAFf*CAFf/MW/1e6
    dCalbRFf   <- koff_alb*CRFb/MW/1e6  - kon_alb*CalbRFf*CRFf/MW/1e6
    dCalbLuFf  <- koff_alb*CLuFb/MW/1e6 - kon_alb*CalbLuFf*CLuFf/MW/1e6
    dCalbSPFf  <- koff_alb*CSPFb/MW/1e6 - kon_alb*CalbSPFf*CSPFf/MW/1e6
    dCalbHFf   <- koff_alb*CHFb/MW/1e6  - kon_alb*CalbHFf*CHFf/MW/1e6
    dCalbBrFf  <- koff_alb*CBrFb/MW/1e6 - kon_alb*CalbBrFf*CBrFf/MW/1e6
    dCalbGoFf  <- koff_alb*CGoFb/MW/1e6 - kon_alb*CalbGoFf*CGoFf/MW/1e6
    dCalbSKFf  <- koff_alb*CSKFb/MW/1e6 - kon_alb*CalbSKFf*CSKFf/MW/1e6
    dCalbBoFf  <- koff_alb*CBoFb/MW/1e6 - kon_alb*CalbBoFf*CBoFf/MW/1e6

    # FABP binding ODEs (tissue)
    dCFabpPTCf <- koff_fabp*CPTCb/MW/1e6 - kon_fabp*CFabpPTCf*CPTCf/MW/1e6
    dCFabpLHCf <- koff_fabp*CLHCb/MW/1e6 - kon_fabp*CFabpLHCf*CLHCf/MW/1e6
    dCFabpDTCf <- koff_fabp*CDTCb/MW/1e6 - kon_fabp*CFabpDTCf*CDTCf/MW/1e6
    dCFabpCDCf <- koff_fabp*CCDCb/MW/1e6 - kon_fabp*CFabpCDCf*CCDCf/MW/1e6
    dCFabpLTf  <- koff_fabp*CLTb/MW/1e6  - kon_fabp*CFabpLTf*CLTf/MW/1e6
    dCalbLuAFf <- koff_alb*CLuAFb/MW/1e6 - kon_alb*CalbLuAFf*CLuAFf/MW/1e6

    # Bound mass ODEs
    dMVenb  <- kon_alb*CalbVenf*CVenf*VVen   - koff_alb*CVenb*VVen
    dMArtb  <- kon_alb*CalbArtf*CArtf*VArt   - koff_alb*CArtb*VArt
    dMKBb   <- kon_alb*CalbKBf*CKBf*VKiB_cap - koff_alb*CKBb*VKiB_cap
    dMLBb   <- kon_alb*CalbLBf*CLBf*VLiB_cap - koff_alb*CLBb*VLiB_cap
    dMSTBb  <- kon_alb*CalbSTBf*CSTBf*VStB_cap - koff_alb*CSTBb*VStB_cap
    dMINBb  <- kon_alb*CalbINBf*CINBf*VInB_cap - koff_alb*CINBb*VInB_cap
    dMMBb   <- kon_alb*CalbMBf*CMBf*VMuB_cap - koff_alb*CMBb*VMuB_cap
    dMABb   <- kon_alb*CalbABf*CABf*VAdB_cap - koff_alb*CABb*VAdB_cap
    dMRBb   <- kon_alb*CalbRBf*CRBf*VReB_cap - koff_alb*CRBb*VReB_cap
    dMLuBb  <- kon_alb*CalbLuBf*CLuBf*VLuB_cap - koff_alb*CLuBb*VLuB_cap
    dMSPBb  <- kon_alb*CalbSPBf*CSPBf*VSpB_cap - koff_alb*CSPBb*VSpB_cap
    dMHBb   <- kon_alb*CalbHBf*CHBf*VHtB_cap - koff_alb*CHBb*VHtB_cap
    dMBrBb  <- kon_alb*CalbBrBf*CBrBf*VBrB_cap - koff_alb*CBrBb*VBrB_cap
    dMGoBb  <- kon_alb*CalbGoBf*CGoBf*VGoB_cap - koff_alb*CGoBb*VGoB_cap
    dMSKBb  <- kon_alb*CalbSKBf*CSKBf*VSkB_cap - koff_alb*CSKBb*VSkB_cap
    dMBoBb  <- kon_alb*CalbBoBf*CBoBf*VBoB_cap - koff_alb*CBoBb*VBoB_cap
    dMKFb   <- kon_alb*CalbKFf*CKFf*VKiF - koff_alb*CKFb*VKiF
    dMLFb   <- kon_alb*CalbLFf*CLFf*VLiF - koff_alb*CLFb*VLiF
    dMSTFb  <- kon_alb*CalbSTFf*CSTFf*VStF - koff_alb*CSTFb*VStF
    dMINFb  <- kon_alb*CalbINFf*CINFf*VInF - koff_alb*CINFb*VInF
    dMMFb   <- kon_alb*CalbMFf*CMFf*VMuF - koff_alb*CMFb*VMuF
    dMAFb   <- kon_alb*CalbAFf*CAFf*VAdF - koff_alb*CAFb*VAdF
    dMRFb   <- kon_alb*CalbRFf*CRFf*VReF - koff_alb*CRFb*VReF
    dMLuFb  <- kon_alb*CalbLuFf*CLuFf*VLuF - koff_alb*CLuFb*VLuF
    dMSPFb  <- kon_alb*CalbSPFf*CSPFf*VSpF - koff_alb*CSPFb*VSpF
    dMHFb   <- kon_alb*CalbHFf*CHFf*VHtF - koff_alb*CHFb*VHtF
    dMBrFb  <- kon_alb*CalbBrFf*CBrFf*VBrF - koff_alb*CBrFb*VBrF
    dMGoFb  <- kon_alb*CalbGoFf*CGoFf*VGoF - koff_alb*CGoFb*VGoF
    dMSKFb  <- kon_alb*CalbSKFf*CSKFf*VSkF - koff_alb*CSKFb*VSkF
    dMBoFb  <- kon_alb*CalbBoFf*CBoFf*VBoF - koff_alb*CBoFb*VBoF
    dMPTCb  <- kon_fabp*CFabpPTCf*CPTCf*VPTC - koff_fabp*CPTCb*VPTC
    dMLHCb  <- kon_fabp*CFabpLHCf*CLHCf*VLHC - koff_fabp*CLHCb*VLHC
    dMDTCb  <- kon_fabp*CFabpDTCf*CDTCf*VDTC - koff_fabp*CDTCb*VDTC
    dMCDCb  <- kon_fabp*CFabpCDCf*CCDCf*VCDC - koff_fabp*CCDCb*VCDC
    dMLTb   <- kon_fabp*CFabpLTf*CLTf*VLiT   - koff_fabp*CLTb*VLiT
    dMLuAFb <- kon_alb*CalbLuAFf*CLuAFf*VLuAF - koff_alb*CLuAFb*VLuAF

    # Free mass ODEs
    dMArtf <- QBLu*CLuBf - CArtf*(QBK+QBL+QBM+QBA+QBR+QBSP+QBH+QBBr+
                                    QBST+QBIN+QBGo+QBSK+QBBo) - QGFR*CArtf +
      koff_alb*CArtb*VArt - kon_alb*CalbArtf*CArtf*VArt

    dMVenf <- kUAB*CUA - CVenf*QBLu + QBK*CKBf + QBLtot*CLBf + QBM*CMBf + QBA*CABf +
      QBR*CRBf + QBH*CHBf + QBBr*CBrBf + QBGo*CGoBf + QBSK*CSKBf + QBBo*CBoBf +
      koff_alb*CVenb*VVen - kon_alb*CalbVenf*CVenf*VVen

    dMKBf  <- QBK*CArtf - QBK*CKBf - (Ptrans_diff_K+PparaKi)*A_peritubular*(CKBf-CKFf) +
      (VmK_baso*CPTCf/(KmK_baso+CPTCf)) +
      koff_alb*CKBb*VKiB_cap - kon_alb*CalbKBf*CKBf*VKiB_cap

    dMKFf  <- (Ptrans_diff_K+PparaKi)*A_peritubular*(CKBf-CKFf) -
      kPtcF*(CKFf-CPTCf) - kDalcF*(CKFf-CLHCf) -
      kDtcF*(CKFf-CDTCf) - kCdcF*(CKFf-CCDCf) +
      koff_alb*CKFb*VKiF - kon_alb*CalbKFf*CKFf*VKiF -
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) - (VmK_Oat3*CKFf/(KmK_Oat3+CKFf))

    dMPTCf <- kPtcF*(CKFf-CPTCf) - kPtcTu*(CPTCf-CPT) +
      (VmK_Oatp*CPT/(KmK_Oatp+CPT)) + (VmK_Urat*CPT/(KmK_Urat+CPT)) +
      (VmK_Oat1*CKFf/(KmK_Oat1+CKFf)) + (VmK_Oat3*CKFf/(KmK_Oat3+CKFf)) -
      (VmK_baso*CPTCf/(KmK_baso+CPTCf)) - (VmK_api*CPTCf/(KmK_api+CPTCf)) -
      (kon_fabp*CFabpPTCf*CPTCf*VPTC - koff_fabp*CPTCb*VPTC)

    dMLHCf <- kDalcF*(CKFf-CLHCf) - kDalcTu*(CLHCf-CLH) -
      (kon_fabp*CFabpLHCf*CLHCf*VLHC - koff_fabp*CLHCb*VLHC)

    dMDTCf <- kDtcF*(CKFf-CDTCf) - kDtcTu*(CDTCf-CDT) -
      (kon_fabp*CFabpDTCf*CDTCf*VDTC - koff_fabp*CDTCb*VDTC)

    dMCDCf <- kCdcF*(CKFf-CCDCf) - kCdcTu*(CCDCf-CCD) -
      (kon_fabp*CFabpCDCf*CCDCf*VCDC - koff_fabp*CCDCb*VCDC)

    dMPT   <- QGFR*CArtf + kPtcTu*(CPTCf-CPT) -
      (VmK_Oatp*CPT/(KmK_Oatp+CPT)) - (VmK_Urat*CPT/(KmK_Urat+CPT)) +
      (VmK_api*CPTCf/(KmK_api+CPTCf)) - QLH_tot*CPT

    dMLH   <- QLH_tot*CPT + kDalcTu*(CLHCf-CLH) - QDT_tot*CLH
    dMDT   <- QDT_tot*CLH + kDtcTu*(CDTCf-CDT)  - QCD_tot*CDT
    dMCD   <- QCD_tot*CDT + kCdcTu*(CCDCf-CCD)  - Qurine*CCD
    dMBladder <- Qurine*CCD - Qelim*CBladder

    dMLBf  <- QBL*CArtf + QBSP*CSPBf + QBIN*CINBf + QBST*CSTBf -
      QBLtot*CLBf - (Ptrans_diff_L+PparaLi)*AL*(CLBf-CLFf) +
      koff_alb*CLBb*VLiB_cap - kon_alb*CalbLBf*CLBf*VLiB_cap

    dMLFf  <- (Ptrans_diff_L+PparaLi)*AL*(CLBf-CLFf) - kLFLT*(CLFf-CLTf) -
      (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) - (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf)) -
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) +
      koff_alb*CLFb*VLiF - kon_alb*CalbLFf*CLFf*VLiF

    dMLTf  <- kLFLT*(CLFf-CLTf) +
      (VmL_Oatp*CLFf/(KmL_Oatp+CLFf)) + (VmL_Oatp2*CLFf/(KmL_Oatp2+CLFf)) +
      (VmL_Ntcp*CLFf/(KmL_Ntcp+CLFf)) +
      koff_fabp*CLTb*VLiT - kon_fabp*CFabpLTf*CLTf*VLiT - CL_hepatobiliary*CLTf

    dMBile <- CL_hepatobiliary*CLTf - CBile*Qbile

    dMSTBf <- QBST*CArtf - QBST*CSTBf - (Ptrans_diff_ST+PparaSt)*AST*(CSTBf-CSTFf) +
      koff_alb*CSTBb*VStB_cap - kon_alb*CalbSTBf*CSTBf*VStB_cap
    dMSTFf <- (Ptrans_diff_ST+PparaSt)*AST*(CSTBf-CSTFf) - kSTFSTT*(CSTFf-CSTTf) +
      koff_alb*CSTFb*VStF - kon_alb*CalbSTFf*CSTFf*VStF
    dMSTTf <- kSTFSTT*(CSTFf-CSTTf) + kabST*CSTL
    dMSTL  <- -QGE*MSTL - kabST*CSTL + CLEua*CUA

    dMINBf <- QBIN*CArtf - QBIN*CINBf - (Ptrans_diff_IN+PparaIn)*AIN*(CINBf-CINFf) +
      koff_alb*CINBb*VInB_cap - kon_alb*CalbINBf*CINBf*VInB_cap
    dMINFf <- (Ptrans_diff_IN+PparaIn)*AIN*(CINBf-CINFf) - kINFINT*(CINFf-CINTf) +
      koff_alb*CINFb*VInF - kon_alb*CalbINFf*CINFf*VInF
    dMINTf <- kINFINT*(CINFf-CINTf) + k_gut_in*CINL - k_gut_out*CINTf +
      (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))
    dMINL  <- QGE*MSTL - (CLfeces*CINL) - k_gut_in*CINL + k_gut_out*CINTf +
      CBile*Qbile - (VmIn_Oatp2*CINL/(KmIn_Oatp2+CINL))

    dMMBf  <- QBM*CArtf - QBM*CMBf - (Ptrans_diff_M+PparaMu)*AM*(CMBf-CMFf) +
      koff_alb*CMBb*VMuB_cap - kon_alb*CalbMBf*CMBf*VMuB_cap
    dMMFf  <- (Ptrans_diff_M+PparaMu)*AM*(CMBf-CMFf) - kMFMT*(CMFf-CMTf) +
      koff_alb*CMFb*VMuF - kon_alb*CalbMFf*CMFf*VMuF
    dMMTf  <- kMFMT*(CMFf-CMTf)

    dMABf  <- QBA*CArtf - QBA*CABf - (Ptrans_diff_A+PparaAd)*AA*(CABf-CAFf) +
      koff_alb*CABb*VAdB_cap - kon_alb*CalbABf*CABf*VAdB_cap
    dMAFf  <- (Ptrans_diff_A+PparaAd)*AA*(CABf-CAFf) - kAFAT*(CAFf-CATf) +
      koff_alb*CAFb*VAdF - kon_alb*CalbAFf*CAFf*VAdF
    dMATf  <- kAFAT*(CAFf-CATf)

    dMRBf  <- QBR*CArtf - QBR*CRBf - (Ptrans_diff_R+PparaRe)*AR*(CRBf-CRFf) +
      koff_alb*CRBb*VReB_cap - kon_alb*CalbRBf*CRBf*VReB_cap
    dMRFf  <- (Ptrans_diff_R+PparaRe)*AR*(CRBf-CRFf) - kRFRT*(CRFf-CRTf) +
      koff_alb*CRFb*VReF - kon_alb*CalbRFf*CRFf*VReF
    dMRTf  <- kRFRT*(CRFf-CRTf)

    dMUA   <- -kUAB*CUA - CLEua*CUA

    dMLuBf <- CVenf*QBLu - QBLu*CLuBf - (Ptrans_diff_Lu+PparaLu)*ALu*(CLuBf-CLuFf) +
      koff_alb*CLuBb*VLuB_cap - kon_alb*CalbLuBf*CLuBf*VLuB_cap
    dMLuFf <- (Ptrans_diff_Lu+PparaLu)*ALu*(CLuBf-CLuFf) + kLuTLuF*(CLuTf-CLuFf) +
      koff_alb*CLuFb*VLuF - kon_alb*CalbLuFf*CLuFf*VLuF -
      (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf))
    dMLuTf <- -kLuTLuF*(CLuTf-CLuFf) - kLuTLuAF*(CLuTf-CLuAFf) +
      (VmLu_Oatp_bas*CLuFf/(KmLu_Oatp_bas+CLuFf)) +
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf))
    dMLuAFf <- k_desorption_fast*MLuAFdust_fast + k_desorption_slow*MLuAFdust_slow +
      kLuTLuAF*(CLuTf-CLuAFf) +
      koff_alb*CLuAFb*VLuAF - kon_alb*CalbLuAFf*CLuAFf*VLuAF -
      (VmLu_Oatp_ap*CLuAFf/(KmLu_Oatp_ap+CLuAFf))
    dMLuAFdust_fast <- -k_desorption_fast*MLuAFdust_fast
    dMLuAFdust_slow <- -k_desorption_slow*MLuAFdust_slow

    dMSPBf <- QBSP*CArtf - QBSP*CSPBf - (Ptrans_diff_SP+PparaSp)*ASP*(CSPBf-CSPFf) +
      koff_alb*CSPBb*VSpB_cap - kon_alb*CalbSPBf*CSPBf*VSpB_cap
    dMSPFf <- (Ptrans_diff_SP+PparaSp)*ASP*(CSPBf-CSPFf) - kSPFSPT*(CSPFf-CSPTf) +
      koff_alb*CSPFb*VSpF - kon_alb*CalbSPFf*CSPFf*VSpF
    dMSPTf <- kSPFSPT*(CSPFf-CSPTf)

    dMHBf  <- QBH*CArtf - QBH*CHBf - (Ptrans_diff_H+PparaHt)*AH*(CHBf-CHFf) +
      koff_alb*CHBb*VHtB_cap - kon_alb*CalbHBf*CHBf*VHtB_cap
    dMHFf  <- (Ptrans_diff_H+PparaHt)*AH*(CHBf-CHFf) - kHFHT*(CHFf-CHTf) +
      koff_alb*CHFb*VHtF - kon_alb*CalbHFf*CHFf*VHtF
    dMHTf  <- kHFHT*(CHFf-CHTf)

    dMBrBf <- QBBr*CArtf - QBBr*CBrBf - (Ptrans_diff_Br+PparaBr)*ABr*(CBrBf-CBrFf) +
      koff_alb*CBrBb*VBrB_cap - kon_alb*CalbBrBf*CBrBf*VBrB_cap
    dMBrFf <- (Ptrans_diff_Br+PparaBr)*ABr*(CBrBf-CBrFf) - kBrFBrT*(CBrFf-CBrTf) +
      koff_alb*CBrFb*VBrF - kon_alb*CalbBrFf*CBrFf*VBrF
    dMBrTf <- kBrFBrT*(CBrFf-CBrTf)

    dMGoBf <- QBGo*CArtf - QBGo*CGoBf - (Ptrans_diff_Go+PparaGo)*AGo*(CGoBf-CGoFf) +
      koff_alb*CGoBb*VGoB_cap - kon_alb*CalbGoBf*CGoBf*VGoB_cap
    dMGoFf <- (Ptrans_diff_Go+PparaGo)*AGo*(CGoBf-CGoFf) - kGoFGoT*(CGoFf-CGoTf) +
      koff_alb*CGoFb*VGoF - kon_alb*CalbGoFf*CGoFf*VGoF
    dMGoTf <- kGoFGoT*(CGoFf-CGoTf)

    dMSKBf <- QBSK*CArtf - QBSK*CSKBf - (Ptrans_diff_SK+PparaSk)*ASK*(CSKBf-CSKFf) +
      koff_alb*CSKBb*VSkB_cap - kon_alb*CalbSKBf*CSKBf*VSkB_cap
    dMSKFf <- (Ptrans_diff_SK+PparaSk)*ASK*(CSKBf-CSKFf) - kSKFSKT*(CSKFf-CSKTf) +
      koff_alb*CSKFb*VSkF - kon_alb*CalbSKFf*CSKFf*VSkF
    dMSKTf <- kSKFSKT*(CSKFf-CSKTf)

    dMBoBf <- QBBo*CArtf - QBBo*CBoBf - (Ptrans_diff_Bo+PparaBo)*ABo*(CBoBf-CBoFf) +
      koff_alb*CBoBb*VBoB_cap - kon_alb*CalbBoBf*CBoBf*VBoB_cap
    dMBoFf <- (Ptrans_diff_Bo+PparaBo)*ABo*(CBoBf-CBoFf) - kBoFBoT*(CBoFf-CBoTf) +
      koff_alb*CBoFb*VBoF - kon_alb*CalbBoFf*CBoFf*VBoF
    dMBoTf <- kBoFBoT*(CBoFf-CBoTf)

    dMfeces <- CLfeces*CINL
    dMurine <- Qelim*CBladder
    dVurine <- Qurine
    dVfeces <- Qfeces

    # Derived concentrations
    Cblood    <- (MVen + MArt) / (VVen + VArt)
    Cplasma   <- Cblood / (1-Hct)
    background_plasma  <- 1.48
    Cblood_background  <- background_plasma * (1 - Hct)
    Mblood    <- MVen + MArt
    Mkidney   <- MKB + MKF + MKT + Mfil
    Mliver    <- MLB + MLF + MLT + MBile
    Mbrain    <- MBrB + MBrF + MBrT
    Mdust     <- MLuAFdust_slow + MLuAFdust_fast
    Ckidney   <- (MKB + MKF + MKT + Mfil) / VKi_tot
    Cliver    <- (MLB + MLF + MLT + MBile) / VLi_tot
    Cintestine <- (MINB + MINF + MINT + MINL) / VIn_tot
    Cstomach  <- (MSTB + MSTF + MSTT + MSTL) / VSt_tot
    Cmuscle   <- (MMB + MMF + MMT) / VMu_tot
    Cadipose  <- (MAB + MAF + MAT) / VAd_tot
    Clungs    <- (MLuB + MLuF + MLuT + MLuAF + Mdust) / VLu_tot
    Clungtissue <- (MLuB + MLuF + MLuT + Mdust*fdust_retained) / (VLuB_cap+VLuF+VLuT)
    CUpperair <- MUA / VUA
    CalveolarLF <- (MLuAF + MLuAFdust_slow + MLuAFdust_fast) / VLuAF
    VBALF_Gustaffson <- 0.005
    CBALF     <- (MLuAF + Mdust*(1-fdust_retained)) / VBALF_Gustaffson
    Cspleen   <- (MSPB + MSPF + MSPT) / VSp_tot
    Cheart    <- (MHB  + MHF  + MHT)  / VHt_tot
    Cbrain    <- (MBrB + MBrF + MBrT) / VBr_tot
    Cgonads   <- (MGoB + MGoF + MGoT) / VGo_tot
    Cskin     <- (MSKB + MSKF + MSKT) / VSk_tot
    Cbones    <- (MBoB + MBoF + MBoT) / VBo_tot
    Crest     <- (MRB  + MRF  + MRT)  / VRe_tot
    Ccarcass  <- (MMB+MMF+MMT+MAB+MAF+MAT+MRB+MRF+MRT+MBoB+MBoF+MBoT+MSKB+MSKF+MSKT) /
                 (VMu_tot+VAd_tot+VRe_tot+VBo_tot+VSk_tot)
    Cfeces    <- Mfeces / (Vfeces * feces_density)
    Curine    <- Murine / Vurine

    list(c(
      'dCalbVenf'=dCalbVenf, 'dCalbArtf'=dCalbArtf,
      'dCalbKBf'=dCalbKBf, 'dCalbLBf'=dCalbLBf, 'dCalbSTBf'=dCalbSTBf,
      'dCalbINBf'=dCalbINBf, 'dCalbMBf'=dCalbMBf, 'dCalbABf'=dCalbABf,
      'dCalbRBf'=dCalbRBf, 'dCalbLuBf'=dCalbLuBf, 'dCalbSPBf'=dCalbSPBf,
      'dCalbHBf'=dCalbHBf, 'dCalbBrBf'=dCalbBrBf, 'dCalbGoBf'=dCalbGoBf,
      'dCalbSKBf'=dCalbSKBf, 'dCalbBoBf'=dCalbBoBf,
      'dCalbKFf'=dCalbKFf, 'dCalbLFf'=dCalbLFf, 'dCalbSTFf'=dCalbSTFf,
      'dCalbINFf'=dCalbINFf, 'dCalbMFf'=dCalbMFf, 'dCalbAFf'=dCalbAFf,
      'dCalbRFf'=dCalbRFf, 'dCalbLuFf'=dCalbLuFf, 'dCalbSPFf'=dCalbSPFf,
      'dCalbHFf'=dCalbHFf, 'dCalbBrFf'=dCalbBrFf, 'dCalbGoFf'=dCalbGoFf,
      'dCalbSKFf'=dCalbSKFf, 'dCalbBoFf'=dCalbBoFf,
      'dCFabpPTCf'=dCFabpPTCf, 'dCFabpLHCf'=dCFabpLHCf,
      'dCFabpDTCf'=dCFabpDTCf, 'dCFabpCDCf'=dCFabpCDCf,
      'dCFabpLTf'=dCFabpLTf, 'dCalbLuAFf'=dCalbLuAFf,
      'dMVenb'=dMVenb, 'dMArtb'=dMArtb, 'dMKBb'=dMKBb,
      'dMLBb'=dMLBb, 'dMSTBb'=dMSTBb, 'dMINBb'=dMINBb, 'dMMBb'=dMMBb,
      'dMABb'=dMABb, 'dMRBb'=dMRBb, 'dMLuBb'=dMLuBb, 'dMSPBb'=dMSPBb,
      'dMHBb'=dMHBb, 'dMBrBb'=dMBrBb, 'dMGoBb'=dMGoBb,
      'dMSKBb'=dMSKBb, 'dMBoBb'=dMBoBb,
      'dMKFb'=dMKFb, 'dMLFb'=dMLFb, 'dMSTFb'=dMSTFb, 'dMINFb'=dMINFb,
      'dMMFb'=dMMFb, 'dMAFb'=dMAFb, 'dMRFb'=dMRFb, 'dMLuFb'=dMLuFb,
      'dMSPFb'=dMSPFb, 'dMHFb'=dMHFb, 'dMBrFb'=dMBrFb,
      'dMGoFb'=dMGoFb, 'dMSKFb'=dMSKFb, 'dMBoFb'=dMBoFb,
      'dMPTCb'=dMPTCb, 'dMLHCb'=dMLHCb, 'dMDTCb'=dMDTCb,
      'dMCDCb'=dMCDCb, 'dMLTb'=dMLTb, 'dMLuAFb'=dMLuAFb,
      'dMArtf'=dMArtf, 'dMVenf'=dMVenf, 'dMKBf'=dMKBf,
      'dMKFf'=dMKFf, 'dMPTCf'=dMPTCf, 'dMLHCf'=dMLHCf,
      'dMDTCf'=dMDTCf, 'dMCDCf'=dMCDCf,
      'dMPT'=dMPT, 'dMLH'=dMLH, 'dMDT'=dMDT, 'dMCD'=dMCD,
      'dMBladder'=dMBladder, 'dMLBf'=dMLBf,
      'dMLFf'=dMLFf, 'dMLTf'=dMLTf, 'dMBile'=dMBile,
      'dMSTBf'=dMSTBf, 'dMSTFf'=dMSTFf, 'dMSTTf'=dMSTTf, 'dMSTL'=dMSTL,
      'dMINBf'=dMINBf, 'dMINFf'=dMINFf, 'dMINTf'=dMINTf, 'dMINL'=dMINL,
      'dMMBf'=dMMBf, 'dMMFf'=dMMFf, 'dMMTf'=dMMTf,
      'dMABf'=dMABf, 'dMAFf'=dMAFf, 'dMATf'=dMATf,
      'dMRBf'=dMRBf, 'dMRFf'=dMRFf, 'dMRTf'=dMRTf,
      'dMUA'=dMUA, 'dMLuBf'=dMLuBf, 'dMLuFf'=dMLuFf, 'dMLuTf'=dMLuTf,
      'dMLuAFf'=dMLuAFf, 'dMLuAFdust_fast'=dMLuAFdust_fast, 'dMLuAFdust_slow'=dMLuAFdust_slow,
      'dMSPBf'=dMSPBf, 'dMSPFf'=dMSPFf, 'dMSPTf'=dMSPTf,
      'dMHBf'=dMHBf, 'dMHFf'=dMHFf, 'dMHTf'=dMHTf,
      'dMBrBf'=dMBrBf, 'dMBrFf'=dMBrFf, 'dMBrTf'=dMBrTf,
      'dMGoBf'=dMGoBf, 'dMGoFf'=dMGoFf, 'dMGoTf'=dMGoTf,
      'dMSKBf'=dMSKBf, 'dMSKFf'=dMSKFf, 'dMSKTf'=dMSKTf,
      'dMBoBf'=dMBoBf, 'dMBoFf'=dMBoFf, 'dMBoTf'=dMBoTf,
      'dMfeces'=dMfeces, 'dMurine'=dMurine, 'dVfeces'=dVfeces, 'dVurine'=dVurine
    ),
    'Cblood'=Cblood, 'Cplasma'=Cplasma, 'Cblood_background'=Cblood_background,
    'Mblood'=Mblood, 'Mkidney'=Mkidney, 'Mliver'=Mliver, 'Mbrain'=Mbrain,
    'Ckidney'=Ckidney, 'Cliver'=Cliver, 'Cintestine'=Cintestine,
    'Cstomach'=Cstomach, 'Cmuscle'=Cmuscle, 'Cadipose'=Cadipose,
    'Clungs'=Clungs, 'Clungtissue'=Clungtissue, 'CUpperair'=CUpperair,
    'CalveolarLF'=CalveolarLF, 'CBALF'=CBALF,
    'Cspleen'=Cspleen, 'Cheart'=Cheart, 'Cbrain'=Cbrain,
    'Cgonads'=Cgonads, 'Cskin'=Cskin, 'Cbones'=Cbones, 'Crest'=Crest,
    'Ccarcass'=Ccarcass, 'CBile'=CBile,
    'CBladder'=CBladder,
    'Cfeces'=Cfeces, 'Curine'=Curine
    )
  })
}

# =============================================================================
# 6. User input (default: Abraham et al. 2024 oral exposure)
# =============================================================================

user.input <- list(
  'BW'         = 82,
  'sex'        = 'M',
  'admin.dose' = 3.96,
  'admin.time' = 0,
  'admin.type' = 'oral'
)

params <- create.params(user.input)
inits  <- create.inits(params)
events <- create.events(params)

sim.start <- 0
sim.end   <- 450 * 24
sim.step  <- 6

sample_time <- seq(sim.start, sim.end, sim.step)

solution <- as.data.frame(
  deSolve::ode(
    times       = sample_time,
    func        = ode.func,
    y           = inits,
    parms       = params,
    events      = events,
    method      = "lsodes", rtol = 1e-07, atol = 1e-07,
    custom.func = custom.func
  )
)

# =============================================================================
# 7. Upload to Jaqpot
# =============================================================================

predicted.feats <- c(
  'Cblood', 'Cplasma', 'Cblood_background',
  'Ckidney', 'Cliver', 'Cintestine', 'Cstomach', 'Cmuscle', 'Cadipose',
  'Clungs', 'Clungtissue', 'Cspleen', 'Cheart', 'Cbrain',
  'Cgonads', 'Cskin', 'Cbones', 'Crest', 'Ccarcass',
  'Murine', 'Mfeces', 'CBladder', 'Vurine', 'Vfeces',
  'CBile', 'CalveolarLF', 'CUpperair', 'CBALF',
  'Cfeces', 'Curine',
  'Mblood', 'Mkidney', 'Mliver', 'Mbrain'
)

# jaqpotr::login.api(login.url = "https://api.jaqpot.org")

jaqpotr::deploy.pbpk(
  user.input    = user.input,
  out.vars      = predicted.feats,
  create.params = create.params,
  create.inits  = create.inits,
  create.events = create.events,
  custom.func   = custom.func,
  ode.fun       = ode.func,
  envFile       = "/Users/vassilis/Desktop/jaqpot.env"
)

