library(deSolve)
library(tidyverse)

#=========================
# 1. Parameters of the model
#=========================

#'
#' @export
#'
.create.params <- function(user_input) {
  with(as.list(user_input), {
    variables <- data.frame(
      PFHpA   = c(364.06,   0.0045,   3300, 263,   4500, 60, 0.032437163, 0.000215127, 0.001693061, 3.8,   3.076585024),
      PFOA    = c(414.07,   0.00245,  2500, 137.5, 4500, 47, 0.008479245, 0.005774641, 0.003513515, 3.96,  9.38649045),
      PFNA    = c(464.08,   0.00187,  1000, 42,    8500, 58, 0,           0.001289201, 0.014197676, 4.02,  1.75136989),
      PFDA    = c(514.09,   0.000618, 2500, 137.5, 6000, 39, 0.002874815, 0.000581054, 0.01291839,  4.01,  0.050889826),
      PFBS    = c(300.1,    0.0365,   2500, 137.5, 4500, 47, 0.000247027, 0.000554781, 0,           18.42, 7.378217792),
      PFHxS   = c(403.11,   0.000695, 5800, 778,   7300, 92, 0.728885662, 0.001637872, 0,           3.76,  0.60841209),
      PFOS    = c(500.13,   0.000753, 2500, 137.5, 6000, 39, 0.004273373, 0.002022908, 0,           3.78,  1.076224703),
      DONA    = c(378.07,   0.0193,   2500, 137.5, 4500, 39, 0.000500224, 8.22E-05,    0.0168116,   3.71,  23.44333),
      HFPO_DA = c(330.0489, 0.0491,   2500, 137.5, 4500, 47, 2.77E-01,    3.21E-06,    5.85E-02,    19.58, 2.45138),
      PFBA    = c(214.04,   0.229,    2500, 137.5, 4500, 47, 0,           2.261705e-04, 0,           4.01,  3.618430e+01),
      PFHxA   = c(314.05,   0.038,    2500, 137.5, 4500, 47, 0.004975478, 0,           0,           3.99,  0)
    ) %>% select(all_of(chemical))

    rownames(variables) <- c(
      'MW', 'Free', 'Vmax_baso_invitro', 'Km_baso', 'Vmax_apical_invitro',
      'Km_apical', 'RAFbaso', 'RAFapi', 'kbilec', 'admin_dose_bolus', 'keffluxc'
    )

    PC <- data.frame(
      PFBA    = c(0.062419681, 0.966507353, 0.030790426, 0.268705326, 0.244400571, 0.314135756,
                  0.345043937, 0.349385162, 0.583201906, 0.075402896, 1, 0.294299611, 0.365478649),
      PFHxA   = c(0.071697159, 0.983691794, 0.053489289, 0.477477294, 0.210055902, 0.535045747,
                  0.505758374, 0.259633774, 0.801059944, 0.128366805, 1, 0.383746016, 0.581766954),
      PFHpA   = c(0.072448616, 0.985020595, 0.013529708, 0.520070803, 0.45599169,  0.574362599,
                  0.554298406, 0.575452047, 0.829094812, 0.137900625, 1, 0.391210366, 0.614488429),
      PFOA    = c(0.062389547, 0.977705241, 0.038844366, 0.397833309, 0.350056979, 0.450153313,
                  0.646887159, 0.228481278, 0.736163418, 0.107021188, 1, 0.303240549, 0.481657533),
      PFNA    = c(0.162230058, 0.953539331, 0.634659348, 0.631468305, 0.624946924, 0.610807936,
                  0.627663989, 0.152975454, 0.709308751, 0.153546363, 1, 0.635624712, 0.628976366),
      PFDA    = c(0.229902,    0.917519484, 0.912564874, 0.821027382, 0.825975253, 0.746101202,
                  0.826144701, 0.1955254,   0.747495821, 0.1913061,   1, 0.868316352, 0.771226649),
      HFPO_DA = c(0.063405738, 0.954105857, 0.159566058, 0.320076723, 0.286706334, 0.346606884,
                  1.353802165, 0.710753658, 0.613353634, 0.083633611, 1, 0.284121859, 0.387329912),
      DONA    = c(0.061417029, 0.98308624,  0.091837315, 0.446218328, 0.387380788, 0.507789131,
                  0.483562817, 0.127921824, 0.789313478, 0.120786751, 1, 0.327689386, 0.544782562),
      PFBS    = c(0.09553098,  0.965214629, 0.338484978, 0.414526508, 0.394411629, 0.434619355,
                  0.438104109, 0.442709614, 0.65005241,  0.1063728,   1, 0.395423443, 0.461523875),
      PFHxS   = c(0.062790197, 0.988369127, 0.039344291, 0.425509501, 0.5611681,   0.498282656,
                  0.448001176, 0.241904006, 0.786332564, 0.117538859, 1, 0.314848005, 0.52383342),
      PFOS    = c(0.289029898, 0.926440135, 0.288439657, 1.038510649, 0.262099455, 0.948121846,
                  1.044334769, 0.496366704, 0.889318484, 0.243172931, 1, 1.09828647,  0.973073472)
    ) %>% select(all_of(chemical))

    rownames(PC) <- c(
      'Adipose', 'Blood', 'Brain', 'Gonads', 'Gut', 'Heart', 'Kidney', 'Liver',
      'Lung', 'Muscle', 'Plasma', 'Skin', 'Spleen'
    )

    if (missing(time_scale) || is.na(time_scale)) {
      time_scale <- 1/24
    } else if (time_scale == "minutes") {
      time_scale <- 60
    } else if (time_scale == "hours") {
      time_scale <- 1
    } else if (time_scale == "days") {
      time_scale <- 1/24
    } else if (time_scale == "weeks") {
      time_scale <- (1/24) / 7
    } else if (time_scale == "months") {
      time_scale <- (1/24) / 30
    } else if (time_scale == "years") {
      time_scale <- (1/24) / 365
    }
    inv_time_scale <- 1 / time_scale

    # --- Cardiac Output and Blood Flow (as fraction of cardiac output) ---
    QCC   <- 12.5 * inv_time_scale  # cardiac output in L/day/kg^0.75; Brown 1997
    QLC   <- 0.065                  # fraction blood flow to liver; Brown 1997
    QKC   <- 0.175                  # fraction blood flow to kidney; Brown 1997
    QAdiC <- 0.05                   # fraction blood flow to adipose tissue; Brown 1997
    QBraC <- 0.12                   # fraction blood flow to brain; Brown 1997
    QGonC <- 0.0005                 # fraction blood flow to gonads; ICRP 2002
    QHeaC <- 0.04                   # fraction blood flow to heart; Brown 1997
    QLunC <- 0.025                  # fraction blood flow to lung; Brown 1997
    QMusC <- 0.17                   # fraction blood flow to muscle; Brown 1997
    QSkiC <- 0.05                   # fraction blood flow to skin; Brown 1997
    QSplC <- 0.03                   # fraction blood flow to spleen; ICRP 2002
    QPanC <- 0.01                   # fraction blood flow to pancreas; ICRP 2002
    QGIC  <- 0.15                   # fraction blood flow to GI tract; ICRP 2002

    Htc <- 0.42  # hematocrit

    # --- Tissue Volumes ---
    VplasC <- 0.0428   # fraction vol. of plasma (L/kg BW); Davies 1993
    VLC    <- 0.026    # fraction vol. of liver (L/kg BW); Brown 1997
    VKC    <- 0.004    # fraction vol. of kidney (L/kg BW); Brown 1997
    VAdiC  <- 0.214    # fraction vol. of adipose tissue (L/kg BW); Brown 1997
    VBraC  <- 0.02     # fraction vol. of brain (L/kg BW); Brown 1997
    VGonC  <- 0.0005   # fraction vol. of gonads (L/kg BW)
    VHeaC  <- 0.005    # fraction vol. of heart (L/kg BW); Brown 1997
    VLunC  <- 0.008    # fraction vol. of lung (L/kg BW); Brown 1997
    VMusC  <- 0.4      # fraction vol. of muscle (L/kg BW); Brown 1997
    VSkiC  <- 0.037    # fraction vol. of skin (L/kg BW); Brown 1997
    VSplC  <- 0.002    # fraction vol. of spleen (L/kg BW); ICRP 2002
    VPanC  <- 0.002    # fraction vol. of pancreas (L/kg BW); ICRP 2002
    VGIC   <- 0.014    # fraction vol. of GI tract (L/kg BW); Brown 1997

    VfilC  <- 4e-4     # fraction vol. of filtrate (L/kg BW)
    VPTCC  <- 1.35e-4  # vol. of proximal tubule cells (L/g kidney)

    # --- Chemical Specific Parameters ---
    MW   <- variables["MW", ]    # PFAS molecular mass (g/mol)
    Free <- variables["Free", ]  # free fraction in plasma (Smeltz 2023); fitted to model

    # --- Kidney Transport Parameters ---
    Vmax_baso_invitro   <- variables["Vmax_baso_invitro", ]
    Km_baso             <- variables["Km_baso", ] * variables["MW", ]
    Vmax_apical_invitro <- variables["Vmax_apical_invitro", ]
    Km_apical           <- variables["Km_apical", ] * variables["MW", ]
    RAFbaso             <- variables["RAFbaso", ]
    RAFapi              <- variables["RAFapi", ]
    protein             <- 2.0e-6   # amount of protein in proximal tubule cells (mg protein/cell)
    GFRC                <- 24.19 * inv_time_scale  # glomerular filtration rate (L/day/kg kidney); Corley 2005

    # --- Partition Coefficients (from Allendorf 2021) ---
    PR   <- 0.01              # rest of body:blood; fitted to model
    PAdi <- PC["Adipose", ]
    PBra <- PC["Brain", ]
    PGon <- PC["Gonads", ]
    PGI  <- PC["Gut", ]
    PHea <- PC["Heart", ]
    PL   <- PC["Liver", ]
    PLun <- PC["Lung", ]
    PMus <- PC["Muscle", ]
    PSki <- PC["Skin", ]
    PSpl <- PC["Spleen", ]
    PPan <- (PGI + PSpl) / 2  # estimated as average of GI tract and spleen

    # --- Rate Constants ---
    kdif     <- 0.001 * inv_time_scale       # diffusion rate from proximal tubule cells (L/day)
    kabsc    <- 2.12 * inv_time_scale        # rate of absorption from small intestine (1/(day*BW^-0.25))
    kunabsc  <- 7.06e-5 * inv_time_scale    # rate of unabsorbed dose to feces; fitted to model
    keffluxc <- variables["keffluxc", ] * inv_time_scale
    kbilec   <- variables["kbilec", ] * inv_time_scale
    kurinec  <- 0.063 * inv_time_scale

    # --- Water Consumption ---
    water_consumption <- 1.36  # L/day

    # === SCALED PARAMETERS ===

    # Cardiac output and blood flows
    QC   <- QCC * (BW^0.75) * (1 - Htc)
    QK   <- QKC  * QC
    QL   <- QLC  * QC
    QAdi <- QAdiC * QC
    QBra <- QBraC * QC
    QGon <- QGonC * QC
    QHea <- QHeaC * QC
    QLun <- QLunC * QC
    QMus <- QMusC * QC
    QSki <- QSkiC * QC
    QSpl <- QSplC * QC
    QPan <- QPanC * QC
    QGI  <- QGIC  * QC
    QR   <- QC - QK - QL - QAdi - QBra - QGon - QHea - QLun -
            QMus - QSki - QSpl - QPan - QGI

    QBal <- QC - (QK + QL + QR + QAdi + QBra + QGon + QHea + QLun +
            QMus + QSki + QSpl + QPan + QGI)  # should equal zero

    # Tissue Volumes
    VPlas <- VplasC * BW
    VK    <- VKC   * BW
    MK    <- VK * 1.0 * 1000
    VKb   <- VK * 0.16
    Vfil  <- VfilC * BW
    VL    <- VLC   * BW
    ML    <- VL * 1.05 * 1000
    VAdi  <- VAdiC * BW
    VBra  <- VBraC * BW
    VGon  <- VGonC * BW
    VHea  <- VHeaC * BW
    VLun  <- VLunC * BW
    VMus  <- VMusC * BW
    VSki  <- VSkiC * BW
    VSpl  <- VSplC * BW
    VPan  <- VPanC * BW
    VGI   <- VGIC  * BW

    # Kidney Parameters
    PTC  <- VKC * 1000 * 6e7
    VPTC <- VK * 1000 * VPTCC
    MPTC <- VPTC * 1000
    VR   <- (0.93 * BW) - VPlas - VPTC - Vfil - VL - VAdi - VBra -
            VGon - VHea - VLun - VMus - VSki - VSpl - VPan - VGI
    VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas + VAdi +
            VBra + VGon + VHea + VLun + VMus + VSki + VSpl + VPan + VGI)  # should equal zero

    Vmax_basoC   <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
    Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
    Vmax_baso    <- Vmax_basoC * BW^0.75
    Vmax_apical  <- Vmax_apicalC * BW^0.75
    kbile        <- kbilec * BW^(-0.25)
    kurine       <- kurinec * BW^(-0.25)
    kefflux      <- keffluxc * BW^(-0.25)
    GFR          <- GFRC * VK

    # GI Tract Parameters
    kabs   <- kabsc   * BW^(-0.25)
    kunabs <- kunabsc * BW^(-0.25)

    return(list(
      "BW" = BW, "QBAL" = QBal, "VBAL" = VBal,
      "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR,
      "QAdi" = QAdi, "QBra" = QBra, "QGon" = QGon, "QHea" = QHea,
      "QLun" = QLun, "QMus" = QMus, "QSki" = QSki, "QSpl" = QSpl,
      "QPan" = QPan, "QGI" = QGI,
      "VPlas" = VPlas, "VKb" = VKb, "Vfil" = Vfil, "VL" = VL,
      "VR" = VR, "ML" = ML,
      "VAdi" = VAdi, "VBra" = VBra, "VGon" = VGon, "VHea" = VHea,
      "VLun" = VLun, "VMus" = VMus, "VSki" = VSki, "VSpl" = VSpl, "VPan" = VPan,
      "VGI" = VGI, "MK" = MK,
      "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
      "kdif" = kdif, "Km_baso" = Km_baso, "Km_apical" = Km_apical,
      "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
      "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs,
      "PR" = PR, "PAdi" = PAdi, "PBra" = PBra, "PGon" = PGon, "PGI" = PGI,
      "PHea" = PHea, "PL" = PL, "PLun" = PLun, "PMus" = PMus, "PSki" = PSki,
      "PSpl" = PSpl, "PPan" = PPan, "water_consumption" = water_consumption,
      "admin_type" = admin_type, "admin_dose" = admin_dose,
      "admin_time" = admin_time, "ingestion" = ingestion,
      "ingestion_time" = ingestion_time, "duration" = duration,
      "time_scale" = time_scale, "exp_type" = exp_type
    ))
  })
}

#===============================================
# 2. Function to create initial values for ODEs
#===============================================

#'
#' @export
#'
.create.inits <- function(parameters) {
  with(as.list(parameters), {
    AR = 0; AAdi = 0; ABra = 0; AGon = 0
    AHea = 0; ALun = 0; AMus = 0; ASki = 0; ASpl = 0; APan = 0
    Adif = 0; A_baso = 0; AKb = 0
    ACl = 0; Aefflux = 0
    A_apical = 0; APTC = 0; Afil = 0
    Aurine = 0; ALumen = 0; AGI = 0
    AabsLumen = 0; Afeces = 0
    AL = 0; Abile = 0; Aplas_free = 0
    ingestion = 0

    return(c(
      "AR" = AR, "AAdi" = AAdi, "ABra" = ABra, "AGon" = AGon,
      "AHea" = AHea, "ALun" = ALun, "AMus" = AMus, "ASki" = ASki,
      "ASpl" = ASpl, "APan" = APan,
      "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
      "ACl" = ACl, "Aefflux" = Aefflux,
      "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
      "Aurine" = Aurine, "ALumen" = ALumen, "AGI" = AGI,
      "AabsLumen" = AabsLumen, "Afeces" = Afeces,
      "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free,
      "ingestion" = ingestion
    ))
  })
}

#===================
# 3. Events function
#===================

#'
#' @export
#'
.create.events <- function(parameters) {
  with(as.list(parameters), {
    if (tolower(admin_type) == "iv") {
      ldose  <- length(admin_dose)
      ltimes <- length(admin_time)
      if (ltimes != ldose) {
        stop("The times of administration should be equal in number to the doses")
      }
      events <- list(data = data.frame(
        var    = "Aplas_free",
        time   = admin_time,
        value  = admin_dose,
        method = "add"
      ))

    } else if (admin_type == "bolus") {
      ldose  <- length(admin_dose)
      ltimes <- length(admin_time)
      if (ltimes != ldose) {
        stop("The times of administration should be equal in number to the doses")
      }
      events <- list(data = data.frame(
        var    = "ALumen",
        time   = admin_time,
        value  = admin_dose,
        method = "add"
      ))

    } else if (admin_type == "oral") {
      lingest      <- length(ingestion)
      lingesttimes <- length(ingestion_time)
      if (lingest != lingesttimes) {
        # stop("The times of ingestion rate change should be equal in vector of ingestion")
      }
      if (exp_type == "pharmacokinetics") {
        events <- list(data = data.frame(
          var    = c("ALumen"),
          time   = ingestion_time,
          value  = ingestion,
          method = c("add")
        ))
      } else if (exp_type == "continuous") {
        events <- list(data = data.frame(
          var    = c("ingestion"),
          time   = ingestion_time,
          value  = ingestion,
          method = c("rep")
        ))
      } else {
        stop("admin_type should be either 'iv', 'bolus', or 'oral'")
      }
      return(events)
    }
  })
}

#==================
# 4. Custom functions
#==================

#'
#' @export
#'
.custom.func <- function() {
  return()
}

#==============
# 5. ODEs System
#==============

#'
#' @export
#'
.ode.func <- function(time, inits, params, custom.func) {
  with(as.list(c(inits, params)), {

    # Concentrations in various compartments
    CR      <- AR / VR
    CVR     <- CR / PR
    CR_free <- CR * Free

    CAdi      <- AAdi / VAdi
    CVAdi     <- CAdi / PAdi
    CAdi_free <- CAdi * Free

    CBra      <- ABra / VBra
    CVBra     <- CBra / PBra
    CBra_free <- CBra * Free

    CGon      <- AGon / VGon
    CVGon     <- CGon / PGon
    CGon_free <- CGon * Free

    CHea      <- AHea / VHea
    CVHea     <- CHea / PHea
    CHea_free <- CHea * Free

    CLun      <- ALun / VLun
    CVLun     <- CLun / PLun
    CLun_free <- CLun * Free

    CMus      <- AMus / VMus
    CVMus     <- CMus / PMus
    CMus_free <- CMus * Free

    CSki      <- ASki / VSki
    CVSki     <- CSki / PSki
    CSki_free <- CSki * Free

    CSpl      <- ASpl / VSpl
    CVSpl     <- CSpl / PSpl
    CSpl_free <- CSpl * Free

    Cpan      <- APan / VPan
    CVPan     <- Cpan / PPan
    CPan_free <- Cpan * Free

    CGI      <- AGI / VGI
    CVGI     <- CGI / PGI
    CGI_free <- CGI * Free

    CKb  <- AKb / VKb
    CVK  <- CKb
    CPTC <- APTC / VPTC
    Cfil <- Afil / Vfil

    CL      <- AL / VL
    CLiver  <- AL / ML
    CVL     <- CL / PL
    CL_free <- CL * Free

    CA_free <- Aplas_free / VPlas
    CA      <- CA_free / Free

    # Rest of Body
    dAR <- QR * (CA - CVR) * Free

    # Adipose Tissue
    dAAdi <- QAdi * (CA - CVAdi) * Free

    # Brain Tissue
    dABra <- QBra * (CA - CVBra) * Free

    # Gonads Tissue
    dAGon <- QGon * (CA - CVGon) * Free

    # Heart Tissue
    dAHea <- QHea * (CA - CVHea) * Free

    # Lung Tissue
    dALun <- QLun * (CA - CVLun) * Free

    # Muscle Tissue
    dAMus <- QMus * (CA - CVMus) * Free

    # Skin Tissue
    dASki <- QSki * (CA - CVSki) * Free

    # Spleen Tissue
    dASpl <- QSpl * (CA - CVSpl) * Free

    # Pancreas Tissue
    dAPan <- QPan * (CA - CVPan) * Free

    # Kidney Blood
    dAdif   <- kdif * (CKb - CPTC)
    dA_baso <- (Vmax_baso * CKb) / (Km_baso + CKb)
    dAKb    <- QK * (CA - CVK) * Free - CA * GFR * Free - dAdif - dA_baso
    dACl    <- CA * GFR * Free

    # Proximal Tubule Cells
    dAefflux  <- kefflux * APTC
    dA_apical <- (Vmax_apical * Cfil) / (Km_apical + Cfil)
    dAPTC     <- dAdif + dA_apical + dA_baso - dAefflux

    # Filtrate
    dAfil <- CA * GFR * Free - dA_apical - Afil * kurine

    # Urinary elimination
    dAurine <- kurine * Afil

    # GI Lumen
    dALumen <- ingestion - kabs * ALumen - kunabs * ALumen

    # GI Tract
    dAGI       <- kabs * ALumen + QGI * (CA - CVGI) * Free
    dAabsLumen <- kabs * ALumen

    # Feces
    dAfeces <- kunabs * ALumen + kbile * AL

    # Liver
    dAL <- QL  * (CA   - CVL)  * Free - kbile * AL +
           QGI * (CVGI - CVL)  * Free +
           QSpl * (CVSpl - CVL) * Free +
           QPan * (CVPan - CVL) * Free

    dAbile <- kbile * AL
    amount_per_gram_liver <- CLiver

    # Plasma
    dAplas_free <- (QR   * CVR   * Free) + (QK  * CVK  * Free) + (QL  * CVL  * Free) +
                  (QAdi * CVAdi * Free) + (QBra * CVBra * Free) + (QGon * CVGon * Free) +
                  (QHea * CVHea * Free) + (QLun * CVLun * Free) + (QMus * CVMus * Free) +
                  (QSki * CVSki * Free) + (QSpl * CVL   * Free) + (QPan * CVL   * Free) +
                  (QGI  * CVL   * Free) - (QC   * CA    * Free) + dAefflux

    dingestion <- 0

    # Mass Balance
    Atissue <- Aplas_free + AR + AAdi + ABra + AGon + AHea + ALun +
               AMus + ASki + ASpl + APan + AKb + Afil + ALumen + AGI + APTC + AL
    Aloss   <- Aurine + Afeces
    Atotal  <- Atissue + Aloss

    list(c(
        "dAR" = dAR, "dAAdi" = dAAdi, "dABra" = dABra, "dAGon" = dAGon,
        "dAHea" = dAHea, "dALun" = dALun, "dAMus" = dAMus, "dASki" = dASki,
        "dASpl" = dASpl, "dAPan" = dAPan,
        "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
        "dACl" = dACl, "dAefflux" = dAefflux,
        "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
        "dAurine" = dAurine, "dALumen" = dALumen, "dAGI" = dAGI,
        "dAabsLumen" = dAabsLumen, "dAfeces" = dAfeces,
        "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
        "dingestion" = dingestion
      ),
      "amount_per_gram_liver" = amount_per_gram_liver,
      "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal,
      "CR" = CR, "CR_free" = CR_free, "CVR" = CVR,
      "CAdi" = CAdi, "CAdi_free" = CAdi_free, "CVAdi" = CVAdi,
      "CBra" = CBra, "CBra_free" = CBra_free, "CVBra" = CVBra,
      "CGon" = CGon, "CGon_free" = CGon_free, "CVGon" = CVGon,
      "CHea" = CHea, "CVHea" = CVHea, "CHea_free" = CHea_free,
      "CLun" = CLun, "CVLun" = CVLun, "CLun_free" = CLun_free,
      "CMus" = CMus, "CVMus" = CVMus, "CMus_free" = CMus_free,
      "CSki" = CSki, "CVSki" = CVSki, "CSki_free" = CSki_free,
      "CSpl" = CSpl, "CVSpl" = CVSpl, "CSpl_free" = CSpl_free,
      "Cpan" = Cpan, "CVPan" = CVPan, "Cpan_free" = CPan_free,
      "CKb" = CKb, "CVK" = CVK, "CPTC" = CPTC,
      "Cfil" = Cfil,
      "CL" = CL, "CVL" = CVL, "CL_free" = CL_free,
      "CGI" = CGI, "CVGI" = CVGI, "CGI_free" = CGI_free,
      "CA_free" = CA_free, "CA" = CA
    )
  })
}
