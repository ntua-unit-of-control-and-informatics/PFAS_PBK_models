library(deSolve)
library(tidyverse)
library(ggplot2)

#=========================
# 1. Parameters of the model
#=========================

variables<-read.csv("estimated_parameters.csv",row.names="Parameters")

PC<-read.csv("PCs.csv",row.names="Organs")
# NOTE: Manually adjust parameter values below as needed
create.params <- function(user_input) {
  with(as.list(user_input), {
  # === USER-ADJUSTABLE PARAMETERS ===
  # Change these default values as needed
  BW <- 82 
  # --- Cardiac Output and Blood Flow (as fraction of cardiac output) ---
  QCC <- 12.5 * 24  # cardiac output in L/day/kg^0.75; Brown 1997

  QLC <- 0.065  # fraction blood flow to liver; Brown 1997
  QKC <- 0.175  # fraction blood flow to kidney; Brown 1997
  QAdiC <-0.05  # fraction blood flow to adipose tissue; Brown 1997
  QBraC <-0.12  # fraction blood flow to brain; Brown 1997
  QGonC <-0.0005  # fraction blood flow to gonads; ICRP 2002
  QHeaC <-0.04  # fraction blood flow to heart; Brown 1997
  QLunC <-0.025  # fraction blood flow to lung; Brown 1997
  QMusC <-0.17  # fraction blood flow to muscle; Brown 1997
  QSkiC <-0.05  # fraction blood flow to skin; Brown 1997
  QSplC <-0.03  # fraction blood flow to spleen; ICRP 2002
  QPanC <-0.01  # fraction blood flow to pancreas; ICRP 2002
  QGIC<-0.15  # fraction blood flow to GI tract; ICRP 2002

  Htc <- 0.467  # hematocrit; Abraham's personal value

  # --- Tissue Volumes ---
  VplasC <- 0.0428  # fraction vol. of plasma (L/kg BW); Davies 1993
  VLC <- 0.026  # fraction vol. of liver (L/kg BW); Brown 1997
  VKC <- 0.004  # fraction vol. of kidney (L/kg BW); Brown 1997
  VAdiC <- 0.214  # fraction vol. of adipose tissue (L/kg BW); Brown 1997
  VBraC <- 0.02  # fraction vol. of brain (L/kg BW); Brown 1997
  VGonC <- 0.0005  # fraction vol. of gonads (L/kg BW); 
  VHeaC <- 0.005  # fraction vol. of heart (L/kg BW); Brown 1997
  VLunC <- 0.008  # fraction vol. of lung (L/kg BW); Brown 1997
  VMusC <- 0.4  # fraction vol. of muscle (L/kg BW); Brown 1997
  VSkiC <- 0.037  # fraction vol. of skin (L/kg BW); Brown 1997
  VSplC <- 0.002  # fraction vol. of spleen (L/kg BW); ICRP 2002
  VPanC <- 0.002 # fraction vol. of pancreas (L/kg BW); ICRP 2002
  VGIC <- 0.014  # fraction vol. of GI tract (L/kg BW); Brown 1997

  VfilC <- 4e-4  # fraction vol. of filtrate (L/kg BW)
  VPTCC <- 1.35e-4  # vol. of proximal tubule cells (L/g kidney)

  # --- Chemical Specific Parameters ---
  MW <- variables["MW", "PFOS"]  # PFOS molecular mass (g/mol)
  Free <- variables["Free", "PFOS"]  # free fraction in plasma (Smeltz 2023); fitted to model

# --- Kidney Transport Parameters ---
  protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
  GFRC <- 24.19 * 24  # glomerular filtration rate (L/day/kg kidney); Corley 2005

  
    # --- Water Consumption ---
  water_consumption <- 1.36  # L/day

  # === SCALED PARAMETERS (calculated from above) ===

  # Cardiac output and blood flows
  QC <- QCC * (BW^0.75) * (1 - Htc)  # cardiac output in L/day; adjusted for plasma
  QK <- (QKC * QC)  # plasma flow to kidney (L/day)
  QL <- (QLC * QC)  # plasma flow to liver (L/day)
  QAdi <- (QAdiC * QC)  # plasma flow to adipose tissue (L/day)
  QBra <- (QBraC * QC)  # plasma flow to brain (L/day)
  QGon <- (QGonC * QC)  # plasma flow to gonads (L/day)
  QHea <- (QHeaC * QC)  # plasma flow to heart (L/day)
  QLun <- (QLunC * QC)  # plasma flow to lung (L/day)
  QMus <- (QMusC * QC)  # plasma flow to muscle (L/day)
  QSki <- (QSkiC * QC)  # plasma flow to skin (L/day)
  QSpl <- (QSplC * QC)  # plasma flow to spleen (L/day)
  QPan <- (QPanC * QC)  # plasma flow to pancreas (L/day)
  QGI <- (QGIC * QC)  # plasma flow to GI tract (L/day)
  QR <- QC - QK - QL - QAdi - QBra - QGon - QHea - QLun -
  QMus - QSki - QSpl - QPan - QGI  # plasma flow to rest of body (L/day)

  QBal <- QC - (QK + QL + QR + QAdi + QBra + QGon + QHea + QLun +
  QMus + QSki + QSpl + QPan + QGI)  # Balance check; should equal zero

  # Tissue Volumes
  VPlas <- VplasC * BW  # volume of plasma (L)
  VK <- VKC * BW  # volume of kidney (L)
  MK <- VK * 1.0 * 1000  # mass of the kidney (g)
  VKb <- VK * 0.16  # volume of blood in the kidney (L); Brown 1997
  Vfil <- VfilC * BW  # volume of filtrate (L)
  VL <- VLC * BW  # volume of liver (L)
  ML <- VL * 1.05 * 1000  # mass of the liver (g)
  VAdi <-17.5 #Abraham's individual value  #VAdiC * BW  # volume of adipose tissue (L)
  VBra <- VBraC * BW  # volume of brain (L)
  VGon <- VGonC * BW  # volume of gonads (L)
  VHea <- VHeaC * BW  # volume of heart (L) 
  VLun <- VLunC * BW  # volume of lung (L)
  VMus <- 32.8 #Abraham's individual value #VMusC * BW  # volume of muscle (L)
  VSki <- VSkiC * BW  # volume of skin (L)
  VSpl <- VSplC * BW  # volume of spleen (L)
  VPan <- VPanC * BW  # volume of pancreas (L)
  VGI <- VGIC * BW  # volume of GI tract (L)


  # Kidney Parameters
  PTC <- VKC * 1000 * 6e7  # number of PTC (cells/kg BW)
  VPTC <- VK * 1000 * VPTCC  # volume of proximal tubule cells (L)
  MPTC <- VPTC * 1000  # mass of the proximal tubule cells (g)
  VR <- (0.93 * BW) - VPlas- VPTC - Vfil - VL -VAdi - VBra - VGon - VHea - VLun - VMus - VSki - VSpl - VPan - VGI  # volume of rest of body (L)
  VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas + VAdi + VBra + VGon + VHea + VLun + VMus + VSki + VSpl + VPan + VGI)  # Balance check; should equal zero

  Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
  Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
  Vmax_baso <- Vmax_basoC * BW^0.75  # (ug/day)
  Vmax_apical <- Vmax_apicalC * BW^0.75  # (ug/day)
  kbile <- kbilec * BW^(-0.25)  # biliary elimination; liver to feces storage (/day)
  kurine <- kurinec * BW^(-0.25)  # urinary elimination, from filtrate (/day)
  kefflux <- keffluxc * BW^(-0.25)  # efflux from PTC to blood (/day)
  GFR <- 163.65  # glomerular filtration rate, (L/day)Abraham's personal value

  # GI Tract Parameters
  kabs <- kabsc * BW^(-0.25)  # rate of absorption from small intestine (/day)
  kunabs <- kunabsc * BW^(-0.25)  # rate of unabsorbed dose to feces (/day)
 

  return(list(
    "BW"=BW, "QBAL" = QBal, "VBAL" = VBal,
    "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR,
    "QAdi" = QAdi, "QBra" = QBra, "QGon" = QGon, "QHea" = QHea,
    "QLun" = QLun, "QMus" = QMus, "QSki" = QSki, "QSpl" = QSpl,
    "QPan" = QPan, "QGI" = QGI,
    "VPlas" = VPlas, "VKb" = VKb, "Vfil" = Vfil, "VL" = VL, "VR" = VR, "ML" = ML,
    "VAdi" = VAdi, "VBra" = VBra, "VGon" = VGon, "VHea" = VHea,
    "VLun" = VLun, "VMus" = VMus, "VSki" = VSki, "VSpl" = VSpl, "VPan" = VPan,
    "VGI" = VGI, "MK" = MK,
    "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
    "kdif" = kdif, "Km_baso" = Km_baso, "Km_apical" = Km_apical,
    "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
    "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs, 
    "PR" = PR, "PAdi" = PAdi, "PBra" = PBra, "PGon" = PGon, "PGI" = PGI,
    "PHea" = PHea, "PL" = PL, "PLun" = PLun, "PMus" = PMus, "PSki" = PSki,
    "PSpl" = PSpl, "PPan" = PPan,"water_consumption" = water_consumption))

  })
}
#===============================================
#2. Function to create initial values for ODEs 
#===============================================
create.inits <- function(parameters){
  with( as.list(parameters),{
    "AR" = 0; "AAdi"=0; "ABra"=0; "AGon"=0;
    "AHea"=0; "ALun"=0; "AMus"=0; "ASki"=0; "ASpl"=0; "APan"=0;
    "Adif" = 0; "A_baso" = 0; "AKb" = 0;
    "ACl" = 0; "Aefflux" = 0;
    "A_apical" = 0; "APTC" = 0; "Afil" = 0;
    "Aurine" = 0; "ALumen" = 0; "AGI" = 0;
    "AabsLumen" = 0; "Afeces" = 0;
    "AL" = 0; "Abile" = 0; "Aplas_free" = 0;
    "ingestion" = 0; "Cwater" = 0

    return(c("AR" = AR, "AAdi"=AAdi, "ABra"=ABra, "AGon"=AGon,
             "AHea"=AHea, "ALun"=ALun, "AMus"=AMus, "ASki"=ASki,
             "ASpl"=ASpl, "APan"=APan,
             "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
             "ACl" = ACl, "Aefflux" = Aefflux,
             "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
             "Aurine" = Aurine, "ALumen" = ALumen, "AGI" = AGI,
             "AabsLumen" = AabsLumen, "Afeces" = Afeces, 
             "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free,
              "Cwater" = Cwater,"ingestion" = ingestion))
  })
}

#===================
# 3. Events function
#===================
create.events <- function(user_input) {

  if (admin_type == "iv") {
    ldose <- length(admin_dose_iv)
    ltimes <- length(admin_time_iv)
    if (ltimes != ldose) {
      stop("The times of administration should be equal in number to the doses")
    }
    events <- list(data = data.frame(
      var = "Aplas_free", time = admin_time_iv,
      value = admin_dose_iv, method = "add"
    ))
  } else if (admin_type == "oral") {
    lcwater <- length(Cwater)
    lcwatertimes <- length(Cwater_time)
    lingest <- length(ingestion)
    lingesttimes <- length(ingestion_time)
    if (lcwater != lcwatertimes) {
      stop("The times of water concentration change should be equal in vector of Cwater")
    } else if (lingest != lingesttimes) {
      stop("The times of ingestion rate change should be equal in vector of ingestion")
    }
    events <- list(data = rbind(
      data.frame(var = "Cwater", time = Cwater_time,
                 value = Cwater, method = "rep"),
      data.frame(var = "ingestion", time = ingestion_time,
                 value = ingestion, method = "rep")
    ))
  } else if (admin_type == "bolus") {
    ldose <- length(admin_dose_bolus)
    ltimes <- length(admin_time_bolus)
    if (ltimes != ldose) {
      stop("The times of administration should be equal in number to the doses")
    }
    events <- list(data = data.frame(
      var = "ALumen", time = admin_time_bolus,
      value = admin_dose_bolus, method = "add"
    ))
  }

  return(events)
}

#==================
# 4. Custom functions
#==================
cumulative_exp_data <- function(df, time_col, concentration_col, multiply_col) {
  working_df <- df[, c(time_col, concentration_col, multiply_col)]
  complete_cases <- complete.cases(working_df)
  clean_df <- working_df[complete_cases, ]

  if (nrow(clean_df) == 0) {
    stop("No complete cases found after removing NA values")
  }

  multiplied_concentration <- clean_df[[concentration_col]] * clean_df[[multiply_col]]

  result_df <- data.frame(
    time = clean_df[[time_col]],
    cumulative_mass = cumsum(multiplied_concentration)
  )

  return(result_df)
}

AUC <- function(x, y){
  individual_auc <- c()
  for (i in 1:(length(x)-1)){
    individual_auc[i] <- (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return(sum(individual_auc))
}



#==============
# 5. ODEs System
#==============
ode.func <- function(time, inits, params) {
  with(as.list(c(inits, params)), {

    # Concentrations in various compartments
    
    CR <- AR / VR  # concentration in rest of body (ug/L)
    CVR <- CR / PR  # concentration in venous blood leaving rest of body (ug/L)
    
    CAdi <- AAdi / VAdi  # concentration in adipose tissue (ug/L)
    CVAdi <- CAdi / PAdi  # concentration in venous blood leaving adipose
    
    CBra <- ABra / VBra  # concentration in brain tissue (ug/L)
    CVBra <- CBra / PBra  # concentration in venous blood leaving brain
    
    CGon <- AGon / VGon  # concentration in gonads tissue (ug/L)
    CVGon <- CGon / PGon  # concentration in venous blood leaving gonads
    
    CHea <- AHea / VHea  # concentration in heart tissue (ug/L)
    CVHea <- CHea / PHea  # concentration in venous blood
    
    CLun <- ALun / VLun  # concentration in lung tissue (ug/L)
    CVLun <- CLun / PLun  # concentration in venous blood leaving lung
    
    CMus <- AMus / VMus  # concentration in muscle tissue (ug/L)
    CVMus <- CMus / PMus  # concentration in venous blood leaving muscle
    
    CSki <- ASki / VSki  # concentration in skin tissue (ug/L)
    CVSki <- CSki / PSki  # concentration in venous blood leaving skin
    
    CSpl <- ASpl / VSpl  # concentration in spleen tissue (ug/L)
    CVSpl <- CSpl / PSpl  # concentration in venous blood leaving spleen

    Cpan <- APan / VPan  # concentration in pancreas tissue (ug/L)
    CVPan <- Cpan / PPan  # concentration in venous blood leaving pancreas
    
    CGI <- AGI / VGI  # concentration in GI tract (ug/L)
    CVGI <- CGI / PGI  # concentration in venous blood leaving GI tract
    
    CKb <- AKb / VKb  # concentration in kidney blood (ug/L)
    CVK <- CKb  # concentration in venous blood leaving kidney (ug/L)
    CPTC <- APTC / VPTC  # concentration in PTC (ug/L)
    Cfil <- Afil / Vfil  # concentration in filtrate (ug/L)
    
    CL <- AL / VL  # concentration in the liver (ug/L)
    CLiver <- AL / ML  # concentration in the liver (ug/g)
    CVL <- CL / PL  # concentration in the venous blood leaving the liver (ug/L)
    
    CA_free <- Aplas_free / VPlas  # free concentration in plasma (ug/L)
    CA <- CA_free / Free  # concentration of total PFOS in plasma (ug/L)
    
    # Rest of Body (Tis)
    dAR <- QR * (CA - CVR) * Free

    # Adipose Tissue (Adi)
    dAAdi <- QAdi * (CA - CVAdi) * Free

    # Brain Tissue (Bra)
    dABra <- QBra * (CA - CVBra) * Free

    # Gonads Tissue (Gon)
    dAGon <- QGon * (CA - CVGon) * Free

    # Heart Tissue (Hea)
    dAHea <- QHea * (CA - CVHea) * Free

    # Lung Tissue (Lun)
    dALun <- QLun * (CA - CVLun) * Free
    
    # Muscle Tissue (Mus)
    dAMus <- QMus * (CA - CVMus) * Free

    # Skin Tissue (Ski)
    dASki <- QSki * (CA - CVSki) * Free

    # Spleen Tissue (Spl)
    dASpl <- QSpl * (CA - CVSpl) * Free

    # Pancreas Tissue (Pan)
    dAPan <- QPan * (CA - CVPan) * Free
    
    # Kidney
    # Kidney Blood (Kb)
    dAdif <- kdif * (CKb - CPTC)
    dA_baso <- (Vmax_baso * CKb) / (Km_baso + CKb)
    dAKb <- QK * (CA - CVK) * Free - CA * GFR * Free - dAdif - dA_baso
    dACl <- CA * GFR * Free

    # Proximal Tubule Cells (PTC)
    dAefflux <- kefflux * APTC
    dA_apical <- (Vmax_apical * Cfil) / (Km_apical + Cfil)
    dAPTC <- dAdif + dA_apical + dA_baso - dAefflux

    # Filtrate (Fil)
    dAfil <- CA * GFR * Free - dA_apical - Afil * kurine

    # Urinary elimination
    dAurine <- kurine * Afil

    # GI Tract (Absorption site of oral dose)
    # Stomach_lumen
    dALumen <- ingestion + Cwater * water_consumption - kabs * ALumen- kunabs * ALumen 


    # GI Tract - Stomach,Small Intestine,Large Intestine 
    dAGI <- kabs * ALumen + QGI*(CA - CVGI) * Free  
    dAabsLumen <- kabs * ALumen


    # Feces compartment
    dAfeces <-  kunabs * ALumen + kbile * AL 

    # Liver
    dAL <- QL * (CA - CVL) * Free - kbile * AL + QGI*(CVGI - CVL) * Free +
           QSpl*(CVSpl-CVL)*Free +QPan*(CVPan-CVL)*Free
    
    dAbile <- kbile * AL
    amount_per_gram_liver <- CLiver

    # Plasma compartment
    dAplas_free <- (QR * CVR * Free) + (QK * CVK * Free) + (QL * CVL * Free) +
    (QAdi*CVAdi*Free) + (QBra*CVBra*Free) + (QGon*CVGon*Free) + (QHea*CVHea*Free) +
    (QLun*CVLun*Free) +(QMus*CVMus*Free) +(QSki*CVSki*Free)+(QSpl*CVL*Free)+(QPan*CVL*Free)+
    (QGI*CVL*Free)- (QC * CA * Free) + dAefflux

    dCwater <- 0
    dingestion <- 0

    # Mass Balance Check
    Atissue <- Aplas_free + AR + AAdi + ABra + AGon + AHea+ ALun + AMus + ASki + ASpl+ APan + AKb + Afil + ALumen+ AGI + APTC + AL 
    Aloss <- Aurine + Afeces
    Atotal <- Atissue + Aloss

    list(
      c(
        "dAR" = dAR, "dAAdi" = dAAdi, "dABra" = dABra, "dAGon" = dAGon,
        "dAHea" = dAHea, "dALun" = dALun, "dAMus" = dAMus, "dASki" = dASki,
        "dASpl" = dASpl, "dAPan"=dAPan,
        "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
        "dACl" = dACl, "dAefflux" = dAefflux,
        "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
        "dAurine" = dAurine, "dALumen" = dALumen, "dAGI" = dAGI,
        "dAabsLumen" = dAabsLumen, "dAfeces" = dAfeces,
        "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
        "dCwater" = dCwater, "dingestion" = dingestion
      ),
      "amount_per_gram_liver" = amount_per_gram_liver,
      "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal,
      "CR" = CR, "CVR" = CVR, "CAdi"=CAdi, "CVAdi"=CVAdi,
      "CBra"=CBra, "CVBra"=CVBra, "Gon"=CGon, "CVGon"=CVGon,
      "CHea"=CHea, "CVHea"=CVHea, "CLun"=CLun, "CVLun"=CVLun,
      "CMus"=CMus, "CVMus"=CVMus, "CSki"=CSki, "CVSki"=CVSki,
      "CSpl"=CSpl, "CVSpl"=CVSpl, "Cpan"=Cpan, "CVPan"=CVPan,
      "CKb" = CKb, "CGI" = CGI,
      "CVGI" = CVGI, "CVK" = CVK, "CPTC" = CPTC,
      "Cfil" = Cfil, "CL" = CL, "CVL" = CVL,
      "CA_free" = CA_free, "CA" = CA
    )
  })
}


# ================================================================================
# 6. SIMULATION SETUP
# ================================================================================
# --- Load initial parameters and partition coefficients ---

# --- Dosing and Administration Settings ---

admin_type <- "bolus"  # administration type: "iv", "oral", or "bolus"
admin_dose_bolus <- c(variables["admin_dose_bolus","PFOS"])  # administered dose through bolus (ug)
admin_time_bolus <- c(0)  # time when bolus doses are administered (days)
admin_dose_iv <- 0  # administered dose through IV (ug)
admin_time_iv <- 0  # time when IV doses are administered (days)
Cwater <- 0.00  # concentration in water (ug/L)
Cwater_time <- 0  # time of water concentration change
ingestion <- 0  # ingestion rate (ug/day)
ingestion_time <- c(0)  # time of ingestion rate change



RAFbaso	= variables["RAFbaso","PFOS"]
RAFapi = variables["RAFapi","PFOS"]
PR = 0.01 # rest of body:blood#; fitted to model
PAdi = PC["Adipose","PFOS"] # adipose tissue:blood;
PBra = PC["Brain","PFOS"]  # brain:blood;
PGon = PC["Gonads","PFOS"]  # gonads:blood;
PGI = PC["Gut","PFOS"]  # GI tract:blood;
PHea = PC["Heart","PFOS"]  # heart:blood;
PL = PC["Liver","PFOS"] # liver:blood;
PLun = PC["Lung","PFOS"]  # lung:blood;
PMus = PC["Muscle","PFOS"]  # muscle:blood;
PSki = PC["Skin","PFOS"] # skin:blood;
PSpl = PC["Spleen","PFOS"]  # spleen:blood;
PPan = (PGI+PSpl)/2  # pancreas:blood; estimated as average of GI tract and spleen
MW = variables["MW", "PFOS"]  # PFOS molecular mass (g/mol)
Free = variables["Free", "PFOS"]  # free fraction in plasma (Smeltz 2023); fitted to model
kdif = 0.001*24
kabsc = 2.12*24
kunabsc = 7.06e-5 * 24
kbilec <- variables["kbilec", "PFOS"]   # biliary elimination rate (1/(day*BW^-0.25)); fitted to model 
kurinec <- 0.063 * 24  # urinary elimination rate (1/(day*BW^-0.25))
Vmax_baso_invitro = variables["Vmax_baso_invitro", "PFOS"] #Vmax of basolateral transporter (pmol/mg protein/min); averaged in vitro value of OAT1 and OAT3 from Nakagawa, 2007
Km_baso = variables["Km_baso", "PFOS"] #Km of basolateral transporter (ug/L) Average of OAT1 and OAT3 from Nakagawa et. al, 2007
Vmax_apical_invitro = variables["Vmax_apical_invitro", "PFOS"] #Vmax of apical transporter (pmol/mg protein/min); invitro value for OAT4 from Yang et al, 2010
Km_apical = variables["Km_apical", "PFOS"] #Km of apical transporter (ug/L), in vitro value for OAT4 and URAT1 from Yang et al, 2010.
keffluxc =variables["keffluxc", "PFOS"] # rate of efflux from PTC to blood (1/(day*BW^-0.25)) 


user_input<-list(
"Free"  = Free,
"RAFbaso"	= RAFbaso,
"RAFapi" = RAFapi,
"PR" = PR, # rest of body:plasma#; fitted to model
"PAdi" = PAdi, # adipose tissue:plasma;
"PBra" = PBra,  # brain:plasma;
"PGon" = PGon,  # gonads:plasma;
"PGI" = PGI,  # GI tract:plasma;
"PHea" = PHea,  # heart:plasma;
"PL" = PL,  # liver:plasma;
"PLun" = PLun,  # lung:plasma;
"PMus" = PMus,  # muscle:plasma;
"PSki" = PSki,  # skin:plasma;
"PSpl" = PSpl,  # spleen:plasma;
"PPan" = PPan,  # pancreas:plasma; estimated as average of GI tract and spleen
"MW" = MW,  # PFOS molecular mass (g/mol)
"keffluxc" = keffluxc,  
"kdif" = kdif,
"kabsc" = kabsc,
"kunabsc" = kunabsc, 
"kbilec" =  kbilec,
"kurinec" =kurinec,
"Vmax_baso_invitro" = Vmax_baso_invitro, #Vmax of basolateral transporter (pmol/mg protein/min); averaged in vitro value of OAT1 and OAT3 from Nakagawa, 2007
"Km_baso" = Km_baso, #Km of basolateral transporter (ug/L) Average of OAT1 and OAT3 from Nakagawa et. al, 2007
"Vmax_apical_invitro" = Vmax_apical_invitro, #Vmax of apical transporter (pmol/mg protein/min); invitro value for OAT4 from Yang et al, 2010
"Km_apical" = Km_apical, #Km of apical transporter (ug/L), in vitro value for OAT4 and URAT1 from Yang et al, 2010.
"admin_type" = admin_type,
"admin_dose_bolus" = admin_dose_bolus,
"admin_time_bolus" = admin_time_bolus,
"admin_dose_iv" = admin_dose_iv,
"admin_time_iv" = admin_time_iv,
"Cwater" = Cwater,
"Cwater_time" = Cwater_time,
"ingestion" = ingestion,  # ingestion rate (ug/day)
"ingestion_time" = ingestion_time
)
#=================================

obj.func<-function(x,dp){
 
if (is.na(x)) {
    # Baseline: use full user_input
    params <- create.params(user_input)
  } else {
    # Perturbed: modify a copy of user_input
    pert_input <- user_input
    pert_input[[x]] <- as.numeric(user_input[[x]]) * (1 + dp)
    params <- create.params(pert_input)
  }
  
  inits <- create.inits(params)
  events <- create.events(params)
  
  sample_time <- seq(0, 6, 0.01)
  solution <- as.data.frame(
    ode(
      times = sample_time,
      func = ode.func,
      y = inits,
      parms = params,
      events = events,
      method = "bdf",
      rtol = 1e-05,
      atol = 1e-05
    )
  )
  
  result <- c(
    Cmax = max(solution[, "CA"]),
    AUC = AUC(sample_time, solution[, "CA"]),
    C_6h = solution[which.min(abs(sample_time - 0.25)), "CA"],
    C_6d = solution[which.min(abs(sample_time - 6)), "CA"],
    Feces_Excreted = tail(solution[, "Afeces"], 1),   # Use final value, not sum!
    Urine_Excreted = tail(solution[, "Aurine"], 1)    # Same here
  )
  
  return(result)
}

#=================================
#8. Sensitivity Analysis Parameter
#=================================
thetas<-list(
"RAFbaso"	,
"RAFapi" ,
"PR" , 
"PAdi" , 
"PBra" ,  
"PGon" ,  
"PGI" , 
"PHea", 
"PL" , 
"PLun" ,  
"PMus" ,  
"PSki" ,  
"PSpl" ,  
"PPan" ,  
"kdif" ,
"kabsc",
"kunabsc",
"keffluxc",
"kbilec",
"kurinec",
"Vmax_baso_invitro", 
"Km_baso", 
"Vmax_apical_invitro",
"Km_apical",
"Free"  )

# =======================
# 9. Sensitivity Analysis (Positive & Negative Perturbations)
# =======================

dp <- 0.5  # Use positive magnitude; we'll apply +dp and -dp

# Get baseline result
baseline <- obj.func(NA, dp = 0)

# Initialize long-format data frame for tidy plotting
sensitivity_df <- tibble(
  Parameter = character(),
  Direction = factor(levels = c("Negative", "Positive")),
  Output = character(),
  RelativeChange = numeric()
)

# Loop over each parameter
for (par in thetas) {
  cat("Processing:", par, "\n")
  
  # Positive perturbation (+dp)
  res_pos <- obj.func(par, dp = dp)
  rel_pos <- ((res_pos - baseline) / baseline)/ dp
  
  # Negative perturbation (-dp)
  res_neg <- obj.func(par, dp = -dp)
  rel_neg <- ((res_neg - baseline) / baseline)/(dp)
  
  # Append to data frame (one row per output metric)
  n_outputs <- length(baseline)
  output_names <- names(baseline)
  
  sensitivity_df <- bind_rows(sensitivity_df,
    tibble(
      Parameter = rep(par, n_outputs),
      Direction = "Positive",
      Output = output_names,
      RelativeChange = rel_pos
    ),
    tibble(
      Parameter = rep(par, n_outputs),
      Direction = "Negative",
      Output = output_names,
      RelativeChange = rel_neg
    )
  )
}

# Convert Parameter to factor for consistent ordering
sensitivity_df$Parameter <- fct_inorder(sensitivity_df$Parameter)

# Generate one plot per output variable
plot_list <- list()
output_vars <- unique(sensitivity_df$Output)

for (j in seq_along(output_vars)) {
  out_var <- output_vars[j]
  
  plot_data <- sensitivity_df %>%
    filter(Output == out_var) %>%
    mutate(Parameter = fct_reorder(Parameter, RelativeChange, .desc = FALSE))
  
  p <- ggplot(plot_data, aes(x = Parameter, y = RelativeChange, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Negative" = "tomato", "Positive" = "steelblue")) +
    scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))  # 5% below, 15% above (adjust as needed)
  ) +

    labs(
      title = paste("PFOS Sensitivity Analysis:", out_var),
      x = "Parameters",
      y = "Relative Change in Output",
      fill = "Perturbation"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      legend.position = "top"
    ) +
    geom_text(
      aes(label = format(RelativeChange, scientific = TRUE, digits = 2)),
      position = position_dodge(width = 0.8),
      hjust = -0.2,
      size = 3
    )
  
  plot_list[[out_var]] <- p
}

# =========================
# Save Plots to Files
# =========================a
output_dir <- "PFOS_posterior_sensitivity_plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (out_var in names(plot_list)) {
  png_file <- file.path(output_dir, paste0("PFOS_", out_var, "_dual_sensitivity.png"))
  ggsave(
    filename = png_file,
    plot = plot_list[[out_var]],
    width = 10,
    height = 8,
    dpi = 300
  )
  cat("Saved:", png_file, "\n")
}