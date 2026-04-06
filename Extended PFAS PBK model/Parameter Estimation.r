library(deSolve)
library(tidyverse)
library(nloptr)



# --- Load initial parameters, partition coefficients and experimental data ---
variables<-read.csv("initial_parameters.csv",row.names="Parameters")
PC<-read.csv("PCs.csv",row.names="Organs")




#=========================
# 1. Parameters of the model
#=========================

  
create.params <- function(user_input,variables,PC) {
     with( as.list(user_input),{
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

  Htc <- 0.467  # hematocrit

  # --- Tissue Volumes ---
  VplasC <- 0.0428  # fraction vol. of plasma (L/kg BW); Davies 1993
  VLC <- 0.026  # fraction vol. of liver (L/kg BW); Brown 1997
  VKC <- 0.004  # fraction vol. of kidney (L/kg BW); Brown 1997
  VAdiC <- 0.214  # fraction vol. of adipose tissue (L/kg BW); Brown 1997
  VBraC <- 0.02  # fraction vol. of brain (L/kg BW); Brown 1997
  VGonC <- 0.0005   # fraction vol. of gonads (L/kg BW); 
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
  MW <- variables["MW","PFOS"]  # PFOS molecular mass (g/mol)
  Free <- variables["Free","PFOS"]  # free fraction in plasma (Smeltz 2023); fitted to model

  # --- Kidney Transport Parameters ---
  Vmax_baso_invitro <- variables["Vmax_baso_invitro","PFOS"]  # Vmax of basolateral transporter (pmol/mg protein/min)
  Km_baso <- variables["Km_baso","PFOS"] * variables["MW","PFOS"]  # Km of basolateral transporter (ug/L)
  Vmax_apical_invitro <- variables["Vmax_apical_invitro", "PFOS"]  # Vmax of apical transporter (pmol/mg protein/min)
  Km_apical <- variables["Km_apical", "PFOS"] * variables["MW", "PFOS"]  # Km of apical transporter (ug/L)
  protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
  GFRC <- 24.19 * 24  # glomerular filtration rate (L/day/kg kidney); Corley 2005

  # --- Partition Coefficients (from Allendorf 2021) ---
  
  PAdi <- PC["Adipose","PFOS"] # adipose tissue:plasma;
  PBra <- PC["Brain","PFOS"]  # brain:plasma;
  PGon <- PC["Gonads","PFOS"]  # gonads:plasma;
  PGI <- PC["Gut","PFOS"]   # GI tract:plasma;
  PHea <- PC["Heart","PFOS"]  # heart:plasma;
  PL<- PC["Liver","PFOS"] # liver:plasma;
  PLun <- PC["Lung","PFOS"]  # lung:plasma;
  PMus <- PC["Muscle","PFOS"]  # muscle:plasma;
  PSki <- PC["Skin","PFOS"]  # skin:plasma;
  PSpl <- PC["Spleen","PFOS"]  # spleen:plasma;
  PPan <- (PGI+PSpl) # pancreas:plasma; estimated as average of GI tract and spleen
  PR <- 0.01 # rest of body:blood

  # --- Rate Constants ---
  kdif <- 0.001 * 24  # diffusion rate from proximal tubule cells (L/day)
  kabsc <-  2.12 *24  # rate of absorption from small intestine (1/(day*BW^-0.25))
  kunabsc <- 7.06e-5 * 24  # rate of unabsorbed dose to feces (1/(day*BW^-0.25)); fitted to model 
  kurinec <- 0.063 * 24  # urinary elimination rate (1/(day*BW^-0.25))
  

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

  Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (variables["MW", "PFOS"] / 1e12) * 1e6) * 24
  Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (variables["MW", "PFOS"] / 1e12) * 1e6) * 24
  Vmax_baso <- Vmax_basoC * BW^0.75  # (ug/day)
  Vmax_apical <- Vmax_apicalC * BW^0.75  # (ug/day)
  kbile <- kbilec * BW^(-0.25)  # biliary elimination; liver to feces storage (/day)
  kurine <- kurinec * BW^(-0.25)  # urinary elimination, from filtrate (/day)
  kefflux <- keffluxc * BW^(-0.25)  # efflux clearance rate, from PTC to blood (/day)
  GFR <- 163.65 # glomerular filtration rate, (L/day)Abraham's personal value

  # GI Tract Parameters
  kabs <- kabsc * BW^(-0.25)  # rate of absorption from small intestine (/day)
  kunabs <- kunabsc * BW^(-0.25)  # rate of unabsorbed dose to feces (/day)
   
  return(list(
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
    "PSpl" = PSpl, "PPan" = PPan,
    "water_consumption" = water_consumption))
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
create.events <- function(parameters) {
   with( as.list(parameters),{
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
 })
}


#==================
#4. Custom function 
#==================
mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  valid_indices <- which(y_obs > 0 & y_pred > 0 & 
                         !is.na(y_obs) & !is.na(y_pred))
  
  if(length(valid_indices) == 0) {
    warning("No valid observations for AAFE calculation (all zeros or NAs)")
    return(NA)
  }
  
  y_obs <- y_obs[valid_indices]
  y_pred <- y_pred[valid_indices]
  
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}
#Cumulative mass function aggragating experimental data
#Used for urine and feces samples
cumulative_exp_data<-function(df,time_col,concentration_col, multiply_col){
  
  working_df <- df[, c(time_col, concentration_col, multiply_col)]
  complete_cases <- complete.cases(working_df)
  clean_df <- working_df[complete_cases, ]
  
  if (nrow(clean_df) == 0) {
    stop("No complete cases found after removing NA values")
  }
  # Multiply and calculate cumulative sum
  multiplied_concentration <- clean_df[[concentration_col]] * clean_df[[multiply_col]]
  
  # Create new dataframe with time and cumulative concentration only
  result_df <- data.frame(
    time = clean_df[[time_col]],
    cumulative_mass = cumsum(multiplied_concentration)
  )
  
  return(result_df)
}

find_nearest <- function(times, target) {
    sapply(target, function(t) which.min(abs(times - t)))}


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
    (QGI*CVL*Free) - (QC * CA * Free) + dAefflux

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



#============================
#7. Parameters for estimation 
#============================
RAFbaso <- variables["RAFbaso", "PFOS"]  # relative activity factor, basolateral transporters
RAFapi <- variables["RAFapi", "PFOS"] # relative activity factor, apical transporter
keffluxc <- variables["keffluxc", "PFOS"] # rate of efflux from PTC to blood (1/(day*BW^-0.25))
kbilec <-variables["kbilec", "PFOS"]  # biliary elimination rate (1/(day*BW^-0.25))



#=====================
#8. Objective function 
#=====================
obj.func<-function(x){
  user_input <- list( "admin_type" = admin_type,
                      "admin_dose_bolus"=admin_dose_bolus,
                      "admin_time_bolus"=admin_time_bolus,
                      "admin_dose_iv" = admin_dose_iv, 
                      "admin_time_iv" = admin_time_iv,
                      "Cwater" = Cwater, 
                      "Cwater_time" = Cwater_time, "ingestion" = ingestion,
                      "ingestion_time" = ingestion_time,
                      "RAFbaso" =x[1],
                      "RAFapi" = x[2],
                      "kbilec"= kbilec, #x[3],
                      "keffluxc"=x[3]
                      )

  params <- create.params(user_input,variables,PC)
  inits <- create.inits(user_input)
  events <- create.events(user_input)


  sample_time=sort(unique(c(seq(0, 450, 0.1), plasma_exp$time#, urine_exp$time#,feces_exp$time
  )))
  
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                      events = events,
                                      method="lsodes",rtol = 1e-05, atol = 1e-05))
  
  #Plasma Score 
  plasma_idx <- find_nearest(solution$time, plasma_exp$time)
  preds_plasma <- solution[plasma_idx, "CA"]
  
  #Urine score
  #urine_idx <- find_nearest(solution$time, urine_exp$time)
  #preds_urine <- solution[urine_idx, "Aurine"]
  
  #Feces score
  #feces_idx <- find_nearest(solution$time, feces_exp$time)
  #preds_feces <- solution[feces_idx, "Afeces"]
  
  # Check for NA values in predictions
  if(any(is.na(preds_plasma)))
  {
    cat("NA values in predictions for parameters:", x, "\n")
    return(1e6)  # Return high error if NAs present
  }
  
  # Check for zero predictions (all zeros would give NA in AAFE)
  if(all(preds_plasma == 0) 
  ) {
    cat("Zero predictions for parameters:", x, "\n")
    return(1e6)
  }
  
  # Calculate scores
  score_plasma <- AAFE(preds_plasma, plasma_exp$PFOS)
  #score_urine <- AAFE(preds_urine, urine_exp$cumulative_mass)
  #score_feces <- AAFE(preds_feces, feces_exp$cumulative_mass)
  
  # Check if any score is NA
  if(is.na(score_plasma) #|| is.na(score_urine) || is.na(score_feces)
  ) {
    cat("NA score for parameters:", x, "\n")
    cat("Scores:", score_plasma
    #, score_urine
    #, score_feces
    , "\n")
    return(1e6)
  }
 
  # Return mean of scores - FIXED: use c() to create a vector
  return(mean(c(score_plasma
  #,score_urine
  #,score_feces
  )))
}
#---------------------------------#
# Set up the Optimization process #
#---------------------------------#



#Continuous amount in feces is at the 6 hour time point
#exp_data_feces <- read.csv("exp_data_feces.csv")
#feces_exp<-cumulative_exp_data(exp_data_feces,"time","PFOS","feces.weight") %>% mutate(cumulative_mass= cumulative_mass/ 1000) %>% filter(time <= 6)

#Creates a table with time and cummulative amount in urine

#Continuous amount in urine is at the 6 hour time point
#exp_data_urine <- read.csv("exp_data_urine.csv")
#urine_exp<-cumulative_exp_data(exp_data_urine,"time","PFOS","urine.volume")%>% mutate(time=time/24)%>% filter(time <= 6)


#Creates a teble with time and plasma concentration
plasma_exp <- read.csv("exp_data_plasma.csv")%>% select("time","PFOS")

x0 <- c( "RAFbaso" = RAFbaso ,
         "RAFapi" = RAFapi,
         #"kbilec"=kbilec,
         "keffluxc"=keffluxc
         )

N_iter <- 2000


# Extra options for the optimization algorithm
opts <- list( "algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NEWUOA",  #"NLOPT_LN_SBPLX" ,
              "xtol_rel" = 1e-05,
              "ftol_rel" = 1e-05,
              "ftol_abs" = 0.0,
              "xtol_abs" = 0.0 ,
              "maxeval" = N_iter,
              "print_level" = 1 )

lower_bounds <- c( "RAFbaso" = 0,
                   "RAFapi" =  0,
                   #"kbilec"=0,
                   "keffluxc"=0
                   )

upper_bounds <- c("RAFbaso"=10,
                  "RAFapi"=10,
                  #"kbilec"=10,
                  "keffluxc"=100
                  )    
#Call the optimization algorithm and provide him with the input data
optimization <- nloptr::nloptr(x0 = x0,
                              eval_f = obj.func,
                              lb	= lower_bounds ,
                              ub = upper_bounds,
                              opts = opts)



# The minimized value of the objective function
optimization$objective

# The values of the optimized params 
x_opt <- optimization$solution


#---------------------------------------------#
# Plot predictions over the experimental data #
#---------------------------------------------#

# Step 1: Solve the ODEs using the optimized values of parameters
user_input <- list( "admin_type" = admin_type,
                    "admin_dose_bolus" = admin_dose_bolus, 
                    "admin_time_bolus" = admin_time_bolus,
                    "Cwater" = Cwater, 
                    "Cwater_time" = Cwater_time, "ingestion" = ingestion,
                    "ingestion_time" = ingestion_time,
                    "RAFbaso" = x_opt[1],
                    "RAFapi" = x_opt[2],
                    "kbilec"=kbilec,#x_opt[3],
                    "keffluxc"=x_opt[3]
                    )

params <- create.params(user_input,variables,PC)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,450,0.1)
#sort(unique(c(0,plasma_exp$time,urine_exp$time,feces_exp$time)))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

compartments <- c('CA','Aurine','Afeces')
color_codes <- scales::hue_pal()(length(compartments))

plot1 <- ggplot()+
 geom_line(data = solution, aes(x = time, y = Aurine, color='Aurine'), size=1.3)+
 geom_line(data = solution, aes(x = time, y = Afeces, color='Afeces'), size=1.3)+
 geom_point(data = urine_exp, aes(x = time, y = cumulative_mass, color='Aurine'), size=5)+
 geom_point(data = feces_exp, aes(x = time, y = cumulative_mass, color='Afeces'), size=5)+

 labs(title = 'Predicted vs Observed Values',
      y = 'Mass (ug)' , x = "Time (days)")+
 xlim(0,6)+
 #ylim(0,0.025)+
 theme(plot.title = element_text(hjust = 0.5,size=30),
       axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
       axis.text.y=element_text(size=22),
       axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
       axis.text.x=element_text(size=22),
       legend.title=element_text(hjust = 0.5,size=25),
       legend.text=element_text(size=22),
       panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
 scale_color_manual("Compartments", values=color_codes)+
 theme(legend.key.size = unit(1.5, 'cm'),
       legend.title = element_text(size=14),
       legend.text = element_text(size=14),
       axis.text = element_text(size = 14))
print(plot1)

plot2 <- ggplot()+
  geom_line(data = solution, aes(x = time, y = CA, color='CA'), size=1.3)+
  geom_point(data = plasma_exp, aes(x = time, y = PFOS, color='CA'), size=5)+
  labs(title = "Predicted vs Observed Values",
       y = "Mass (ug/L)" , x = "Time (days)")+
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25),
        legend.text=element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
        xlim(0,6)+
  scale_color_manual("Compartments", values=color_codes)+
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))
print(plot2)


