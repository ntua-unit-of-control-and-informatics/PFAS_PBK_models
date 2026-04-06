  library(deSolve)
  library(tidyverse)

  #=========================
  # 1. Parameters of the model
  #=========================
  # --- Load initial parameters and partition coefficients ---
  variables<-read.csv("estimated_parameters.csv",row.names="Parameters")
  PC<-read.csv("PCs.csv",row.names="Organs")

  create.params <- function(variables,PC,BW) {
    # === USER-ADJUSTABLE PARAMETERS ===
    # Change these default values as needed
    
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
    VGonC <- 0.0005 # fraction vol. of gonads (L/kg BW); 
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
    MW <- variables["MW", "PFOS"]  # PFOSmolecular mass (g/mol)
    Free <- variables["Free", "PFOS"]  # free fraction in plasma (Smeltz 2023); fitted to model

    # --- Kidney Transport Parameters ---
    Vmax_baso_invitro <- variables["Vmax_baso_invitro", "PFOS"]  # Vmax of basolateral transporter (pmol/mg protein/min)
    Km_baso <- variables["Km_baso", "PFOS"] * variables["MW", "PFOS"]  # Km of basolateral transporter (ug/L)
    Vmax_apical_invitro <- variables["Vmax_apical_invitro", "PFOS"]  # Vmax of apical transporter (pmol/mg protein/min)
    Km_apical <- variables["Km_apical", "PFOS"] * variables["MW", "PFOS"]  # Km of apical transporter (ug/L)
    RAFbaso <- variables["RAFbaso", "PFOS"]  # relative activity factor, basolateral transporters (male)
    RAFapi <- variables["RAFapi", "PFOS"] # relative activity factor, apical transporters (male); fitted to model
    protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
    GFRC <- 24.19 * 24  # glomerular filtration rate (L/day/kg kidney); Corley 2005

    # --- Partition Coefficients (from Allendorf 2021) ---
    PR <- 0.01 # rest of body:blood#; fitted to model
    PAdi <- PC["Adipose","PFOS"] # adipose tissue:blood;
    PBra <- PC["Brain","PFOS"]  # brain:blood;
    PGon <- PC["Gonads","PFOS"]  # gonads:blood;
    PGI <- PC["Gut","PFOS"]# GI tract:blood;
    PHea <- PC["Heart","PFOS"]  # heart:blood;
    PL<- PC["Liver","PFOS"]  # liver:blood;
    PLun <- PC["Lung","PFOS"]  # lung:blood;
    PMus <- PC["Muscle","PFOS"] # muscle:blood;
    PSki <- PC["Skin","PFOS"]  # skin:blood;
    PSpl <- PC["Spleen","PFOS"]  # spleen:blood;
    PPan <- (PGI+PSpl)/2  # pancreas:blood; estimated as average of GI tract and spleen

    # --- Rate Constants ---
    kdif <- 0.001 * 24  # diffusion rate from proximal tubule cells (L/day)
    kabsc <-  2.12 * 24  # rate of absorption from small intestine (1/(day*BW^-0.25))
    kunabsc <- 7.06e-5 * 24  # rate of unabsorbed dose to feces (1/(day*BW^-0.25)); fitted to model 
    keffluxc <-variables["keffluxc", "PFOS"]  # rate of efflux from PTC to blood (1/(day*BW^-0.25))
    kbilec <- variables["kbilec", "PFOS"]  # biliary elimination rate (1/(day*BW^-0.25)); fitted to model 
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

    Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
    Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
    Vmax_baso <- Vmax_basoC * BW^0.75  # (ug/day)
    Vmax_apical <- Vmax_apicalC * BW^0.75  # (ug/day)
    kbile <- kbilec * BW^(-0.25)  # biliary elimination; liver to feces storage (/day)
    kurine <- kurinec * BW^(-0.25)  # urinary elimination, from filtrate (/day)
    kefflux <- keffluxc * BW^(-0.25)  # efflux clearance rate, from PTC to blood (/day)
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
  create.events <- function(admin_type, admin_dose_bolus, admin_time_bolus,
                          admin_dose_iv, admin_time_iv,
                          Cwater, Cwater_time, ingestion, ingestion_time) {

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
      print("No complete cases found after removing NA values")
      result_df <- data.frame(
        time = numeric(0),
        cumulative_mass = numeric(0)
      )
    }

    multiplied_concentration <- clean_df[[concentration_col]] * clean_df[[multiply_col]]

    result_df <- data.frame(
      time = clean_df[[time_col]],
      cumulative_mass = cumsum(multiplied_concentration)
    )

    return(result_df)
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
      
      CA_free <- tids / VPlas  # free concentration in plasma (ug/L)
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
  BW <- 82 # body weight (kg); Abraham's individual value


  # --- Simulation Time Settings ---
  simulation_time <- 450  # total simulation time (days)
  time_step <- 0.1  # time step for ODE solver (days)

  cat("Setting up simulation...\n")

  # Generate parameters
  params<-create.params(variables,PC,BW)
  inits <- create.inits(params)
  events <- create.events(admin_type, admin_dose_bolus, admin_time_bolus,
                          admin_dose_iv, admin_time_iv,
                          Cwater, Cwater_time, ingestion, ingestion_time)

  cat("Solving ODEs...\n")

  # Solve ODEs
  sample_time <- seq(0, simulation_time, time_step)
  solution <- data.frame(deSolve::ode(
    times = sample_time, func = ode.func, y = inits, parms = params,
    events = events,
    method = "lsodes", rtol = 1e-05, atol = 1e-05
  ))

  cat("ODEs solved successfully!\n")

  # ================================================================================
  # 7. LOAD EXPERIMENTAL DATA
  # ================================================================================

  cat("Loading experimental data...\n")

  # Feces experimental data
  exp_data_feces <- read.csv("C:\\Users\\fotis\\Documents\\GitHub\\PFAS_PBK_models\\Extended PFAS PBK model\\exp_data_feces.csv")
  feces_exp <- cumulative_exp_data(exp_data_feces, "time", "PFOS", "feces.weight") %>%
    mutate(cumulative_mass = cumulative_mass / 1000) %>%
    filter(time <= 6)

  # Urine experimental data
  exp_data_urine <- read.csv("C:\\Users\\fotis\\Documents\\GitHub\\PFAS_PBK_models\\Extended PFAS PBK model\\exp_data_urine.csv")
  urine_exp <- cumulative_exp_data(exp_data_urine, "time", "PFOS", "urine.volume") %>%
    mutate(time = time / 24) %>%
    filter(time <= 6)

  # Plasma experimental data
  plasma_exp <- read.csv("C:\\Users\\fotis\\Documents\\GitHub\\PFAS_PBK_models\\Extended PFAS PBK model\\exp_data_plasma.csv") %>%
    select("time", "PFOS")

  cat("Experimental data loaded!\n")

  # ================================================================================
  # 8. PLOTTING
  # ================================================================================

  cat("Creating plots...\n")

  # --- Plot 1: Mass in Liver, Rest of Body, and Plasma ---
  plot_mass <- ggplot(solution) +
    geom_line(aes(x = time, y = AL, color = "Liver"), linewidth = 1.3) +
    geom_line(aes(x = time, y = AR, color = "Rest of Body"), linewidth = 1.3) +
    geom_line(aes(x = time, y = Aplas_free, color = "Plasma"), linewidth = 1.3) +
    geom_line(aes(x = time, y = ASpl, color = "Spleen"), linewidth = 1.3) +
    geom_line(aes(x = time, y = APan, color = "Pancreas"), linewidth = 1.3) +
    geom_line(aes(x = time, y = AGI, color = "GI Tract"), linewidth = 1.3) +
    geom_line(aes(x = time, y = AAdi, color = "Adipose"), linewidth = 1.3) +
    geom_line(aes(x = time, y = ABra, color = "Brain"), linewidth = 1.3) +
    geom_line(aes(x = time, y = AGon, color = "Gonads"), linewidth = 1.3) +
    geom_line(aes(x = time, y = AHea, color = "Heart"), linewidth = 1.3) +
    geom_line(aes(x = time, y = ALun, color = "Lung"), linewidth = 1.3) +
    geom_line(aes(x = time, y = AMus, color = "Muscle"), linewidth = 1.3) +
    geom_line(aes(x = time, y = ASki, color = "Skin"), linewidth = 1.3) +
    labs(
      title = "Mass in Different Compartments",
      x = "Time (days)",
      y = "Mass (ug)",
      color = "Compartment"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0)
    )

  print(plot_mass)

  # --- Plot 2: Concentration in Liver, Rest of Body, and Plasma ---
  plot_concentration <- ggplot(solution) +
    geom_line(aes(x = time, y = CL, color = "Liver"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CR, color = "Rest of Body"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CA, color = "Plasma"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CSpl, color = "Spleen"), linewidth = 1.3) +
    geom_line(aes(x = time, y = Cpan, color = "Pancreas"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CGI, color = "GI Tract"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CAdi, color = "Adipose"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CBra, color = "Brain"), linewidth = 1.3) +
    geom_line(aes(x = time, y = Gon, color = "Gonads"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CHea, color = "Heart"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CLun, color = "Lung"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CMus, color = "Muscle"), linewidth = 1.3) +
    geom_line(aes(x = time, y = CSki, color = "Skin"), linewidth = 1.3) +
    labs(
      title = "Concentration in Different Compartments",
      x = "Time (days)",
      y = "Concentration (ug/L)",
      color = "Compartment"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0)
    )

  print(plot_concentration)

  # --- Plot 3: Feces and Urine - Predictions vs Data ---
  plot_excretion <- ggplot() +
    geom_line(data = solution, aes(x = time, y = Aurine, color = "Urine (Model)"),
              linewidth = 1.3) +
    geom_line(data = solution, aes(x = time, y = Afeces, color = "Feces (Model)"),
              linewidth = 1.3) +
    geom_point(data = urine_exp, aes(x = time, y = cumulative_mass, color = "Urine (Data)"),
              size = 4) +
    geom_point(data = feces_exp, aes(x = time, y = cumulative_mass, color = "Feces (Data)"),
              size = 4) +
    labs(
      title = "Feces and Urine: Model Predictions vs Experimental Data",
      x = "Time (days)",
      y = "Cumulative Mass (ug)",
      color = "Legend"
    ) +
    xlim(0, 6) +
    #ylim(0, 0.05) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0)
    )

  print(plot_excretion)

  # --- Plot 4: Plasma - Predictions vs Data ---
  plot_plasma <- ggplot() +
    geom_line(data = solution, aes(x = time, y = CA, color = "Model Prediction"),
              linewidth = 1.3) +
    geom_point(data = plasma_exp, aes(x = time, y = PFOS, color = "Experimental Data"),
              size = 4) +
    labs(
      title = "Plasma Concentration: Model Predictions vs Experimental Data",
      x = "Time (days)",
      y = "Concentration (ug/L)",
      color = "Legend"
    ) +
    #xlim(0, 6) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0)
    )

  print(plot_plasma)

  cat("\nSimulation complete! All plots generated.\n")

  # ================================================================================
  # 9. SUMMARY OUTPUT
  # ================================================================================

  cat("\n=== SIMULATION SUMMARY ===\n")
  cat("Body Weight:", params$BW, "kg\n")
  cat("Administration Type:", params$admin_type, "\n")
  if (admin_type == "bolus") {
    cat("Bolus Dose:", admin_dose_bolus, "ug at time", admin_time_bolus, "days\n")
  } else if (admin_type == "iv") {
    cat("IV Dose:", admin_dose_iv, "ug at time",admin_time_iv, "days\n")
  }
  cat("Simulation Time:", simulation_time, "days\n")
  cat("Final Mass in Plasma:", round(tail(solution$Aplas_free, 1), 4), "ug\n")
  cat("Final Mass in Liver:", round(tail(solution$AL, 1), 4), "ug\n")
  cat("Final Mass in Rest of Body:", round(tail(solution$AR, 1), 4), "ug\n")
  cat("Total Mass Excreted in Urine:", round(tail(solution$Aurine, 1), 4), "ug\n")
  cat("Total Mass Excreted in Feces:", round(tail(solution$Afeces, 1), 4), "ug\n")
  cat("==========================\n")

  output_dir <- "PFOS_simulation_plots"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  png_file <- file.path(output_dir, paste0("PFOS_MASS_IN_BODY.png"))
    ggsave(
      filename = png_file,
      plot = plot_mass,
      width = 10,
      height = 8,
      dpi = 300
    )

  png_file <- file.path(output_dir, paste0("PFOS_CONCENTRATION_IN_BODY.png"))
    ggsave(
      filename = png_file,
      plot = plot_concentration,
      width = 10,
      height = 8,
      dpi = 300
    )   

  png_file <- file.path(output_dir, paste0("PFOS_URINE_FECES_COMPARISON.png"))
    ggsave(
      filename = png_file,
      plot = plot_excretion,
      width = 10,
      height = 8,
      dpi = 300
    )
  png_file <- file.path(output_dir, paste0("PFOS_PLASMA_COMPARISON.png"))
    ggsave( 
      filename = png_file,
      plot = plot_plasma,
      width = 10,
      height = 8,
      dpi = 300
    )

