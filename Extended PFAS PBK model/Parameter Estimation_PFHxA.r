
  library(tidyverse)
  library(nloptr)
  library(deSolve)

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
    MW <- variables["MW","PFHxA"]  # PFHxA molecular mass (g/mol)
    Free <- variables["Free","PFHxA"]  # free fraction in plasma (Smeltz 2023); fitted to model

    # --- Kidney Transport Parameters ---
    Vmax_baso_invitro <- variables["Vmax_baso_invitro","PFHxA"]  # Vmax of basolateral transporter (pmol/mg protein/min)
    Km_baso <- variables["Km_baso","PFHxA"] * variables["MW","PFHxA"]  # Km of basolateral transporter (ug/L)
    Vmax_apical_invitro <- variables["Vmax_apical_invitro", "PFHxA"]  # Vmax of apical transporter (pmol/mg protein/min)
    Km_apical <- variables["Km_apical", "PFHxA"] * variables["MW", "PFHxA"]  # Km of apical transporter (ug/L)
    protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
    GFRC <- 24.19 * 24  # glomerular filtration rate (L/day/kg kidney); Corley 2005

    # --- Partition Coefficients (from Allendorf 2021) ---
    PAdi <- (1-Htc)*PC["Adipose","PFHxA"] # adipose tissue:plasma;
    PBra <- (1-Htc)*PC["Brain","PFHxA"] # brain:plasma;
    PGon <- (1-Htc)*PC["Gonads","PFHxA"] # gonads:plasma;
    PGI <-  (1-Htc)*PC["Gut","PFHxA"] # GI tract:plasma;
    PHea <- (1-Htc)*PC["Heart","PFHxA"] # heart:plasma;
    PL <-   (1-Htc)*PC["Liver","PFHxA"] # liver:plasma;
    PLun <- (1-Htc)*PC["Lung","PFHxA"]  # lung:plasma;
    PMus <- (1-Htc)*PC["Muscle","PFHxA"] # muscle:plasma;
    PSki <- (1-Htc)*PC["Skin","PFHxA"]   # skin:plasma;
    PSpl <- (1-Htc)*PC["Spleen","PFHxA"]  # spleen:plasma;
    PPan <- (PGI+PSpl)/2 # pancreas:plasma; estimated as average of GI tract and spleen
    PR <- 0.1 # rest of body:blood

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

    Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (variables["MW", "PFHxA"] / 1e12) * 1e6) * 24
    Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (variables["MW", "PFHxA"] / 1e12) * 1e6) * 24
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
  
  
  rmsd <- function(predictions, observations, weight=2,threshold=10, times=NULL){
    y_obs <- unlist(observations)
    y_pred <- unlist(predictions)
    # Total number of observations
    N<- length(y_obs)
    summation <- 0
    
    if(!is.null(times)) {
      for ( i in 1:N){
        if(times[i]>threshold) {
          summation <- summation + (weight*(y_obs[i]-y_pred[i]))^2
        } else {
          summation <- summation + (y_obs[i]-y_pred[i])^2
        }
      }
    }else{
      for ( i in 1:N){
        summation <- summation + (y_obs[i]-y_pred[i])^2
      }
    }

    rmsd <- sqrt(summation/N)

    return(rmsd)
  }
  
  
  AAFE <- function(predictions, observations,weight=2,threshold=10, times=NULL){
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
   
       if(!is.null(times)) {
      times <- times[valid_indices]
       for ( i in 1:N){
        if(times[i]>threshold) {
          log_ratio[i] <- weight*abs(log((y_pred[i]/y_obs[i]), base = 10)) 
        } else {
      log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
    }
    }
    }else{
      for ( i in 1:N){
        log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
      }
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
      CA <- CA_free / Free  # concentration of total PFHxA in plasma (ug/L)
      
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
  admin_dose_bolus <- c(variables["admin_dose_bolus","PFHxA"])  # administered dose through bolus (ug)
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
  RAFbaso <-variables["RAFbaso", "PFHxA"]  # relative activity factor, basolateral transporters
  RAFapi <- variables["RAFapi", "PFHxA"] # relative activity factor, apical transporter
  keffluxc <-variables["keffluxc", "PFHxA"] # rate of efflux from PTC to blood (1/(day*BW^-0.25))
  kbilec <-variables["kbilec", "PFHxA"]  # biliary elimination rate (1/(day*BW^-0.25))


  #=====================
  #8. Objective function 
  #=====================
  obj.func<-function(x, weight, threshold, metric){
    user_input <- list( "admin_type" = admin_type,
                        "admin_dose_bolus"=admin_dose_bolus,
                        "admin_time_bolus"=admin_time_bolus,
                        "admin_dose_iv" = admin_dose_iv, 
                        "admin_time_iv" = admin_time_iv,
                        "Cwater" = Cwater, 
                        "Cwater_time" = Cwater_time, "ingestion" = ingestion,
                        "ingestion_time" = ingestion_time,
                        "RAFbaso" = x[1],
                        "RAFapi" = x[2],
                        "kbilec"= x[3],
                        "keffluxc"=x[4]
                        )

    params <- create.params(user_input,variables,PC)
    inits <- create.inits(user_input)
    events <- create.events(user_input)


    sample_time=sort(unique(c(seq(0, 100, 1), plasma_exp$time#, urine_exp$time,feces_exp$time
    )))
    
    solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                        events = events,
                                        method="lsoda",rtol = 1e-05, atol = 1e-05))
    
    #Plasma Score 
    plasma_idx <- find_nearest(solution$time, plasma_exp$time)
    preds_plasma <- solution[plasma_idx, "CA"]
    time_plasma <- solution[plasma_idx, "time"]
    
    #Urine score
    urine_idx <- find_nearest(solution$time, urine_exp$time)
    preds_urine <- solution[urine_idx, "Aurine"]

    
    # #Feces score
    feces_idx <- find_nearest(solution$time, feces_exp$time)
    preds_feces <- solution[feces_idx, "Afeces"]

    
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
    
    if (metric == "AAFE"){
    # Calculate scores AAFE
    score_plasma <- AAFE(preds_plasma, plasma_exp$PFHxA
    ,weight = weight,threshold = threshold ,times = time_plasma)
    }else if(metric == "rmsd"){
      score_plasma <- rmsd(preds_plasma, plasma_exp$PFHxA
                           ,weight = weight,threshold = threshold ,times = time_plasma)
    }
    score_urine <- AAFE(preds_urine, urine_exp$cumulative_mass)
     score_feces <- AAFE(preds_feces, feces_exp$cumulative_mass)
    
    # #Calculate scores RMSD
    # score_plasma <- rmsd(preds_plasma, plasma_exp$PFHxA)
    # score_urine <- rmsd(preds_urine, urine_exp$cumulative_mass)
    # score_feces <- rmsd(preds_feces, feces_exp$cumulative_mass)

    # Check if any score is NA
    if(is.na(score_plasma) #|| is.na(score_urine) || is.na(score_feces)
    ) {
      cat("NA score for parameters:", x, "\n")
      cat("Scores:", score_plasma
      , "\n")
      return(1e6)
    }
  
    #Return mean of scores - FIXED: use c() to create a vector

    return(mean(c(score_plasma
    ,score_urine
     ,score_feces
    )))

    # return((score_plasma*length(plasma_exp$PFHxA) + score_urine*length(exp_data_urine) +
    #  score_feces*length(exp_data_feces))
    # /(length(plasma_exp$PFHxA) + length(exp_data_urine) + length(exp_data_feces)))
  }
  #---------------------------------#
  # Set up the Optimization process #
  #---------------------------------#



  #Continuous amount in feces is at the 6 hour time point
  # ==========================================
  # 1. FECES EXTENSION (DENSE DT)
  # ==========================================
   exp_data_feces <- read.csv("exp_data_feces.csv")

  # # Original data up to 6 hours
  feces_exp <- cumulative_exp_data(exp_data_feces, "time", "PFHxA", "feces.weight") %>% 
    mutate(cumulative_mass = cumulative_mass / 1000) %>% 
    filter(time <= 6) %>% 
    arrange(time)

  #  # Extension data beyond 6 hours
  # feces_ext <- exp_data_feces %>% 
  #   filter(time > 6, PFHxA>0) %>% 
  #   select(time, PFHxA) %>% 
  #   mutate(PFHxA = PFHxA / 1000) %>%  # Convert from ng/g to ug/g
  #   arrange(time)

  # fecal_excretion_rate <- 128 # g/day
  # step_feces <- 5       

  # # Create dense time grid
  # t_dense_feces <- seq(min(feces_ext$time), max(feces_ext$time), by = step_feces)


  # # Interpolate PFHxA onto dense grid & compute cumulative mass
  # feces_cum_dense <- approx(feces_ext$time, feces_ext$PFHxA, xout = t_dense_feces)$y %>% 
  #   { cumsum(fecal_excretion_rate * . * step_feces) + tail(feces_exp$cumulative_mass, 1) }

  # df_feces_dense <- data.frame(time = t_dense_feces, cumulative_mass = feces_cum_dense)

  # # Bind & sort
  # feces_exp <- bind_rows(feces_exp, df_feces_dense) %>% arrange(time)


  # # ==========================================
  # # 2. URINE EXTENSION (DENSE DT)
  # # ==========================================
  exp_data_urine <- read.csv("exp_data_urine.csv")

  # # Fixed: merged duplicate mutate lines & applied time/24 consistently
  urine_exp <- cumulative_exp_data(exp_data_urine, "time", "PFHxA", "urine.volume") %>% 
    mutate(cumulative_mass = cumulative_mass , time = time / 24) %>% 
    filter(time <= 6) %>% 
    arrange(time)

  # urine_ext <- exp_data_urine%>% 
  #   mutate(time = time / 24)%>%
  #   filter(time > 6, PFHxA>0) %>% 
  #   select(time, PFHxA)
  

  # urine_excretion_rate <- 1.4 # L/day
  # step_urine <- 5       #DENSITY STEP

  # # Create dense time grid
  # t_dense_urine <- seq(min(urine_ext$time), max(urine_ext$time), by = step_urine)

  # # Interpolate & compute cumulative mass
  # urine_cum_dense <- approx(urine_ext$time, urine_ext$PFHxA, xout = t_dense_urine)$y %>% 
  #   { cumsum(urine_excretion_rate * . * step_urine) + tail(urine_exp$cumulative_mass, 1) }

  # # Map back to original time points
  # df_urine_dense <- data.frame(time = t_dense_urine, cumulative_mass = urine_cum_dense)

  # # #Bind & sort
  # urine_exp <- rbind(urine_exp, df_urine_dense) 

  #Creates a teble with time and plasma concentration
  plasma_exp <- read.csv("exp_data_plasma.csv")%>% select("time","PFHxA") #%>%filter(time<=28)#%>%
  #slice(c(seq(1,12,1),seq(13,55,4),seq(56,n(),1))) 
  # # 1. Extend time points

  # early_time <- seq(0.1, 12, by = 0.05)
  # #late_time<- seq(10, 2500, by = 1)

  # total_time <- sort(c( early_time#, late_time
  # ))


  # C_el= 1.365*exp(-0.246*total_time)+0.607*exp(-7.749*total_time)

  # plasma_art<- data.frame(time = c(total_time), PFHxA = C_el)



  x0 <- c("RAFbaso" = RAFbaso,
          "RAFapi" = RAFapi,
          "kbilec"=kbilec,
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

  lower_bounds <- c("RAFbaso" =1e-05,
                    "RAFapi"= 1e-07,
                    "kbilec"= 1e-05,
                    "keffluxc"= 1e-05
                                        )

  upper_bounds <- c("RAFbaso"=1e4,
                    "RAFapi"= 1e4,
                    "kbilec"= 1e4,
                    "keffluxc"= 1e2
                    )    
  #Call the optimization algorithm and provide him with the input data
  optimization <- nloptr::nloptr(x0 = x0,
                                eval_f = obj.func,
                                lb	= lower_bounds ,
                                ub = upper_bounds,
                                opts = opts,
                                weight = 5, 
                                threshold = 0.005,
                                metric = "AAFE")



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
                      "kbilec"= x_opt[3],
                      "keffluxc"=x_opt[4]
                      
                      )
  
  params <- create.params(user_input,variables,PC)
  inits <- create.inits(params)
  events <- create.events(params)

  # sample_time: a vector of time points to solve the ODEs
  sample_time=unique(sort(c(seq(0,450,1),plasma_exp$time)))


  # ode(): The solver of the ODEs
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                      events = events,
                                      method="lsoda",rtol = 1e-05, atol = 1e-05))

  compartments <- c('CA','Aurine','Afeces')
  color_codes <- scales::hue_pal()(length(compartments))

  plot1 <- ggplot()+
  geom_line(data = solution, aes(x = time, y = Aurine, color='Aurine'), size=1.3)+
  geom_line(data = solution, aes(x = time, y = Afeces, color='Afeces'), size=1.3)+
  geom_point(data = urine_exp, aes(x = time, y = cumulative_mass, color='Aurine'), size=5)+
  geom_point(data = feces_exp, aes(x = time, y = cumulative_mass, color='Afeces'), size=5)+
  
  labs(title = 'Predicted vs Observed Values',
        y = 'Mass (ug)' , x = "Time (days)")+
  xlim(0,10)+
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
    geom_line(data = solution, aes(x = time/365, y = CA, color='CA'), size=1.3)+
    geom_point(data = plasma_exp, aes(x = time/365, y = PFHxA, color='CA'), size=5)+
    labs(title = "Predicted vs Observed Values",
        y = "Mass (ug/L)" , x = "Time (years)")+
    theme(plot.title = element_text(hjust = 0.5,size=30),
          axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.y=element_text(size=22),
          axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
          axis.text.x=element_text(size=22),
          legend.title=element_text(hjust = 0.5,size=25),
          legend.text=element_text(size=22),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
          xlim(0,12/365)+
    scale_color_manual("Compartments", values=color_codes)+
    theme(legend.key.size = unit(1.5, 'cm'),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          axis.text = element_text(size = 14))
  print(plot2)

  x_opt

# ===================================================================
# NON-COMPARTMENTAL ANALYSIS (NCA) using PKNCA
# ===================================================================
library(PKNCA)

# --- Step 1: Prepare concentration data ---
conc_exp <- data.frame(Subject = "Experimental", Time = plasma_exp$time, Concentration = plasma_exp$PFHxA)
conc_pred <- data.frame(Subject = "Predicted", Time = solution$time, Concentration = solution$CA)
conc_df <- rbind(conc_exp, conc_pred)

# FIX: Add explicit time=0, concentration=0 rows for both groups
# This allows PKNCA to calculate AUC starting from t=0
zero_rows <- data.frame(
  Subject = c("Experimental", "Predicted"),
  Time = 0,
  Concentration = 0
)
conc_df <- rbind(conc_df, zero_rows)

# Remove exact duplicates (in case solution already contains Time=0) & sort
conc_df <- conc_df[!duplicated(conc_df[, c("Subject", "Time")]), ]
conc_df <- conc_df[order(conc_df$Subject, conc_df$Time), ]

# --- Step 2: Prepare dose data ---
dose_df <- data.frame(
  Subject = c("Experimental", "Predicted"),
  Time = 0,
  Dose = admin_dose_bolus,
  Route = "extravascular"  # Change to "intravascular" if using IV
)

# --- Step 3: Create PKNCA objects ---
my_conc <- PKNCAconc(conc_df, Concentration ~ Time | Subject)
my_dose <- PKNCAdose(dose_df, Dose ~ Time | Subject)

# --- Step 4: Define NCA intervals ---
max_time <- max(conc_df$Time, na.rm = TRUE)

my_data <- PKNCAdata(
  my_conc,
  my_dose,
  intervals = data.frame(
    start = 0,  # Now valid because we added Time=0, Conc=0
    end = max_time,
    cmax = TRUE,          # Maximum concentration
    tmax = TRUE,          # Time of maximum concentration
    auclast = TRUE,       # AUC from 0 to last observed time
    aucinf.obs = TRUE,    # AUC from 0 to infinity (observed terminal slope)
    aucinf.pred = TRUE,   # AUC from 0 to infinity (predicted terminal slope)
    half.life = TRUE,     # Terminal elimination half-life
    lambda.z = TRUE,      # Terminal elimination rate constant
    r.squared = TRUE,     # R² of terminal slope fit
    stringsAsFactors = FALSE
  )
)

# --- Step 5: Perform NCA ---
my_nca <- pk.nca(my_data)

# --- Step 6: View & Export Results ---
cat("\n========================================\n")
cat("       NCA SUMMARY RESULTS\n")
cat("========================================\n")
print(summary(my_nca))

nca_results_df <- as.data.frame(my_nca)
View(nca_results_df)


# # Optional: Save to CSV
 write.csv(nca_results_df, "NCA_results_PFHxA.csv", row.names = FALSE)
 write.csv(c(params,x_opt), "Optimized_Parameters_PFHxA.csv", row.names = TRUE)
