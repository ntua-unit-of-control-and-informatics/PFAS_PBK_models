
  # ================================================================================
  # Generalized parameter estimation for all PFAS (Extended PFAS PBK model)
  #
  # Replaces the per-compound "Parameter Estimation_<PFAS>.r" scripts. The model,
  # ODEs and optimisation set-up are shared. One PFAS is estimated per run; the PFAS,
  # the plasma AAFE weight and the threshold are set in USER INPUTS below (or on the
  # command line). Bounds are common to all PFAS (section 7). Partition
  # coefficients are taken unscaled from PCs_original.csv.
  #
  # If a compound has no quantified sample in exp_data_feces.csv, kbilec is set to 0
  # and removed from the estimated parameters.
  #
  # Experimental data (exp_data_plasma.csv, exp_data_urine.csv, exp_data_feces.csv)
  # are screened per compound: samples that are NA/empty, reported as <LOQ / <LOD
  # (or ND, BLQ, BDL, "<value"), or reported as 0 are flagged and EXCLUDED from the
  # AAFE calculation. A matrix with no quantified sample for a compound is dropped
  # from the objective function automatically.
  # ================================================================================

  library(tidyverse)
  library(nloptr)
  library(deSolve)
  library(PKNCA)

  # ================================================================================
  # USER INPUTS
  # ================================================================================
  compound  <- "PFOA"  # PFAS to estimate: PFBA, PFHxA, PFHpA, PFOA, PFNA, PFDA,
                       # PFBS, PFHxS, PFOS, HFPO_DA or DONA
  weight    <- 10      # weight applied to the plasma AAFE of samples with time > threshold
  threshold <- 20      # time (days) after which plasma samples are weighted; both unused
                       # when metric = "halflife" (below)

  # metric: what the optimizer minimizes.
  #   "AAFE"     - mean of the (weighted) plasma/urine/feces AAFE, as before.
  #   "rmsd"     - same, with RMSD on plasma instead of AAFE.
  #   "halflife" - the fold-discrepancy between this candidate's predicted terminal half-life
  #                (from a lean NCA fit on the simulated profile) and the experimental
  #                half-life, on the same 10^|log10(ratio)| scale as AAFE (so < 2 means within
  #                2-fold). Plasma/urine/feces AAFE is not part of the search target at all in
  #                this mode, only reported afterwards for information.
  metric <- "AAFE"

  # success_aafe_threshold: the AAFE bar a fit must clear (on all matrices present) to count
  # as "successful" (see aafe_success()). Half-life is still reported for every run but does
  # not gate success.
  success_aafe_threshold <- 2

  N_iter <- 2000  # max evaluations of the objective

  # Running from the command line overrides the values above:
  #   Rscript "Parameter Estimation_all_PFAS.r" <PFAS> <weight> <threshold> [<N_iter>]
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    if (!length(args) %in% c(3, 4)) {
      stop("Usage: Rscript \"Parameter Estimation_all_PFAS.r\" <PFAS> <weight> <threshold> [<N_iter>]")
    }
    compound  <- args[1]
    weight    <- as.numeric(args[2])
    threshold <- as.numeric(args[3])
    if (length(args) == 4) N_iter <- as.integer(args[4])
  }

  # --- Run settings ---
  print_level <- 0                # nloptr verbosity (1 = print every iteration)
  max_time_excreta <- 6           # urine/feces samples used up to this time (days)
  zero_is_censored <- TRUE        # reported 0 values are treated as <LOD
  censored_mass_contribution <- 0 # mass (ug) added to cumulative excreta for an excluded sample
  output_dir <- paste0("Optimized_parameters_generalized_", compound)  # one folder per PFAS

  # Data files live in "Extended PFAS PBK model/"; allow running from there or
  # from Parameter_estimation_codes/
  data_dir <- if (file.exists("initial_parameters.csv")) "." else ".."
  output_dir <- file.path(data_dir, output_dir)

  # --- Load initial parameters ---
  variables <- read.csv(file.path(data_dir, "initial_parameters.csv"),
                        row.names = "Parameters", fileEncoding = "UTF-8-BOM")

  # --- Load partition coefficients (used as is, no scaling) ---
  PC <- read.csv(file.path(data_dir, "PCs_original.csv"), row.names = "Organs",
                 fileEncoding = "UTF-8-BOM")

  # --- Check the user inputs ---
  if (!compound %in% names(variables) || !compound %in% names(PC)) {
    stop("Unknown PFAS '", compound, "'. Available: ",
         paste(intersect(names(variables), names(PC)), collapse = ", "))
  }
  if (length(weight) != 1 || is.na(weight) || weight <= 0) stop("weight must be a positive number")
  if (length(threshold) != 1 || is.na(threshold) || threshold < 0) {
    stop("threshold must be a non-negative number (days)")
  }
  if (length(N_iter) != 1 || is.na(N_iter) || N_iter < 1) stop("N_iter must be a positive integer")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Results will be written to:", normalizePath(output_dir), "\n")


  #=========================
  # 1. Parameters of the model
  #=========================

  create.params <- function(user_input, variables, PC, compound) {
    with( as.list(user_input),{
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
    MW <- variables["MW",compound]  # molecular mass (g/mol)
    Free <- variables["Free",compound]  # free fraction in plasma (Ryu 2024)

    # --- Kidney Transport Parameters ---
    Vmax_baso_invitro <- variables["Vmax_baso_invitro",compound]  # Vmax of basolateral transporter (pmol/mg protein/min)
    Km_baso <- variables["Km_baso",compound] * variables["MW",compound]  # Km of basolateral transporter (ug/L)
    Vmax_apical_invitro <- variables["Vmax_apical_invitro", compound]  # Vmax of apical transporter (pmol/mg protein/min)
    Km_apical <- variables["Km_apical", compound] * variables["MW", compound]  # Km of apical transporter (ug/L)
    protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
    GFRC <- 24.19 * 24  # glomerular filtration rate (L/day/kg kidney); Corley 2005

    # --- Partition Coefficients (from Allendorf 2021, PCs_original.csv) ---
    PAdi <- (1-Htc)* PC["Adipose",compound] # adipose tissue:plasma;
    PBra <- (1-Htc)* PC["Brain",compound] # brain:plasma;
    PGon <- (1-Htc)* PC["Gonads",compound] # gonads:plasma;
    PGI <-  (1-Htc)* PC["Gut",compound] # GI tract:plasma;
    PHea <- (1-Htc)* PC["Heart",compound] # heart:plasma;
    PL <-   (1-Htc)* PC["Liver",compound] # liver:plasma;
    PLun <- (1-Htc)* PC["Lung",compound]  # lung:plasma;
    PMus <- (1-Htc)* PC["Muscle",compound] # muscle:plasma;
    PSki <- (1-Htc)* PC["Skin",compound]   # skin:plasma;
    PSpl <- (1-Htc)* PC["Spleen",compound]  # spleen:plasma;
    PPan <- (PGI+PSpl)/2 # pancreas:plasma; estimated as average of GI tract and spleen
    PR <- 0.1 # rest of body:blood

    # --- Rate Constants ---
    kdif <- 0.001 * 24  # diffusion rate from proximal tubule cells (L/day)
    kabsc <-  2.12 *24  # rate of absorption from small intestine (1/(day*BW^-0.25))
    kunabsc <- 7.06e-5 * 24  # rate of unabsorbed dose to feces (1/(day*BW^-0.25)); fitted to model
    # kurinec: urinary elimination rate (1/(day*BW^-0.25)); estimated (see fit_params), so it
    # comes from user_input, not hardcoded here (initial value 0.063*24 = 1.512 in
    # initial_parameters.csv, matching the constant this replaced).


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

    Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (variables["MW", compound] / 1e12) * 1e6) * 24
    Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (variables["MW", compound] / 1e12) * 1e6) * 24
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
  #4. Custom functions
  #==================

  # --- Screening of experimental measurements ---
  # Classifies each raw entry as:
  #   "quantified" : numeric value > 0
  #   "missing"    : NA, empty cell, "NA", "-", ...
  #   "<LOQ"       : entries mentioning LOQ/LLOQ/BLQ/BQL/NQ
  #   "<LOD"       : entries mentioning LOD/BDL/ND, or reported as 0 (if zero_is_censored)
  #   "<value"     : entries such as "<0.01" (censored at a reported limit)
  #   "invalid"    : any other non-numeric text, or a negative value
  # Only "quantified" samples enter the AAFE.
  classify_measurements <- function(raw, zero_is_censored = TRUE) {
    txt <- toupper(trimws(as.character(raw)))
    value <- suppressWarnings(as.numeric(txt))
    status <- rep("quantified", length(txt))

    is_missing <- is.na(txt) | txt %in% c("", "NA", "N/A", "NAN", "-", "--", ".")
    is_loq <- !is_missing & grepl("LOQ|BLQ|BQL|^NQ$", txt)
    is_lod <- !is_missing & !is_loq & grepl("LOD|BDL|^N\\.?D\\.?$", txt)
    is_lt  <- !is_missing & !is_loq & !is_lod & grepl("^<", txt)

    status[is.na(value)] <- "invalid"
    status[!is.na(value) & value < 0] <- "invalid"
    if (zero_is_censored) status[!is.na(value) & value == 0] <- "<LOD"
    status[is_lt] <- "<value"
    status[is_lod] <- "<LOD"
    status[is_loq] <- "<LOQ"
    status[is_missing] <- "missing"

    bad <- status == "invalid"
    if (any(bad)) warning("Unrecognised entries treated as excluded: ",
                          paste(unique(raw[bad]), collapse = ", "))

    value[status != "quantified"] <- NA
    data.frame(value = value, status = status, stringsAsFactors = FALSE)
  }

  # Reads an experimental data file keeping every cell as text, so that entries like
  # "<LOQ" survive reading and can be classified
  read_exp_data <- function(file) {
    df <- read.csv(file, colClasses = "character", na.strings = character(0),
                   check.names = FALSE, fileEncoding = "UTF-8-BOM")
    names(df) <- make.names(trimws(names(df)))
    df
  }

  AAFE <- function(predictions, observations, weight = 1, threshold = Inf, times = NULL,
                   valid = NULL){
    y_obs <- as.numeric(unlist(observations))
    y_pred <- as.numeric(unlist(predictions))
    if (length(y_obs) != length(y_pred)) {
      stop("AAFE: predictions and observations have different lengths")
    }
    # valid: logical vector from the data screening; FALSE for samples that are
    # NA, <LOQ, <LOD or 0 in the experimental data files
    if (is.null(valid)) valid <- rep(TRUE, length(y_obs))
    valid_indices <- which(valid & !is.na(valid) &
                           !is.na(y_obs) & !is.na(y_pred) &
                           y_obs > 0 & y_pred > 0)

    if(length(valid_indices) == 0) {
      warning("No valid observations for AAFE calculation (all NA, <LOQ, <LOD or zero)")
      return(NA)
    }
    y_obs <- y_obs[valid_indices]
    y_pred <- y_pred[valid_indices]

    log_ratio <- abs(log10(y_pred / y_obs))
    if(!is.null(times)) {
      times <- times[valid_indices]
      log_ratio <- ifelse(times > threshold, weight * log_ratio, log_ratio)
    }
    aafe <- 10^(sum(log_ratio)/length(log_ratio))
    return(aafe)
  }

  rmsd <- function(predictions, observations, weight = 1, threshold = Inf, times = NULL,
                   valid = NULL){
    y_obs <- as.numeric(unlist(observations))
    y_pred <- as.numeric(unlist(predictions))
    if (is.null(valid)) valid <- rep(TRUE, length(y_obs))
    keep <- which(valid & !is.na(y_obs) & !is.na(y_pred))
    if (length(keep) == 0) return(NA)
    w <- if (is.null(times)) 1 else ifelse(times[keep] > threshold, weight, 1)
    sqrt(sum((w * (y_obs[keep] - y_pred[keep]))^2) / length(keep))
  }

  # Cumulative excreted mass for urine and feces. Excluded samples add
  # `censored_contribution` to the running total and are flagged valid = FALSE,
  # so the cumulative curve stays continuous but those points are not scored.
  # Samples without urine volume / feces weight cannot be converted to mass and are dropped.
  cumulative_exp_data <- function(time, conc, status, amount, censored_contribution = 0){
    keep <- !is.na(time) & !is.na(amount)
    time <- time[keep]; conc <- conc[keep]; status <- status[keep]; amount <- amount[keep]
    ord <- order(time)
    time <- time[ord]; conc <- conc[ord]; status <- status[ord]; amount <- amount[ord]

    valid <- status == "quantified"
    mass <- ifelse(valid, conc * amount, censored_contribution)
    data.frame(time = time,
               cumulative_mass = cumsum(mass),
               status = status,
               valid = valid)
  }

  find_nearest <- function(times, target) {
      sapply(target, function(t) which.min(abs(times - t)))}

  # Builds the plasma / urine / feces observation tables for one compound.
  # A matrix without any quantified sample is returned as NULL.
  prepare_exp_data <- function(compound, data_dir) {
    col <- make.names(compound)
    out <- list()

    plasma_raw <- read_exp_data(file.path(data_dir, "exp_data_plasma.csv"))
    urine_raw  <- read_exp_data(file.path(data_dir, "exp_data_urine.csv"))
    feces_raw  <- read_exp_data(file.path(data_dir, "exp_data_feces.csv"))

    # Plasma (ug/L, time in days)
    if (col %in% names(plasma_raw)) {
      cl <- classify_measurements(plasma_raw[[col]], zero_is_censored)
      out$plasma <- data.frame(time = as.numeric(plasma_raw$time),
                               obs = cl$value, status = cl$status,
                               valid = cl$status == "quantified")
    }

    # Urine: concentration (ug/L) x volume (L) -> ug; time h -> days
    if (col %in% names(urine_raw)) {
      cl <- classify_measurements(urine_raw[[col]], zero_is_censored)
      out$urine <- cumulative_exp_data(as.numeric(urine_raw$time) / 24, cl$value, cl$status,
                                       as.numeric(urine_raw$urine.volume),
                                       censored_mass_contribution) %>%
        filter(time <= max_time_excreta) %>%
        rename(obs = cumulative_mass)
    }

    # Feces: concentration (ng/g) x weight (g) / 1000 -> ug; time in days
    if (col %in% names(feces_raw)) {
      cl <- classify_measurements(feces_raw[[col]], zero_is_censored)
      out$feces <- cumulative_exp_data(as.numeric(feces_raw$time), cl$value / 1000, cl$status,
                                       as.numeric(feces_raw$feces.weight),
                                       censored_mass_contribution) %>%
        filter(time <= max_time_excreta) %>%
        rename(obs = cumulative_mass)
    }

    # Summary of what is used / excluded per matrix
    summary <- bind_rows(lapply(names(out), function(m) {
      d <- out[[m]]
      data.frame(compound = compound, matrix = m,
                 n_total = nrow(d), n_used = sum(d$valid),
                 n_missing = sum(d$status == "missing"),
                 n_LOQ = sum(d$status == "<LOQ"),
                 n_LOD = sum(d$status %in% c("<LOD", "<value")),
                 n_invalid = sum(d$status == "invalid"))
    }))

    usable <- out[vapply(out, function(d) sum(d$valid) > 0, logical(1))]
    for (m in setdiff(names(out), names(usable))) {
      message(sprintf("  %s: no quantified %s samples -> %s excluded from the objective",
                      compound, m, m))
    }
    if (is.null(usable$plasma)) stop("No quantified plasma samples for ", compound)

    # all: every screened sample (for reporting); data: matrices usable in the objective
    list(all = out, data = usable, summary = summary)
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
      CA <- CA_free / Free  # concentration of total PFAS in plasma (ug/L)

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

  # --- Dosing and Administration Settings (dose is compound specific) ---
  admin_type <- "bolus"  # administration type: "iv", "oral", or "bolus"
  admin_time_bolus <- c(0)  # time when bolus doses are administered (days)
  admin_dose_iv <- 0  # administered dose through IV (ug)
  admin_time_iv <- 0  # time when IV doses are administered (days)
  Cwater <- 0.00  # concentration in water (ug/L)
  Cwater_time <- 0  # time of water concentration change
  ingestion <- 0  # ingestion rate (ug/day)
  ingestion_time <- c(0)  # time of ingestion rate change

  build_user_input <- function(theta, compound, variables) {
    list("admin_type" = admin_type,
         "admin_dose_bolus" = c(variables["admin_dose_bolus", compound]),
         "admin_time_bolus" = admin_time_bolus,
         "admin_dose_iv" = admin_dose_iv,
         "admin_time_iv" = admin_time_iv,
         "Cwater" = Cwater,
         "Cwater_time" = Cwater_time, "ingestion" = ingestion,
         "ingestion_time" = ingestion_time,
         "RAFbaso" = theta[["RAFbaso"]],
         "RAFapi" = theta[["RAFapi"]],
         "kbilec" = theta[["kbilec"]],
         "keffluxc" = theta[["keffluxc"]],
         "kurinec" = theta[["kurinec"]])
  }

  # ctx$variables / ctx$PC carry the point-estimate parameter tables (globals) for this run.
  run_model <- function(theta, ctx, sample_time) {
    user_input <- build_user_input(theta, ctx$compound, ctx$variables)
    params <- create.params(user_input, ctx$variables, ctx$PC, ctx$compound)
    inits <- create.inits(user_input)
    events <- create.events(user_input)
    solution <- data.frame(deSolve::ode(times = sample_time, func = ode.func, y = inits,
                                        parms = params, events = events,
                                        method = "lsoda", rtol = 1e-05, atol = 1e-05))
    list(solution = solution, params = params)
  }


  #===================================
  #7. Parameters to estimate and bounds
  #===================================
  # Initial values are taken from initial_parameters.csv. kbilec is fixed to 0 and
  # not estimated when the compound has no quantified feces sample.
  # Bounds were collected from the individual "Parameter Estimation_<PFAS>.r" scripts:
  # lower = smallest non-zero lower bound, upper = largest upper bound per parameter.
  # kurinec: urinary elimination rate constant, added to the fit because RAFbaso/RAFapi/
  # keffluxc only govern a kidney-blood <-> PTC recycling loop (uptake, reabsorption, efflux
  # back to blood) with no path out of the body; kurinec (along with kbilec) is one of only
  # two constants that actually set the terminal elimination rate, and it was previously
  # fixed at 0.063*24 = 1.512 for every compound.
  #
  # fit_params: every kinetic parameter theta0/the bounds tables know about (so any parameter
  # not in estimated_params still gets its initial_parameters.csv point estimate).
  # estimated_params: the subset actually optimized.
  fit_params <- c("RAFbaso", "RAFapi", "kbilec", "keffluxc", "kurinec")
  # estimated_params: the subset actually optimized; kbilec/keffluxc/kurinec stay fixed at
  # their initial_parameters.csv point estimate.
  estimated_params <- c("RAFbaso", "RAFapi")
  lower_bounds_all <- c("RAFbaso" = 1e-07,   # PFOS
                        "RAFapi" = 1e-07,    # HFPO_DA, PFHxA, PFOS
                        "kbilec" = 1e-05,    # DONA, HFPO_DA, PFDA, PFHpA, PFHxA, PFNA
                        "keffluxc" = 1e-05,  # PFOS
                        "kurinec" = 1e-05)
  upper_bounds_all <- c("RAFbaso" = 1e05,    # PFOA
                        "RAFapi" = 1e04,
                        "kbilec" = 1e04,
                        "keffluxc" = 1e04,
                        "kurinec" = 1e02)
  # all_matrices: matrices allowed into the objective (cfg$matrices in estimate_pfas).
  # Restricted to plasma only per current request; a matrix would still need quantified
  # samples to enter (prepare_exp_data() drops any that don't, regardless of this list).
  all_matrices <- c("plasma")


  #=====================
  #8. Objective function
  #=====================
  # Plasma is scored with the time-weighted AAFE; cumulative urine and feces with the
  # unweighted AAFE. Samples flagged invalid by the data screening are excluded.
  score_fit <- function(solution, ctx, weight, threshold, metric = "AAFE") {
    scores <- c()
    for (m in names(ctx$data)) {
      d <- ctx$data[[m]]
      state <- switch(m, plasma = "CA", urine = "Aurine", feces = "Afeces")
      idx <- find_nearest(solution$time, d$time)
      preds <- solution[idx, state]
      if (any(is.na(preds))) return(NA)
      if (m == "plasma") {
        f <- if (metric == "rmsd") rmsd else AAFE
        scores[m] <- f(preds, d$obs, weight = weight, threshold = threshold,
                       times = solution$time[idx], valid = d$valid)
      } else {
        scores[m] <- AAFE(preds, d$obs, valid = d$valid)
      }
    }
    scores
  }

  make_objective <- function(ctx, metric = "AAFE") {
    obs_times <- unlist(lapply(ctx$data, function(d) d$time))
    sample_time <- sort(unique(c(seq(0, ctx$cfg$sim_end, 1), obs_times)))

    function(x) {
      theta <- ctx$theta0
      theta[ctx$cfg$fit] <- x
      solution <- tryCatch(suppressWarnings(run_model(theta, ctx, sample_time)$solution),
                           error = function(e) NULL)
      if (is.null(solution) || nrow(solution) < length(sample_time)) {
        cat("ODE failure for parameters:", x, "\n")
        return(1e6)
      }
      if (all(solution$CA == 0)) {
        cat("Zero predictions for parameters:", x, "\n")
        return(1e6)
      }

      if (metric == "halflife") {
        m <- suppressWarnings(tryCatch(
          nca_metrics(solution$time, solution$CA, ctx$dose),
          error = function(e) c(cmax = NA_real_, half.life = NA_real_, auclast = NA_real_)))
        pred_hl <- unname(m["half.life"])
        if (is.na(pred_hl) || pred_hl <= 0) {
          cat("NCA half-life unavailable for parameters:", x, "\n")
          return(1e6)
        }
        return(10^abs(log10(pred_hl / ctx$exp_half_life)))
      }

      scores <- suppressWarnings(score_fit(solution, ctx, ctx$cfg$weight,
                                           ctx$cfg$threshold, metric))
      if (length(scores) == 0 || any(is.na(scores))) {
        cat("NA score for parameters:", x, "\n")
        return(1e6)
      }
      mean(scores)
    }
  }


  #=====================================
  #9. Estimation, plots, NCA per compound
  #=====================================
  run_nca <- function(ctx, solution, dose) {
    p <- ctx$data$plasma
    conc_exp <- data.frame(Subject = "Experimental", Time = p$time[p$valid],
                           Concentration = p$obs[p$valid])
    conc_pred <- data.frame(Subject = "Predicted", Time = solution$time,
                            Concentration = solution$CA)
    # Explicit time = 0, concentration = 0 rows so AUC starts from t = 0
    zero_rows <- data.frame(Subject = c("Experimental", "Predicted"), Time = 0,
                            Concentration = 0)
    conc_df <- rbind(conc_exp, conc_pred, zero_rows)
    conc_df <- conc_df[!duplicated(conc_df[, c("Subject", "Time")]), ]
    conc_df <- conc_df[order(conc_df$Subject, conc_df$Time), ]

    dose_df <- data.frame(Subject = c("Experimental", "Predicted"), Time = 0,
                          Dose = dose, Route = "extravascular")

    my_conc <- PKNCAconc(conc_df, Concentration ~ Time | Subject)
    my_dose <- PKNCAdose(dose_df, Dose ~ Time | Subject)
    my_data <- PKNCAdata(
      my_conc, my_dose,
      intervals = data.frame(
        start = 0, end = max(conc_df$Time, na.rm = TRUE),
        cmax = TRUE, tmax = TRUE, auclast = TRUE,
        aucinf.obs = TRUE, aucinf.pred = TRUE,
        half.life = TRUE, lambda.z = TRUE, r.squared = TRUE,
        stringsAsFactors = FALSE))
    as.data.frame(pk.nca(my_data))
  }

  # Lean single-subject NCA: Cmax, half-life and AUClast for one concentration-time series.
  # Used to compute the experimental NCA once and to check the "halflife" objective metric
  # without rebuilding the two-subject table run_nca() returns.
  nca_metrics <- function(time, conc, dose) {
    df <- data.frame(Subject = "S", Time = time, Concentration = conc)
    df <- rbind(df, data.frame(Subject = "S", Time = 0, Concentration = 0))
    df <- df[!duplicated(df[, c("Subject", "Time")]), ]
    df <- df[order(df$Time), ]
    dose_df <- data.frame(Subject = "S", Time = 0, Dose = dose, Route = "extravascular")
    my_conc <- PKNCAconc(df, Concentration ~ Time | Subject)
    my_dose <- PKNCAdose(dose_df, Dose ~ Time | Subject)
    my_data <- PKNCAdata(
      my_conc, my_dose,
      intervals = data.frame(start = 0, end = max(df$Time, na.rm = TRUE),
                             cmax = TRUE, auclast = TRUE, half.life = TRUE, lambda.z = TRUE,
                             stringsAsFactors = FALSE))
    res <- as.data.frame(pk.nca(my_data))
    vals <- setNames(res$PPORRES, res$PPTESTCD)
    c(cmax = unname(vals["cmax"]), half.life = unname(vals["half.life"]),
      auclast = unname(vals["auclast"]))
  }

  # A predicted metric counts as "close" to its experimental counterpart within a 2-fold
  # window, matching the AAFE < 2 acceptance threshold used elsewhere.
  is_close <- function(pred, exp, fold = 2) {
    !is.na(pred) & !is.na(exp) & pred > 0 & exp > 0 &
      (pred / exp) >= 1 / fold & (pred / exp) <= fold
  }

  # AAFE-only success flag: AAFE_plasma < threshold, and AAFE_urine/AAFE_feces < threshold
  # where computed (NA when that matrix wasn't in the objective, e.g. no quantified samples).
  aafe_success <- function(row, threshold = success_aafe_threshold) {
    row$AAFE_plasma < threshold &
      (is.na(row$AAFE_urine) | row$AAFE_urine < threshold) &
      (is.na(row$AAFE_feces) | row$AAFE_feces < threshold)
  }

  # One predicted-vs-observed plot per matrix present in ctx$data (plasma always; urine/feces
  # when quantified samples exist), returned as a named list keyed by matrix name.
  plot_fit <- function(ctx, solution) {
    compound <- ctx$compound
    plots <- list()

    if (!is.null(ctx$data$plasma)) {
      p <- ctx$data$plasma %>% mutate(used = ifelse(valid, "used", "excluded"))
      plots$plasma <- ggplot() +
        geom_line(data = filter(solution, time > 0), aes(x = time, y = CA), linewidth = 1.1) +
        geom_point(data = filter(p, !is.na(obs)), aes(x = time, y = obs, shape = used),
                   size = 3, colour = "red") +
        scale_y_log10() +
        scale_shape_manual(values = c(used = 16, excluded = 1), drop = FALSE) +
        labs(title = paste(compound, "- plasma"), x = "Time (days)",
             y = "Concentration (ug/L)", shape = "Observation") +
        theme_bw()
    }

    excreta_state <- c(urine = "Aurine", feces = "Afeces")
    for (m in intersect(names(excreta_state), names(ctx$data))) {
      d <- ctx$data[[m]] %>% mutate(used = ifelse(valid, "used", "excluded"))
      sim <- data.frame(time = solution$time, mass = solution[[excreta_state[[m]]]]) %>%
        filter(time <= max_time_excreta + 4)
      plots[[m]] <- ggplot() +
        geom_line(data = sim, aes(x = time, y = mass), linewidth = 1.1) +
        geom_point(data = d, aes(x = time, y = obs, shape = used), size = 3, colour = "red") +
        scale_shape_manual(values = c(used = 16, excluded = 1), drop = FALSE) +
        labs(title = paste(compound, "-", m, "(cumulative)"), x = "Time (days)",
             y = "Mass (ug)", shape = "Observation") +
        theme_bw()
    }
    plots
  }

  estimate_pfas <- function(compound, weight, threshold, metric = "AAFE") {
    cat("\n==================== ", compound, " ====================\n")
    cat("weight =", weight, "| threshold =", threshold, "days\n")
    # sim_end: end of the simulation grid used in the objective function (days)
    cfg <- list(fit = estimated_params, weight = weight, threshold = threshold,
                sim_end = 100, matrices = all_matrices)

    exp <- prepare_exp_data(compound, data_dir)
    data <- exp$data[intersect(cfg$matrices, names(exp$data))]
    print(exp$summary, row.names = FALSE)
    cat("Matrices in objective:", paste(names(data), collapse = ", "), "\n")

    # Initial values from initial_parameters.csv
    theta0 <- setNames(sapply(fit_params, function(p) variables[p, compound]), fit_params)

    # No quantified feces sample -> no information on biliary excretion:
    # kbilec is fixed to 0 and not estimated
    if (is.null(exp$data$feces)) {
      theta0[["kbilec"]] <- 0
      cfg$fit <- setdiff(cfg$fit, "kbilec")
      cat("No feces data: kbilec fixed to 0 and excluded from estimation\n")
    }
    cat("Estimated parameters:", paste(cfg$fit, collapse = ", "), "\n")

    lower <- lower_bounds_all[cfg$fit]
    upper <- upper_bounds_all[cfg$fit]
    x0 <- pmin(pmax(theta0[cfg$fit], lower), upper)
    print(data.frame(initial = x0, lower = lower, upper = upper))

    ctx <- list(compound = compound, cfg = cfg, data = data, theta0 = theta0,
                variables = variables, PC = PC)

    if (metric == "halflife") {
      p <- data$plasma
      ctx$dose <- variables["admin_dose_bolus", compound]
      ctx$exp_half_life <- tryCatch(
        unname(nca_metrics(p$time[p$valid], p$obs[p$valid], ctx$dose)["half.life"]),
        error = function(e) NA_real_)
      if (is.na(ctx$exp_half_life) || ctx$exp_half_life <= 0) {
        warning("Experimental half-life could not be computed for ", compound,
                "; falling back to metric = 'AAFE'")
        metric <- "AAFE"
      } else {
        cat("Objective = half-life fold-discrepancy (predicted vs experimental); ",
            "experimental half-life =", ctx$exp_half_life, "days\n")
      }
    }

    opts <- list("algorithm" = "NLOPT_LN_SBPLX",
                 "xtol_rel" = 1e-05,
                 "ftol_rel" = 1e-05,
                 "ftol_abs" = 0.0,
                 "xtol_abs" = 0.0,
                 "maxeval" = N_iter,
                 "print_level" = print_level)

    # The fitted parameters span many orders of magnitude (e.g. RAFapi's bounds are 1e-07 to
    # 1e04). NLOPT_LN_SBPLX is a local pattern-search method that takes steps sized relative
    # to the starting simplex in the space it searches; in linear space that means it can get
    # stuck exploring only within a small multiple of x0 and never reach a solution many
    # orders of magnitude away. Searching in log10 space instead lets an order-of-magnitude
    # move cost the same as a fractional one, which linear-space search cannot offer.
    obj <- make_objective(ctx, metric)
    log_objective <- function(log_x) obj(10^log_x)
    optimization <- nloptr::nloptr(x0 = log10(x0), eval_f = log_objective,
                                   lb = log10(lower), ub = log10(upper), opts = opts)
    optimization$solution <- 10^optimization$solution  # back to linear scale for reporting

    theta_opt <- theta0
    theta_opt[cfg$fit] <- optimization$solution
    cat("Objective:", optimization$objective, "\n")
    print(theta_opt)

    # Simulation with the optimized parameters
    obs_times <- unlist(lapply(data, function(d) d$time))
    sample_time <- sort(unique(c(seq(0, 450, 1), obs_times)))
    fit <- run_model(theta_opt, ctx, sample_time)
    solution <- fit$solution
    scores <- score_fit(solution, ctx, 1, Inf)  # unweighted AAFE per matrix

    plots <- plot_fit(ctx, solution)
    for (m in names(plots)) {
      ggsave(file.path(output_dir, paste0("Fit_", m, "_", compound, ".png")), plots[[m]],
             width = 8, height = 6)
    }

    nca <- tryCatch(run_nca(ctx, solution, variables["admin_dose_bolus", compound]),
                    error = function(e) { cat("NCA failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(nca)) {
      write.csv(nca, file.path(output_dir, paste0("NCA_results_", compound, ".csv")),
                row.names = FALSE)
    }
    write.csv(c(fit$params, as.list(theta_opt)),
              file.path(output_dir, paste0("Optimized_Parameters_", compound, ".csv")),
              row.names = TRUE)

    # Sample-level record of what entered the AAFE
    screening <- bind_rows(lapply(names(exp$all), function(m) {
      d <- exp$all[[m]]
      data.frame(compound = compound, matrix = m, time = d$time, obs = d$obs,
                 status = d$status, used_in_AAFE = d$valid & m %in% names(data))
    }))

    list(compound = compound,
         fitted = cfg$fit,
         weight = weight,
         threshold = threshold,
         optimization = optimization,
         theta = theta_opt,
         scores = scores,
         screening = screening,
         data_summary = exp$summary,
         solution = solution,
         plots = plots,
         nca = nca)
  }


  #=====================
  #10. Run the estimation
  #=====================
  summary_row <- function(res) {
    data.frame(
      compound = res$compound,
      weight = res$weight,
      threshold = res$threshold,
      objective = res$optimization$objective,
      status = res$optimization$status,
      iterations = res$optimization$iterations,
      fitted = paste(res$fitted, collapse = ", "),
      RAFbaso = res$theta[["RAFbaso"]], RAFapi = res$theta[["RAFapi"]],
      kbilec = res$theta[["kbilec"]], keffluxc = res$theta[["keffluxc"]],
      kurinec = res$theta[["kurinec"]],
      AAFE_plasma = unname(res$scores["plasma"]),
      AAFE_urine = unname(res$scores["urine"]),
      AAFE_feces = unname(res$scores["feces"]))
  }

  result <- estimate_pfas(compound, weight, threshold, metric = metric)

  estimation_summary <- summary_row(result)
  print(estimation_summary)

  cat(sprintf("\nAAFE success (< %g on plasma, and on urine/feces where present): %s\n",
             success_aafe_threshold, aafe_success(estimation_summary)))

  if (!is.null(result$nca)) {
    get_metric <- function(subject, param) {
      row <- result$nca[result$nca$Subject == subject & result$nca$PPTESTCD == param, ]
      if (nrow(row) == 1) row$PPORRES else NA_real_
    }
    hl_pred <- get_metric("Predicted", "half.life")
    hl_exp <- get_metric("Experimental", "half.life")
    cat(sprintf("Half-life (informational): predicted = %.1f days, experimental = %.1f days",
               hl_pred, hl_exp))
    cat(sprintf(" -> %s (fold = %.2f)\n",
               if (is_close(hl_pred, hl_exp)) "within 2-fold" else "not within 2-fold",
               hl_pred / hl_exp))
  }

  write.csv(estimation_summary,
            file.path(output_dir, paste0("Estimation_summary_", compound, ".csv")),
            row.names = FALSE)
  write.csv(result$data_summary,
            file.path(output_dir, paste0("Data_screening_summary_", compound, ".csv")),
            row.names = FALSE)
  write.csv(result$screening,
            file.path(output_dir, paste0("Data_screening_samples_", compound, ".csv")),
            row.names = FALSE)

  for (m in names(result$plots)) print(result$plots[[m]])
