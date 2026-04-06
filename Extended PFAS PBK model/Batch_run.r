library(deSolve)
library(tidyverse)

#=========================
# PFAS Compound List
#=========================
pfas_compounds <- c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS", 
                    "PFHxS", "PFOS", "DONA", "HFPO_DA", "PFBA", "PFHxA")

#=========================
# PK PARAMETER CALCULATION FUNCTION
#=========================
calculate_pk_parameters <- function(time, conc, route = c("IV", "Oral"), 
                                    plot = TRUE, min_points = 3, 
                                    output_dir = NULL, pfas_name = NULL,
                                    data_type = "Model") {
  # --- 1. Input Validation & Cleaning ---
  route <- match.arg(route)
  
  if (length(time) != length(conc)) {
    stop("Time and Concentration vectors must be the same length.")
  }
  
  # Create data frame and remove NA or Zero/Negative concentrations for log calculations
  df <- data.frame(time = time, conc = conc)
  df_clean <- df[!is.na(df$conc) & df$conc > 0, ]
  df_clean <- df_clean[order(df_clean$time), ]
  
  if (nrow(df_clean) < min_points) {
    warning(paste("Not enough valid data points. Need at least", min_points))
    return(list(
      half_life = NA, lambda_z = NA, cmax = NA, tmax = NA, 
      auc_0_last = NA, auc_0_inf = NA, reason = "Insufficient data points",
      data_type = data_type
    ))
  }
  
  # --- 2. Cmax and Tmax Calculation ---
  cmax <- max(df$conc, na.rm = TRUE)
  tmax <- df$time[which.max(df$conc)]
  
  # --- 3. AUC Calculation (Trapezoidal Rule) ---
  # AUC from time 0 to last measured time point
  auc_0_last <- 0
  for (i in 2:nrow(df)) {
    dt <- df$time[i] - df$time[i-1]
    avg_conc <- (df$conc[i] + df$conc[i-1]) / 2
    auc_0_last <- auc_0_last + (dt * avg_conc)
  }
  
  # --- 4. Route Specific Filtering for Half-Life ---
  if (route == "Oral") {
    max_conc_idx <- which.max(df_clean$conc)
    t_max <- df_clean$time[max_conc_idx]
    df_elim <- df_clean[df_clean$time >= t_max, ]
    
    if (nrow(df_elim) < min_points) {
      warning("Oral route selected, but not enough points after Tmax to calculate half-life.")
      return(list(
        half_life = NA, lambda_z = NA, cmax = cmax, tmax = tmax,
        auc_0_last = auc_0_last, auc_0_inf = NA,
        reason = "Insufficient post-Tmax data", data_type = data_type
      ))
    }
  } else {
    df_elim <- df_clean
  }
  
  # --- 5. Automated Terminal Phase Selection ---
  n <- nrow(df_elim)
  best_model <- NULL
  best_adj_r2 <- -Inf
  best_range <- NULL
  
  for (i in 1:(n - min_points + 1)) {
    subset_data <- df_elim[i:n, ]
    model <- lm(log(conc) ~ time, data = subset_data)
    slope <- coef(model)[2]
    
    if (slope < 0) {
      adj_r2 <- summary(model)$adj.r.squared
      if (adj_r2 > best_adj_r2) {
        best_adj_r2 <- adj_r2
        best_model <- model
        best_range <- subset_data
      }
    }
  }
  
  # --- 6. Half-Life Calculation ---
  if (is.null(best_model)) {
    warning("Could not identify a terminal elimination phase (no negative slope found).")
    return(list(
      half_life = NA, lambda_z = NA, cmax = cmax, tmax = tmax,
      auc_0_last = auc_0_last, auc_0_inf = NA,
      reason = "No negative slope found", data_type = data_type
    ))
  }
  
  lambda_z <- -coef(best_model)[2]
  half_life <- log(2) / lambda_z
  
  # AUC from 0 to infinity (AUC_0-inf = AUC_0-last + Clast/lambda_z)
  clast <- best_range$conc[nrow(best_range)]
  tlast <- best_range$time[nrow(best_range)]
  auc_extra <- clast / lambda_z
  auc_0_inf <- auc_0_last + auc_extra
  
  
  # --- 8. Return Results ---
  return(list(
    half_life = half_life,
    lambda_z = lambda_z,
    cmax = cmax,
    tmax = tmax,
    auc_0_last = auc_0_last,
    auc_0_inf = auc_0_inf,
    clast = clast,
    tlast = tlast,
    route = route,
    points_used = nrow(best_range),
    adj_r_squared = best_adj_r2,
    terminal_data = best_range,
    data_type = data_type
  ))
}

#=========================
# 1. Parameters of the model
#=========================
variables <- read.csv("estimated_parameters.csv", row.names = "Parameters")
PC <- read.csv("PCs.csv", row.names = "Organs")

create.params <- function(variables, PC, BW, pfas_name) {
  QCC <- 12.5 * 24
  QLC <- 0.065
  QKC <- 0.175
  QAdiC <- 0.05
  QBraC <- 0.12
  QGonC <- 0.0005
  QHeaC <- 0.04
  QLunC <- 0.025
  QMusC <- 0.17
  QSkiC <- 0.05
  QSplC <- 0.03
  QPanC <- 0.01
  QGIC <- 0.15
  Htc <- 0.467

  VplasC <- 0.0428
  VLC <- 0.026
  VKC <- 0.004
  VAdiC <- 0.214
  VBraC <- 0.02
  VGonC <- 0.0005
  VHeaC <- 0.005
  VLunC <- 0.008
  VMusC <- 0.4
  VSkiC <- 0.037
  VSplC <- 0.002
  VPanC <- 0.002
  VGIC <- 0.014
  VfilC <- 4e-4
  VPTCC <- 1.35e-4

  MW <- variables["MW", pfas_name]
  Free <- variables["Free", pfas_name]

  Vmax_baso_invitro <- variables["Vmax_baso_invitro", pfas_name]
  Km_baso <- variables["Km_baso", pfas_name] * variables["MW", pfas_name]
  Vmax_apical_invitro <- variables["Vmax_apical_invitro", pfas_name]
  Km_apical <- variables["Km_apical", pfas_name] * variables["MW", pfas_name]
  RAFbaso <- variables["RAFbaso", pfas_name]
  RAFapi <- variables["RAFapi", pfas_name]
  protein <- 2.0e-6
  GFRC <- 24.19 * 24

  PR <- 0.01
  PAdi <- PC["Adipose", pfas_name]
  PBra <- PC["Brain", pfas_name]
  PGon <- PC["Gonads", pfas_name]
  PGI <- PC["Gut", pfas_name]
  PHea <- PC["Heart", pfas_name]
  PL <- PC["Liver", pfas_name]
  PLun <- PC["Lung", pfas_name]
  PMus <- PC["Muscle", pfas_name]
  PSki <- PC["Skin", pfas_name]
  PSpl <- PC["Spleen", pfas_name]
  PPan <- (PGI + PSpl) / 2

  kdif <- 0.001 * 24
  kabsc <- 2.12 * 24
  kunabsc <- 7.06e-5 * 24
  keffluxc <- variables["keffluxc", pfas_name]
  kbilec <- variables["kbilec", pfas_name]
  kurinec <- 0.063 * 24
  water_consumption <- 1.36

  QC <- QCC * (BW^0.75) * (1 - Htc)
  QK <- (QKC * QC)
  QL <- (QLC * QC)
  QAdi <- (QAdiC * QC)
  QBra <- (QBraC * QC)
  QGon <- (QGonC * QC)
  QHea <- (QHeaC * QC)
  QLun <- (QLunC * QC)
  QMus <- (QMusC * QC)
  QSki <- (QSkiC * QC)
  QSpl <- (QSplC * QC)
  QPan <- (QPanC * QC)
  QGI <- (QGIC * QC)
  QR <- QC - QK - QL - QAdi - QBra - QGon - QHea - QLun - QMus - QSki - QSpl - QPan - QGI
  QBal <- QC - (QK + QL + QR + QAdi + QBra + QGon + QHea + QLun + QMus + QSki + QSpl + QPan + QGI)

  VPlas <- VplasC * BW
  VK <- VKC * BW
  MK <- VK * 1.0 * 1000
  VKb <- VK * 0.16
  Vfil <- VfilC * BW
  VL <- VLC * BW
  ML <- VL * 1.05 * 1000
  VAdi <- 17.5
  VBra <- VBraC * BW
  VGon <- VGonC * BW
  VHea <- VHeaC * BW
  VLun <- VLunC * BW
  VMus <- 32.8
  VSki <- VSkiC * BW
  VSpl <- VSplC * BW
  VPan <- VPanC * BW
  VGI <- VGIC * BW

  PTC <- VKC * 1000 * 6e7
  VPTC <- VK * 1000 * VPTCC
  MPTC <- VPTC * 1000
  VR <- (0.93 * BW) - VPlas - VPTC - Vfil - VL - VAdi - VBra - VGon - VHea - VLun - VMus - VSki - VSpl - VPan - VGI
  VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas + VAdi + VBra + VGon + VHea + VLun + VMus + VSki + VSpl + VPan + VGI)

  Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
  Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
  Vmax_baso <- Vmax_basoC * BW^0.75
  Vmax_apical <- Vmax_apicalC * BW^0.75
  kbile <- kbilec * BW^(-0.25)
  kurine <- kurinec * BW^(-0.25)
  kefflux <- keffluxc * BW^(-0.25)
  GFR <- 163.65

  kabs <- kabsc * BW^(-0.25)
  kunabs <- kunabsc * BW^(-0.25)

  return(list(
    "BW" = BW, "QBAL" = QBal, "VBAL" = VBal,
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
    "PSpl" = PSpl, "PPan" = PPan, "water_consumption" = water_consumption,
    "pfas_name" = pfas_name))
}

create.inits <- function(parameters){
  with(as.list(parameters), {
    "AR" = 0; "AAdi" = 0; "ABra" = 0; "AGon" = 0;
    "AHea" = 0; "ALun" = 0; "AMus" = 0; "ASki" = 0; "ASpl" = 0; "APan" = 0;
    "Adif" = 0; "A_baso" = 0; "AKb" = 0;
    "ACl" = 0; "Aefflux" = 0;
    "A_apical" = 0; "APTC" = 0; "Afil" = 0;
    "Aurine" = 0; "ALumen" = 0; "AGI" = 0;
    "AabsLumen" = 0; "Afeces" = 0;
    "AL" = 0; "Abile" = 0; "Aplas_free" = 0;
    "ingestion" = 0; "Cwater" = 0

    return(c("AR" = AR, "AAdi" = AAdi, "ABra" = ABra, "AGon" = AGon,
            "AHea" = AHea, "ALun" = ALun, "AMus" = AMus, "ASki" = ASki,
            "ASpl" = ASpl, "APan" = APan,
            "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
            "ACl" = ACl, "Aefflux" = Aefflux,
            "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
            "Aurine" = Aurine, "ALumen" = ALumen, "AGI" = AGI,
            "AabsLumen" = AabsLumen, "Afeces" = Afeces, 
            "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free,
            "Cwater" = Cwater, "ingestion" = ingestion))
  })
}

create.events <- function(admin_type, admin_dose_bolus, admin_time_bolus,
                          admin_dose_iv, admin_time_iv,
                          Cwater, Cwater_time, ingestion, ingestion_time) {
  if (admin_type == "iv") {
    events <- list(data = data.frame(
      var = "Aplas_free", time = admin_time_iv,
      value = admin_dose_iv, method = "add"
    ))
  } else if (admin_type == "oral") {
    events <- list(data = rbind(
      data.frame(var = "Cwater", time = Cwater_time, value = Cwater, method = "rep"),
      data.frame(var = "ingestion", time = ingestion_time, value = ingestion, method = "rep")
    ))
  } else if (admin_type == "bolus") {
    events <- list(data = data.frame(
      var = "ALumen", time = admin_time_bolus,
      value = admin_dose_bolus, method = "add"
    ))
  }
  return(events)
}

cumulative_exp_data <- function(df, time_col, concentration_col, multiply_col) {
  working_df <- df[, c(time_col, concentration_col, multiply_col)]
  complete_cases <- complete.cases(working_df)
  clean_df <- working_df[complete_cases, ]
  if (nrow(clean_df) == 0) {
    result_df <- data.frame(time = numeric(0), cumulative_mass = numeric(0))
  } else {
    multiplied_concentration <- clean_df[[concentration_col]] * clean_df[[multiply_col]]
    result_df <- data.frame(
      time = clean_df[[time_col]],
      cumulative_mass = cumsum(multiplied_concentration)
    )
  }
  return(result_df)
}

ode.func <- function(time, inits, params) {
  with(as.list(c(inits, params)), {
    CR <- AR / VR; CVR <- CR / PR
    CAdi <- AAdi / VAdi; CVAdi <- CAdi / PAdi
    CBra <- ABra / VBra; CVBra <- CBra / PBra
    CGon <- AGon / VGon; CVGon <- CGon / PGon
    CHea <- AHea / VHea; CVHea <- CHea / PHea
    CLun <- ALun / VLun; CVLun <- CLun / PLun
    CMus <- AMus / VMus; CVMus <- CMus / PMus
    CSki <- ASki / VSki; CVSki <- CSki / PSki
    CSpl <- ASpl / VSpl; CVSpl <- CSpl / PSpl
    Cpan <- APan / VPan; CVPan <- Cpan / PPan
    CGI <- AGI / VGI; CVGI <- CGI / PGI
    CKb <- AKb / VKb; CVK <- CKb
    CPTC <- APTC / VPTC; Cfil <- Afil / Vfil
    CL <- AL / VL; CLiver <- AL / ML; CVL <- CL / PL
    CA_free <- Aplas_free / VPlas; CA <- CA_free / Free

    dAR <- QR * (CA - CVR) * Free
    dAAdi <- QAdi * (CA - CVAdi) * Free
    dABra <- QBra * (CA - CVBra) * Free
    dAGon <- QGon * (CA - CVGon) * Free
    dAHea <- QHea * (CA - CVHea) * Free
    dALun <- QLun * (CA - CVLun) * Free
    dAMus <- QMus * (CA - CVMus) * Free
    dASki <- QSki * (CA - CVSki) * Free
    dASpl <- QSpl * (CA - CVSpl) * Free
    dAPan <- QPan * (CA - CVPan) * Free
    
    dAdif <- kdif * (CKb - CPTC)
    dA_baso <- (Vmax_baso * CKb) / (Km_baso + CKb)
    dAKb <- QK * (CA - CVK) * Free - CA * GFR * Free - dAdif - dA_baso
    dACl <- CA * GFR * Free
    dAefflux <- kefflux * APTC
    dA_apical <- (Vmax_apical * Cfil) / (Km_apical + Cfil)
    dAPTC <- dAdif + dA_apical + dA_baso - dAefflux
    dAfil <- CA * GFR * Free - dA_apical - Afil * kurine
    dAurine <- kurine * Afil
    dALumen <- ingestion + Cwater * water_consumption - kabs * ALumen - kunabs * ALumen
    dAGI <- kabs * ALumen + QGI * (CA - CVGI) * Free
    dAabsLumen <- kabs * ALumen
    dAfeces <- kunabs * ALumen + kbile * AL
    dAL <- QL * (CA - CVL) * Free - kbile * AL + QGI * (CVGI - CVL) * Free +
           QSpl * (CVSpl - CVL) * Free + QPan * (CVPan - CVL) * Free
    dAbile <- kbile * AL
    amount_per_gram_liver <- CLiver
    dAplas_free <- (QR * CVR * Free) + (QK * CVK * Free) + (QL * CVL * Free) +
      (QAdi * CVAdi * Free) + (QBra * CVBra * Free) + (QGon * CVGon * Free) + 
      (QHea * CVHea * Free) + (QLun * CVLun * Free) + (QMus * CVMus * Free) + 
      (QSki * CVSki * Free) + (QSpl * CVL * Free) + (QPan * CVL * Free) +
      (QGI * CVL * Free) - (QC * CA * Free) + dAefflux
    dCwater <- 0; dingestion <- 0

    Atissue <- Aplas_free + AR + AAdi + ABra + AGon + AHea + ALun + AMus + ASki + 
               ASpl + APan + AKb + Afil + ALumen + AGI + APTC + AL
    Aloss <- Aurine + Afeces
    Atotal <- Atissue + Aloss

    list(
      c("dAR" = dAR, "dAAdi" = dAAdi, "dABra" = dABra, "dAGon" = dAGon,
        "dAHea" = dAHea, "dALun" = dALun, "dAMus" = dAMus, "dASki" = dASki,
        "dASpl" = dASpl, "dAPan" = dAPan,
        "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
        "dACl" = dACl, "dAefflux" = dAefflux,
        "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
        "dAurine" = dAurine, "dALumen" = dALumen, "dAGI" = dAGI,
        "dAabsLumen" = dAabsLumen, "dAfeces" = dAfeces,
        "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
        "dCwater" = dCwater, "dingestion" = dingestion),
      "amount_per_gram_liver" = amount_per_gram_liver,
      "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal,
      "CR" = CR, "CVR" = CVR, "CAdi" = CAdi, "CVAdi" = CVAdi,
      "CBra" = CBra, "CVBra" = CVBra, "CGon" = CGon, "CVGon" = CVGon,
      "CHea" = CHea, "CVHea" = CVHea, "CLun" = CLun, "CVLun" = CVLun,
      "CMus" = CMus, "CVMus" = CVMus, "CSki" = CSki, "CVSki" = CVSki,
      "CSpl" = CSpl, "CVSpl" = CVSpl, "Cpan" = Cpan, "CVPan" = CVPan,
      "CKb" = CKb, "CGI" = CGI, "CVGI" = CVGI, "CVK" = CVK, 
      "CPTC" = CPTC, "Cfil" = Cfil, "CL" = CL, "CVL" = CVL,
      "CA_free" = CA_free, "CA" = CA
    )
  })
}

# ================================================================================
# 6. SIMULATION LOOP FOR ALL PFAS COMPOUNDS
# ================================================================================

all_results <- list()
all_summaries <- data.frame()
all_pk_parameters <- data.frame()

for (pfas in pfas_compounds) {
  cat("\n========================================\n")
  cat(sprintf("Running simulation for: %s\n", pfas))
  cat("========================================\n")
  
  admin_type <- "bolus"
  admin_dose_bolus <- c(variables["admin_dose_bolus", pfas])
  admin_time_bolus <- c(0)
  admin_dose_iv <- 0
  admin_time_iv <- 0
  Cwater <- 0.00
  Cwater_time <- 0
  ingestion <- 0
  ingestion_time <- c(0)
  BW <- 82

  simulation_time <- 450
  time_step <- 0.1

  cat("Setting up simulation...\n")

  params <- create.params(variables, PC, BW, pfas)
  inits <- create.inits(params)
  events <- create.events(admin_type, admin_dose_bolus, admin_time_bolus,
                          admin_dose_iv, admin_time_iv,
                          Cwater, Cwater_time, ingestion, ingestion_time)

  cat("Solving ODEs...\n")

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

  plasma_exp <- read.csv("C:\\Users\\fotis\\Documents\\GitHub\\PFAS_PBK_models\\Extended PFAS PBK model\\exp_data_plasma.csv") %>%
    select("time", all_of(pfas)) %>%
    filter(!is.na(!!sym(pfas)))

  cat("Experimental data loaded!\n")

cat("Available columns:", paste(colnames(plasma_exp), collapse=", "), "\n")
cat(pfas, "in columns?", pfas %in% colnames(plasma_exp), "\n")

if (pfas %in% colnames(plasma_exp)) {
  plasma_exp <- plasma_exp %>% select("time", all_of(pfas))
  cat("Total rows:", nrow(plasma_exp), "\n")
  cat("Non-NA values:", sum(!is.na(plasma_exp[[pfas]])), "\n")
  cat("Positive values:", sum(plasma_exp[[pfas]] > 0, na.rm=TRUE), "\n")
  cat("Time range:", range(plasma_exp$time, na.rm=TRUE), "\n")
  cat("Conc range:", range(plasma_exp[[pfas]], na.rm=TRUE), "\n")
} else {
  cat("ERROR: Column", pfas, "NOT FOUND!\n")
}  


cat("\n=== Checking observed data for", pfas, "===\n")
cat("Number of rows in plasma_exp:", nrow(plasma_exp), "\n")
cat("Columns in plasma_exp:", paste(colnames(plasma_exp), collapse=", "), "\n")
if (pfas %in% colnames(plasma_exp)) {
  cat("Non-NA concentrations:", sum(!is.na(plasma_exp[[pfas]])), "\n")
  cat("Positive concentrations:", sum(plasma_exp[[pfas]] > 0, na.rm=TRUE), "\n")
  cat("Concentration range:", range(plasma_exp[[pfas]], na.rm=TRUE), "\n")
} else {
  cat("WARNING: Column '", pfas, "' not found in plasma_exp!\n")
}

  # ================================================================================
  # 8. PK PARAMETER CALCULATIONS
  # ================================================================================
  cat("Calculating PK parameters...\n")

  # Model predicted plasma PK parameters
  model_pk <- calculate_pk_parameters(
    time = solution$time,
    conc = solution$CA,
    route = "Oral",
    plot = TRUE,
    min_points = 5,
    output_dir = output_dir,
    pfas_name = pfas,
    data_type = "Model"
  )

  # Observed plasma PK parameters
  if (nrow(plasma_exp) > 0 && pfas %in% colnames(plasma_exp)) {
    obs_pk <- calculate_pk_parameters(
      time = plasma_exp$time,
      conc = plasma_exp[[pfas]],
      route = "Oral",
      plot = TRUE,
      min_points = 3,
      output_dir = output_dir,
      pfas_name = pfas,
      data_type = "Observed"
    )
  } else {
    obs_pk <- list(
      half_life = NA, lambda_z = NA, cmax = NA, tmax = NA,
      auc_0_last = NA, auc_0_inf = NA, reason = "No observed data",
      data_type = "Observed"
    )
  }

  cat("PK parameter calculations complete!\n")



  # ================================================================================
  # 10. SUMMARY OUTPUT
  # ================================================================================
  cat("\n=== SIMULATION SUMMARY ===\n")
  cat("PFAS Compound:", pfas, "\n")
  cat("Body Weight:", params$BW, "kg\n")
  cat("Administration Type:", admin_type, "\n")
  cat("Bolus Dose:", admin_dose_bolus, "ug at time", admin_time_bolus, "days\n")
  cat("Simulation Time:", simulation_time, "days\n")
  cat("\n=== PK PARAMETERS ===\n")
  cat("Model Cmax:", round(model_pk$cmax, 4), "ug/L\n")
  cat("Model Tmax:", round(model_pk$tmax, 2), "days\n")
  cat("Model AUC(0-last):", round(model_pk$auc_0_last, 4), "ug·day/L\n")
  cat("Model AUC(0-inf):", round(model_pk$auc_0_inf, 4), "ug·day/L\n")
  cat("Model Half-Life:", round(model_pk$half_life, 2), "days\n")
  cat("\nObserved Cmax:", round(obs_pk$cmax, 4), "ug/L\n")
  cat("Observed Tmax:", round(obs_pk$tmax, 2), "days\n")
  cat("Observed AUC(0-last):", round(obs_pk$auc_0_last, 4), "ug·day/L\n")
  cat("Observed AUC(0-inf):", round(obs_pk$auc_0_inf, 4), "ug·day/L\n")
  cat("Observed Half-Life:", round(obs_pk$half_life, 2), "days\n")
  cat("==========================\n")

  all_results[[pfas]] <- solution
  
  summary_row <- data.frame(
    PFAS = pfas,
    BW = params$BW,
    Admin_Type = admin_type,
    Dose = admin_dose_bolus,
    Sim_Time = simulation_time,
    Final_Plasma = round(tail(solution$Aplas_free, 1), 4),
    Final_Liver = round(tail(solution$AL, 1), 4),
    Final_Rest = round(tail(solution$AR, 1), 4),
    Urine_Excreted = round(tail(solution$Aurine, 1), 4),
    Feces_Excreted = round(tail(solution$Afeces, 1), 4)
  )
  all_summaries <- rbind(all_summaries, summary_row)

  pk_row <- data.frame(
    PFAS = pfas,
    # Model PK Parameters
    Model_Cmax_ug_L = round(model_pk$cmax, 4),
    Model_Tmax_days = round(model_pk$tmax, 2),
    Model_AUC_0_last = round(model_pk$auc_0_last, 4),
    Model_AUC_0_inf = round(model_pk$auc_0_inf, 4),
    Model_Half_Life_days = round(model_pk$half_life, 2),
    Model_Lambda_z = round(model_pk$lambda_z, 4),
    Model_Points_Used = model_pk$points_used,
    Model_Adj_R2 = round(model_pk$adj_r_squared, 3),
    # Observed PK Parameters
    Observed_Cmax_ug_L = round(obs_pk$cmax, 4),
    Observed_Tmax_days = round(obs_pk$tmax, 2),
    Observed_AUC_0_last = round(obs_pk$auc_0_last, 4),
    Observed_AUC_0_inf = round(obs_pk$auc_0_inf, 4),
    Observed_Half_Life_days = round(obs_pk$half_life, 2),
    Observed_Lambda_z = round(obs_pk$lambda_z, 4),
    Observed_Points_Used = obs_pk$points_used,
    Observed_Adj_R2 = round(obs_pk$adj_r_squared, 3),
    # Ratios (Model/Observed)
    Cmax_Ratio = round(model_pk$cmax / obs_pk$cmax, 2),
    AUC_Ratio = round(model_pk$auc_0_last / obs_pk$auc_0_last, 2),
    Half_Life_Ratio = round(model_pk$half_life / obs_pk$half_life, 2)
  )
  all_pk_parameters <- rbind(all_pk_parameters, pk_row)
}

write.csv(all_summaries, "PFAS_simulation_summary.csv", row.names = FALSE)
write.csv(all_pk_parameters, "PFAS_PK_parameters_comparison.csv", row.names = FALSE)

cat("\n\n========================================\n")
cat("ALL SIMULATIONS COMPLETE!\n")
cat("========================================\n")
cat("Summary saved to: PFAS_simulation_summary.csv\n")
cat("PK parameters saved to: PFAS_PK_parameters_comparison.csv\n")
cat("========================================\n")

print(all_pk_parameters)
