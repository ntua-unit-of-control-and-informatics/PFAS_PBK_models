library(deSolve)
library(tidyverse)
library(ggplot2)

#=========================
# 1. Parameters of the model
#=========================

variables<-read.csv("estimated_parameters.csv",row.names="Parameters")
PC<-read.csv("PCs.csv",row.names="Organs")

# NOTE: Manually adjust parameter values below as needed
create.params <- function(user_input, pfas_name) {
  with(as.list(user_input), {
  # === USER-ADJUSTABLE PARAMETERS ===
  BW <- 82 
  QCC <- 12.5 * 24
  QLC <- 0.065
  QKC <- 0.175
  QAdiC <-0.05
  QBraC <-0.12
  QGonC <-0.0005
  QHeaC <-0.04
  QLunC <-0.025
  QMusC <-0.17
  QSkiC <-0.05
  QSplC <-0.03
  QPanC <-0.01
  QGIC<-0.15
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
  protein <- 2.0e-6
  GFRC <- 24.19 * 24
  water_consumption <- 1.36

  # === SCALED PARAMETERS ===
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
  VAdi <-17.5
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
  VR <- (0.93 * BW) - VPlas- VPTC - Vfil - VL -VAdi - VBra - VGon - VHea - VLun - VMus - VSki - VSpl - VPan - VGI
  VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas + VAdi + VBra + VGon + VHea - VLun + VMus + VSki + VSpl + VPan + VGI)

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
    CR <- AR / VR
    CVR <- CR / PR
    CAdi <- AAdi / VAdi
    CVAdi <- CAdi / PAdi
    CBra <- ABra / VBra
    CVBra <- CBra / PBra
    CGon <- AGon / VGon
    CVGon <- CGon / PGon
    CHea <- AHea / VHea
    CVHea <- CHea / PHea
    CLun <- ALun / VLun
    CVLun <- CLun / PLun
    CMus <- AMus / VMus
    CVMus <- CMus / PMus
    CSki <- ASki / VSki
    CVSki <- CSki / PSki
    CSpl <- ASpl / VSpl
    CVSpl <- CSpl / PSpl
    Cpan <- APan / VPan
    CVPan <- Cpan / PPan
    CGI <- AGI / VGI
    CVGI <- CGI / PGI
    CKb <- AKb / VKb
    CVK <- CKb
    CPTC <- APTC / VPTC
    Cfil <- Afil / Vfil
    CL <- AL / VL
    CLiver <- AL / ML
    CVL <- CL / PL
    CA_free <- Aplas_free / VPlas
    CA <- CA_free / Free
    
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
    dALumen <- ingestion + Cwater * water_consumption - kabs * ALumen- kunabs * ALumen
    dAGI <- kabs * ALumen + QGI*(CA - CVGI) * Free
    dAabsLumen <- kabs * ALumen
    dAfeces <-  kunabs * ALumen + kbile * AL
    dAL <- QL * (CA - CVL) * Free - kbile * AL + QGI*(CVGI - CVL) * Free +
           QSpl*(CVSpl-CVL)*Free +QPan*(CVPan-CVL)*Free
    dAbile <- kbile * AL
    amount_per_gram_liver <- CLiver
    
    dAplas_free <- (QR * CVR * Free) + (QK * CVK * Free) + (QL * CVL * Free) +
    (QAdi*CVAdi*Free) + (QBra*CVBra*Free) + (QGon*CVGon*Free) + (QHea*CVHea*Free) +
    (QLun*CVLun*Free) +(QMus*CVMus*Free) +(QSki*CVSki*Free)+(QSpl*CVL*Free)+(QPan*CVL*Free)+
    (QGI*CVL*Free)- (QC * CA * Free) + dAefflux

    dCwater <- 0
    dingestion <- 0

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

admin_type <- "bolus"
admin_dose_bolus <- c(variables["admin_dose_bolus",pfas])
admin_time_bolus <- c(0)
admin_dose_iv <- 0
admin_time_iv <- 0
Cwater <- 0.00
Cwater_time <- 0
ingestion <- 0
ingestion_time <- c(0)

kdif = 0.001*24
kabsc = 2.12*24
kunabsc = 7.06e-5 * 24

thetas<-list(
"RAFbaso", "RAFapi", "PR", "PAdi", "PBra", "PGon", "PGI", "PHea", "PL", 
"PLun", "PMus", "PSki", "PSpl", "PPan", "kdif", "kabsc", "kunabsc",
"keffluxc", "kbilec", "kurinec", "Vmax_baso_invitro", "Km_baso", 
"Vmax_apical_invitro", "Km_apical", "Free")

# ================================================================================
# 7. PFAS LIST
# ================================================================================

pfas_list <- c("PFHpA", "PFOA", "PFNA", "PFDA", "PFBS", "PFHxS", "PFOS", 
               "DONA", "HFPO_DA", "PFBA", "PFHxA")

# ================================================================================
# 8. OBJECTIVE FUNCTION (modified to accept PFAS name)
# ================================================================================

obj.func <- function(x, dp, pfas_name, user_input){
  if (is.na(x)) {
    params <- create.params(user_input, pfas_name)
  } else {
    pert_input <- user_input
    pert_input[[x]] <- as.numeric(user_input[[x]]) * (1 + dp)
    params <- create.params(pert_input, pfas_name)
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
    #C_6h = solution[which.min(abs(sample_time - 0.25)), "CA"],
    #C_6d = solution[which.min(abs(sample_time - 6)), "CA"],
    Feces_Excreted = tail(solution[, "Afeces"], 1),
    Urine_Excreted = tail(solution[, "Aurine"], 1)
  )
  
  return(result)
}

# ================================================================================
# 9. SENSITIVITY ANALYSIS FOR ALL PFAS
# ================================================================================

dp <- 0.5

# Initialize storage for all results
all_sensitivity_df <- tibble(
  PFAS = character(),
  Parameter = character(),
  Direction = factor(levels = c("Negative", "Positive")),
  Output = character(),
  RelativeChange = numeric()
)

for (pfas in pfas_list) {
  cat("\n========== Processing:", pfas, "==========\n")
  
  # Set PFAS-specific parameters
  PR = 0.01
  PAdi = PC["Adipose", pfas]
  PBra = PC["Brain", pfas]
  PGon = PC["Gonads", pfas]
  PGI = PC["Gut", pfas]
  PHea = PC["Heart", pfas]
  PL = PC["Liver", pfas]
  PLun = PC["Lung", pfas]
  PMus = PC["Muscle", pfas]
  PSki = PC["Skin", pfas]
  PSpl = PC["Spleen", pfas]
  PPan = (PGI + PSpl) / 2
  MW = variables["MW", pfas]
  Free = variables["Free", pfas]
  kbilec <- variables["kbilec", pfas]
  kurinec <- 0.063 * 24
  Vmax_baso_invitro = variables["Vmax_baso_invitro", pfas]
  Km_baso = variables["Km_baso", pfas]
  Vmax_apical_invitro = variables["Vmax_apical_invitro", pfas]
  Km_apical = variables["Km_apical", pfas]
  keffluxc = variables["keffluxc", pfas]
  RAFbaso = variables["RAFbaso", pfas]
  RAFapi = variables["RAFapi", pfas]
  
  user_input <- list(
    "Free" = Free, "RAFbaso" = RAFbaso, "RAFapi" = RAFapi,
    "PR" = PR, "PAdi" = PAdi, "PBra" = PBra, "PGon" = PGon,
    "PGI" = PGI, "PHea" = PHea, "PL" = PL, "PLun" = PLun,
    "PMus" = PMus, "PSki" = PSki, "PSpl" = PSpl, "PPan" = PPan,
    "MW" = MW, "keffluxc" = keffluxc, "kdif" = kdif,
    "kabsc" = kabsc, "kunabsc" = kunabsc, "kbilec" = kbilec,
    "kurinec" = kurinec, "Vmax_baso_invitro" = Vmax_baso_invitro,
    "Km_baso" = Km_baso, "Vmax_apical_invitro" = Vmax_apical_invitro,
    "Km_apical" = Km_apical, "admin_type" = admin_type,
    "admin_dose_bolus" = admin_dose_bolus, "admin_time_bolus" = admin_time_bolus,
    "admin_dose_iv" = admin_dose_iv, "admin_time_iv" = admin_time_iv,
    "Cwater" = Cwater, "Cwater_time" = Cwater_time,
    "ingestion" = ingestion, "ingestion_time" = ingestion_time
  )
  
  # Get baseline
  baseline <- obj.func(NA, dp = 0, pfas_name = pfas, user_input = user_input)
  
  # Loop over each parameter
  for (par in thetas) {
    cat("  Parameter:", par, "\n")
    
    # Positive perturbation
    res_pos <- obj.func(par, dp = dp, pfas_name = pfas, user_input = user_input)
    rel_pos <- ((res_pos - baseline) / baseline) / dp
    
    # Negative perturbation
    res_neg <- obj.func(par, dp = -dp, pfas_name = pfas, user_input = user_input)
    rel_neg <- ((res_neg - baseline) / baseline) / (dp)
    
    n_outputs <- length(baseline)
    output_names <- names(baseline)
    
    all_sensitivity_df <- bind_rows(all_sensitivity_df,
      tibble(
        PFAS = rep(pfas, n_outputs),
        Parameter = rep(par, n_outputs),
        Direction = "Positive",
        Output = output_names,
        RelativeChange = rel_pos
      ),
      tibble(
        PFAS = rep(pfas, n_outputs),
        Parameter = rep(par, n_outputs),
        Direction = "Negative",
        Output = output_names,
        RelativeChange = rel_neg
      )
    )
  }
}
# ================================================================================
# 10. CREATE HEATMAPS FOR EACH OUTPUT VARIABLE
# ================================================================================

# Set PFAS and Parameter as factors for consistent ordering
all_sensitivity_df$PFAS <- factor(all_sensitivity_df$PFAS, levels = pfas_list)
all_sensitivity_df$Parameter <- factor(all_sensitivity_df$Parameter, levels = thetas)

output_vars <- c(#"C_6h", "C_6d",
 "Cmax", "AUC", "Feces_Excreted", "Urine_Excreted")
output_dir <- "Extended PFAS PBK model/PFAS_sensitivity_heatmaps"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

scale_limits <- c(min(all_sensitivity_df$RelativeChange, na.rm = TRUE),
                  max(all_sensitivity_df$RelativeChange, na.rm = TRUE))

# Elsevier & LaTeX optimized theme (VERTICAL X-AXIS LABELS)
elsevier_theme <- theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 1, face = "bold", family = "Helvetica"),
    axis.title = element_text(size = 12, family = "Helvetica"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 11, family = "Helvetica", color = "black"),
    axis.text.y = element_text(size = 11, family = "Helvetica", color = "black"),
    axis.ticks = element_line(color = "black", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    legend.title = element_text(size = 11, family = "Helvetica"),
    legend.text = element_text(size = 10, family = "Helvetica"),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  )

fill_scale <- scale_fill_gradient2(
  low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
  name = "Sensitivity Index", limits = scale_limits,
  guide = guide_colorbar(barwidth = 0.8, barheight = 6)
)
# %.0f ensures it prints as a whole number without decimals
pert_pct <- sprintf("%.0f", dp * 100)

for (output_var in output_vars) {
  cat("\nCreating heatmap for:", output_var, "\n")
  
  positive_df <- all_sensitivity_df %>% filter(Direction == "Positive", Output == output_var)
  negative_df <- all_sensitivity_df %>% filter(Direction == "Negative", Output == output_var)
  
  # =========================
  # POSITIVE PERTURBATION HEATMAP
  # =========================
  p_positive <- ggplot(positive_df, aes(x = PFAS, y = Parameter, fill = RelativeChange)) +
    geom_tile(color = "white", size = 0.3) + 
    fill_scale +
    labs(
      title = paste0("Sensitivity Index - ", output_var, "\n(Positive Perturbation +", pert_pct, "%)"),
      x = "PFAS Compounds", 
      y = "Parameters"
    ) + 
    elsevier_theme
  
  ggsave(filename = file.path(output_dir, paste0("heatmap_", output_var, "_positive.png")),
         plot = p_positive, width = 7, height = 5.5, dpi = 600)
  
  # =========================
  # NEGATIVE PERTURBATION HEATMAP
  # =========================
  p_negative <- ggplot(negative_df, aes(x = PFAS, y = Parameter, fill = RelativeChange)) +
    geom_tile(color = "white", size = 0.3) + 
    fill_scale +
    labs(
      # Fixed syntax: properly concatenates the negative sign and percentage
      title = paste0("Sensitivity Index - ", output_var, "\n(Negative Perturbation -", pert_pct, "%)"),
      x = "PFAS Compounds", 
      y = "Parameters"
    ) + 
    elsevier_theme
  
  ggsave(filename = file.path(output_dir, paste0("heatmap_", output_var, "_negative.png")),
         plot = p_negative, width = 7, height = 5.5, dpi = 600)
  
  cat("  Saved: heatmap_", output_var, "_positive.png & _negative.png\n")
}
# ================================================================================
# OPTIONAL: COMBINED HEATMAP (2-COLUMN + BOTTOM LEGEND + SPACED Y-AXIS)
# ================================================================================
library(patchwork)

plot_list <- list()

# Theme for individual subplots with SPACED y-axis labels
combined_theme <- theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Helvetica"),
    axis.title = element_text(size = 12, family = "Helvetica"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, 
                              family = "Helvetica", color = "black"),
    # === KEY FIX: Space out y-axis parameter labels ===
    axis.text.y = element_text(
      size = 14, 
      family = "Helvetica", 
      color = "black",
      lineheight = 2,              # ← Increase vertical spacing between labels
      margin = margin(t = 6, b = 6)  # ← Add small top/bottom padding per label
    ),
    axis.ticks = element_line(color = "black", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    legend.position = "right",  # Keep as "right" for patchwork collection
    legend.key.height = unit(1, "cm"),
    plot.margin = margin(8, 8, 8, 8, "pt")
  )

# Shared fill scale with HORIZONTAL colorbar
shared_fill_scale <- scale_fill_gradient2(
  low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
  name = "Sensitivity Index", limits = scale_limits,
  guide = guide_colorbar(
    barwidth = 6,
    barheight = 0.8,
    title.position = "top",
    title.theme = element_text(size = 11, family = "Helvetica", face = "bold"),
    label.theme = element_text(size = 10, family = "Helvetica"),
    direction = "horizontal"
  )
)

# Build plot list
for (output_var in output_vars) {
  positive_df <- all_sensitivity_df %>% 
    filter(Direction == "Positive", Output == output_var)
  
  plot_list[[output_var]] <- ggplot(positive_df, aes(x = PFAS, y = Parameter, fill = RelativeChange)) +
    geom_tile(color = "white", size = 0.3) + 
    shared_fill_scale +
    labs(title = output_var, x = "", y = "") + 
    combined_theme
}

# Assemble: 2 columns, collect legends, position at bottom
combined_plot <- patchwork::wrap_plots(
  plot_list, 
  ncol = 2,
  plot_spacing = unit(0.6, "cm")
) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "PFAS Sensitivity Analysis - All Output Variables",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Helvetica"),
      plot.margin = margin(10, 15, 20, 15, "pt")  # Extra bottom margin for legend
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, family = "Helvetica", face = "bold"),
    legend.text = element_text(size = 11, family = "Helvetica"),
    legend.key.width = unit(1.5, "cm"),
    legend.margin = margin(10, 0, 0, 0)
  )

# === KEY FIX: Increase height significantly to accommodate spaced y-labels ===
ggsave(filename = file.path(output_dir, "heatmap_all_outputs_combined.png"),
       plot = combined_plot, 
       width = 13,    
       height = 16,   # ← Increased from 12 to 16 for taller rows + spaced labels
       dpi = 600)
cat("Saved: heatmap_all_outputs_combined.png\n")

output_dir <- "Extended PFAS PBK model/PFAS_sensitivity_heatmaps"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pert_pct <- sprintf("%.0f", dp * 100)

for (output_var in output_vars) {
  cat("\nCreating Elsevier-ready heatmap for:", output_var, "\n")
  
  positive_df <- all_sensitivity_df %>% 
    filter(Direction == "Positive", Output == output_var)
  negative_df <- all_sensitivity_df %>% 
    filter(Direction == "Negative", Output == output_var)
  
  # === POSITIVE PERTURBATION ===
  p_positive <- ggplot(positive_df, aes(x = PFAS, y = Parameter, fill = RelativeChange)) +
    geom_tile(color = "white", size = 0.25) + 
    fill_scale_elsevier +
    labs(
      title = paste0(output_var, " (+", pert_pct, "%)"),
      x = "PFAS Compound", 
      y = "Model Parameter"
    ) + 
    elsevier_latex_theme
  
  # Export as vector PDF (single-column)
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_", output_var, "_positive.pdf")),
    plot = p_positive,
    width = 3.4, height = 4.2, units = "in",
    device = cairo_pdf,  # Ensures font embedding & transparency support
    dpi = 300  # fallback for raster elements
  )
  
  # === NEGATIVE PERTURBATION ===
  p_negative <- ggplot(negative_df, aes(x = PFAS, y = Parameter, fill = RelativeChange)) +
    geom_tile(color = "white", size = 0.25) + 
    fill_scale_elsevier +
    labs(
      title = paste0(output_var, " (−", pert_pct, "%)"),  # Use proper minus sign (U+2212)
      x = "PFAS Compound", 
      y = "Model Parameter"
    ) + 
    elsevier_latex_theme
  
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_", output_var, "_negative.pdf")),
    plot = p_negative,
    width = 3.4, height = 4.2, units = "in",
    device = cairo_pdf,
    dpi = 300
  )
  
  cat("  ✓ Saved PDFs for", output_var, "\n")
}

write.csv2(all_sensitivity_df, file.path(output_dir, "PFAS_sensitivity_results.csv"), row.names = FALSE)

# ================================================================================
# REQUIRED PACKAGES
# ================================================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(tikzDevice)

# ================================================================================
# SETUP (Overleaf-friendly: NO SPACES IN PATHS)
# ================================================================================
output_dir <- "PFAS_Sensitivity_Heatmaps"  # ← Overleaf/LaTeX breaks on spaces
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pert_pct <- sprintf("%.0f", dp * 100)

# ================================================================================
# THEME & SCALE DEFINITIONS
# ================================================================================
combined_theme <- theme_bw(base_family = "Helvetica") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "Helvetica"),
    axis.title = element_text(size = 12, family = "Helvetica"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 14, 
                              family = "Helvetica", color = "black"),
    axis.text.y = element_text(
      size = 14, family = "Helvetica", color = "black",
      lineheight = 2, margin = margin(t = 6, b = 6)
    ),
    axis.ticks = element_line(color = "black", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1, "cm"),
    plot.margin = margin(8, 8, 8, 8, "pt")
  )

shared_fill_scale <- scale_fill_gradient2(
  low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
  name = "Sensitivity Index", limits = scale_limits,
  guide = guide_colorbar(
    barwidth = 6, barheight = 0.8,
    title.position = "top",
    title.theme = element_text(size = 11, family = "Helvetica", face = "bold"),
    label.theme = element_text(size = 10, family = "Helvetica"),
    direction = "horizontal"
  )
)

# ================================================================================
# BUILD CUMULATIVE PLOT
# ================================================================================
plot_list <- list()

for (output_var in output_vars) {
  positive_df <- all_sensitivity_df %>% 
    filter(Direction == "Positive", Output == output_var)
  
  plot_list[[output_var]] <- ggplot(positive_df, aes(x = PFAS, y = Parameter, fill = RelativeChange)) +
    geom_tile(color = "white", size = 0.3) + 
    shared_fill_scale +
    labs(title = output_var, x = "", y = "") + 
    combined_theme
}

combined_plot <- patchwork::wrap_plots(
  plot_list, 
  ncol = 2,
  plot_spacing = unit(0.6, "cm")
) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "PFAS Sensitivity Analysis - All Output Variables",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Helvetica"),
      plot.margin = margin(10, 15, 30, 15, "pt")
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, family = "Helvetica", face = "bold"),
    legend.text = element_text(size = 11, family = "Helvetica"),
    legend.key.width = unit(1.5, "cm"),
    legend.margin = margin(10, 0, 0, 0)
  )

# ================================================================================
# OVERLEAF-READY TIKZ EXPORT
# ================================================================================
tex_filename <- file.path(output_dir, "heatmap_cumulative.tex")

tikzDevice::tikz(
  filename = tex_filename,
  width = 7,          # inches (adjust to match your column/text width)
  height = 9,         # inches
  standAlone = FALSE, # Ready for \input{} in Overleaf
  sanitize = TRUE,    # Escapes LaTeX special characters automatically
  encoding = "utf-8"  # Ensures Unicode compatibility in Overleaf
)

print(combined_plot)
dev.off()

cat("\n✓ Overleaf-ready TikZ file exported:", tex_filename, "\n")