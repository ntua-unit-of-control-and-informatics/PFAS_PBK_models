# ================================================================================
# PFAS PBK Model - Manual Simulation Script
# Based on Loccisano et al. (2012) PBK model
# ================================================================================
# This script simulates PFAS distribution with manually adjustable parameters
# All parameters are defined in create.params() function with default values

library(deSolve)
library(tidyverse)

#=========================
# 1. Parameters of the model
#=========================
# NOTE: Manually adjust parameter values below as needed
create.params <- function(BW) {

  # === USER-ADJUSTABLE PARAMETERS ===
  # Change these default values as needed

  # --- Cardiac Output and Blood Flow (as fraction of cardiac output) ---
  QCC <- 12.5 * 24  # cardiac output in L/day/kg^0.75; Brown 1997
  QLC <- 0.25  # fraction blood flow to liver; Brown 1997
  QKC <- 0.175  # fraction blood flow to kidney; Brown 1997
  Htc <- 0.467  # hematocrit for the rat; Davies 1993

  # --- Tissue Volumes ---
  VplasC <- 0.0428  # fraction vol. of plasma (L/kg BW); Davies 1993
  VLC <- 0.026  # fraction vol. of liver (L/kg BW); Brown 1997
  VKC <- 0.004  # fraction vol. of kidney (L/kg BW); Brown 1997
  VfilC <- 4e-4  # fraction vol. of filtrate (L/kg BW)
  VPTCC <- 1.35e-4  # vol. of proximal tubule cells (L/g kidney)

  # --- Chemical Specific Parameters ---
  MW <- 414.07  # PFOA molecular mass (g/mol)
  Free <- 0.001  # free fraction in plasma (Smeltz 2023)

  # --- Kidney Transport Parameters ---
  Vmax_baso_invitro <- 439.2  # Vmax of basolateral transporter (pmol/mg protein/min)
  Km_baso <- 20100  # Km of basolateral transporter (ug/L)
  Vmax_apical_invitro <- 37400  # Vmax of apical transporter (pmol/mg protein/min)
  Km_apical <- 77500  # Km of apical transporter (ug/L)
  RAFbaso <- 1  # relative activity factor, basolateral transporters (male)
  RAFapi <- 0.0007  # relative activity factor, apical transporters (male)
  protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
  GFRC <- 24.19 * 24  # glomerular filtration rate (L/day/kg kidney); Corley 2005

  # --- Partition Coefficients (from Allendorf 2021) ---
  PL <- 0.434698291544763  # liver:blood
  PK <- 0.413283707125888  # kidney:blood
  PR <- 0.08#0.534206103089267  # rest of body:blood

  # --- Rate Constants ---
  kdif <- 0.001 * 24  # diffusion rate from proximal tubule cells (L/day)
  kabsc <- 2.12 * 24  # rate of absorption from small intestine (1/(day*BW^-0.25))
  kunabsc <- 7.06e-5 * 24  # rate of unabsorbed dose to feces (1/(day*BW^-0.25))
  GEC <- 3.5 * 24  # gastric emptying time (1/(day*BW^-0.25)); Yang 2013
  k0C <- 1.0 * 24  # rate of uptake from stomach to liver (1/(day*BW^-0.25))
  keffluxc <- 0.1 * 24  # rate of efflux from PTC to blood (1/(day*BW^-0.25))
  kbilec <- 0.0001 * 24  # biliary elimination rate (1/(day*BW^-0.25))
  kurinec <- 0.063 * 24  # urinary elimination rate (1/(day*BW^-0.25))
  kvoid <- 0.06974 * 24  # daily urine volume rate (L/day); Van Haarst 2004

  # --- Water Consumption ---
  water_consumption <- 1.36  # L/day

  # === SCALED PARAMETERS (calculated from above) ===

  # Cardiac output and blood flows
  QC <- QCC * (BW^0.75) * (1 - Htc)  # cardiac output in L/day; adjusted for plasma
  QK <- (QKC * QC)  # plasma flow to kidney (L/day)
  QL <- (QLC * QC)  # plasma flow to liver (L/day)
  QR <- QC - QK - QL  # plasma flow to rest of body (L/day)
  QBal <- QC - (QK + QL + QR)  # Balance check; should equal zero

  # Tissue Volumes
  VPlas <- VplasC * BW  # volume of plasma (L)
  VK <- VKC * BW  # volume of kidney (L)
  MK <- VK * 1.0 * 1000  # mass of the kidney (g)
  VKb <- VK * 0.16  # volume of blood in the kidney (L); Brown 1997
  Vfil <- VfilC * BW  # volume of filtrate (L)
  VL <- VLC * BW  # volume of liver (L)
  ML <- VL * 1.05 * 1000  # mass of the liver (g)

  # Kidney Parameters
  PTC <- VKC * 1000 * 6e7  # number of PTC (cells/kg BW)
  VPTC <- VK * 1000 * VPTCC  # volume of proximal tubule cells (L)
  MPTC <- VPTC * 1000  # mass of the proximal tubule cells (g)
  VR <- (0.93 * BW) - VPlas - VPTC - Vfil - VL  # volume of remaining tissue (L)
  VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas)  # Balance check; should equal zero

  Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
  Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24
  Vmax_baso <- Vmax_basoC * BW^0.75  # (ug/day)
  Vmax_apical <- Vmax_apicalC * BW^0.75  # (ug/day)
  kbile <- kbilec * BW^(-0.25)  # biliary elimination; liver to feces storage (/day)
  kurine <- kurinec * BW^(-0.25)  # urinary elimination, from filtrate (/day)
  kefflux <- keffluxc * BW^(-0.25)  # efflux clearance rate, from PTC to blood (/day)
  GFR <- 163.65  # glomerular filtration rate, (L/day)

  # GI Tract Parameters
  kabs <- kabsc * BW^(-0.25)  # rate of absorption from small intestine (/day)
  kunabs <- kunabsc * BW^(-0.25)  # rate of unabsorbed dose to feces (/day)
  GE <- GEC * BW^(-0.25)  # gastric emptying time (/day)
  k0 <- k0C * BW^(-0.25)  # rate of uptake from stomach to liver (/day)

  return(list(
    "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR,
    "VPlas" = VPlas, "VKb" = VKb, "Vfil" = Vfil, "VL" = VL, "VR" = VR, "ML" = ML,
    "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
    "kdif" = kdif, "Km_baso" = Km_baso, "Km_apical" = Km_apical,
    "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
    "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs, "GE" = GE, "k0" = k0,
    "PL" = PL, "PK" = PK, "PR" = PR, "kvoid" = kvoid,
    "water_consumption" = water_consumption
  ))
}

#===============================================
# 2. Function to create initial values for ODEs
#===============================================
create.inits <- function(parameters) {
  with(as.list(parameters), {
    return(c(
      "AR" = 0, "Adif" = 0, "A_baso" = 0, "AKb" = 0,
      "ACl" = 0, "Aefflux" = 0, "A_apical" = 0, "APTC" = 0, "Afil" = 0,
      "Aurine" = 0, "AST" = 0, "AabsST" = 0, "ASI" = 0, "AabsSI" = 0,
      "Afeces" = 0, "AL" = 0, "Abile" = 0, "Aplas_free" = 0,
      "ingestion" = 0, "Cwater" = 0
    ))
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
      var = "AST", time = admin_time_bolus,
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

#==============
# 5. ODEs System
#==============
ode.func <- function(time, inits, params) {
  with(as.list(c(inits, params)), {

    CR <- AR / VR  # concentration in rest of body (ug/L)
    CVR <- CR / PR  # concentration in venous blood leaving rest of body (ug/L)
    CKb <- AKb / VKb  # concentration in kidney blood (ug/L)
    CVK <- CKb  # concentration in venous blood leaving kidney (ug/L)
    CPTC <- APTC / VPTC  # concentration in PTC (ug/L)
    Cfil <- Afil / Vfil  # concentration in filtrate (ug/L)
    CL <- AL / VL  # concentration in the liver (ug/L)
    CLiver <- AL / ML  # concentration in the liver (ug/g)
    CVL <- CL / PL  # concentration in the venous blood leaving the liver (ug/L)
    CA_free <- Aplas_free / VPlas  # free concentration in plasma (ug/L)
    CA <- CA_free / Free  # concentration of total PFOA in plasma (ug/L)
    Curine <- Aurine / kvoid

    # Rest of Body (Tis)
    dAR <- QR * (CA - CVR) * Free

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
    # Stomach
    dAST <- ingestion + Cwater * water_consumption - k0 * AST - GE * AST
    dAabsST <- k0 * AST

    # Small Intestine
    dASI <- GE * AST - kabs * ASI - kunabs * ASI
    dAabsSI <- kabs * ASI

    total_oral_uptake <- AabsSI + AabsST

    # Feces compartment
    dAfeces <- kbile * AL + kunabs * ASI

    # Liver
    dAL <- QL * (CA - CVL) * Free - kbile * AL + kabs * ASI + k0 * AST
    dAbile <- kbile * AL
    amount_per_gram_liver <- CLiver

    # Plasma compartment
    dAplas_free <- (QR * CVR * Free) + (QK * CVK * Free) + (QL * CVL * Free) -
      (QC * CA * Free) + dAefflux

    dCwater <- 0
    dingestion <- 0

    # Mass Balance Check
    Atissue <- Aplas_free + AR + AKb + Afil + APTC + AL + AST + ASI
    Aloss <- Aurine + Afeces
    Atotal <- Atissue + Aloss

    list(
      c(
        "dAR" = dAR, "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
        "dACl" = dACl, "dAefflux" = dAefflux,
        "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
        "dAurine" = dAurine, "dAST" = dAST,
        "dAabsST" = dAabsST, "dASI" = dASI, "dAabsSI" = dAabsSI, "dAfeces" = dAfeces,
        "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
        "dCwater" = dCwater, "dingestion" = dingestion
      ),
      "total_oral_uptake" = total_oral_uptake,
      "amount_per_gram_liver" = amount_per_gram_liver,
      "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal,
      "CR" = CR, "CVR" = CVR, "CKb" = CKb,
      "CVK" = CVK, "CPTC" = CPTC,
      "Cfil" = Cfil, "CL" = CL, "CVL" = CVL,
      "CA_free" = CA_free, "CA" = CA
    )
  })
}

# ================================================================================
# 6. SIMULATION SETUP
# ================================================================================

# --- Dosing and Administration Settings ---
BW <- 82  # Body weight (kg)
admin_type <- "bolus"  # administration type: "iv", "oral", or "bolus"
admin_dose_bolus <- c(3.96)  # administered dose through bolus (ug)
admin_time_bolus <- c(0)  # time when bolus doses are administered (days)
admin_dose_iv <- 0  # administered dose through IV (ug)
admin_time_iv <- 0  # time when IV doses are administered (days)

Cwater <- 0.00  # concentration in water (ug/L)
Cwater_time <- 0  # time of water concentration change
ingestion <- 0  # ingestion rate (ug/day)
ingestion_time <- c(0)  # time of ingestion rate change

# --- Simulation Time Settings ---
simulation_time <- 450  # total simulation time (days)
time_step <- 0.1  # time step for ODE solver (days)

cat("Setting up simulation...\n")

# Generate parameters
params <- create.params(BW)
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
exp_data_feces <- read.csv("Human_PBK_PFAS/exp_data_feces.csv")
feces_exp <- cumulative_exp_data(exp_data_feces, "time", "PFOA", "feces.weight") %>%
  mutate(cumulative_mass = cumulative_mass / 1000) %>%
  filter(time <= 6)

# Urine experimental data
exp_data_urine <- read.csv("Human_PBK_PFAS/exp_data_urine.csv")
urine_exp <- cumulative_exp_data(exp_data_urine, "time", "PFOA", "urine.volume") %>%
  mutate(time = time / 24) %>%
  filter(time <= 6)

# Plasma experimental data
plasma_exp <- read.csv("Human_PBK_PFAS/exp_data_plasma.csv") %>%
  select("time", "PFOA")

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
  ylim(0, 0.01) +
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
  geom_point(data = plasma_exp, aes(x = time, y = PFOA, color = "Experimental Data"),
             size = 4) +
  labs(
    title = "Plasma Concentration: Model Predictions vs Experimental Data",
    x = "Time (days)",
    y = "Concentration (ug/L)",
    color = "Legend"
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

print(plot_plasma)

cat("\nSimulation complete! All plots generated.\n")

# ================================================================================
# 9. SUMMARY OUTPUT
# ================================================================================

cat("\n=== SIMULATION SUMMARY ===\n")
cat("Body Weight:", BW, "kg\n")
cat("Administration Type:", admin_type, "\n")
if (admin_type == "bolus") {
  cat("Bolus Dose:", admin_dose_bolus, "ug at time", admin_time_bolus, "days\n")
} else if (admin_type == "iv") {
  cat("IV Dose:", admin_dose_iv, "ug at time", admin_time_iv, "days\n")
}
cat("Simulation Time:", simulation_time, "days\n")
cat("Final Mass in Plasma:", round(tail(solution$Aplas_free, 1), 4), "ug\n")
cat("Final Mass in Liver:", round(tail(solution$AL, 1), 4), "ug\n")
cat("Final Mass in Rest of Body:", round(tail(solution$AR, 1), 4), "ug\n")
cat("Total Mass Excreted in Urine:", round(tail(solution$Aurine, 1), 4), "ug\n")
cat("Total Mass Excreted in Feces:", round(tail(solution$Afeces, 1), 4), "ug\n")
cat("==========================\n")
