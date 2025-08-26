BW_shah <- 280  # reference body weight for Shah rat (g)

shah_rat <- data.frame(
  Compartment = c("Heart", "Lung", "Muscle", "Skin", "Adipose", "Bone", "Brain",
                  "Kidney", "Liver", "S_Intestine", "L_Intestine", "Pancreas", "Thymus",
                  "Spleen", "Ly_Node", "Other"),
  Total = c(1.02, 1.40, 122, 49.9, 33.1, 21.0, 2.28, 2.41, 15.7, 4.99, 2.87, 1.00, 0.096, 2.77, 1.15, 6.09),
  Plasma = c(0.0394, 0.231, 2.68, 1.87, 0.364, 0.462, 0.0502, 0.132, 1.34, 0.0795, 0.0458, 0.0547, 0.00528, 0.335, NA, 0.256),
  Blood_cells = c(0.0323, 0.189, 2.19, 1.53, 0.298, 0.378, 0.0410, 0.108, 1.10, 0.0651, 0.0375, 0.0448, 0.00432, 0.274, NA, 0.209),
  Interstitial = c(0.146, 0.263, 15.8, 16.5, 5.63, 3.90, 0.410, 0.361, 2.56, 0.867, 0.5, 0.173, 0.0163, 0.554, NA, 1.04),
  Cellular = c(0.801, 0.710, 100, 29.8, 26.7, 16.1, 1.77, 1.79, 10.7, 3.95, 2.28, 0.717, 0.0696, 1.59, 1.15, 4.55)
)

# Function to return a table with all values for each compartment
get_shah_volumes_table <- function(bw_user) {
  scale_factor <- bw_user / BW_shah
  df <- shah_rat
  df$Total_scaled <- df$Total * scale_factor
  df$Capillary_Blood_scaled <- (df$Plasma + df$Blood_cells) * scale_factor
  df$Interstitial_scaled <- df$Interstitial * scale_factor
  df$Intracellular_scaled <- df$Cellular * scale_factor
  result <- df[, c("Compartment", "Total_scaled", "Capillary_Blood_scaled", "Interstitial_scaled", "Intracellular_scaled")]
  colnames(result) <- c("Compartment", "Total_mL", "Capillary_Blood_mL", "Interstitial_mL", "Intracellular_mL")
  return(result)
}

# Example usage:
bw_user <- 250
shah_table <- get_shah_volumes_table(bw_user)
View(shah_table)


BW_kawai <- 250  # reference body weight for Kawai model (g)
organ_data_kawai <- data.frame(
  Organ = c("Lung", "Brain", "Heart", "Kidneys", "Bone", "Thymus", "Muscle",
            "Stomach", "Spleen", "Liver", "Gut", "Pancreas", "Skin", "Adipose"),
  Mass_g = c(1.0, 1.7, 0.8, 2.3, 15.8, 0.7, 122, 1.1, 0.6, 10.3, 10, 1.3, 40, 10),
  Fwa = c(0.175, 0.014, 0.061, 0.046, 0.019, 0.009, 0.004, 0.011, 0.321, 0.057, 0.010, 0.032, 0.002, 0.005),
  Fwb = c(0.262, 0.037, 0.262, 0.105, 0.041, 0.030, 0.026, 0.032, 0.282, 0.115, 0.024, 0.180, 0.019, 0.010),
  Fvi = c(0.188, 0.004, 0.100, 0.200, 0.100, 0.150, 0.120, 0.100, 0.150, 0.163, 0.094, 0.120, 0.302, 0.135)
)

# Add residual blood fraction (large vessel blood) as Fwb - Fwa
organ_data_kawai$f_residual <- organ_data_kawai$Fwb - organ_data_kawai$Fwa

get_kawai_volumes_table <- function(bw_user) {
  scale_factor <- bw_user / BW_kawai
  df <- organ_data_kawai
  df$Mass_scaled <- df$Mass_g * scale_factor
  df$Capillary_Blood_mL <- df$Mass_scaled * df$Fwa
  df$Large_vessels_mL <- df$Mass_scaled * (df$Fwb - df$Fwa)
  df$Residual_Blood_mL <- df$Mass_scaled * df$Fwb
  df$Interstitial_mL <- df$Mass_scaled * df$Fvi
  df$Intracellular_mL <- df$Mass_scaled * (1 - (df$Fwb + df$Fvi))
  result <- df[, c("Organ", "Mass_scaled", "Capillary_Blood_mL", "Large_vessels_mL",
                   "Residual_Blood_mL", "Interstitial_mL", "Intracellular_mL")]
  colnames(result)[2] <- "Total_mL"
  return(result)
}

# Example usage:
bw_user <- 250
kawai_table <- get_kawai_volumes_table(bw_user)
View(kawai_table)


# Example: Organ names, volumes (ml), and tissue fractions from PK-Sim screenshots

organ_pksim <- data.frame(
  Organ = c("Venous Blood", "Arterial Blood", "Bone", "Brain", "Fat", "Gonads", "Heart", "Kidney",
            "Large Intestine", "Liver", "Lung", "Muscle", "Pancreas", "Portal Vein", "Skin",
            "Small Intestine", "Spleen", "Stomach"),
  Volume_mL = c(7.42, 3.23, 17.34, 1.87, 10.98, 2.74, 0.88, 2.52, 2.39, 11.3, 1.10, 133.90, 1.43, 1650, 43.90, 5.49, 0.66, 1.21),
          f_vascular = c(1, 1, 0.04, 0.04, 0.01, 0.14, 0.26, 0.11, 0.03, 0.12, 0.63, 0.03, 0.18, 1, 0.02, 0.02, 0.28, 0.03),
  f_intracellular = c(NA, NA, 0.86, 0.96, 0.86, 0.79, 0.64, 0.70, 0.87, 0.72, 0.19, 0.85, 0.70, NA, 0.68, 0.88, 0.57, 0.87),
  f_interstitial = c(NA, NA, 0.10, 0.004, 0.14, 0.07, 0.10, 0.20, 0.10, 0.16, 0.19, 0.12, 0.12, NA, 0.30, 0.09, 0.15, 0.10)
  
)

# Compute actual space volumes (ml)
organ_pksim$Capillary_Blood_mL   <- organ_pksim$Volume_mL * organ_pksim$f_vascular
organ_pksim$Interstitial_mL      <- organ_pksim$Volume_mL * organ_pksim$f_interstitial
organ_pksim$Intracellular_mL     <- organ_pksim$Volume_mL * organ_pksim$f_intracellular

# Preview result
print(organ_pksim[, c("Organ", "Volume_mL", "Capillary_Blood_mL", "Interstitial_mL", "Intracellular_mL")])


