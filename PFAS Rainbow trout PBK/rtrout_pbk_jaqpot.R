library(deSolve)

# =============================================================================
# 1. Parameters
# =============================================================================

create.params <- function(user.input) {
  with(as.list(user.input), {

    # --- Fitted parameters (common across all substances, from optimization) ---
    Ku       <- 1.466892e+00
    CLU_coef <- 5.719031e-04
    Cl_feces <- 1.306465e+00

    # --- Substance-specific fitted partition coefficients (from optimization) ---
    fitted_P_params <- list(
      PFOS  = c(P_liver=1.5685318,  P_muscle=0.11316576, P_kidney=0.4398762,
                P_skin=0.2715996,   P_gills=0.2291869,   P_carcass=0.1074253,
                P_viscera=3.699111e+00),
      PFOA  = c(P_liver=2.0035690,  P_muscle=0.03687850, P_kidney=0.8511917,
                P_skin=0.3187692,   P_gills=0.3428033,   P_carcass=0.1799704,
                P_viscera=5.637376e-01),
      PFBS  = c(P_liver=1.7414537,  P_muscle=0.13872068, P_kidney=0.7630533,
                P_skin=0.2242364,   P_gills=0.2771880,   P_carcass=0.1160583,
                P_viscera=1.076958e-05),
      PFHxS = c(P_liver=1.6979372,  P_muscle=0.04341057, P_kidney=0.3607912,
                P_skin=0.2934588,   P_gills=0.1537518,   P_carcass=0.0411729,
                P_viscera=8.846102e-06),
      PFNA  = c(P_liver=0.8032868,  P_muscle=0.06493063, P_kidney=0.2460815,
                P_skin=0.2335112,   P_gills=0.2204153,   P_carcass=0.1135412,
                P_viscera=1.274696e+00)
    )

    P         <- fitted_P_params[[substance]]
    P_liver   <- P[["P_liver"]]
    P_muscle  <- P[["P_muscle"]]
    P_kidney  <- P[["P_kidney"]]
    P_skin    <- P[["P_skin"]]
    P_gills   <- P[["P_gills"]]
    P_carcass <- P[["P_carcass"]]
    P_viscera <- P[["P_viscera"]]

    # --- Temperature correction ---
    Texp_K         <- 273 + Texp
    Tref           <- 273 + c(6, 12, 18)
    keep_ref_value <- which.min(abs(Tref - Texp_K))

    # Cardiac output and BW reference values at 6, 12, 18 °C (Barron et al. 1987)
    F_card_ref <- c(1.188, 2.322, 3.75)[keep_ref_value]   # ml/h/g
    BW_ref     <- c(270.1, 296.4, 414.5)[keep_ref_value]  # g

    TA <- 6930  # Arrhenius temperature (K), Grech et al. 2018
    Tr <- Tref[keep_ref_value]
    KT <- exp(TA/Tr - TA/Texp_K)

    # --- Body weight: clamp user-supplied BW to the calibrated range [314, 808] g ---
    BW_min <- 314  # minimum BW in calibration dataset (Falk et al. 2015)
    BW_max <- 808  # maximum BW in calibration dataset (Falk et al. 2015)
    BW_eff <- max(min(BW, BW_max), BW_min)

    # --- Tissue weight fractions (Vidal et al. 2019) ---
    fw_Liver   <- 0.0120
    fw_Blood   <- 0.0450
    fw_Skin    <- 0.0640
    fw_Muscle  <- 0.566
    fw_Gills   <- 0.0200
    fw_Kidney  <- 0.0160
    fw_Viscera <- 0.051
    fw_lumen   <- 0.012

    # --- Blood flow fractions (Vidal et al. 2019) ---
    fb_Liver   <- 0.0035
    fb_Skin    <- 0.0728
    fb_Muscle  <- 0.655
    fb_Gills   <- 0.0021
    fb_Kidney  <- 0.071
    fb_Viscera <- 0.069

    # --- Substance-specific chemical parameters ---
    if (substance == 'PFOA') {
      a          <- 0.138   # Sun et al. 2022, Goeritz et al. 2013
      f_reab_hep <- 0.30    # Cao et al. 2022
      K_urine    <- 2.08
      Cl_urine   <- 0.029 * 3600  # 1/h, Sun et al. 2022
    } else if (substance == 'PFNA') {
      a          <- 0.522
      f_reab_hep <- 0.34
      K_urine    <- 1.35
      Cl_urine   <- 0.050 * 3600
    } else if (substance == 'PFBS') {
      a          <- 0.0598
      f_reab_hep <- 0.23
      K_urine    <- 5.88
      Cl_urine   <- 0.023 * 3600
    } else if (substance == 'PFHxS') {
      a          <- 0.558
      f_reab_hep <- 0.30
      K_urine    <- 5.88
      Cl_urine   <- 0.023 * 3600
    } else if (substance == 'PFOS') {
      a          <- 0.721
      f_reab_hep <- 0.42
      K_urine    <- 1.35
      Cl_urine   <- 0.050 * 3600
    }

    # --- Physiological flow/volume constants ---
    Q_bile_coef  <- 7.5e-05   # ml/g BW/h  (Grosell et al. 2000)
    Q_urine_coef <- 2.755e-03 # ml/h/g BW  (Curtis et al. 1981)
    V_urine_coef <- 2.2e-03   # ml/g BW    (Curtis et al. 1981)

    a_skin   <- 0.9  # fraction of skin venous blood routed to kidney (Nichols et al. 1996)
    a_muscle <- 0.6  # fraction of muscle venous blood routed to kidney (Nichols et al. 1996)
    plasma   <- 0.7

    return(list(
      'BW'         = BW_eff,
      'F_card_ref' = F_card_ref, 'BW_ref' = BW_ref, 'KT' = KT,

      'admin.dose' = admin.dose,
      'admin.time' = admin.time,

      'fw_Liver'   = fw_Liver,   'fw_Blood'   = fw_Blood,  'fw_Skin'    = fw_Skin,
      'fw_Muscle'  = fw_Muscle,  'fw_Gills'   = fw_Gills,  'fw_Kidney'  = fw_Kidney,
      'fw_Viscera' = fw_Viscera, 'fw_lumen'   = fw_lumen,

      'fb_Liver'   = fb_Liver,   'fb_Skin'    = fb_Skin,   'fb_Muscle'  = fb_Muscle,
      'fb_Gills'   = fb_Gills,   'fb_Kidney'  = fb_Kidney, 'fb_Viscera' = fb_Viscera,

      'a_skin' = a_skin, 'a_muscle' = a_muscle,
      'Q_bile_coef'  = Q_bile_coef,
      'Q_urine_coef' = Q_urine_coef, 'V_urine_coef' = V_urine_coef,
      'K_urine'    = K_urine,   'Cl_urine'   = Cl_urine,
      'f_reab_hep' = f_reab_hep, 'plasma' = plasma,
      'Free' = 1, 'a' = a,

      'Ku'       = Ku,
      'CLU_coef' = CLU_coef,
      'Cl_feces' = Cl_feces,

      'P_liver'   = P_liver,   'P_muscle'  = P_muscle,  'P_kidney'  = P_kidney,
      'P_skin'    = P_skin,    'P_gills'   = P_gills,   'P_carcass' = P_carcass,
      'P_viscera' = P_viscera
    ))
  })
}

# =============================================================================
# 2. Initial conditions
# =============================================================================

create.inits <- function(parameters) {
  with(as.list(parameters), {
    M_art <- 0; M_venous <- 0
    M_gills <- 0; M_lumen <- 0; M_lumen_2 <- 0; M_viscera <- 0
    M_liver <- 0; M_kidney <- 0; M_muscle <- 0; M_skin <- 0
    M_carcass <- 0; M_storage <- 0; M_urine <- 0; M_feces <- 0; M_input <- 0

    return(c(
      'M_art'     = M_art,     'M_venous'  = M_venous,
      'M_gills'   = M_gills,   'M_lumen'   = M_lumen,   'M_lumen_2' = M_lumen_2,
      'M_viscera' = M_viscera, 'M_liver'   = M_liver,   'M_kidney'  = M_kidney,
      'M_muscle'  = M_muscle,  'M_skin'    = M_skin,    'M_carcass' = M_carcass,
      'M_storage' = M_storage, 'M_urine'   = M_urine,   'M_feces'   = M_feces,
      'M_input'   = M_input
    ))
  })
}

# =============================================================================
# 3. Dosing events
# =============================================================================

create.events <- function(parameters) {
  with(as.list(parameters), {
    ldose  <- length(admin.dose)
    ltimes <- length(admin.time)

    if (ltimes != ldose) {
      stop("admin.time and admin.dose must have the same length")
    }

    events <- data.frame(
      var    = rep(c('M_lumen', 'M_input'), ldose),
      time   = sort(rep(admin.time, 2)),
      value  = rep(admin.dose, each = 2),
      method = 'add'
    )
    return(list(data = events))
  })
}

# =============================================================================
# 4. Custom function (placeholder required by deploy.pbpk)
# =============================================================================

custom.func <- function() {
  return()
}

# =============================================================================
# 5. ODE system
# =============================================================================

ode.func <- function(time, inits, params, custom.func) {
  with(as.list(c(inits, params)), {

    # BW is constant (already clamped to calibration range in create.params)

    # Cardiac output (plasma flow, ml/h)
    Q_total <- F_card_ref * KT * (BW / BW_ref)^(-0.1) * BW * plasma

    # Tissue masses (g)
    w_blood   <- fw_Blood   * BW * plasma
    w_liver   <- fw_Liver   * BW
    w_skin    <- fw_Skin    * BW
    w_muscle  <- fw_Muscle  * BW
    w_gills   <- fw_Gills   * BW
    w_kidney  <- fw_Kidney  * BW
    w_viscera <- fw_Viscera * BW
    w_lumen   <- fw_lumen   * BW
    w_art     <- (1/3) * w_blood
    w_venous  <- (2/3) * w_blood
    w_carcass <- BW - (w_blood / plasma + w_liver + w_skin + w_muscle +
                         w_gills + w_kidney + w_viscera + w_lumen)

    # Regional blood flows (ml/h)
    Q_liver   <- fb_Liver   * Q_total
    Q_skin    <- fb_Skin    * Q_total
    Q_muscle  <- fb_Muscle  * Q_total
    Q_gills   <- Q_total
    Q_kidney  <- fb_Kidney  * Q_total
    Q_viscera <- fb_Viscera * Q_total
    Q_carcass <- Q_total - (Q_liver + Q_skin + Q_muscle + Q_kidney + Q_viscera)

    # Auxiliary flows
    Q_bile  <- Q_bile_coef  * BW  # ml/h
    Q_urine <- Q_urine_coef * BW  # ml/h
    v_urine <- V_urine_coef * BW  # ml

    # Urinary reabsorption rate (1/h)
    f_reab_urine <- Cl_urine * CLU_coef / K_urine

    # Tissue concentrations (ug/g or ug/ml)
    C_gills   <- M_gills   / w_gills
    C_viscera <- M_viscera / w_viscera
    C_liver   <- M_liver   / w_liver
    C_kidney  <- M_kidney  / w_kidney
    C_muscle  <- M_muscle  / w_muscle
    C_skin    <- M_skin    / w_skin
    C_carcass <- M_carcass / w_carcass
    C_lumen   <- (M_lumen + M_lumen_2) / w_lumen
    C_art     <- M_art     / w_art
    C_venous  <- M_venous  / w_venous
    C_blood   <- (M_art + M_venous) / w_blood
    C_storage <- M_storage / v_urine

    # ODEs
    dM_art <- Free * Q_gills * C_gills / P_gills -
      (Q_viscera + Q_liver + Q_kidney + Q_muscle + Q_skin + Q_carcass) * Free * C_art

    dM_venous <- -Free * Q_total * C_venous +
      ((Q_liver + Q_viscera) * C_liver / P_liver +
         (Q_kidney + a_muscle * Q_muscle + a_skin * Q_skin) * C_kidney / P_kidney +
         (1 - a_muscle) * Q_muscle * C_muscle / P_muscle +
         (1 - a_skin)   * Q_skin   * C_skin   / P_skin +
         Q_carcass * C_carcass / P_carcass) * Free

    dM_gills <- Q_gills * Free * (C_venous - C_gills / P_gills)

    dM_input <- 0

    # Lumen – fraction a is available for absorption
    dM_lumen   <- -Ku * a * M_lumen - Cl_feces * (1 - a) * M_lumen

    # Lumen_2 – bile contents, only eliminable via feces
    dM_lumen_2 <- (1 - f_reab_hep) * Q_bile * C_liver - Cl_feces * M_lumen_2

    dM_viscera <- Q_viscera * Free * (C_art - C_viscera / P_viscera) +
      Ku * a * M_lumen + f_reab_hep * Q_bile * C_liver

    dM_Liver <- Q_liver * Free * C_art + Q_viscera * Free * C_viscera / P_viscera -
      (Q_liver + Q_viscera) * Free * C_liver / P_liver - Q_bile * C_liver

    dM_kidney <- Q_kidney * Free * C_art -
      (Q_kidney + a_muscle * Q_muscle + a_skin * Q_skin) * Free * C_kidney / P_kidney +
      a_muscle * Q_muscle * Free * C_muscle / P_muscle +
      a_skin   * Q_skin   * Free * C_skin   / P_skin -
      Cl_urine * CLU_coef * M_kidney + f_reab_urine * M_storage

    dM_muscle  <- Q_muscle  * Free * (C_art - C_muscle  / P_muscle)
    dM_skin    <- Q_skin    * Free * (C_art - C_skin    / P_skin)
    dM_carcass <- Q_carcass * Free * (C_art - C_carcass / P_carcass)

    dM_storage <- Cl_urine * CLU_coef * M_kidney - f_reab_urine * M_storage -
      Q_urine * C_storage

    dM_urine <- Q_urine * C_storage

    dM_feces <- Cl_feces * ((1 - a) * M_lumen + M_lumen_2)

    Mass_balance <- M_input - (M_art + M_venous + M_gills + M_lumen + M_lumen_2 +
                                 M_viscera + M_liver + M_kidney + M_muscle +
                                 M_skin + M_carcass + M_storage + M_urine + M_feces)

    return(list(
      c(
        'dM_art'     = dM_art,     'dM_venous'  = dM_venous,
        'dM_gills'   = dM_gills,   'dM_lumen'   = dM_lumen,   'dM_lumen_2' = dM_lumen_2,
        'dM_viscera' = dM_viscera, 'dM_Liver'   = dM_Liver,   'dM_kidney'  = dM_kidney,
        'dM_muscle'  = dM_muscle,  'dM_skin'    = dM_skin,    'dM_carcass' = dM_carcass,
        'dM_storage' = dM_storage, 'dM_urine'   = dM_urine,   'dM_feces'   = dM_feces,
        'dM_input'   = dM_input
      ),
      'C_Gills'      = C_gills,
      'C_Viscera'    = C_viscera,
      'C_Liver'      = C_liver,
      'C_Kidney'     = C_kidney,
      'C_Muscle'     = C_muscle,
      'C_Skin'       = C_skin,
      'C_Carcass'    = C_carcass,
      'C_Lumen'      = C_lumen,
      'C_Blood'      = C_blood * plasma,
      'Mass_balance' = Mass_balance,
      'BW'           = BW
    ))
  })
}

# =============================================================================
# 6. User input
# =============================================================================
# Default values reproduce the Falk et al. 2015 PFOS dietary exposure experiment:
#   - 500 ug PFAS/kg food, 2.6% BW/day feeding rate, 28-day uptake, 56-day total
#   - Using constant BW = 314 g (initial fish weight from the experiment)

user.input <- list(
  'substance'  = 'PFOS',
  'Texp'       = 15,                         # experimental temperature (°C)
  'BW'         = 314,                        # fish body weight (g); clamped to [314, 808]
  'admin.dose' = rep(314*0.026*500/1000, 28),# dose per feeding event (ug)
  'admin.time' = seq(0, 27*24, 24)          # time of each feeding event (hours)
)


params <- create.params(user.input)
inits  <- create.inits(params)
events <- create.events(params)

sim.start  = 0                          # simulation start time (hours)
sim.end    = 56 * 24                    # simulation end time (hours)
sim.step   = 1                           # output time step (hours)

sample_time <- seq(sim.start, sim.end, sim.step)

solution <- as.data.frame(
  deSolve::ode(
    times  = sample_time,
    func   = ode.func,
    y      = inits,
    parms  = params,
    events = events,
    method = "lsodes", rtol = 1e-03, atol = 1e-03,
    custom.func = custom.func
  )
)

# =============================================================================
# 7. Upload to Jaqpot
# =============================================================================

predicted.feats <- c(
  'M_art', 'M_venous', 'M_gills', 'M_lumen', 'M_lumen_2',
  'M_viscera', 'M_liver', 'M_kidney', 'M_muscle', 'M_skin',
  'M_carcass', 'M_storage', 'M_urine', 'M_feces',
  'C_Gills', 'C_Viscera', 'C_Liver', 'C_Kidney',
  'C_Muscle', 'C_Skin', 'C_Carcass', 'C_Lumen', 'C_Blood',
  'Mass_balance', 'BW'
)

# jaqpotr::login.api(login.url = "https://api.jaqpot.org")

jaqpotr::deploy.pbpk(
  user.input    = user.input,
  out.vars      = predicted.feats,
  create.params = create.params,
  create.inits  = create.inits,
  create.events = create.events,
  custom.func   = custom.func,
  ode.fun       = ode.func,
  envFile       = ".env"
)
