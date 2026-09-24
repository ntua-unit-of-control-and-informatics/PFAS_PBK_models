# =============================================================================
#  popgen_port.R
#  Base-R port of the PopGen virtual-population algorithm, following the
#  PyPopGenBE implementation (PyPI package `pypopgenbe` 1.5.5, MIT licence,
#  https://pypopgen.github.io/), which is itself the maintained Python port of
#  the PopGen MATLAB backend described in:
#    McNally K, Cotton R, Hogg A, Loizou G (2014) PopGen: A virtual human
#    population generator. Toxicology 315:70-85. doi:10.1016/j.tox.2013.07.009
#
#  SCOPE: reproduces the "ICRP" and "P3M" adult datasets (the two used in this
#  project). HSE and NDNS, and the paediatric ("Child") P3M/HSE branches, are
#  not ported; calling popgen_r() with age <= 16 under P3M stops with an error
#  rather than silently using the wrong reference table.
#
#  IMPORTANT - what "same results as the PyPopGen tool" can mean here:
#  R and Python/NumPy use different Mersenne-Twister sampling implementations
#  and draw random numbers in a different order even for an identical
#  algorithm, so INDIVIDUAL-level output cannot be made bit-identical to
#  pypopgenbe's, seed or no seed. What this port reproduces is the ALGORITHM
#  (every equation and rejection rule below is taken from the pypopgenbe
#  source, not re-derived) and therefore the population's DISTRIBUTION:
#  means, GM/GSD, percentiles and correlations match pypopgenbe's output for
#  large N, as shown in the validation block at the end of this file and
#  in popgen_port_validation.csv / .png.
#
#  Requires: popgen_consts.R (constants extracted from pypopgenbe's
#  popgenconsts.pkl) in the same directory. No other dependency (base R only).
# =============================================================================

# Run source("popgen_consts.R") yourself first if it isn't already loaded
# (e.g. from a different working directory); this only covers the common
# case of both files sitting together in the current directory.
if (!exists("POPGEN_CONSTS")) {
  if (file.exists("popgen_consts.R")) source("popgen_consts.R") else
    stop("popgen_port.R requires popgen_consts.R -- source() it first, or run from the folder containing both files")
}
CO <- POPGEN_CONSTS

# -----------------------------------------------------------------------------
#  Low-level samplers (impl/agecdf.py, ageinv.py, agernd.py, truncatednormrnd.py,
#  decreasingsquarelyrnd.py/decreasingsquarelyinv.py, normrnd0.py, lognrnd0.py)
# -----------------------------------------------------------------------------

# Piecewise-linear age CDF: rises steeply to age 45, more gently to age 100
# (this is what makes .pg_age_rnd() draw a realistic, not uniform, age).
.pg_age_cdf <- function(x, mid_age = 45, max_age = 100) {
  h <- 1 / (mid_age + 0.5 * (max_age - mid_age))
  p <- rep(NA_real_, length(x))
  p[x <= 0] <- 0; p[x >= max_age] <- 1
  young <- x > 0 & x <= mid_age
  p[young] <- h * x[young]
  old <- x > mid_age & x <= max_age
  # np.polyval([-0.5, 100, -1012.5], x) = -0.5*x^2 + 100*x - 1012.5
  p[old] <- (h / (max_age - mid_age)) * (-0.5 * x[old]^2 + 100 * x[old] - 1012.5)
  p
}
.pg_age_inv <- function(p, mid_age = 45, max_age = 100) {
  if (p == 0) return(0)
  if (p == 1) return(max_age)
  h <- 1 / (mid_age + 0.5 * (max_age - mid_age)); m <- h * mid_age
  if (p <= m) return(p / h)
  max_age - sqrt(7975 * (1 - p))
}
.pg_age_rnd <- function(lower, upper) {
  bounds <- .pg_age_cdf(c(lower, upper))
  q <- runif(1, bounds[1], bounds[2])
  .pg_age_inv(q)
}
.pg_assign_age <- function(population_type, age_range) {
  a <- if (population_type == "Realistic") .pg_age_rnd(age_range[1], age_range[2])
       else runif(1, age_range[1], age_range[2])
  max(sqrt(.Machine$double.eps), a)
}

.pg_assign_sex <- function(prob_of_male) if (runif(1) > prob_of_male) 2L else 1L  # 1=Male, 2=Female

.pg_create_ethnicity_breaks <- function(p) { stopifnot(length(p) == 3, all(p >= 0)); c(cumsum(p), 1 + .Machine$double.eps) }
.pg_assign_ethnicity <- function(breaks) sum(breaks < runif(1)) + 1L  # matches np.searchsorted(breaks, u) + 1

.pg_truncated_norm_rnd <- function(mu, sigma, lower, upper) {
  bl <- pnorm(lower, mu, sigma); bu <- pnorm(upper, mu, sigma)
  qnorm(runif(1, bl, bu), mu, sigma)
}
# "Decreasing squarely" inverse-CDF (cubic root of a 1/x^2-type density), used
# only for the HighVariation body-weight sampler.
.pg_decreasing_squarely_inv <- function(p, l, u) {
  a <- 1; b <- -3 * (l + u); c <- 3 * (l + u)^2; d <- -(l^3 + 3 * l^2 * u + 3 * l * u^2) - (u^3 - l^3) * p
  r <- polyroot(c(d, c, b, a))                       # polyroot() takes coefs low-to-high degree
  re <- Re(r)[abs(Im(r)) < 1e-6 & Re(r) >= l - 1e-6 & Re(r) <= u + 1e-6]
  if (length(re)) re[1] else Re(r)[which.min(abs(Im(r)))]  # fall back to the most-real root
}
.pg_decreasing_squarely_rnd <- function(l, u) .pg_decreasing_squarely_inv(runif(1), l, u)

.pg_normrnd0  <- function(mean, cv) (rnorm(length(mean)) * cv + 1) * mean
.pg_lognrnd0  <- function(mean, cv) {
  sigma <- sqrt(log(cv^2 + 1)); mu <- log(mean) - sigma^2 / 2
  exp(rnorm(length(mean)) * sigma + mu)
}
.pg_add_stochastic_variation <- function(target, cv, is_lognormal) {
  out <- target
  out[!is_lognormal] <- .pg_normrnd0(target[!is_lognormal], cv[!is_lognormal])
  out[is_lognormal]  <- .pg_lognrnd0(target[is_lognormal],  cv[is_lognormal])
  out
}

.pg_assign_target_height <- function(population_type, height_range, mean_height, cv) {
  if (population_type == "Realistic") rnorm(1, mean_height, mean_height * cv)
  else runif(1, height_range[1], height_range[2])
}
.pg_calculate_mass <- function(bmi, height_cm) bmi * (0.01 * height_cm)^2
.pg_assign_target_body_weight <- function(population_type, bw_range, mean_bw, cv) {
  if (population_type == "Realistic") .pg_truncated_norm_rnd(mean_bw, mean_bw * cv, bw_range[1], bw_range[2])
  else .pg_decreasing_squarely_rnd(bw_range[1], bw_range[2])
}

# -----------------------------------------------------------------------------
#  Reference height/body-mass at age (impl/generatepop.py dataset branches)
# -----------------------------------------------------------------------------

.pg_interp_lin_extrap <- function(age, ages, vals) {
  # equivalent to scipy interp1d(..., fill_value="extrapolate"): linear
  # interpolation within range, linear extrapolation from the two nearest
  # end points outside it (used for ICRP height/weight and the bone-density
  # decline curve; the ICRP cardiac-output table is looked up without
  # extrapolation in the source, see .pg_reference_height_weight())
  n <- length(ages)
  if (age < ages[1]) {
    slope <- (vals[2] - vals[1]) / (ages[2] - ages[1])
    return(vals[1] + slope * (age - ages[1]))
  }
  if (age > ages[n]) {
    slope <- (vals[n] - vals[n - 1]) / (ages[n] - ages[n - 1])
    return(vals[n] + slope * (age - ages[n]))
  }
  approx(ages, vals, xout = age)$y
}
.pg_poly_eval <- function(coefs, age) { # coefs highest-degree first (numpy Polynomial convention already applied)
  n <- length(coefs); s <- 0
  for (i in seq_len(n)) s <- s + coefs[i] * age^(n - i)
  s
}

.pg_reference_height_weight <- function(dataset, sex_name, ethnicity_name, age) {
  if (dataset == "ICRP") {
    h_vals <- if (sex_name == "Male") CO$icrp_height_male else CO$icrp_height_female
    w_vals <- if (sex_name == "Male") CO$icrp_bw_male     else CO$icrp_bw_female
    mean_height <- .pg_interp_lin_extrap(age, CO$icrp_height_ages, h_vals)
    mean_bw     <- .pg_interp_lin_extrap(age, CO$icrp_bw_ages,     w_vals)
    age_at_maturity <- min(age, if (sex_name == "Male") 20 else 16)
    mean_bw_maturity <- .pg_interp_lin_extrap(age_at_maturity, CO$icrp_bw_ages, w_vals)
  } else if (dataset == "P3M") {
    if (age <= 16) stop("popgen_r(): P3M 'Child' branch (age <= 16) is not ported; use ICRP or restrict age_range to > 16")
    key <- paste0(tolower(sex_name), ".", gsub("[- ]", "", tolower(ethnicity_name)))
    hc <- CO$p3m_height_poly[[key]]; wc <- CO$p3m_bw_poly[[key]]
    if (is.null(hc)) stop("popgen_r(): unknown P3M sex/ethnicity combination '", key, "'")
    mean_height <- .pg_poly_eval(hc, age)
    mean_bw     <- .pg_poly_eval(wc, age)
    age_at_maturity <- min(age, if (sex_name == "Male") 20 else 16)
    mean_bw_maturity <- .pg_poly_eval(wc, age_at_maturity)
  } else stop("popgen_r(): dataset must be 'ICRP' or 'P3M' in this port")
  list(mean_height = mean_height, mean_bw = mean_bw, mean_bw_maturity = mean_bw_maturity)
}

# -----------------------------------------------------------------------------
#  Organ mass targets (impl/calculatetargetorganmass.py + the brain/muscle/
#  bone/skin overrides it calls)
# -----------------------------------------------------------------------------

.pg_brain_mass <- function(age, sex_name) { # impl/calculatebrainmass.py, cites Bosgra et al. 2012
  B <- if (sex_name == "Male") 0.405 else 0.373
  B * (3.68 - 2.68 * exp(-age / 0.89)) * exp(-age / 629)
}
.pg_muscle_adj <- function(age, sex_name) { # impl/calculatemusclemassadjustmentfactor.py, cites Janssen et al. 2000
  if (age <= 45) return(1)
  c <- if (sex_name == "Male") -0.006181912 else -0.005522342
  1 + c * (age - 45)
}
.pg_skin_mass <- function(body_mass_kg, sex_name) { # impl/calculateskinmass.py
  costeff <- (4 * body_mass_kg + 7) / (body_mass_kg + 90)               # surface area, m^2
  thickness <- if (sex_name == "Male") 0.42 * 0.225 + 0.29 * 0.135 + 0.29 * 0.140
               else                    0.42 * 0.180 + 0.29 * 0.110 + 0.29 * 0.115
  costeff * thickness * 10
}
# impl/calculatebonemassadjustmentfactor.py
.pg_bone_density_tables <- list(
  ages = c(25, 35, 45, 55, 65, 75, 85),
  male = c(1.203, 1.214, 1.200, 1.186, 1.161, 1.140, 1.118),
  male_black = c(1.297, 1.294, 1.258, 1.244, 1.245, 1.180, 1.177),
  male_nbh   = c(1.170, 1.156, 1.149, 1.146, 1.135, 1.118, 1.074),
  female = c(1.104, 1.129, 1.118, 1.090, 1.032, 0.982, 0.917),
  female_black = c(1.191, 1.194, 1.191, 1.133, 1.093, 1.023, 0.978),
  female_nbh   = c(1.104, 1.115, 1.104, 1.104, 0.989, 0.939, 0.884)
)
.pg_bone_mass_adj <- function(age, sex_name, ethnicity_name) {
  # Exact port of impl/calculatebonemassadjustmentfactor.py: a "young" ramp
  # (peak bone accrual, ages 20-24 M / 18-22 F), an ethnicity offset (fixed,
  # White/Other = 0), and an "old" age-related decline read off a normalised
  # bone-density-vs-age curve (ethnicity-specific for Black / Non-black
  # Hispanic, common curve otherwise), each folded through
  # 1 + (adjustment) * 0.44 (bone_mineral_fraction) and multiplied together.
  bd <- .pg_bone_density_tables
  calc_af <- function(p) 1 + p * 0.44
  eth <- tolower(gsub("[- ]", "", ethnicity_name))
  is_black <- eth == "black"; is_nbh <- eth == "nonblackhispanic"
  if (sex_name == "Male") {
    af_young <- if (age >= 20 && age <= 24) calc_af(0.05 - 0.01 * (24 - age))
                else if (age >= 24) calc_af(0.05) else 1
    af_eth <- if (is_black) calc_af(0.08) else if (is_nbh) calc_af(-0.03) else 1
    curve <- if (is_black) bd$male_black / bd$male_black[1] else
             if (is_nbh)   bd$male_nbh   / bd$male_nbh[1]   else bd$male / bd$male[1]
  } else {
    af_young <- if (age >= 18 && age <= 22) calc_af(0.05 - 0.01 * (22 - age))
                else if (age >= 22) calc_af(0.05) else 1
    af_eth <- if (is_black) calc_af(0.055) else if (is_nbh) calc_af(-0.01) else 1
    curve <- if (is_black) bd$female_black / bd$female_black[1] else
             if (is_nbh)   bd$female_nbh   / bd$female_nbh[1]   else bd$female / bd$female[1]
  }
  af_old <- if (age >= 25) calc_af(.pg_interp_lin_extrap(age, bd$ages, curve) - 1) else 1
  af_young * af_old * af_eth
}

.pg_calculate_target_organ_mass <- function(age, sex_name, ethnicity_name,
                                             mean_bw_maturity, mean_height,
                                             target_bw, target_height, organ_mass_mean) {
  idx <- setNames(seq_along(CO$organ_names), CO$organ_names)
  scaled_height <- (target_height / mean_height)^0.75
  m <- scaled_height * organ_mass_mean * mean_bw_maturity
  m[idx["Brain"]]  <- .pg_brain_mass(age, sex_name)
  m[idx["Muscle"]] <- m[idx["Muscle"]] * .pg_muscle_adj(age, sex_name)
  m[idx["Bone"]]   <- m[idx["Bone"]] * .pg_bone_mass_adj(age, sex_name, ethnicity_name)
  m[idx["Skin"]]   <- .pg_skin_mass(target_bw, sex_name)
  m[idx["Adipose"]] <- target_bw - (sum(m) - m[idx["Adipose"]])
  list(scaled_height = scaled_height, mass = m)
}

# -----------------------------------------------------------------------------
#  Main generator: one virtual individual per call to .pg_one(); popgen_r()
#  loops it up to population_size, discarding and resampling exactly as
#  pypopgenbe's generatepop._generate_pop() does.
# -----------------------------------------------------------------------------

.pg_one <- function(dataset, age_range, bmi_range, height_range, prob_of_male,
                     ethnicity_probs, population_type) {
  idx <- setNames(seq_along(CO$organ_names), CO$organ_names)
  sex <- .pg_assign_sex(prob_of_male); sex_name <- c("Male", "Female")[sex]
  age <- .pg_assign_age(population_type, age_range)

  eth_breaks <- .pg_create_ethnicity_breaks(ethnicity_probs)
  eth_code <- .pg_assign_ethnicity(eth_breaks)
  ethnicity_name <- c("White", "Black", "NonBlackHispanic", "Other")[eth_code]

  ref <- .pg_reference_height_weight(dataset, sex_name, ethnicity_name, age)

  # ICRP cardiac-output-vs-age table: plain linear interpolation, no
  # extrapolation in the source (safe: age is always within [0, 80] <= 120)
  mean_co <- approx(CO$co_age_groups,
                     if (sex == 1) CO$co_values["Male", ] else CO$co_values["Female", ],
                     xout = age)$y

  bw_range <- c(.pg_calculate_mass(bmi_range[1], height_range[1]),
                .pg_calculate_mass(bmi_range[2], height_range[2]))
  target_height <- .pg_assign_target_height(population_type, height_range, ref$mean_height, CO$height_cv[sex])
  target_bw     <- .pg_assign_target_body_weight(population_type, bw_range, ref$mean_bw, CO$bw_cv[sex])

  om <- .pg_calculate_target_organ_mass(age, sex_name, ethnicity_name, ref$mean_bw_maturity,
                                         ref$mean_height, target_bw, target_height,
                                         if (sex == 1) CO$mass_mean["Male", ] else CO$mass_mean["Female", ])
  # suppressWarnings: a not-yet-rebased adipose estimate can be <= 0 here,
  # producing a log()-of-nonpositive NaN that the source explicitly silences
  # (np.seterr(invalid='ignore')) and lets the discard-on-NaN check catch below
  mass <- suppressWarnings(.pg_add_stochastic_variation(om$mass,
                                        if (sex == 1) CO$mass_cv["Male", ] else CO$mass_cv["Female", ],
                                        if (sex == 1) CO$mass_is_lognormal["Male", ] else CO$mass_is_lognormal["Female", ]))

  adipose_mean <- max(target_bw - (sum(mass) - mass[idx["Adipose"]]), 0.01)
  mass[idx["Adipose"]] <- suppressWarnings(.pg_lognrnd0(adipose_mean, 0.42))

  # A NaN can appear here (log() of a non-positive intermediate adipose
  # estimate inside the lognormal draw a few lines up, before the mean is
  # re-based on the actual target body weight) -- the source suppresses the
  # numpy warning and lets the ensuing "> 0" test fail (NaN > 0 is False in
  # numpy); R's NA > 0 is NA, not FALSE, so it is treated as a failure here too
  if (anyNA(mass) || !all(mass > 0)) return(NULL)                    # discard: negative/NaN tissue mass
  target_bw <- sum(mass)
  if (mass[idx["Adipose"]] / target_bw <= CO$min_adipose_fraction) return(NULL)  # discard: adipose <= 10%
  bmi <- target_bw / (0.01 * target_height)^2
  if (!(bmi >= bmi_range[1] && bmi <= bmi_range[2])) return(NULL)     # discard: BMI out of range

  flow_mean <- if (sex == 1) CO$flow_mean["Male", ] else CO$flow_mean["Female", ]
  flow_cv   <- if (sex == 1) CO$flow_cv["Male", ]   else CO$flow_cv["Female", ]
  target_co <- om$scaled_height * mean_co
  flow <- .pg_add_stochastic_variation(flow_mean * target_co, flow_cv, rep(FALSE, length(flow_cv)))
  target_co <- sum(flow) - flow[idx["Lung"]]                          # re-close CO after variation
  flow[idx["Lung"]] <- target_co                                      # Lung "flow" column := CO (pulmonary)

  list(age = age, sex = sex_name, ethnicity = ethnicity_name, height = target_height,
       body_mass = target_bw, cardiac_output = target_co, mass = mass, flow = flow)
}

#' Generate a PopGen-equivalent virtual population in R
#'
#' Algorithmic port of pypopgenbe::generate_pop() (ICRP and P3M adult datasets
#' only -- see the header of popgen_port.R). Returns one row per individual,
#' with mass in kg and flow in the units of `flow_units`. This is a
#' STATISTICAL, not a bit-for-bit, reproduction of the PyPopGen tool's output
#' (see the header note on RNG streams); validate any downstream use against
#' pypopgenbe for your own N before relying on it.
#'
#' @param n            population size (accepted individuals)
#' @param dataset      "ICRP" or "P3M"
#' @param age_range    c(lower, upper) years
#' @param bmi_range    c(lower, upper) kg/m2
#' @param height_range c(lower, upper) cm
#' @param prob_male    probability an individual is male
#' @param ethnicity_probs c(White, Black, NonBlackHispanic); must sum <= 1
#' @param population_type "Realistic" (default) or "HighVariation"
#' @param flow_units   "LitresPerHour" (default) or "MilliLitresPerMinute"
#' @param seed         optional RNG seed (R's own stream; not comparable to a
#'                      Python seed -- see header)
#' @param lump_GI      if TRUE (default), also add a "GI tract" mass/flow
#'                      column = Stomach + Small intestine + Large intestine
#' @return list(pop = data.frame, n_discarded = integer)
popgen_r <- function(n, dataset = "ICRP", age_range = c(18, 70), bmi_range = c(18.5, 30),
                      height_range = c(150, 200), prob_male = 0.5,
                      ethnicity_probs = c(1, 0, 0), population_type = "Realistic",
                      flow_units = c("LitresPerHour", "MilliLitresPerMinute"),
                      seed = NULL, lump_GI = TRUE) {
  flow_units <- match.arg(flow_units)
  if (!is.null(seed)) set.seed(seed)
  # underscores, not the raw organ names (which contain spaces in "Small
  # intestine"/"Large intestine") -- rbind.data.frame() silently mangles
  # spaces to dots in column names, which would desync every downstream
  # column reference below, so sanitise once, up front
  organ_names <- CO$organ_names
  organ_col   <- gsub(" ", "_", organ_names)
  rows <- vector("list", n); got <- 0L; discarded <- 0L
  while (got < n) {
    ind <- .pg_one(dataset, age_range, bmi_range, height_range, prob_male, ethnicity_probs, population_type)
    if (is.null(ind)) { discarded <- discarded + 1L; next }
    got <- got + 1L
    row <- c(list(id = got, age = ind$age, sex = ind$sex, ethnicity = ind$ethnicity,
                  height_cm = ind$height, body_mass_kg = ind$body_mass,
                  cardiac_output = ind$cardiac_output),
             setNames(as.list(ind$mass), paste0(organ_col, "_mass_kg")),
             setNames(as.list(ind$flow), paste0(organ_col, "_flow")))
    rows[[got]] <- row
  }
  pop <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
  flow_cols <- paste0(organ_col, "_flow")
  to_Lph <- if (flow_units == "LitresPerHour") 0.06 else 1              # PyPopGenBE flows are generated in mL/min
  pop[flow_cols] <- pop[flow_cols] * to_Lph
  pop$cardiac_output <- pop$cardiac_output * to_Lph
  if (lump_GI) {
    pop$GI_tract_mass_kg <- pop$Stomach_mass_kg + pop$Small_intestine_mass_kg + pop$Large_intestine_mass_kg
    pop$GI_tract_flow    <- pop$Stomach_flow    + pop$Small_intestine_flow    + pop$Large_intestine_flow
  }
  # impl/generatepop.py post-processing, reproduced for parity with the
  # pypopgenbe CSV columns used elsewhere in this project (popgen_import_csv.R):
  #  - "Liver Total" = sum of the flows feeding the liver (portal system),
  #    used for the Sec. 2.1.3 portal-vein balance check
  #  - "Lung Bronchial" flow = the actual pulmonary-tissue perfusion (2.5% of
  #    the Lung/CO placeholder); this, not the "Lung" column, is what
  #    popgen_import_csv.R maps onto the PBK model's lung compartment
  pop$Liver_Total_flow  <- pop$Liver_flow + pop$Pancreas_flow + pop$Spleen_flow + pop$Stomach_flow +
                            pop$Small_intestine_flow + pop$Large_intestine_flow
  pop$Lung_Bronchial_flow <- pop$Lung_flow * CO$bronchial_flow_fraction
  attr(pop, "flow_units") <- flow_units
  list(pop = pop, n_discarded = discarded)
}

#' Reshape a popgen_r() population into the popgen_read_csv() / popgen_population()
#' layout (id, sex M/F, age, height in m, BM, BMI, M_<organ>, Q_<organ>, CO,
#' CO_popgen, ethnicity), so it can be dropped into run_population_mc.R via a
#' new `pop_source` branch without touching the rest of that pipeline. Q_lung
#' is the bronchial (pulmonary-tissue) flow and CO is the sum of all 15 flows
#' with bronchial standing in for Lung -- the same "[A]" convention documented
#' in popgen_import_csv.R, not pypopgenbe's own Cardiac Output (kept as
#' CO_popgen). Requires popgen_r(..., flow_units = "LitresPerHour").
popgen_r_as_population <- function(pop) {
  if (!identical(attr(pop, "flow_units"), "LitresPerHour"))
    stop("popgen_r_as_population(): call popgen_r(..., flow_units = 'LitresPerHour') first")
  m <- c(lung = "Lung", brain = "Brain", kidney = "Kidneys", liver = "Liver", pancreas = "Pancreas",
         spleen = "Spleen", stomach = "Stomach", int_small = "Small_intestine", int_large = "Large_intestine",
         muscle = "Muscle", heart = "Heart", adipose = "Adipose", skeleton = "Bone", skin = "Skin", gonads = "Gonads")
  M <- as.matrix(pop[paste0(m, "_mass_kg")]); colnames(M) <- paste0("M_", names(m))
  Q <- as.matrix(pop[paste0(m, "_flow")]);    colnames(Q) <- paste0("Q_", names(m))
  Q[, "Q_lung"] <- pop$Lung_Bronchial_flow                          # bronchial, not the CO placeholder
  height_m <- pop$height_cm / 100
  data.frame(id = pop$id, sex = ifelse(pop$sex == "Male", "M", "F"), age = pop$age,
             height = height_m, BM_int = NA_real_, BM = pop$body_mass_kg, BMI = pop$body_mass_kg / height_m^2,
             as.data.frame(M), as.data.frame(Q), CO = rowSums(Q), CO_popgen = pop$cardiac_output,
             ethnicity = pop$ethnicity, stringsAsFactors = FALSE)
}
