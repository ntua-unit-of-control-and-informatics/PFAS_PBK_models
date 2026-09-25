# =============================================================================
#  Reverse_Dosimetry_Stochastic.r
#  Stochastic (population) counterpart to shiny_app/R/reverse_dosimetry.r's
#  .reverse_dosimetry().
#
#  .reverse_dosimetry(chemical, BW, duration, time_scale, POD, compartment,
#  admin_type, admin_time, optim_tol, solver, rtol, atol) takes a single
#  fixed body weight and numerically searches (nloptr, 1-D, log-dose space)
#  for the constant continuous exposure level (or, for admin_type != "oral",
#  the single bolus dose) that reproduces a target internal concentration
#  (POD) in a chosen compartment at the end of the simulation, for that one
#  fixed BW. NOTE: the roxygen block above .reverse_dosimetry_set_up() in
#  the source file documents a much richer reverse_dosimetry() -- multi-row
#  bmd_df, Bayesian model averaging, isStochastic/isParallel/Nsamples --
#  but that richer function is NOT actually implemented in
#  reverse_dosimetry.r; only the simpler single-POD .reverse_dosimetry() is
#  real code, and that is what this file mirrors. Its "isStochastic" also
#  means something different from what's built here: BMA posterior-model-
#  weight/log-normal error uncertainty, not popgen inter-individual
#  variability.
#
#  reverse_dosimetry_stochastic() below takes the SAME dosing/target
#  arguments MINUS admin_type/admin_dose/admin_time, which no longer exist
#  anywhere in Extended_Model_Population.r (see that file's .create.events():
#  dosing is now selected purely by exp_type -- "pharmacokinetics" for a
#  one-off dose, "continuous" for a repeating intake rate -- with no
#  separate bolus/oral/iv switch). It also replaces the single fixed BW with
#  a popgen virtual population -- the same population-generation inputs
#  Extended_Model_Population.r's .generate_population() requires (n,
#  dataset, age_range, bmi_range, height_range, prob_male, ethnicity_probs,
#  population_type, optional seed). The SAME 1-D nloptr search
#  .reverse_dosimetry() runs once is run once PER INDIVIDUAL, using that
#  individual's own popgen physiology (organ mass/flow, sex-specific
#  hematocrit, age-corrected GFR) and own BW in the per-kg-BW ->
#  absolute-amount conversion (Forward_Dosimetry_Stochastic.r's/
#  forward_dosimetry()'s `exp(x) * BW / 1000` formula) -- so the result is a
#  POPULATION DISTRIBUTION of "what exposure reproduces this POD" instead of
#  one point estimate.
#
#  Design note: since .create.events() dropped the admin_type distinction,
#  the optimizer's search variable exp(x) is ALWAYS a per-kg-BW rate here,
#  BW-scaled per individual -- for BOTH exp_type = "pharmacokinetics" (a
#  one-off dose search) and exp_type = "continuous" (a steady-rate search).
#  This mirrors Forward_Dosimetry_Stochastic.r's own design choice (there,
#  `ingestion` is always per-kg-BW regardless of exp_type). It differs from
#  the OLD non-oral/bolus reverse-dosimetry path, which searched for a flat
#  ABSOLUTE dose identical for every individual regardless of their BW; that
#  option no longer exists.
#
#  POD uncertainty note: `POD` itself is normally a measured/estimated
#  quantity with its own uncertainty (replicate assay measurements, a
#  reported CI, inter-study spread, ...), not a known-exact constant. When
#  `POD_sd` > 0, EACH individual's search target is instead a fresh
#  lognormal draw from POD_mean = `POD`, sd = `POD_sd` (.sample_POD(),
#  section 1b below) -- the same mean/CV-parameterised lognormal form
#  Extended_Model_Population_PC_MC.r's .sample_PC_mc() already uses for
#  partition-coefficient uncertainty, for the same reason: POD is a
#  positive-valued, typically right-skewed measured concentration, so a
#  lognormal can't draw a non-physical POD <= 0 the way a plain normal
#  could. This pairs POD measurement uncertainty with each individual's own
#  physiology draw, so the returned exposure distribution reflects BOTH
#  inter-individual variability AND POD uncertainty at once (a combined,
#  not separated, 1-D Monte Carlo -- same caveat as the file's own
#  variability/uncertainty note above). `POD_sd = 0` (the default) keeps the
#  old exact-POD behaviour.
#
#  Performance note: this runs a full nloptr search (up to ~1000 ODE solves
#  each, though a smooth 1-D search usually converges in far fewer) for
#  EVERY individual, so runtime scales with n roughly the same way
#  .run_population_extended()'s forward runs do, but multiplied by the
#  optimizer's iteration count per individual. Start with a small n (10-30)
#  to gauge runtime before scaling up.
#
#  Dependencies: Extended_Model_Population.r (and, transitively,
#  popgen_consts.R / popgen_port.R, estimated_parameters.csv), and the
#  nloptr package, must be available.
# =============================================================================

library(deSolve)
library(tidyverse)

#=========================
# 0. Load Extended_Model_Population.r's function definitions
#=========================
# Source everything up to (not including) its own "# --- Example usage"
# block, so we reuse .get_chemical_vars()/.generate_population()/
# .map_popgen_individual()/.gfr_age_factor()/.sample_Htc()/
# .create.params.population()/.create.inits()/.create.events()/.ode.func()/
# .check_required() without re-running its own 100-individual demo
# simulation as a side effect of loading this file.
if (!exists(".create.params.population")) {
  if (!file.exists("Extended_Model_Population.r"))
    stop("Reverse_Dosimetry_Stochastic.r requires Extended_Model_Population.r in the working directory")
  .rds_lines <- readLines("Extended_Model_Population.r")
  .rds_cut <- grep("^# --- Example usage", .rds_lines)[1]
  if (is.na(.rds_cut))
    stop("Could not find the '# --- Example usage' marker in Extended_Model_Population.r -- has the file structure changed?")
  .rds_tmp <- tempfile(fileext = ".R")
  writeLines(.rds_lines[seq_len(.rds_cut - 1)], .rds_tmp)
  source(.rds_tmp)
  rm(.rds_lines, .rds_cut, .rds_tmp)
}

#=========================
# 1. Per-individual reverse-dosimetry objective (copied logic from
#    reverse_dosimetry.r's .reverse_obj_func_models(), retargeted at one
#    popgen individual instead of a single fixed BW)
#=========================
#' Relative error between a candidate exposure and the target POD, for one
#' popgen individual
#'
#' Same log-dose parameterisation as reverse_dosimetry.r's
#' .reverse_obj_func_models(): the optimizer searches over `x` with the
#' actual exposure = exp(x) (keeps the search unconstrained while the dose
#' itself stays positive). exp(x) is always a per-kg-BW rate, converted to
#' an absolute amount using THIS individual's own BW
#' (`ind[["body_mass_kg"]]`) -- the same formula
#' Forward_Dosimetry_Stochastic.r's forward_dosimetry_stochastic() uses.
#' `user_input$exp_type` (fixed for the whole search, set by the caller)
#' decides how .create.events() applies the resulting `ingestion`: a one-off
#' dose ("pharmacokinetics") or a steady rate ("continuous") -- see the file
#' header's design note.
#'
#' @param x Log-dose (optimizer's search variable).
#' @param POD Target internal concentration.
#' @param compartment One of the names in `col_map` below (case-insensitive).
#' @param user_input Base user_input list (chemical/dosing/time fields,
#'   including `exp_type`); `ingestion` is overwritten here per candidate `x`.
#' @param ind One popgen individual (row of .generate_population()$pop).
#' @param chem list(variables, PC) from .get_chemical_vars(user_input$chemical).
#' @keywords internal
#'
.reverse_obj_func_population <- function(x, POD, compartment, user_input, ind, chem,
                                          solver = "lsodes", rtol = 1e-4, atol = 1e-4) {

  BW         <- ind[["body_mass_kg"]]
  time_scale <- user_input$time_scale

  # The day/time_scale rescaling factor (x365 for years, x30 for months, ...)
  # only makes sense for exp_type = "continuous": it converts a per-DAY rate
  # into "amount per one repeat interval" so the periodic .create.events()
  # "rep" re-application matches whatever time_scale the simulation runs in.
  # For exp_type = "pharmacokinetics", exp(x) is a single administered dose,
  # not a rate -- the same ng/kg dose must convert to the same absolute ug
  # amount regardless of which time_scale unit the simulation uses, so no
  # rescaling factor applies there (see Forward_Dosimetry_Stochastic.r's
  # .ingestion_for_BW() for the same fix).
  user_input$ingestion <- if (identical(user_input$exp_type, "pharmacokinetics")) {
    exp(x) * BW / 1000
  } else {
    switch(time_scale,
      "minutes" = (1/24/60) * exp(x) * BW / 1000,
      "hours"   = (1/24)    * exp(x) * BW / 1000,
      "days"    =             exp(x) * BW / 1000,
      "weeks"   = 7         * exp(x) * BW / 1000,
      "months"  = 30        * exp(x) * BW / 1000,
      "years"   = 365       * exp(x) * BW / 1000
    )
  }

  params <- .create.params.population(user_input, ind, chem)
  inits  <- .create.inits(params)
  events <- .create.events(params)

  solution <- as.data.frame(deSolve::ode(
    times = user_input$exposure_time, func = .ode.func, y = inits,
    parms = params, events = events, method = solver, rtol = rtol, atol = atol
  ))

  col_map <- c(serum = "CA", liver = "CL", adipose = "CAdi", brain = "CBra",
               gonads = "CGon", gut = "CGI", heart = "CHea", lung = "CLun",
               muscle = "CMus", skin = "CSki", kidney = "CKb")

  col <- col_map[tolower(compartment)]
  if (is.na(col)) stop("Invalid compartment: ", compartment)

  final_concentration <- tail(solution, 1)[[col]]
  rel_error <- abs((final_concentration - POD) / POD)
  return(rel_error)
}

#' Run the 1-D nloptr search for one popgen individual
#' (copied structure from reverse_dosimetry.r's .reverse_dosimetry_set_up(),
#' retargeted at .reverse_obj_func_population())
#'
#' @keywords internal
#'
.reverse_dosimetry_set_up_population <- function(user_input, ind, chem, POD, compartment,
                                                  optim_tol = 1e-3,
                                                  solver = "lsodes", rtol = 1e-4, atol = 1e-4) {

  opts <- list(
    algorithm   = "NLOPT_LN_SBPLX",
    xtol_rel    = optim_tol,
    xtol_abs    = optim_tol,
    ftol_rel    = optim_tol,
    ftol_abs    = optim_tol,
    maxeval     = 1000,
    print_level = 0
  )

  optimizer <- nloptr::nloptr(
    x0          = 1,
    eval_f      = .reverse_obj_func_population,
    lb          = -8,
    ub          = 8,
    opts        = opts,
    POD         = POD,
    compartment = compartment,
    user_input  = user_input,
    ind         = ind,
    chem        = chem,
    solver      = solver,
    rtol        = rtol,
    atol        = atol
  )

  list(
    id        = ind[["id"]],
    exposure  = exp(optimizer$solution),
    rel_error = optimizer$objective
  )
}

#' Draw one random POD from measurement uncertainty around a mean/sd
#'
#' Lognormal, parameterised by mean and CV (sd/mean) -- the same form
#' Extended_Model_Population_PC_MC.r's .sample_PC_mc() uses for partition-
#' coefficient uncertainty (see the file header's POD uncertainty note for
#' why lognormal, not normal). If you have raw replicate POD measurements
#' rather than a single sd, pass `POD_sd = sd(your_replicates)` and
#' `POD = mean(your_replicates)` to reverse_dosimetry_stochastic() --
#' this function itself only takes the summarised mean/sd, matching how
#' PC uncertainty is already consumed elsewhere in this project.
#'
#' @param POD_mean Point estimate (mean) of the POD.
#' @param POD_sd Standard deviation of the POD measurement/estimate.
#' @return One random POD draw (same units as POD_mean/POD_sd).
#' @keywords internal
#'
.sample_POD <- function(POD_mean, POD_sd) {
  cv <- POD_sd / POD_mean
  sigma <- sqrt(log(cv^2 + 1))
  mu <- log(POD_mean) - sigma^2 / 2
  exp(rnorm(1) * sigma + mu)
}

#=========================
# 2. Population reverse-dosimetry driver
#=========================
#' Population counterpart to reverse_dosimetry.r's .reverse_dosimetry()
#'
#' Same dosing/target arguments as .reverse_dosimetry() (chemical, duration,
#' time_scale, POD, compartment, optim_tol, solver, rtol, atol) -- except
#' there is no single BW argument, and admin_type/admin_dose/admin_time are
#' replaced by `exp_type` (see the file header's design note: every search
#' is now BW-scaled, whether it is a one-off dose or a continuous rate).
#' `n` virtual individuals are generated by popgen (age_range/bmi_range/
#' height_range/prob_male/ethnicity_probs/dataset/population_type, exactly
#' as .generate_population() requires), and the SAME nloptr search
#' .reverse_dosimetry() runs once is run once PER INDIVIDUAL, using that
#' individual's own popgen physiology and own BW.
#'
#' @param n,dataset,age_range,bmi_range,height_range,prob_male,ethnicity_probs,
#'   population_type,seed Population-generation inputs, passed straight to
#'   .generate_population() (see Extended_Model_Population.r for the exact
#'   meaning of each).
#' @param chemical,duration,time_scale,POD,compartment,optim_tol Same names
#'   and meaning as .reverse_dosimetry()'s. `POD` is the target internal
#'   concentration to reproduce at `time = duration`; `compartment` is one
#'   of "serum", "liver", "adipose", "brain", "gonads", "gut", "heart",
#'   "lung", "muscle", "skin", "kidney" (case-insensitive).
#' @param POD_sd Standard deviation of the POD measurement/estimate (default
#'   0 = exact POD, the old behaviour). When > 0, each individual searches
#'   against its own fresh lognormal(POD, POD_sd) draw instead of the exact
#'   POD -- see the file header's POD uncertainty note and .sample_POD().
#' @param exp_type "continuous" (default; searches for the steady per-kg-BW
#'   intake RATE that reproduces the POD) or "pharmacokinetics" (searches
#'   for a one-off per-kg-BW dose, given at `ingestion_time`, instead).
#' @param ingestion_time When the dose/rate search begins (default 0);
#'   passed straight through to .create.events() as `ingestion_time`.
#' @param solver,rtol,atol Passed through to deSolve::ode() inside each
#'   individual's nloptr search.
#' @return list(
#'   population = data.frame, one row per individual (popgen demographics,
#'     organ mass/flow, and the id used to join against `exposure_estimates`),
#'   exposure_estimates = data.frame, one row per individual: id, the
#'     estimated per-kg-BW exposure that reproduces the POD for that
#'     individual (same units .reverse_dosimetry()'s `exposure` is in for
#'     its oral/continuous case), rel_error (how exactly the optimizer hit
#'     the target), and the POD (the individual's own random draw, if
#'     POD_sd > 0 -- not just the mean)/compartment/time, for reference.
#' )
#' @export
#'
reverse_dosimetry_stochastic <- function(n, dataset, age_range, bmi_range, height_range,
                                          prob_male, ethnicity_probs, population_type, seed = NULL,
                                          chemical, duration, time_scale = "years",
                                          POD, POD_sd = 0, compartment, exp_type = "continuous",
                                          ingestion_time = 0,
                                          optim_tol = 1e-3,
                                          solver = "lsodes", rtol = 1e-4, atol = 1e-4) {

  exposure_time <- sort(unique(c(seq(0, floor(duration), 1), duration)))

  user_input <- list(
    n = n, dataset = dataset, age_range = age_range, bmi_range = bmi_range,
    height_range = height_range, prob_male = prob_male, ethnicity_probs = ethnicity_probs,
    population_type = population_type, seed = seed,
    exposure_time  = exposure_time,
    chemical       = chemical,
    ingestion      = 0,
    ingestion_time = ingestion_time,
    exp_type       = exp_type,
    time_scale     = time_scale,
    duration       = duration
  )

  .check_required(user_input, c(.required_population_fields,
                                 "chemical", "exp_type", "duration", "exposure_time"))

  gen <- .generate_population(user_input)
  pop <- gen$pop
  message(sprintf("Generated %d individuals (%d discarded by popgen's rejection rules)",
                   nrow(pop), gen$n_discarded))

  chem <- .get_chemical_vars(chemical)

  estimates_list <- vector("list", nrow(pop))
  for (i in seq_len(nrow(pop))) {
    ind <- pop[i, ]
    # per-individual POD draw when POD_sd > 0 (measurement/estimate
    # uncertainty), else the exact POD every time -- see the file header's
    # POD uncertainty note and .sample_POD()
    POD_i <- if (POD_sd > 0) .sample_POD(POD, POD_sd) else POD
    res <- .reverse_dosimetry_set_up_population(
      user_input = user_input, ind = ind, chem = chem,
      POD = POD_i, compartment = compartment,
      optim_tol = optim_tol,
      solver = solver, rtol = rtol, atol = atol
    )
    estimates_list[[i]] <- data.frame(
      id = res$id, exposure = res$exposure, rel_error = res$rel_error,
      POD = POD_i, compartment = compartment, time = duration
    )
    message(sprintf("  individual %d/%d (id=%s): POD=%.4g, exposure=%.6g, rel_error=%.4g",
                     i, nrow(pop), res$id, POD_i, res$exposure, res$rel_error))
  }

  exposure_estimates <- do.call(rbind, estimates_list)

  list(population = pop, exposure_estimates = exposure_estimates)
}

