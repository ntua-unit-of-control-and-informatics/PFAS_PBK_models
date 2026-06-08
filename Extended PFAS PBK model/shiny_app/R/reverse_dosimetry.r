#' Run extended PBK model (internal)
#'
#' Sets up and solves the extended model via `deSolve::ode()`.
#'
#' @param user_input List with fields required by the model (e.g., `exposure_time`, dosing, BW, etc.).
#' @param solver Character or function passed to `deSolve::ode(method=)`. Defaults taken from calling scope.
#' @param rtol Numeric relative tolerance for the ODE solver.
#' @param atol Numeric absolute tolerance for the ODE solver.
#' @return A `data.frame` containing the time grid and state variables.
#' @seealso create.params, create.inits, create.events, ode.func
#' @keywords internal
#' @export
#' 




.forward_dosimetry <- function(chemical, ingestion, BW, duration, time_scale = "years",
                              solver = "lsodes", rtol = 1e-4, atol = 1e-4){ 
  vars <- list(ingestion = ingestion, BW = BW, duration = duration)
  mins <- c(ingestion = 0, BW = 0, duration = 1)
  
  for (name in names(vars)) {
    x <- vars[[name]]
    if (!is.numeric(x)) stop(name, " should be numeric")
    if (x < mins[name]) {
      stop(name, " should be greater than ", mins[name])
    }
  }
  
  exposure_time <- sort(unique(c(seq(0, floor(duration), 1), duration)))  #in same time units as time scale
  
 if (time_scale == "minutes"){
   ingestion_input <- (1/24/60)*ingestion*BW/1000 #ug/minute
  }else if (time_scale == "hours"){
    ingestion_input <- (1/24)*ingestion*BW/1000 #ug/hour
  }else if (time_scale == "days"){
    ingestion_input <- ingestion*BW/1000 #ug/day
  }else if (time_scale == "weeks"){
    ingestion_input <- 7*ingestion*BW/1000 #ug/week
  }else if (time_scale == "months"){
    ingestion_input <- 30*ingestion*BW/1000 #ug/month
  }else if (time_scale == "years"){
    ingestion_input <- 365*ingestion*BW/1000 #ug/year
  }
  
  user_input <- list(
    'BW' = BW,
    "exposure_time" = exposure_time,
    'chemical' = chemical,
    "ingestion" = ingestion_input,
    "ingestion_time" = 0,
    "admin_dose" = 0,
    "admin_time" = 0,
    "admin_type" = "oral",
    "exp_type" = "continuous",
    "time_scale" = time_scale
  )

  results <- .run_extended(user_input = user_input,solver = solver, rtol = rtol, atol = atol)
  return(results)
}


changing_exp_forward_dosimetry<- function(chemical, ingestion,ingestion_time, BW, duration, time_scale = "years",
                              solver = "lsodes", rtol = 1e-4, atol = 1e-4){ 
   vars <- list(ingestion = ingestion, ingestion_time = ingestion_time, BW = BW, duration = duration)
  mins <- c(ingestion = 0, ingestion_time = 0, BW = 0, duration = 1)
  
  for (name in names(vars)) {
    x <- vars[[name]]
    if (!is.numeric(x)) stop(name, " should be numeric")
    # Use any() to handle both scalars and vectors safely
    if (any(x < mins[name], na.rm = TRUE)) {
      stop(name, " should be greater than ", mins[name])
    }
  }
  
  exposure_time <- sort(unique(c(seq(0, floor(duration), 1), duration)))  #in same time units as time_scale
  
 if (time_scale == "minutes"){
   ingestion_input <- (1/24/60)*ingestion*BW/1000 #ug/minute
  }else if (time_scale == "hours"){
   ingestion_input <- (1/24)*ingestion*BW/1000 #ug/hour
  }else if (time_scale == "days"){
   ingestion_input <- ingestion*BW/1000 #ug/day
  }else if (time_scale == "weeks"){
   ingestion_input <- 7*ingestion*BW/1000 #ug/week
  }else if (time_scale == "months"){
   ingestion_input <- 30*ingestion*BW/1000 #ug/month
  }else if (time_scale == "years"){
   ingestion_input <- 365*ingestion*BW/1000 #ug/year
  }
  
  user_input <- list(
    'BW' = BW,
    "exposure_time" = exposure_time,
    'chemical' = chemical,
    "ingestion" = ingestion_input,
    "ingestion_time" = ingestion_time,
    "admin_dose" = 0,
    "admin_time" = 0,
    "admin_type" = "oral",
    "exp_type" = "continuous",
    "time_scale" = time_scale
  )

  results <- .run_extended(user_input = user_input, solver = solver, rtol = rtol, atol = atol)
  return(results)
}









#' Reverse dosimetry
#'
#' Estimates continuous exposure levels that reproduce user-specified
#' points of departure (PODs) in serum or liver at selected time points,
#' using a physiologically based kinetic (PBK) model with Bayesian
#' model averaging (BMA).
#'
#' @param chemical Character scalar. Which chemical to run:
#'   \code{"PFOA"}, \code{"PFOS"}, \code{"PFBA"}, \code{"PFHxS"}, \code{"PFNA"},
#'   \code{"PFDA"},
#'   \code{"HFPO_DA"}, \code{"PFBS"} or \code{"PFHxS"} (case-insensitive).
#' . Matching is case-insensitive.
#' @param BW Numeric scalar. Body weight (e.g. in kg).
#' @param duration Numeric scalar. Total simulation time in the units specified
#'   by \code{time_scale}. Must satisfy \eqn{\ge 1}. The exposure time grid is
#'   internally generated as \code{seq(0, duration, 1)}.
#' @param time_scale Character string, one of \code{"minutes"}, \code{"hours"},
#'   \code{"days"}, \code{"weeks"}, \code{"months"} or \code{"years"}. Times in
#'   \code{bmd_df$time} are interpreted in these units.
#'   @param free Logical; if \code{TRUE}, free concentrations will be estimated.
#'   Default value is \code{FALSE}
#'   @param bmd_df A \code{data.frame} with exactly four columns named:
#'   \itemize{
#'   \item \code{POD} — numeric target internal concentration
#'   (point of departure).
#'   \item \code{compartment} — character; \code{"serum"} or \code{"liver"},
#'   \code{"adipose"},\code{"brain"}, \code{"gonads"},
#'   \code{"gut"},\code{"heart"},
#'   \code{"lung"},\code{"muscle"},\code{"skin"}
#'   or \code{"kidney"}, or \code{"other"} (case-insensitive).
#'           (case-insensitive).
#'     \item \code{key_event} — character label for the biological key event
#'           (free text, used for tracing).
#'     \item \code{time} — numeric time at which the POD should be met, in
#'           the units defined by \code{time_scale}.
#'   }
#'   Column names must match these exactly (order may differ).
#' @param isStochastic Logical. If \code{FALSE} (default), a single
#'   BMA-based exposure estimate is returned for each row of \code{bmd_df}
#'   (deterministic mode, using mean posterior model weights and no
#'   stochastic error).
#'   If \code{TRUE}, reverse dosimetry is performed under a stochastic BMA
#'   scheme that samples from the posterior model weights and a log-normal
#'   error model. In this mode, \code{isParallel} must be \code{TRUE};
#'   otherwise an error is raised.
#' @param isParallel Logical. If \code{TRUE}, Monte Carlo computations used by
#'   the stochastic reverse-dosimetry procedure are run in parallel (where
#'   implemented). For \code{isStochastic = TRUE}, \code{isParallel} must be
#'   \code{TRUE}. In deterministic mode (\code{isStochastic = FALSE}), this
#'   flag is ignored.
#' @param Nsamples Integer. Number of Monte Carlo / posterior samples used in
#'   the stochastic mode (\code{isStochastic = TRUE}). Ignored in deterministic
#'   mode.
#' @param n_cores Integer or \code{NULL}. Number of CPU cores to use when
#'   \code{isParallel = TRUE}. If \code{NULL}, a sensible default is chosen
#'   internally (typically based on \code{parallel::detectCores()}).
#' @param stat_sum Logical; controls the stochastic output format when
#'   \code{isStochastic = TRUE}.  
#'   If \code{FALSE}, the function returns all \code{Nsamples} reverse-dosimetry
#'   draws (one value per Monte Carlo sample) in the \code{exposure} and
#'   \code{rel_error} fields.  
#'   If \code{TRUE}, the function returns summary statistics of these draws
#'   (mean, sd, and selected quantiles) in \code{exposure} and
#'   \code{rel_error}. In deterministic mode (\code{isStochastic = FALSE}),
#'   \code{stat_sum} is ignored and single values are returned.
#' @param seed Integer. Random seed used when drawing Monte Carlo samples in
#'   stochastic mode. When \code{isParallel = TRUE}, this is also used to
#'   initialise the parallel random number streams.
#' @param optim_tol Numeric scalar. Tolerance passed to the numerical optimiser
#'   for the reverse-dosimetry objective (e.g. \code{1e-3}).
#' @param solver Character scalar. ODE solver to use internally
#'   (e.g. \code{"lsodes"}).
#' @param rtol,atol Numeric scalars. Relative and absolute tolerances for the
#'   numerical ODE solver.
#'
#' @details
#' For each requested chemical, the function:
#' \enumerate{
#'   \item Builds an internal \code{user_input} structure describing continuous
#'         oral exposure (no bolus dosing) over the specified \code{duration}.
#'   \item Simulates the PBK system until the requested times.
#'   \item For each row of \code{bmd_df}, searches for the continuous intake
#'         that yields the specified \code{POD} in the requested
#'         \code{compartment} at the given \code{time}.
#' }
#'
#' In deterministic mode (\code{isStochastic = FALSE}), a single exposure
#' estimate and its associated relative error are produced for each
#' POD/compartment/time combination. These are returned as scalar values in
#' the \code{exposure} and \code{rel_error} fields.
#'
#' In stochastic mode (\code{isStochastic = TRUE}), the reverse-dosimetry step
#' is repeated \code{Nsamples} times using draws from the posterior model
#' weights and the log-normal error distribution, yielding a distribution of
#' exposure estimates per POD/compartment/time combination:
#' \itemize{
#'   \item If \code{stat_sum = FALSE}, \code{exposure} and \code{rel_error}
#'         are numeric vectors containing all \code{Nsamples} draws.
#'   \item If \code{stat_sum = TRUE}, \code{exposure} and \code{rel_error}
#'         are named numeric vectors containing summary statistics (mean, sd,
#'         and quantiles such as p1, p5, p25, p50, p75, p95, p99), analogous
#'         to the forward-dosimetry summaries.
#' }
#'
#' Basic input checks include:
#' \itemize{
#'   \item \code{BW} must be numeric and \eqn{> 0}.
#'   \item \code{duration} must be numeric and \eqn{\ge 1}.
#'   \item \code{bmd_df} must have exactly the four required columns.
#'   \item \code{compartment} must be either \code{"serum"} or \code{"liver"}
#'         (case-insensitive).
#' }
#'
#' @return A named \code{list} with elements \code{PFOA} and \code{PFOS}.
#' Each element is either:
#' \itemize{
#'   \item \code{NULL}, if the corresponding chemical was not requested, or
#'   \item an object containing one reverse-dosimetry result per row of
#'         \code{bmd_df}, including at least:
#'         \itemize{
#'           \item \code{exposure} — estimated continuous intake(s) that achieve
#'                 the requested POD(s); either a single value (deterministic),
#'                 a vector of draws, or a vector of summary statistics,
#'                 depending on \code{isStochastic} and \code{stat_sum}.
#'           \item \code{rel_error} — relative error(s) between the achieved
#'                 internal concentrations and the target POD(s), in the same
#'                 format as \code{exposure}.
#'           \item \code{POD}, \code{compartment}, \code{key_event}, \code{time}
#'                 — the corresponding input POD definition.
#'         }
#' }
#'
#' @examples
#' \dontrun{
#' bmd_df <- data.frame(
#'   POD = c(5, 10),
#'   compartment = c("serum", "liver"),
#'   key_event = c("KE1", "KE2"),
#'   time = c(1, 2) # in time_scale units, e.g. years
#' )
#'
#' # Deterministic reverse dosimetry for PFOA
#' out_pfoa <- reverse_dosimetry(
#'   chemical = "PFOA",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df,
#'   solver = "lsodes", rtol = 1e-4, atol = 1e-4
#' )
#'
#' # Stochastic reverse dosimetry (full draws, requires isParallel = TRUE)
#' out_pfoa_stoch_draws <- reverse_dosimetry(
#'   chemical = "PFOA",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df,
#'   isStochastic = TRUE,
#'   isParallel = TRUE,
#'   Nsamples = 50,
#'   stat_sum = FALSE
#' )
#'
#' # Stochastic reverse dosimetry (summary statistics)
#' out_pfoa_stoch_sum <- reverse_dosimetry(
#'   chemical = "PFOA",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df,
#'   isStochastic = TRUE,
#'   isParallel = TRUE,
#'   Nsamples = 50,
#'   stat_sum = TRUE
#' )
#'
#' # Both chemicals
#' out_both <- reverse_dosimetry(
#'   chemical = "BOTH",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df
#' )
#' }
#'
#' @export



.reverse_obj_func_models <- function(x, POD, compartment, user_input,
                                     admin_type = "oral",
                                     solver = "lsodes", rtol = 1e-4, atol = 1e-4) {

  BW         <- user_input$BW
  time_scale <- user_input$time_scale

  if (admin_type == "oral") {
    user_input$ingestion <- switch(time_scale,
      "minutes" = (1/24/60) * exp(x) * BW / 1000,
      "hours"   = (1/24)    * exp(x) * BW / 1000,
      "days"    =             exp(x) * BW / 1000,
      "weeks"   = 7         * exp(x) * BW / 1000,
      "months"  = 30        * exp(x) * BW / 1000,
      "years"   = 365       * exp(x) * BW / 1000
    )
  } else {
    user_input$admin_dose <- exp(x)
    user_input$ingestion  <- 0
    user_input$admin_type <- "iv"
  }

  solution <- .run_extended(user_input, solver = solver, rtol = rtol, atol = atol)

  col_map <- c(serum = "CA", liver = "CL", adipose = "CAdi", brain = "CBra",
               gonads = "CGon", gut = "CGI", heart = "CHea", lung = "CLun",
               muscle = "CMus", skin = "CSki", kidney = "CKb")

  col <- col_map[tolower(compartment)]
  if (is.na(col)) stop("Invalid compartment: ", compartment)

  final_concentration <- tail(solution, 1)[[col]]
  rel_error <- abs((final_concentration - POD) / POD)
  return(rel_error)
}


.reverse_dosimetry_set_up <- function(user_input, bmd_df_row, admin_type = "oral",
                                      optim_tol = 1e-3, solver = "lsodes",
                                      rtol = 1e-4, atol = 1e-4) {

  POD         <- bmd_df_row$POD
  compartment <- bmd_df_row$compartment

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
    x0         = 1,
    eval_f     = .reverse_obj_func_models,
    lb         = -8,
    ub         = 8,
    opts       = opts,
    POD        = POD,
    compartment = compartment,
    user_input = user_input,
    admin_type = admin_type,
    solver     = solver,
    rtol       = rtol,
    atol       = atol
  )

  return(list(
    exposure    = exp(optimizer$solution),
    rel_error   = optimizer$objective,
    POD         = POD,
    compartment = compartment,
    key_event   = bmd_df_row$key_event,
    time        = bmd_df_row$time
  ))
}

#' @examples
#' \dontrun{
#' bmd_df <- data.frame(
#'   POD = c(5, 10),
#'   compartment = c("serum", "liver"),
#'   key_event = c("KE1", "KE2"),
#'   time = c(1, 2) # in time_scale units, e.g. years
#' )
#'
#' # Deterministic reverse dosimetry for PFOA
#' out_pfoa <- reverse_dosimetry(
#'   chemical = "PFOA",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df,
#'   solver = "lsodes", rtol = 1e-4, atol = 1e-4
#' )
#'
#' # Stochastic reverse dosimetry (full draws, requires isParallel = TRUE)
#' out_pfoa_stoch_draws <- reverse_dosimetry(
#'   chemical = "PFOA",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df,
#'   isStochastic = TRUE,
#'   isParallel = TRUE,
#'   Nsamples = 50,
#'   stat_sum = FALSE
#' )
#'
#' # Stochastic reverse dosimetry (summary statistics)
#' out_pfoa_stoch_sum <- reverse_dosimetry(
#'   chemical = "PFOA",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df,
#'   isStochastic = TRUE,
#'   isParallel = TRUE,
#'   Nsamples = 50,
#'   stat_sum = TRUE
#' )
#'
#' # Both chemicals
#' out_both <- reverse_dosimetry(
#'   chemical = "BOTH",
#'   BW = 70,
#'   duration = 5,
#'   time_scale = "years",
#'   bmd_df = bmd_df
#' )
#' }
#'
#' @export
.reverse_dosimetry <- function(chemical, BW, duration, time_scale = "years",
                               POD, compartment,
                               admin_type = "oral", admin_time = 0,
                               optim_tol = 1e-3,
                               solver = "lsodes", rtol = 1e-4, atol = 1e-4) {

  exposure_time <- sort(unique(c(seq(0, floor(duration), 1), duration)))

  user_input <- list(
    BW             = BW,
    exposure_time  = exposure_time,
    chemical       = chemical,
    ingestion      = 0,
    ingestion_time = 0,
    admin_dose     = 0,
    admin_time     = admin_time,
    admin_type     = admin_type,
    exp_type       = "continuous",
    time_scale     = time_scale
  )

  bmd_df_row <- list(
    POD         = POD,
    compartment = compartment,
    key_event   = "POD",
    time        = duration
  )

  results <- .reverse_dosimetry_set_up(
    user_input = user_input,
    bmd_df_row = bmd_df_row,
    admin_type = admin_type,
    optim_tol  = optim_tol,
    solver     = solver,
    rtol       = rtol,
    atol       = atol
  )

  return(results)
}
