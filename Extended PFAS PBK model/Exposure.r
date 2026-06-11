library(deSolve)
library(nloptr)
library(tidyverse)

source("Extended_Model.r")
.run_extended <- function(user_input, solver = "lsodes", rtol = 1e-5, atol = 1e-5){
  params<- .create.params(user_input)
  inits <- .create.inits(params)
  events <- .create.events(params)
  solution <-  as.data.frame(deSolve::ode(times = user_input$exposure_time,  func = .ode.func, 
                                        y = inits, parms = params,
                                        events = events,
                                       method= solver, rtol = rtol, atol = atol))
  return(solution)
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
  ingestion_input <- ingestion

  
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

.reverse_obj_func_models_exposure <- function(x, bio_con, bio_time, user_input, 
                                     solver = "lsodes",  rtol = 1e-4, atol =1e-4){
  
  BW <- user_input$BW
  exposure_time <- user_input$exposure_time
  chemical <- user_input$chemical
  time_scale <- user_input$time_scale


  # Update user_input with current exposure guess
  if (time_scale == "minutes"){
    user_input$ingestion <- (1/24/60)*exp(x)*BW/1000 #ug/minute
  }else if (time_scale == "hours"){
    user_input$ingestion <- (1/24)*exp(x)*BW/1000 #ug/hour
  }else if (time_scale == "days"){
    user_input$ingestion <- exp(x)*BW/1000 #ug/day
  }else if (time_scale == "weeks"){
    user_input$ingestion <- 7*exp(x)*BW/1000 #ug/week
  }else if (time_scale == "months"){
    user_input$ingestion <- 30*exp(x)*BW/1000 #ug/month
  }else if (time_scale == "years"){
    user_input$ingestion <- 365*exp(x)*BW/1000 #ug/year
  }

  
    solution <- .run_extended(user_input,
                              solver = solver,
                              rtol = rtol,
                              atol = atol)
    
    concentration_points<-slice(solution, match(bio_time, solution$time))$CA #replace Cserum with the relevant compartment

    rel_error <- mean(abs((concentration_points - bio_con) / bio_con))
  
  return(rel_error)
  
}

.reverse_dosimetry_set_up_exposure <- function(user_input, bio_mon,
                                  optim_tol =1e-3 , solver ="lsodes" , rtol=1e-4 , atol =1e-4 ) {
  
  BW <- user_input$BW
  exposure_time <- user_input$exposure_time
  chemical <- user_input$chemical
  bio_con <- bio_mon$bio_con
  bio_time <- bio_mon$bio_time
 
  #Initialise optimiser to NULL for better error handling later
  opts <- list(
    "algorithm"  = "NLOPT_LN_SBPLX",   # or "NLOPT_LN_NEWUOA" for smooth functions
    "xtol_rel"   = optim_tol,               # relative tolerance in x (e.g., stop if Δx/x < 0.1%)
    "xtol_abs"   = optim_tol,               # absolute tolerance in x (e.g., stop if Δx < 0.001)
    "ftol_rel"   = optim_tol,               # relative tolerance in objective (e.g., stop if error stops improving)
    "ftol_abs"   = optim_tol,               # absolute tolerance in objective (e.g., stop if error < 0.01)
    "maxeval"    = 1000,           # number of allowed function evaluations
    "print_level"= 0                  # change to 1 if you want progress output
  )

    optimizer <- nloptr::nloptr(
      x0 = rep(log(1),length(bio_con)),  # initial guess for log(exposure)
      eval_f = .reverse_obj_func_models_exposure,
      lb = rep(log(1e-4), length(bio_con)),  # Set lower bound to log(1e-10) to prevent negative exposures
      ub = rep(log(1e4), length(bio_con)),   # Set upper bound to log(1e10) to prevent excessively large exposures
      opts = opts,
      bio_con = bio_con,
      bio_time = bio_time,
      user_input = user_input,
      solver = solver, atol = atol, rtol = rtol
    )
    
exposure <- exp(optimizer$solution)
rel_error <- optimizer$objective

return(list("exposure" = exposure,
            "rel_error" = rel_error,
            "bio_con" = bio_con,
            "bio_time" = bio_time
            ))
}

.reverse_dosimetry_exposure<- function( chemical, BW, duration, time_scale = "years", bio_mon,
                              optim_tol = 1e-3, solver = "lsodes", rtol = 1e-4, atol = 1e-4){
  vars <- list(ingestion = 0, ingestion_time = 0, BW = BW, duration = duration)
  mins <- c(ingestion = 0, ingestion_time = 0, BW = 0, duration = 1)
  
  for (name in names(vars)) {
    x <- vars[[name]]
    if (!is.numeric(x)) stop(name, " should be numeric")
    if (x < mins[name]) {
      stop(name, " should be greater than ", mins[name])
    }
  }

  required_cols <- c("bio_con","bio_time")
  
  # if (ncol(bio_mon) != length(required_cols)) {
  #   stop("bio_mon should have exactly 3 columns")
  # }
  
  if (!setequal(colnames(bio_mon), required_cols)) {
    stop("The colnames of bio_mon should be exactly: 'bio_con', 'bio_time'")
  }
  

  exposure_time <- sort(unique(c(seq(0, floor(duration), 0.1), duration))) #in same time units as time scale

  user_input <- list(
    'BW' = BW,
    "exposure_time" = exposure_time,
    'chemical' = chemical,
    "ingestion_time" = c(0, bio_mon$bio_time[1:(length(bio_mon$bio_time)-1)]) , 
    "ingestion" = numeric(length(bio_mon$bio_time)), 
    "admin_dose" = 0,
    "admin_time" = 0,
    "admin_type" = "oral",
    "exp_type" = "continuous",
    "time_scale" = time_scale
  )
  # Check here how you can take a case insensitive input
    print("Starting optimization...")
    results <- as.list(.reverse_dosimetry_set_up_exposure(user_input = user_input, bio_mon,
                                          optim_tol = optim_tol, solver = solver, rtol = rtol, atol = atol))
  
  
  return(results)
}

