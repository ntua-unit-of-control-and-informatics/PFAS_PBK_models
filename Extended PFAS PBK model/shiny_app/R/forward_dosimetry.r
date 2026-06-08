.run_extended <- function(user_input, solver = "lsodes", rtol = 1e-6, atol = 1e-6) {
  params   <- .create.params(user_input)
  inits    <- .create.inits(params)
  events   <- .create.events(params)
  solution <- as.data.frame(deSolve::ode(
    times  = user_input$exposure_time,
    func   = .ode.func,
    y      = inits,
    parms  = params,
    events = events,
    method = solver,
    rtol   = rtol,
    atol   = atol
  ))
  return(solution)
}

forward_dosimetry <- function(BW, duration, time_step, chemical, ingestion,
                              ingestion_time, admin_dose, admin_time,
                              admin_type, exp_type, time_scale) {

  exposure_time <- sort(unique(c(seq(0, duration, by = time_step), duration)))

  ingestion_input <- switch(time_scale,
    "minutes" = (1/24/60) * ingestion * BW / 1000,
    "hours"   = (1/24)    * ingestion * BW / 1000,
    "days"    =             ingestion * BW / 1000,
    "weeks"   = 7         * ingestion * BW / 1000,
    "months"  = 30        * ingestion * BW / 1000,
    "years"   = 365       * ingestion * BW / 1000
  )

  user_input <- list(
    BW             = BW,
    exposure_time  = exposure_time,
    chemical       = chemical,
    ingestion      = ingestion_input,
    ingestion_time = ingestion_time,
    admin_dose     = admin_dose,
    admin_time     = admin_time,
    admin_type     = admin_type,
    exp_type       = exp_type,
    time_scale     = time_scale
  )

  return(.run_extended(user_input))
}
