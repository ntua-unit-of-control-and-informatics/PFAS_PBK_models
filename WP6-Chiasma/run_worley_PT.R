user_input <- list(
  'BW' = 70,
  "exposure_time" = seq(0,30,1),
  "ingestion" = 10,
  "ingestion_time" = 0,
  "admin_dose" = 0,
  "admin_time" = 0,
  "admin_type" = "oral",
  "exp_type" = "continuous",
  "time_scale" = "years"
)

run_worley <- function(user_input, solver = "lsodes", rtol = 1e-7, atol = 1e-7){
  params_worley <- create.params.worley(user_input)
  inits_worley <- create.inits.worley(params_worley)
  events_worley <- create.events.worley(params_worley)
  solution_worley <-  as.data.frame(deSolve::ode(times = user_input$exposure_time,  
                                                 func = ode.func.worley, 
                                                 y = inits_worley, parms = params_worley,
                                                 events = events_worley, 
                                                 method= solver, rtol = rtol, atol = atol))
  
  return(solution_worley)
}