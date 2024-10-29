# Set up simulations for the 1st case, i.e. kudo (2007) high dose, tissues
BW <- 0.29  # body weight (kg)
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[1]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time=seq(0,2,0.1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-05, atol = 1e-05))

preds_kudo_high <-  solution[, c("time","Cblood","Cliver","Ckidney", "Ccarcass","Clungs", 
                                 "Cspleen", "Cheart","Cbrain", "Cgonads", "Cstomach", 
                                 "Cintestine")]
###############################################################################
# Set up simulations for the 2nd case, i.e. kudo (2007) low dose, tissues
BW <- 0.29  # body weight (kg)
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[2]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,2,0.1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kudo_low <- solution[, c("time","Cblood","Cliver","Ckidney", "Ccarcass","Clungs", 
                               "Cspleen", "Cheart","Cbrain", "Cgonads", "Cstomach", 
                               "Cintestine")]
############################################################################
# Set up simulations for the 3rd case, i.e. kim (2016) IV male tissues
BW <- 0.25 #kg, from Kim et al. 2018
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[3]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_IV_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

#########################################################################
# Set up simulations for the 4th case, i.e. kim (2016) ORAL male tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[4]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_OR_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

###################################################################
# Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[5]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_IV_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

#################################################################################
# Set up simulations for the 6th case, i.e. kim (2016) IV female tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[6]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

###############################################################################################
# Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[7]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,864,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]

#############################################################################################
# Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
BW <- 0.2  # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[8]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,24,0.1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]

######################################################################################
# Set up simulations for the 9th case, i.e. Kim (2016) ORAL male blood
BW <- 0.25  #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[9]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,288,1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_OR_Mblood <-  solution[, c("time", "Cplasma")]

######################################################################################
# Set up simulations for the 10th case, i.e. Kim (2016) IV male blood
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[10]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs

sample_time=c(0, 5/60, seq(1,288,1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_IV_Mblood <-  solution[, c("time", "Cplasma")]

##########################
#-------------------------
# Kemper 2003 (Worley)
#-------------------------
##########################
# Set up simulations for the 11th case, i.e.Kemper 2003 (Worley) ORAL female LOW

sex <- "F"
BW <- 0.2 #kg
sample_time <- seq(0,192,1)
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[11]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Furine_low <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Ffeces_low <-  solution[, c("time", "Mfeces")]

##########################################################################################
# Set up simulations for the 12th case, i.e.Kemper 2003 (Worley) ORAL female  MEDIUM
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[12]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Furine_med <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Ffeces_med <-  solution[, c("time", "Mfeces")]

##########################################################################################
# Set up simulations for the 13th case, i.e.Kemper 2003 (Worley) ORAL female HIGH
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[13]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Furine_high <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Ffeces_high <-  solution[, c("time", "Mfeces")]

###################################################################################################
# Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL male LOW

sex <- "M"
BW <- 0.3 #kg
sample_time <- seq(0,673,1)
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[14]], variable_params)
inits <- create.inits (params)
events <- create.events(params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Murine_low <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_low <-  solution[, c("time", "Mfeces")]

####################################################################################################
# Set up simulations for the 15th case, i.e.Kemper 2003 (Worley) ORAL male MEDIUM
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[15]], variable_params)
inits <- create.inits (params)
events <- create.events(params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Murine_med <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_med <-  solution[, c("time", "Mfeces")]

######################################################################################################
# Set up simulations for the 16th case, i.e.Kemper 2003 (Worley) ORAL male  HIGH
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[16]], variable_params)
inits <- create.inits (params)
events <- create.events(params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Murine_high <-  solution[, c("time", "Murine")]
preds_Kemp_OR_Mfeces_high <-  solution[, c("time", "Mfeces")]

#####################################################################################################
# Set up simulations for the 17th case, i.e. Dzierlenga 2021, IV male serum
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[17]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.083, 0.25, 0.5, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_IV_Mserum <-  solution[, c("time", "Cplasma")]

##############################################################################################
# Set up simulations for the 18th case, i.e. Dzierlenga 2021, ORAL male serum low
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[18]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Mserum_low <-  solution[, c("time", "Cplasma")]

########################################################################################################
# Set up simulations for the 19th case, i.e. Dzierlenga 2021, ORAL male serum medium
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[19]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Mserum_medium <-  solution[, c("time", "Cplasma")]

##############################################################################################
# Set up simulations for the 20th case, i.e. Dzierlenga 2021, ORAL male serum high
BW <- 0.3   # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[20]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, 1, 3, 6, seq(12, 1200, 4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Mserum_high <-  solution[, c("time", "Cplasma")]

###############################################################################################
# Set up simulations for the 21th case, i.e. Dzierlenga 2021, IV female serum
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[21]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.083, 0.25, seq(0.5, 192, 0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_IV_Fserum <-  solution[, c("time", "Cplasma")]

#########################################################################################
# Set up simulations for the 22th case, i.e. Dzierlenga 2021, ORAL female serum low
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[22]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0, 96, 0.25)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]

################################################################################################
# Set up simulations for the 23th case, i.e. Dzierlenga 2021, ORAL female serum medium
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[23]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 0.25, seq(1, 192, 0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Fserum_medium <-  solution[, c("time", "Cplasma")]

########################################################################################
# Set up simulations for the 24th case, i.e. Dzierlenga 2021, ORAL female serum high
BW <- 0.2    # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "F"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[24]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time <- seq(0, 96, 0.25)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_dzi_OR_Fserum_high <-  solution[, c("time", "Cplasma")]

##########################################################################################
# Set up simulations for the 25th case, i.e. Kim (2016) ORAL female blood
BW <- 0.25  #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[25]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= c(seq(0, 1.2, 0.2), seq(1.5,24,0.5))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_OR_Fblood <-  solution[, c("time", "Cplasma")]

########################################################################################
# Set up simulations for the 26th case, i.e. Kim (2016) IV female blood
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[26]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0, 25, 0.5)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]

###############################################################################################
# Set up simulations for the 27th case, i.e. Gustafsson (2022) oral male blood
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M" 
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[27]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0,48,0.2)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_gus_OR_Mblood <-  solution[, c("time", "Cplasma")]

##########################################################################################
# Set up simulations for the 28st case, i.e. Gustafsson (2022) Inhalation male tissues
BW <- 0.5125  #kg, from Gustafsson et al., 2022
sex <- "M"
variable_params <- create_variable_params(BW,sex,estimated_params)
params <- c(fixed_params[[28]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= seq(0,48,0.2)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_gus_OR_Mtissues <-  solution[, c("time", "CalveolarLF","Cliver", "Clungtissue", "Ckidney")]
