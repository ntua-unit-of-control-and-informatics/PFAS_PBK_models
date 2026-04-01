
# Set up simulations for the 1st case, i.e. kudo (2007) high dose, tissues
BW <- 0.29  # body weight (kg)
sex <- "M" 
variable_params <- create_variable_params(test_pars,BW,sex, fixed_params[[1]])
params <- c(fixed_params[[1]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time=c(0,  1e-5, 1e-2, seq(0.1,2,0.1))
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-4, atol = 1e-4))

preds_kudo_high <-  solution[, c("time","Cblood","Cliver","Ckidney", "Ccarcass","Clungs", 
                                 "Cspleen", "Cheart","Cbrain", "Cgonads", "Cstomach", 
                                 "Cintestine")]
preds_kudo_high[2,3:dim(preds_kudo_high)[2]] <- 0 
###############################################################################
# Set up simulations for the 2nd case, i.e. kudo (2007) low dose, tissues
BW <- 0.29  # body weight (kg)
sex <- "M" 
variable_params <- create_variable_params(test_pars,BW,sex, fixed_params[[2]])
params <- c(fixed_params[[2]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time= c(0, 1e-10, 1e-5, 1e-2,  seq(0.1,2,0.1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kudo_low <- solution[, c("time","Cblood","Cliver","Ckidney", "Ccarcass","Clungs", 
                               "Cspleen", "Cheart","Cbrain", "Cgonads", "Cstomach", 
                               "Cintestine")]
preds_kudo_low[2,3:dim(preds_kudo_low)[2]] <- 0 

############################################################################
# Set up simulations for the 3rd case, i.e. kim (2016) IV male tissues
BW <- 0.25 #kg, from Kim et al. 2018
variable_params <- create_variable_params(test_pars,BW,sex, fixed_params[[3]])
params <- c(fixed_params[[3]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time= c(0, 1e-10, 1e-5, 1e-2, seq(1,288,1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_IV_Mtissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]
preds_kim_IV_Mtissues[2,2:dim(preds_kim_IV_Mtissues)[2]] <- 0 

#########################################################################
# Set up simulations for the 4th case, i.e. kim (2016) ORAL male tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "M" 
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[4]])
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
preds_kim_OR_Mtissues[2,2:dim(preds_kim_OR_Mtissues)[2]] <- 0 

###################################################################
# Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[5]])
params <- c(fixed_params[[5]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time= c(0, 1e-10, 1e-5, 1e-2,  seq(1,24,1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_IV_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]
preds_kim_IV_Ftissues[2,2:dim(preds_kim_IV_Ftissues)[2]] <- 0 

#################################################################################
# Set up simulations for the 6th case, i.e. kim (2016) IV female tissues
BW <- 0.25 #kg, from Kim et al. 2018
sex <- "F" 
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[6]])
params <- c(fixed_params[[6]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time= c(0, 1e-10, 1e-5, 1e-2,  seq(1,24,1))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_kim_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]
preds_kim_OR_Ftissues[2,2:dim(preds_kim_OR_Ftissues)[2]] <- 0 

###############################################################################################
# Set up simulations for the 7th case, i.e. Dzierlenga (2021) ORAL male tissues
BW <- 0.3  # body weight (kg) not reported, based on 8 week male rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
sex <- "M" 
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[7]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[8]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[9]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[10]])
params <- c(fixed_params[[10]], variable_params)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs

sample_time= c(0, 1e-10, 1e-5, 1e-2, 5/60, seq(1,288,1))

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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[11]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[12]])
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
sample_time <- seq(0,672,1)
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[13]])
params <- c(fixed_params[[13]], variable_params)
events <- create.events(params)
inits <- create.inits (params)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

preds_Kemp_OR_Furine_high <-  solution[1:169, c("time", "Murine")]
preds_Kemp_OR_Ffeces_high <-  solution[, c("time", "Mfeces")]

###################################################################################################
# Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL male LOW

sex <- "M"
BW <- 0.3 #kg
sample_time <- seq(0,673,1)
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[14]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[15]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[16]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[17]])
params <- c(fixed_params[[17]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 1e-10, 1e-5, 1e-2,  0.083, 0.25, 0.5, 1, 3, 6, seq(12, 1200, 4))

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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[18]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[19]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[20]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[21]])
params <- c(fixed_params[[21]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- c(0, 1e-10, 1e-5, 1e-2, 0.083, 0.25, seq(0.5, 192, 0.5))

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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[22]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[23]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[24]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[25]])
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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[26]])
params <- c(fixed_params[[26]], variable_params)
inits <- create.inits(params)
events <- create.events(params)
sample_time= c(0, 1e-10, 1e-5, 1e-2,  seq(0.5, 25, 0.5))

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
variable_params <- create_variable_params(test_pars, BW,sex, fixed_params[[27]])
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
preds_gus_OR_Mtissues <-  solution[, c("time", "CBALF","Cliver_bled", "Clungtissue_bled", "Ckidney_bled")]
##########################################################################################



#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets
preds_kudo_high[,2:dim(preds_kudo_high)[2]] <- preds_kudo_high[,2:dim(preds_kudo_high)[2]] /1000 
preds_kudo_high <- preds_kudo_high[2:dim(preds_kudo_high)[1],]
preds_kudo_low[,2:dim(preds_kudo_low)[2]] <- preds_kudo_low[,2:dim(preds_kudo_low)[2]] /1000 
preds_kudo_low <- preds_kudo_low[2:dim(preds_kudo_low)[1],]
preds_kim_IV_Mtissues[,2:dim(preds_kim_IV_Mtissues)[2]] <- preds_kim_IV_Mtissues[,2:dim(preds_kim_IV_Mtissues)[2]] /1000
preds_kim_IV_Mtissues <- preds_kim_IV_Mtissues[2:dim(preds_kim_IV_Mtissues)[1],]
preds_kim_OR_Mtissues[,2:dim(preds_kim_OR_Mtissues)[2]] <- preds_kim_OR_Mtissues[,2:dim(preds_kim_OR_Mtissues)[2]] /1000
preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] <- preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] /1000
preds_kim_IV_Ftissues <- preds_kim_IV_Ftissues[2:dim(preds_kim_IV_Ftissues)[1],]
preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] <- preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] /1000
preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] <- preds_dzi_OR_Mtissues[,2:dim(preds_dzi_OR_Mtissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000
preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] <- preds_kim_OR_Mblood[,2:dim(preds_kim_OR_Mblood)[2]] /1000
preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] <- preds_kim_IV_Mblood[,2:dim(preds_kim_IV_Mblood)[2]] /1000
preds_kim_IV_Mblood <- preds_kim_IV_Mblood[2:dim(preds_kim_IV_Mblood)[1],]
preds_Kemp_OR_Ffeces_high[,2:dim(preds_Kemp_OR_Ffeces_high)[2]] <- preds_Kemp_OR_Ffeces_high[,2:dim(preds_Kemp_OR_Ffeces_high)[2]] /1000
preds_Kemp_OR_Mfeces_high[,2:dim(preds_Kemp_OR_Mfeces_high)[2]] <- preds_Kemp_OR_Mfeces_high[,2:dim(preds_Kemp_OR_Mfeces_high)[2]] /1000
preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] <- preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] /1000
preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] <- preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] /1000
preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] <- preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] /1000
preds_Kemp_OR_Murine_low[,2:dim(preds_Kemp_OR_Murine_low)[2]] <- preds_Kemp_OR_Murine_low[,2:dim(preds_Kemp_OR_Murine_low)[2]] /1000
preds_Kemp_OR_Murine_med[,2:dim(preds_Kemp_OR_Murine_med)[2]] <- preds_Kemp_OR_Murine_med[,2:dim(preds_Kemp_OR_Murine_med)[2]] /1000
preds_Kemp_OR_Murine_high[,2:dim(preds_Kemp_OR_Murine_high)[2]] <- preds_Kemp_OR_Murine_high[,2:dim(preds_Kemp_OR_Murine_high)[2]] /1000
preds_dzi_IV_Mserum[,2:dim(preds_dzi_IV_Mserum)[2]] <- preds_dzi_IV_Mserum[,2:dim(preds_dzi_IV_Mserum)[2]] /1000
preds_dzi_IV_Mserum <- preds_dzi_IV_Mserum[2:dim(preds_dzi_IV_Mserum)[1],]
preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] <- preds_dzi_OR_Mserum_low[,2:dim(preds_dzi_OR_Mserum_low)[2]] /1000
preds_dzi_OR_Mserum_medium[,2:dim(preds_dzi_OR_Mserum_medium)[2]] <- preds_dzi_OR_Mserum_medium[,2:dim(preds_dzi_OR_Mserum_medium)[2]] /1000
preds_dzi_OR_Mserum_high[,2:dim(preds_dzi_OR_Mserum_high)[2]] <- preds_dzi_OR_Mserum_high[,2:dim(preds_dzi_OR_Mserum_high)[2]] /1000
preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] <- preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] /1000
preds_dzi_IV_Fserum <- preds_dzi_IV_Fserum[2:dim(preds_dzi_IV_Fserum)[1],]
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] <- preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] /1000
preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] <- preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] /1000
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000
preds_kim_IV_Fblood <- preds_kim_IV_Fblood[2:dim(preds_kim_IV_Fblood)[1],]
preds_Kemp_OR_Ffeces_med[,2:dim(preds_Kemp_OR_Ffeces_med)[2]] <- preds_Kemp_OR_Ffeces_med[,2:dim(preds_Kemp_OR_Ffeces_med)[2]] /1000
preds_Kemp_OR_Mfeces_med[,2:dim(preds_Kemp_OR_Mfeces_med)[2]] <- preds_Kemp_OR_Mfeces_med[,2:dim(preds_Kemp_OR_Mfeces_med)[2]] /1000
preds_Kemp_OR_Ffeces_low[,2:dim(preds_Kemp_OR_Ffeces_low)[2]] <- preds_Kemp_OR_Ffeces_low[,2:dim(preds_Kemp_OR_Ffeces_low)[2]] /1000
preds_Kemp_OR_Mfeces_low[,2:dim(preds_Kemp_OR_Mfeces_low)[2]] <- preds_Kemp_OR_Mfeces_low[,2:dim(preds_Kemp_OR_Mfeces_low)[2]] /1000
preds_gus_OR_Mblood[,2:dim(preds_gus_OR_Mblood)[2]] <- preds_gus_OR_Mblood[,2:dim(preds_gus_OR_Mblood)[2]] /1000
preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] <- preds_gus_OR_Mtissues[,2:dim(preds_gus_OR_Mtissues)[2]] /1000

# ######################################################################################


# Convert Kudo High dose from long to wide format using reshape
experiment1 <- list()
experiment1[[1]] <- reshape(kudo_high_dose[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment1[[1]]) <- c("Time",kudo_high_dose$Tissue )

# Convert Kudo Low dose from long to wide format using reshape
experiment2 <- list()
experiment2[[1]] <- reshape(kudo_low_dose[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment2[[1]]) <- c("Time",kudo_low_dose$Tissue )

# Convert Kim IV Male tissues from long to wide format using reshape
experiment3 <- list()
experiment3[[1]] <- reshape(kim_IV_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment3[[1]]) <- c("Time",kim_IV_Mtissues$Tissue )

# Convert Kim ORAL Male tissues from long to wide format using reshape
experiment4 <- list()
experiment4[[1]] <- reshape(kim_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment4[[1]]) <- c("Time",kim_OR_Mtissues$Tissue )

# Convert Kim IV female tissues from long to wide format using reshape
experiment5 <- list()
experiment5[[1]] <- reshape(kim_IV_Ftissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment5[[1]]) <- c("Time",kim_IV_Ftissues$Tissue )

# Convert Kim ORAL female tissues from long to wide format using reshape
experiment6 <- list()
experiment6[[1]] <- reshape(kim_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment6[[1]]) <- c("Time",kim_OR_Ftissues$Tissue )

# Convert Dzierlenga ORAL male tissues from long to wide format using reshape
experiment7 <- list()
experiment7[[1]] <- reshape(dzi_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microM")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment7[[1]]) <- c("Time",unique(dzi_OR_Mtissues$Tissue))


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- list()
experiment8[[1]] <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                              "Concentration_microM")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8[[1]]) <- c("Time",unique(dzi_OR_Ftissues$Tissue))

experiment9 <- list()
# Convert Kim ORAL male blood from long to wide format using reshape
experiment9[[1]] <- reshape(kim_OR_Mblood[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9[[1]]) <- c("Time",unique(kim_OR_Mblood$Tissue))


# Convert Kim IV male blood from long to wide format using reshape
experiment9[[2]] <- reshape(kim_IV_Mblood[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment9[[2]]) <- c("Time",unique(kim_IV_Mblood$Tissue))


experiment10 <- list()
# Convert Kemper ORAL female urine low from long to wide format using reshape
experiment10[[1]] <- reshape(Kemp_OR_Furine_low [c("Tissue" ,"Time_h", 
                                                   "Cum_dose_%")], 
                             idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[1]]) <- c("Time",unique(Kemp_OR_Furine_low$Tissue))
experiment10[[1]]$Urine = (experiment10[[1]]$Urine/100)*0.2*1


# Convert Kemper ORAL female feces low from long to wide format using reshape
experiment10[[2]] <- reshape(Kemp_OR_Ffeces_low[c("Tissue" ,"Time_h", 
                                                  "Cum_dose_%")], 
                             idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[2]] ) <- c("Time",unique(Kemp_OR_Ffeces_low$Tissue))
experiment10[[2]]$Feces = (experiment10[[2]] $Feces/100)*0.2*1

# Convert Kemper ORAL female urine med from long to wide format using reshape
experiment10[[3]]  <- reshape(Kemp_OR_Furine_med[c("Tissue" ,"Time_h", 
                                                   "Cum_dose_%")], 
                              idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[3]]  ) <- c("Time",unique(Kemp_OR_Furine_med$Tissue))
experiment10[[3]]$Urine = (experiment10[[3]]  $Urine/100)*0.2*5


# Convert Kemper ORAL female feces medium from long to wide format using reshape
experiment10[[4]]   <- reshape(Kemp_OR_Ffeces_med[c("Tissue" ,"Time_h", 
                                                    "Cum_dose_%")], 
                               idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[4]]  ) <- c("Time",unique(Kemp_OR_Ffeces_med$Tissue))
experiment10[[4]]$Feces = (experiment10[[4]]$Feces/100)*0.2*5

# Convert Kemper ORAL female urine high from long to wide format using reshape
experiment10[[5]] <- reshape(Kemp_OR_Furine_high[c("Tissue" ,"Time_h", 
                                                   "Cum_dose_%")], 
                             idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[5]] ) <- c("Time",unique(Kemp_OR_Furine_high$Tissue))
experiment10[[5]]$Urine = (experiment10[[5]]$Urine/100)*0.2*25

# Convert Kemper ORAL female feces from long to wide format using reshape
experiment10[[6]]  <- reshape(Kemp_OR_Ffeces_high[c("Tissue" ,"Time_h", 
                                                    "Cum_dose_%")], 
                              idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[6]] ) <- c("Time",unique(Kemp_OR_Ffeces_high$Tissue))
experiment10[[6]]$Feces = (experiment10[[6]]$Feces/100)*0.2*25

# Convert Kemper ORAL male urine low from long to wide format using reshape
experiment10[[7]] <- reshape(Kemp_OR_Murine_low [c("Tissue" ,"Time_h", 
                                                   "Cum_dose_%")], 
                             idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[7]]) <- c("Time",unique(Kemp_OR_Murine_low$Tissue))
experiment10[[7]]$Urine = (experiment10[[7]]$Urine/100)*0.3*1


# Convert Kemper ORAL male feces low from long to wide format using reshape
experiment10[[8]] <- reshape(Kemp_OR_Mfeces_low[c("Tissue" ,"Time_h", 
                                                  "Cum_dose_%")], 
                             idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[8]]) <- c("Time",unique(Kemp_OR_Mfeces_low$Tissue))
experiment10[[8]]$Feces = (experiment10[[8]]$Feces/100)*0.3*1

# Convert Kemper ORAL male urine med from long to wide format using reshape
experiment10[[9]] <- reshape(Kemp_OR_Murine_med[c("Tissue" ,"Time_h", 
                                                  "Cum_dose_%")], 
                             idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[9]]) <- c("Time",unique(Kemp_OR_Murine_med$Tissue))
experiment10[[9]]$Urine = (experiment10[[9]]$Urine/100)*0.3*5


# Convert Kemper ORAL male feces medium from long to wide format using reshape
experiment10[[10]] <- reshape(Kemp_OR_Mfeces_med[c("Tissue" ,"Time_h", 
                                                   "Cum_dose_%")], 
                              idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[10]]) <- c("Time",unique(Kemp_OR_Mfeces_med$Tissue))
experiment10[[10]]$Feces = (experiment10[[10]]$Feces/100)*0.3*5


# Convert Kemper ORAL male urine high from long to wide format using reshape
experiment10[[11]] <- reshape(Kemp_OR_Murine_high[c("Tissue" ,"Time_h", 
                                                    "Cum_dose_%")], 
                              idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[11]]) <- c("Time",unique(Kemp_OR_Furine_high$Tissue))
experiment10[[11]]$Urine = (experiment10[[11]]$Urine/100)*0.3*25


# Convert Kemper ORAL male feces from long to wide format using reshape
experiment10[[12]] <- reshape(Kemp_OR_Mfeces_high[c("Tissue" ,"Time_h", 
                                                    "Cum_dose_%")], 
                              idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment10[[12]]) <- c("Time",unique(Kemp_OR_Mfeces_high$Tissue))
experiment10[[12]]$Feces = (experiment10[[12]]$Feces/100)*0.3*25

experiment11 <- list()
# Convert Dzierlenga 2021, IV male serum from long to wide format using reshape
experiment11[[1]] <- reshape(dzi_IV_Mserum[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[1]]) <- c("Time",unique(dzi_IV_Mserum$Tissue))

# Convert Dzierlenga 2021, ORAL male serum low from long to wide format using reshape
experiment11[[2]] <- reshape(dzi_OR_Mserum_low[c("Tissue" ,"Time_hours", 
                                                 "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[2]]) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))

# Convert Dzierlenga 2021, ORAL male serum medium from long to wide format using reshape
experiment11[[3]] <- reshape(dzi_OR_Mserum_medium[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[3]]) <- c("Time",unique(dzi_OR_Mserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL male serum high from long to wide format using reshape
experiment11[[4]] <- reshape(dzi_OR_Mserum_high[c("Tissue" ,"Time_hours", 
                                                  "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[4]]) <- c("Time",unique(dzi_OR_Mserum_high$Tissue))

#Convert Dzierlenga 2021, IV female serum from long to wide format using reshape
experiment11[[5]] <- reshape(dzi_IV_Fserum[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[5]]) <- c("Time",unique(dzi_IV_Fserum$Tissue))

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment11[[6]] <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                                 "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[6]]) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL female serum medium from long to wide format using reshape
experiment11[[7]] <- reshape(dzi_OR_Fserum_medium[c("Tissue" ,"Time_hours", 
                                                    "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[7]]) <- c("Time",unique(dzi_OR_Fserum_medium$Tissue))

#Convert Dzierlenga 2021, ORAL female serum high from long to wide format using reshape
experiment11[[8]] <- reshape(dzi_OR_Fserum_high[c("Tissue" ,"Time_hours", 
                                                  "Concentration_microM")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment11[[8]]) <- c("Time",unique(dzi_OR_Fserum_high$Tissue))

experiment12 <- list()
#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment12[[1]] <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment12[[1]]) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment12[[2]]<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                            "Concentration_microg_per_g_organ")], 
                            idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment12[[2]]) <- c("Time",unique(kim_IV_Fblood$Tissue))


# Convert Gustafsson Oral male blood from long to wide format using reshape
experiment13 <- list()
experiment13[[1]] <- reshape(gus_OR_Mblood[c("Tissue" ,"Time_hours", 
                                             "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment13[[1]]) <- c("Time",unique(gus_OR_Mblood$Tissue))


# Convert Gustafsson Oral male tissues from long to wide format using reshape
experiment13[[2]] <- reshape(gus_OR_Mtissues[c("Tissue" ,"Time_hours", 
                                               "Concentration_microg_per_g_organ")], 
                             idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment13[[2]]) <- c("Time",unique(gus_OR_Mtissues$Tissue))


# Put the experiments in a list
experiments <- list(experiment1 = experiment1, experiment2 = experiment2, experiment3 = experiment3, 
                    experiment4 = experiment4,experiment5 = experiment5, experiment6 = experiment6,
                    experiment7 = experiment7, experiment8 = experiment8,experiment9 = experiment9, 
                    experiment10 = experiment10, experiment11 = experiment11, 
                    experiment12 = experiment12, experiment13 = experiment13)


# Rename predictions so that they share the same name as the names of the experimental data dataframe
colnames(preds_kudo_high) <- c( "Time", "Blood", "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                                "Heart", "Brain", "Gonads", "Stomach", "Intestine")

colnames(preds_kudo_low) <-  colnames(preds_kudo_high) 
colnames(preds_kim_IV_Mtissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")

colnames(preds_kim_OR_Mtissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_IV_Ftissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_OR_Ftissues) <- colnames(preds_kim_IV_Mtissues)

colnames(preds_dzi_OR_Mtissues) <- c("Time","Liver","Kidney","Brain")
colnames(preds_dzi_OR_Ftissues) <- colnames(preds_dzi_OR_Mtissues)

colnames(preds_kim_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_kim_IV_Mblood) <- c ("Time", "Plasma")

colnames(preds_Kemp_OR_Ffeces_high) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_high) <- c ("Time", "Feces")

colnames(preds_Kemp_OR_Furine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_high) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_high) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Ffeces_med) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_med) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Ffeces_low) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_low) <- c ("Time", "Feces")

colnames(preds_dzi_IV_Mserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_high) <- c ("Time", "Serum")
colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")


colnames(preds_gus_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_gus_OR_Mtissues) <- c ("Time", "BALF", "Liver", "Lung", "Kidney")



# Create a list containing the corresponding predictions
simulations <- list(predictions1 = list(preds_kudo_high = preds_kudo_high), 
                    predictions2 = list(preds_kudo_low = preds_kudo_low), 
                    predictions3 = list(preds_kim_IV_Mtissues =  preds_kim_IV_Mtissues), 
                    predictions4 = list(preds_kim_OR_Mtissues = preds_kim_OR_Mtissues), 
                    predictions5 = list(preds_kim_IV_Ftissues = preds_kim_IV_Ftissues),
                    predictions6 = list(preds_kim_OR_Ftissues = preds_kim_OR_Ftissues),
                    predictions7 = list(preds_dzi_OR_Mtissues =  preds_dzi_OR_Mtissues), 
                    predictions8 = list(preds_dzi_OR_Ftissues = preds_dzi_OR_Ftissues), 
                    predictions9 = list(preds_kim_OR_Mblood = preds_kim_OR_Mblood,  
                                        preds_kim_IV_Mblood = preds_kim_IV_Mblood), 
                    predictions10 = list(preds_Kemp_OR_Furine_low = preds_Kemp_OR_Furine_low,  
                                         preds_Kemp_OR_Ffeces_low = preds_Kemp_OR_Ffeces_low,
                                         preds_Kemp_OR_Furine_med = preds_Kemp_OR_Furine_med, 
                                         preds_Kemp_OR_Ffeces_med = preds_Kemp_OR_Ffeces_med,
                                         preds_Kemp_OR_Furine_high = preds_Kemp_OR_Furine_high, 
                                         preds_Kemp_OR_Ffeces_high = preds_Kemp_OR_Ffeces_high,
                                         preds_Kemp_OR_Murine_low = preds_Kemp_OR_Murine_low,  
                                         preds_Kemp_OR_Mfeces_low = preds_Kemp_OR_Mfeces_low, 
                                         preds_Kemp_OR_Murine_med = preds_Kemp_OR_Murine_med, 
                                         preds_Kemp_OR_Mfeces_med = preds_Kemp_OR_Mfeces_med,
                                         preds_Kemp_OR_Murine_high = preds_Kemp_OR_Murine_high,  
                                         preds_Kemp_OR_Mfeces_high = preds_Kemp_OR_Mfeces_high),
                    prediction11 = list(preds_dzi_IV_Mserum = preds_dzi_IV_Mserum,  
                                        preds_dzi_OR_Mserum_low = preds_dzi_OR_Mserum_low,
                                        preds_dzi_OR_Mserum_medium = preds_dzi_OR_Mserum_medium,  
                                        preds_dzi_OR_Mserum_high = preds_dzi_OR_Mserum_high, 
                                        preds_dzi_IV_Fserum = preds_dzi_IV_Fserum,  
                                        preds_dzi_OR_Fserum_low = preds_dzi_OR_Fserum_low,
                                        preds_dzi_OR_Fserum_medium = preds_dzi_OR_Fserum_medium,
                                        preds_dzi_OR_Fserum_high  = preds_dzi_OR_Fserum_high), 
                    predictions12  = list(preds_kim_OR_Fblood = preds_kim_OR_Fblood, 
                                          preds_kim_IV_Fblood = preds_kim_IV_Fblood), 
                    predictions13 = list(preds_gus_OR_Mblood = preds_gus_OR_Mblood, 
                                         preds_gus_OR_Mtissues = preds_gus_OR_Mtissues))



plot_names <- list(predictions1 = c("Blood", "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                                    "Heart", "Brain", "Gonads", "Stomach", "Intestine"), 
                   predictions2 = c("Blood", "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                                    "Heart", "Brain", "Gonads", "Stomach", "Intestine"), 
                   predictions3 = list(preds_kim_IV_Mtissues =  c("Liver",  "Kidney", "Lung",
                                                                  "Spleen", "Heart")), 
                   predictions4 = list(preds_kim_OR_Mtissues = c("Liver",  "Kidney", "Lung",
                                                                 "Spleen", "Heart")), 
                   predictions5 = list(preds_kim_IV_Ftissues = c("Liver",  "Kidney", "Lung",
                                                                 "Spleen", "Heart")),
                   predictions6 = list(preds_kim_OR_Ftissues = c("Liver",  "Kidney", "Lung",
                                                                 "Spleen", "Heart")),
                   predictions7 = list(preds_dzi_OR_Mtissues =   c("Liver","Kidney","Brain")), 
                   predictions8 = list(preds_dzi_OR_Ftissues = c("Liver","Kidney","Brain")), 
                   predictions9 = list(preds_kim_OR_Mblood = "Oral",  
                                       preds_kim_IV_Mblood = "IV"), 
                   predictions10 = list(preds_Kemp_OR_Furine_low = "Urine Female 1 mg/kg",  
                                        preds_Kemp_OR_Ffeces_low = "Feces Female 1 mg/kg",
                                        preds_Kemp_OR_Furine_med = "Urine Female 5 mg/kg", 
                                        preds_Kemp_OR_Ffeces_med = "Feces Female 5 mg/kg",
                                        preds_Kemp_OR_Furine_high = "Urine Female 25 mg/kg", 
                                        preds_Kemp_OR_Ffeces_high = "Feces Female 25 mg/kg",
                                        preds_Kemp_OR_Murine_low = "Urine Male 1 mg/kg",  
                                        preds_Kemp_OR_Mfeces_low = "Feces Male 1 mg/kg", 
                                        preds_Kemp_OR_Murine_med = "Urine Male 5 mg/kg", 
                                        preds_Kemp_OR_Mfeces_med = "Feces Male 5 mg/kg",
                                        preds_Kemp_OR_Murine_high = "Urine Male 25 mg/kg",  
                                        preds_Kemp_OR_Mfeces_high = "Feces Male 25 mg/kg"),
                   prediction11 = list(preds_dzi_IV_Mserum = "IV Male 6 mg/kg",  
                                       preds_dzi_OR_Mserum_low = "Oral Male 6 mg/kg",
                                       preds_dzi_OR_Mserum_medium = "Oral Male 12 mg/kg",  
                                       preds_dzi_OR_Mserum_high = "Oral Male 48 mg/kg", 
                                       preds_dzi_IV_Fserum =  "IV Female 40 mg/kg",  
                                       preds_dzi_OR_Fserum_low = "Oral Female 40 mg/kg",
                                       preds_dzi_OR_Fserum_medium = "Oral Female 80 mg/kg",
                                       preds_dzi_OR_Fserum_high  = "Oral Female 320 mg/kg"), 
                   predictions12  = list(preds_kim_OR_Fblood = "Oral", 
                                         preds_kim_IV_Fblood = "IV"), 
                   predictions13 = list(preds_gus_OR_Mblood = "Plasma", 
                                        preds_gus_OR_Mtissues = c("BALF", "Liver", "Lung", "Kidney")))



colnames(preds_kim_OR_Mtissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_IV_Ftissues) <- colnames(preds_kim_IV_Mtissues)
colnames(preds_kim_OR_Ftissues) <- colnames(preds_kim_IV_Mtissues)

colnames(preds_dzi_OR_Mtissues) <- c("Time","Liver","Kidney","Brain")
colnames(preds_dzi_OR_Ftissues) <- colnames(preds_dzi_OR_Mtissues)

colnames(preds_kim_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_kim_IV_Mblood) <- c ("Time", "Plasma")

colnames(preds_Kemp_OR_Ffeces_high) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_high) <- c ("Time", "Feces")

colnames(preds_Kemp_OR_Furine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_high) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Murine_high) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Ffeces_med) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_med) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Ffeces_low) <- c ("Time", "Feces")
colnames(preds_Kemp_OR_Mfeces_low) <- c ("Time", "Feces")

colnames(preds_dzi_IV_Mserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Mserum_high) <- c ("Time", "Serum")
colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")


colnames(preds_gus_OR_Mblood) <- c ("Time", "Plasma")
colnames(preds_gus_OR_Mtissues) <- c ("Time", "BALF", "Liver", "Lung", "Kidney")


library(ggplot2)
library(patchwork)

create.plots <- function(predictions, observations, compartment, plot_name) {
  cls <- c("Prediction" = "#0072B2", "Observation" = "#D55E00")
  
  if(plot_name %in% c( "ELF",  "Brain",  "Liver",  "Kidney", "Carcass", "Lung",  "Spleen", 
                       "Heart", "Brain", "Gonads", "Stomach", "Intestine") ){
    y_name <-  expression("PFOA concentration (" * mu * "g/g tissue)")
  }else if(plot_name %in% c( "Urine Female 1 mg/kg", "Feces Female 1 mg/kg",  "Urine Female 5 mg/kg", 
                             "Feces Female 5 mg/kg", "Urine Female 25 mg/kg", "Feces Female 25 mg/kg",
                             "Urine Male 1 mg/kg","Feces Male 1 mg/kg",  "Urine Male 5 mg/kg", 
                             "Feces Male 5 mg/kg","Urine Male 25 mg/kg", "Feces Male 25 mg/kg")){
    y_name <- expression("PFOA Mass (" * mu * "g)")
  }else{
    y_name <- expression("PFOA concentration (" * mu * "g/mL)")
  }
  ggplot(data = predictions) +
    geom_line(aes_string(x = "Time", y = rlang::expr(!!compartment), color = '"Prediction"'), 
              size = 1.2, alpha = 0.9) +
    geom_point(data = observations, 
               aes_string(x = "Time", y = rlang::expr(!!compartment), color = '"Observation"'), 
               size = 3.5, shape = 21, stroke = 1, fill = "white") +
    labs(title = plot_name, 
         y = y_name,
         x = "Time (hours)",
         color = "Type") +
    scale_color_manual(values = cls) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
}

# Iterate over experiments
for (i in seq_along(experiments)) {
  all_plots <- list()
  
  for (j in seq_along(experiments[[i]])) {
    observations <- experiments[[i]][[j]]
    predictions <- simulations[[i]][[j]]
    compartments <- names(predictions)[-1]  # Skip "Time"
    
    if(length(compartments) > 1){
      compartment_plots <- lapply(compartments, function(compartment) {
        create.plots(
          predictions = predictions,
          observations = observations,
          compartment = compartment,
          plot_name = compartment
        )
      })
    }else{
      compartment_plots <- lapply(compartments, function(compartment) {
        create.plots(
          predictions = predictions,
          observations = observations,
          compartment = compartment,
          plot_name = plot_names[[i]][[j]]
        )
      })
    }
    
    all_plots <- c(all_plots, compartment_plots)  # Flatten here
  }
  
  # Combine with patchwork
  if (length(all_plots) == 1) {
    final_plot <- all_plots[[1]] + theme(legend.position = "bottom")
  } else {
    final_plot <- wrap_plots(all_plots, ncol = 3) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
  }
  
  ggsave(
    filename = paste0("experiment", i, ".png"),
    plot = final_plot,
    device = "png",
    dpi = 300,
    width = 14,
    height = ceiling(length(all_plots) / 3) * 4,
    units = "in",
    bg = "white"
  )
}