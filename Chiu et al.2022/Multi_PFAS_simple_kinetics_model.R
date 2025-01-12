# The code replicates the Loccisano et al.(2012) PBK model
# and simulates IV or Oral PFOA distribution for female or male rats of constant or varying
# body weight. The Bernstein et al. (2021) code was used for replicating the
# model and the parameter set was expanded to also describe female rats
# (the original code only described male rats). 
# We did not include simultaneous IV and oral because there are no such experiments

library(deSolve)
library(ggplot2)

#=========================
#1. Parameters of the model
#=========================

create.params  <- function(user_input){
  with( as.list(user_input),{
    # Vd: L/kg
    # T1/2:years
    print("Hi")
    
   if (substance == "PFOA"){
     Vd <- 0.43 * BW
     T_half <-  3.14*365
     k <-   log(2)/T_half 
   }else if (substance == "PFOS"){
      Vd <- 0.32* BW
      T_half <- 3.36*365
      k <- log(2)/T_half
   }else if (substance == "PFNA"){
      Vd <- 0.19* BW
      T_half <- 2.35*365
      k <- log(2)/T_half
   }else if (substance == "PFHxS"){
      Vd <- 0.29* BW
      T_half <- 8.30*365
      k <- log(2)/T_half
   }
    
    return(list("BW" = BW, "Vd" = Vd, "k"  = k, "Cserum_init" =Cserum_init,
                "DW_rate" = DW_rate, "Cwater" = Cwater,
                "Cwater_times" = Cwater,
                "non_DW_intake" = non_DW_intake,
                "ingestion_times" = ingestion_times
                 
    ))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
   Cserum <- Cserum_init
   non_DW_intake<- 0
   Cwater <- 0
    return(c("Cserum" = Cserum, "non_DW_intake" = non_DW_intake, "Cwater" = Cwater))
  })
}

#===================
#3. Events function
#===================

create.events <- function(parameters){
  with(as.list(parameters), {
      # Calculate number of administrated doses and corresponding administration time
      lcwater <- length(Cwater)
      lcwatertimes <- length(Cwater_times)
      lingest <- length(non_DW_intake)
      lingesttimes <- length(ingestion_times)
      # If not equal, then stop 
      if (lcwater != lcwatertimes){
        stop("The times of water concentration change should be equal in vector of Cwater")
      }else if (lingest != lingesttimes){
        stop("The times of ingestion rate change should be equal in vector of Cwater")
      }else{
        events <- list(data = rbind(data.frame(var = c("Cwater"),  time = Cwater_times, 
                                               value = Cwater, method = c("rep")),
                                    
                                    data.frame(var = c("non_DW_intake"),  time = ingestion_times, 
                                               value = non_DW_intake, method = c("rep"))))
      }
      
    
    return(events)
  })
}

#==================
#4. Custom function 
#==================
custom.func <- function(){
  return()
}

#==============
#5. ODEs System
#==============

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    # Cserum: in ug/L
    dCserum <- ((non_DW_intake+DW_rate*Cwater)/Vd) - k*Cserum 
    dCwater <- 0
    dnon_DW_intake <- 0
    list(c( "dCserum"  =  dCserum, "dCwater" = dCwater,  "dnon_DW_intake" = dnon_DW_intake))
    
  })
}

#=============
#6. User input 
#=============
# Parameters for reproducing Figure 2 of Worley et al. 2017
#North Alabama exposure
substance <- "PFOA"
BW <- 82.3# average body weight in kg [ATSDR, 2016]
Cserum_init <- #ug/L
Cwater <- 0.04 #ug/L
Cwater_times <- 0
non_DW_intake <- 0.01*24 #ug/day
ingestion_times <- 0
Cserum_init <- 0
DW_rate <- 1.5 #L/day

user_input <- list( substance = substance, "BW" = BW, "DW_rate" = DW_rate,
                    "Cwater" = Cwater,
                    "Cwater_times" = Cwater_times, "non_DW_intake" = non_DW_intake,
                    "ingestion_times" = ingestion_times, "Cserum_init" = Cserum_init)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,50*365,1)  #days

solution <-  as.data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                               events = events, method="bdf",rtol = 1e-05, atol = 1e-05))
ggplot(solution, aes(x = time/365, y = Cserum)) +
  geom_line() +
  theme_minimal()   # Set a minimal theme

predicted.feats <- c("Cserum")
# Log in Jaqpot server

# Deploy the model on the Jaqpot server to create a web service
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func,
                     envFile = ".env")