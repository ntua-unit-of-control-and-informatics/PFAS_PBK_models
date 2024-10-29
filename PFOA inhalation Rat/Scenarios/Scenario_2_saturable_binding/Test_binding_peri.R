#female increase in BW: 3.5 g/d
#Male increase in BW: 5.9 g.d

library(deSolve)

create.params <- function(user.input){
  with(as.list(user.input),{

    # https://doi.org/10.3390/toxics12040253, Table 1
    
    n=1
    Ka <- 7e4
    CalbB_init <- 486*1e-06 #mol/L
    Calb_exp_init <- 600*1e-06 #mol/L
    Cpfoa_init <-  1.92/414.07#5*1e-06 #mol/L
    koff_alb <- 0.0005 #1/s
    kon_alb <- Ka * koff_alb #1/M/s

    return(list("CalbB_init" = CalbB_init, "Cpfoa_init" = Cpfoa_init, "kon_alb" = kon_alb, 
                "koff_alb" = koff_alb,
                "Calb_exp_init" = Calb_exp_init,"n"=n
    ))
    
  })
}  

ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{

    
    dCalb_free = koff_alb*Calb_bound - n*kon_alb*Calb_free*CPFOA_free
    dCPFOA_free = koff_alb*Calb_bound - n*kon_alb*Calb_free*CPFOA_free
    
    dCalb_bound = - koff_alb*Calb_bound + n*kon_alb*Calb_free*CPFOA_free
    dCPFOA_bound = - koff_alb*Calb_bound + n*kon_alb*Calb_free*CPFOA_free
    
   ff = CPFOA_free / (CPFOA_free  + CPFOA_bound )

    
    list(c('dCalb_free' = dCalb_free, 'dCPFOA_free' = dCPFOA_free, 'dCalb_bound' = dCalb_bound,
           'dCPFOA_bound' = dCPFOA_bound), ff = ff)
    
  })
}

#Initial condition for each compartment.

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    Calb_free <- CalbB_init; CPFOA_free <- Cpfoa_init
    Calb_bound<- 0; CPFOA_bound<-0
    
    
    return(c('Calb_free' = Calb_free, 'CPFOA_free' = CPFOA_free,
             'Calb_bound' = Calb_bound, 'CPFOA_bound' = CPFOA_bound))
    
    
  })
}


create.events <- function(parameters){
  with(as.list(parameters), {
    
    
    return()
  })
}


###############################################
###############################################

params <- create.params(NULL)
inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,48,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-07, atol = 1e-07))
plot(solution$CPFOA_free, solution$CPFOA_bound, type = "l")

Bmax = 486*1e-06
Ka <- 5.8e5
kd <- 1/Ka
CalbB_init <- 486*1e-06 #mol/L
Calb_exp_init <- 600*1e-06 #mol/L
Cpfoa_init <-  1.92/414.07#5*1e-06 #mol/L
koff_alb <- 0.01 #1/s
kon_alb <- Ka * koff_alb #1/M/s

Cfree = seq(1e-6,1e-4, 1e-6)
Cbound = n*Bmax*Cfree/(Cfree+kd)
plot(Cfree, Cbound)
