# The code replicates the Fischer et al.(2017) in vitro kinetics model 

library(deSolve)
library(ggplot2)

#=========================
#1. Parameters of the model
#=========================



create.params  <- function(user_input){
  with( as.list(user_input),{
     
    return(list( "C_nominal" = C_nominal,
                 "cells_per_well" =cells_per_well,
                 "cell_water_content"  =cell_water_content,
                 "cell_protein_content" =cell_protein_content, 
                 "cell_lipid_content" = cell_lipid_content,
                 "medium_protein_content" = medium_protein_content,
                 "medium_lipid_content" = medium_lipid_content,
                 "serum_protein_content" = serum_protein_content,
                 "serum_lipid_content" = serum_lipid_content,
                 "serum_fraction" = serum_fraction,
                 "pKa" = pKa,
                 "pH" = pH,
                 "log_D_bsa_w"  = log_D_bsa_w,
                 "log_D_lip_w" = log_D_lip_w,
                 "log_K_bsa_w_neutral" = log_K_bsa_w_neutral,
                 "log_K_bsa_w_ion" = log_K_bsa_w_ion,
                 "log_K_lip_w_neutral" = log_K_lip_w_neutral,
                 "log_K_lip_w_ion"  = log_K_lip_w_ion,
                 "V_medium"  = V_medium
                 
    ))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
    N <- 0
    
    return(c("N" = N))
  })
}

#===================
#3. Events function
#===================

create.events <- function(parameters){
  with(as.list(parameters), {
      events <- list()
      
    
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
      dN <- 0
      
      # Medium volumes
      V_prot_medium <- ((1-serum_fraction)*medium_protein_content + serum_fraction*serum_protein_content)*V_medium/1000
      V_lip_medium <-  ((1-serum_fraction)*medium_lipid_content + serum_fraction*serum_lipid_content)*V_medium/1000; 
      V_water_medium <- V_medium - V_prot_medium - V_lip_medium
      VF_prot_medium <- V_prot_medium/V_medium 
      VF_lip_medium <- V_lip_medium/V_medium 
      VF_water_medium <- V_water_medium/V_medium 
      
      
      # Cell volumes
      V_water_cell <- 1000*cell_water_content*cells_per_well/1e6 
      V_prot_cell <- 1000*cell_protein_content*cells_per_well/1e6  
      V_lip_cell <- 1000*cell_lipid_content*cells_per_well /1e6
      V_cell <- V_water_cell + V_prot_cell + V_lip_cell
      VF_prot_cell <-  V_prot_cell/V_cell  
      VF_lip_cell <- V_lip_cell/V_cell 
      VF_water_cell <- V_water_cell/V_cell
      
      
      V_sys <- V_medium+V_cell
      V_water_sys <- V_water_medium + V_water_cell
      V_prot_sys <- V_prot_medium + V_prot_cell 
      V_lip_sys <- V_lip_medium + V_lip_cell
      VF_prot_sys <- V_prot_sys/V_sys 
      VF_lip_sys <- V_lip_sys/V_sys 
      VF_water_sys <- V_water_sys/V_sys 
      
  
      frac_neutral <- 1/(1+10^(pH-pKa)) 
      frac_ion <- 1-frac_neutral
      
      if(log_D_bsa_w<0 | log_D_lip_w<0 ){
        K_bsa_w_neutral <- 10^(log_K_bsa_w_neutral)
        K_bsa_w_ion <- 10^(log_K_bsa_w_ion)
        K_lip_w_neutral <- 10^(log_K_lip_w_neutral)
        K_lip_w_ion <- 10^(log_K_lip_w_ion)
        D_bsa_w <- frac_neutral*K_bsa_w_neutral + frac_ion*K_bsa_w_ion; 
        D_lip_w <- frac_neutral*K_lip_w_neutral + frac_ion*K_lip_w_ion
      }else{
        D_bsa_w <- 10^(log_D_bsa_w)
        D_lip_w <- 10^(log_D_lip_w)
      }  
      
      D_medium_w <- VF_prot_medium*D_bsa_w+VF_lip_medium*D_lip_w+VF_water_medium
      D_cell_w <- VF_prot_cell*D_bsa_w+VF_lip_cell*D_lip_w+VF_water_cell
      D_medium_cell <- D_medium_w/D_cell_w
      
      f_cell <- (1 + D_medium_cell*V_medium/V_cell)^(-1)
      f_medium <- 1-f_cell
      
      f_free_medium <- (1+D_bsa_w*V_prot_medium/V_water_medium + D_lip_w*V_lip_medium/V_water_medium +
                          D_cell_w*V_cell/V_water_medium)^(-1)
      
      f_mem <- (1 +(1/D_lip_w)*V_water_cell/V_lip_cell + (D_bsa_w/D_lip_w)* V_prot_cell/V_lip_cell+
                  (D_medium_w/D_lip_w)*V_water_medium/V_lip_cell)^(-1)
      
      C_free <- C_nominal * f_free_medium * V_sys/V_water_medium


    
    list(c("dN" = dN), "C_nominal" = C_nominal , "C_free" = C_free, 
         "frac_ion" = frac_ion, "f_cell" = f_cell, "f_medium" = f_medium, 
         "f_free_medium" = f_free_medium)
    
  })
}

#=============
#6. User input 
#=============
C_nominal <- 10 #mol/L
cells_per_well <- 1500
cell_water_content <- 2.65E-03 * 3000/850 #mL per million cells
cell_protein_content <- 2.87E-04 * 3000/850 #mL per million cells
cell_lipid_content <- 9.69E-05 * 3000/850  #mL per million cells
medium_protein_content <- 0 #mL/L
medium_lipid_content<- 0 #mL/L
serum_protein_content<- 0  #mL/L
serum_lipid_content<- 0 # mL/L
serum_fraction<- 0  # %
pKa <-8
pH <- 7.4
log_D_bsa_w <- -1#3.939
log_D_lip_w <- -1#3.52#
log_K_bsa_w_neutral <- 7
log_K_bsa_w_ion <- 2
log_K_lip_w_neutral <- 6
log_K_lip_w_ion <- 2
V_medium <-  100#uL
  
user_input <- list( "C_nominal" = C_nominal,
                    "cells_per_well" =cells_per_well,
                    "cell_water_content"  =cell_water_content,
                    "cell_protein_content" =cell_protein_content, 
                    "cell_lipid_content" = cell_lipid_content,
                    "medium_protein_content" = medium_protein_content,
                    "medium_lipid_content" = medium_lipid_content,
                    "serum_protein_content" = serum_protein_content,
                    "serum_lipid_content" = serum_lipid_content,
                    "serum_fraction" = serum_fraction,
                    "pKa" = pKa,
                    "pH" = pH,
                    "log_D_bsa_w"  = log_D_bsa_w,
                    "log_D_lip_w" = log_D_lip_w,
                    "log_K_bsa_w_neutral" = log_K_bsa_w_neutral,
                    "log_K_bsa_w_ion" = log_K_bsa_w_ion,
                    "log_K_lip_w_neutral" = log_K_lip_w_neutral,
                    "log_K_lip_w_ion"  = log_K_lip_w_ion,
                    "V_medium"  = V_medium)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,1,1)  #days
solution <-  as.data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                               events = events, method="bdf",rtol = 1e-05, atol = 1e-05))
#====================
#7. Upload on Jaqpot 
#===================

# Subset of features to be displayed on the user interface
predicted.feats <- c( "C_nominal", "C_free",  "frac_ion" , "f_cell" , "f_medium", 
                      "f_free_medium")
# Log in Jaqpot server

# Deploy the model on the Jaqpot server to create a web service
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func,
                     ode.fun =  ode.func, envFile = ".env")