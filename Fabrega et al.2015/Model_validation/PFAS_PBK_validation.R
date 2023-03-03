library(deSolve)
wdir <- 'C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/Fabrega et al.2015/Model_validation'
setwd(wdir)

#============================#
#1. Parameters of the model  #
#============================#

create.params <- function(user_input){
  with(as.list(user_input),{
    # Units of parameters
    # Partition Coefficients: (unitless)
    # V_i: Volume of tissue i (ml) 
    # Q_i: Blood flow of tissue i (L/h)
    # Tm: Resoprtion maximum (ng/h)
    # Kt: Affinity constant (ng/h)
    # Free: the unbound fraction of substance in the plasma (unitless)
    
    # Load the physiological parameters (Fabrega et al.(2021): Table 1)
    physiological_params <- openxlsx::read.xlsx('Physiological_parameters.xlsx')
    colnames(physiological_params) <- c('Compartment', 'Volume', 'Cardiac_output')
    
    # N_compertments is the number of compartments of the model (total volume and 
    # cardiac output included)
    N_compartments = dim(physiological_params)[1]
    
    physiological_list <- list()  
    # Create a V_i and Q_i for each compartment and assign the corresponding values
    for (i in 1:N_compartments) {
      physiological_list[i] <- physiological_params$Volume[i]
      names(physiological_list)[i] <- paste0('V_', physiological_params$Compartment[i])
      physiological_list[N_compartments + i] <- physiological_params$Cardiac_output[i]
      names(physiological_list)[N_compartments + i] <- paste0('Q_', physiological_params$Compartment[i])
    }
    physiological_list <- within(physiological_list, rm(Q_Plasma)) # remove Q_Plasma because it does not exist
    
    substance_params <- openxlsx::read.xlsx('Partition_Coefficients.xlsx')
    substance_list <- list()
    keep_params <- substance_params[substance_params$Substance==substance,2:dim(substance_params)[2]]
    for (i in 1:length(keep_params)) {
      substance_list[i] <- keep_params[i]
      if(i < 6){
        names(substance_list)[i] <- paste0('P_', names(keep_params)[i])
      }else{
        names(substance_list)[i] <- names(keep_params)[i]
      }
    }
    # According to Fabrega et al. (2015) the Partition coefficients of
    # Fat, Rob and Gut compartments come from the data on PFOA in rats 
    # from Loccisano et al.(2011)
    constant_list <- list('P_Fat' = 0.04, 'P_Rob' = 0.12, 'P_Gut' = 0.05)
    substance_list <- append(constant_list, substance_list)
    
    # Load the data for the intake (Fabrega et al.(2015) Table 1)
    # The column of total intake is given in ng/day
    intake_data <- openxlsx::read.xlsx('PFAS_Intake.xlsx')
    if(is.na(exposure)){
      Intake_day = intake_data[intake_data$Substance==substance, "Total"]
      exposure <- Intake_day/24 # Intake given in ng/h
    }else{
      exposure <- exposure
    }
    
    substance_list <- append(substance_list, c('exposure'=exposure, 
                                               'exposure_times'=exposure_times))
    params_list <- append(physiological_list, substance_list)  
    return(params_list)
  })
}

#===============================================#
#2. Function to create initial values for ODEs  #
#===============================================#

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    C_Gut <- 0; C_Liver <- 0; C_Kidney <- 0; C_Filtrate <- 0; C_Fat <- 0
    C_Lung <- 0; C_Brain <- 0; C_Rob <- 0; C_Plasma <- 0; Intake <- 0
    
    return(c('C_Gut' = C_Gut, 'C_Liver' = C_Liver, 'C_Kidney' = C_Kidney,
             'C_Filtrate' = C_Filtrate, 'C_Fat' = C_Fat, 
             'C_Lung' = C_Lung, 'C_Brain' = C_Brain, 'C_Rob' = C_Rob,
             'C_Plasma' = C_Plasma, 'Intake'=Intake))
  })
}

#===================#
#3. Events function #
#===================#

create.events <- function(parameters){
  with(as.list(parameters),{
    ltimes <- length(exposure_times)
    lexposure <- length(exposure)
    
    events <- data.frame(var = rep('Intake', ltimes), time = exposure_times,
                         value = exposure, method = rep('rep',ltimes))
    return(list(data=events))
  })
}

#==================
#4. Custom function 
#==================
custom.func <- function(parameters){
  
  return()
}

#===============#
#5. ODEs System #
#===============#

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits, params)),{
    # Units
    # C_i: Concentration in tissue i (ng/L)
    # C_Plasma: Concentration in plasma (ng/L)
    # Intake: hourly ingestion of PFASs (ng/h)
    
    # Concentration in Gut
    dC_Gut <- (Q_Gut*Free*(C_Plasma - C_Gut/P_Gut) + Intake)/V_Gut
    
    # Concentration in Liver
    dC_Liver <- Free*(Q_Liver*C_Plasma + Q_Gut*C_Gut/P_Gut - (Q_Liver+Q_Gut)*C_Liver/P_Liver)/V_Liver
    
    # Concentration in Kidney
    dC_Kidney <- ((Q_Kidney*Free*(C_Plasma-C_Kidney)/P_Kidney) + (Tm*C_Filtrate)/(Kt+C_Filtrate))/V_Kidney
    
    # Concentration in Filtrate
    dC_Filtrate = (Q_Filtrate*(Free*C_Plasma - C_Filtrate) - (Tm*C_Filtrate)/(Kt+C_Filtrate))/V_Filtrate
    
    # Concentration in Fat
    dC_Fat <- Q_Fat*Free*(C_Plasma - C_Fat/P_Fat)/V_Fat
    
    # Concentration in Lung
    dC_Lung <- Q_Lung*Free*(C_Plasma - C_Lung/P_Lung)/V_Lung
    
    # Concentration in Brain
    dC_Brain <- Q_Brain*Free*(C_Plasma - C_Brain/P_Brain)/V_Brain
    
    # Concentration in Rob
    dC_Rob <- Q_Rob*Free*(C_Plasma - C_Rob/P_Rob)/V_Rob
    
    # Concentration in Plasma
    dC_Plasma = Q_Total*(C_Gut/P_Gut + C_Liver/P_Liver + C_Kidney/P_Kidney + C_Filtrate +
                         C_Fat/P_Fat + C_Lung/P_Lung + C_Brain/P_Brain + C_Rob/P_Rob - C_Plasma)
    dIntake <- 0
    
    return(list(c("dC_Gut"=dC_Gut, "dC_Liver"=dC_Liver, "dC_Kidney"=dC_Kidney,
                  "dC_Filtrate"=dC_Filtrate, "dC_Fat"=dC_Fat,
                  "dC_Lung"=dC_Lung, "dC_Brain"=dC_Brain, "dC_Rob"=dC_Rob,
                  "dC_Plasma"=dC_Plasma, "dIntake"=dIntake)))
    
  })
}

################################################################################
substance <- 'PFOS'
# Considered mean lifetime of 80 years
L <- 80*360*24 #hours 
exposure <- NA
exposure_times <- 0 

user_input <- list( "substance" = substance,
                    "exposure"=exposure,
                    "exposure_times"= exposure_times)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
L = 10
sample_time <- seq(0,L,1) #hours

solution <- ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                events = events, 
                method="bdf",rtol = 1e-05, atol = 1e-05)
print(head(solution))




