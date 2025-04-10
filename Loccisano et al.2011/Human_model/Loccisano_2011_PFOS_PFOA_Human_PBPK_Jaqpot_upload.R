library(deSolve)

create.params <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al 1997)
    #fractional blood flows
    QCC = 12.5 		# Cardiac blood output (L/h/kg^0.75)
    QFC = 0.052		# Fraction cardiac output going to fat
    QLC = 0.25 		# Fraction cardiac output going to liver
    QKC = 0.175		# Fraction cardiac output going to kidney
    #QfilC = 0.035		# Fraction cardiac output to the filtrate compartment (20% of kidney blood flow)
    QSkC = 0.058	# Fraction cardiac output going to skin
    QGC = 0.181		# Fraction cardiac output going to gut
    
    #fractional tissue volumes
    VLC = 0.026		# Fraction liver volume
    VFC = 0.214		# Fraction fat volume
    VKC = 0.004		# Fraction kidney volume
    VfilC = 0.0004	# Fraction filtrate compartment volume (10% of kidney volume)
    VGC = 0.0171		# Fraction gut volume
    VPlasC = 0.0428	# Fraction plasma volume (58% of blood)
    Htc = 0.44    # hematocrit
    
    # for dermal exposure
    # Dermal exposure
    Dermconc <- 0.0		# Dermal concentration (mg/mL)
    Dermvol <- 0.0		  # Dermal exposure volume (mL)
    Dermdose <- Dermconc*Dermvol*1000  # (ug)
    Skinarea = 5		#Exposed area on skin (cm^2)
    SkinTarea = 9.1*((BW*1000)**0.666)	  # Total area of skin (cm^2)
    Skinthickness = 0.1	                  # Skin thickness (cm)
    
    # Scaling parameters
    QC = QCC*BW**0.75	# Cardiac output (L/h)
    QCP = QC*(1-Htc)	# adjust for plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/h)
    QF = QFC*QCP			# Plasma flow to fat (L/h)
    QK = QKC*QCP	    # Plasma flow to kidney (L/h)
    Qfil = 0.2*QK		  # Plasma flow to filtrate compartment (L/h)# 20% of QK
    QG = QGC*QCP		  # Plasma flow to gut (L/h)
    QSk <- ifelse(Dermconc > 0, QSkC*QCP*(Skinarea/SkinTarea), 0.0) # plasma flow to skin
    #QSk <- QSkC*QCP*(Skinarea/SkinTarea) # plasma flow to skin
    QR = QCP - QL - QF - QK - Qfil - QG - QSk	# Plasma flow to rest of the body (L/h)
    
    Qbal = QCP - (QL+QF+QK+Qfil+QG+QSk+QR)       # balance check--better be 0
    
    VL = VLC*BW			# Liver volume (L)
    VF = VFC*BW			# Fat volume (L)
    VK = VKC*BW			# Kidney volume (L)
    Vfil = VfilC*BW			# Fitrate compartment volume (L)
    VG = VGC*BW			# Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VSk = (Skinarea*Skinthickness)/1000	                                    # Skin volume (L)
    VR = 0.84*BW - VL - VF - VK - Vfil - VG - VPlas - VSk	          # Rest of the body volume (L)
    Vbal = (0.84*BW)-(VL+VF+VK+Vfil+VG+VPlas+VSk+VR)               # Balance check--better be 0
    
    if(substance == 'PFOS'){
      # Chemical-specific parameters (PFOS)
      Tmc = 3.5			# Maximum resorption rate (mg/h/kg^(0.75))
      Kt = 0.023		# Resorption affinity (mg/L)
      Free = 0.025	# Free fraction of PFOS in plasma
      PL = 3.72			# Liver/blood partition coefficient
      PF = 0.14			# Fat/blood partition coefficient
      PK = 0.8			# Kidney/blood partition coefficient
      PSk = 0.29		# Skin/blood partition coefficient
      PR = 0.2			# Rest of the body/blood partition coefficient
      PG = 0.57     # Gut/blood partition coeff. 
      kurinec = 0.001		# urinary elimination rate constant  (/h/kg^-0.25)# estimated from Harada, et al 2005
    } else if(substance == 'PFOA'){
      # Chemical-specific parameters (PFOA)
      Tmc = 10			# Maximum resorption rate (mg/h/kg^(0.75))
      Kt = 0.055		# Resorption affinity (mg/L)
      Free = 0.02		# Free fraction of PFOA in plasma# same as monkey
      PL = 2.2			# Liver/plasma partition coefficient
      PF = 0.04			# Fat/plasma partition coefficient
      PK = 1.05			# Kidney/plasma partition coefficient
      PSk = 0.1			# Skin/plasma partition coefficient
      PR = 0.12			# Rest of the body/plasma partition coefficient
      PG = 0.05     # Gut/blood plasma coeff. 
      kurinec = 0.0003		# Elimination rate (1/h)# estimated from data of Harada, et al 2005
    }
    
    kurine = kurinec*BW**(-0.25)
    Tm = Tmc*BW**0.75   #transporter maximum
    
    # Free fraction of chemical in tissues
    FreeL = Free/PL  #liver
    FreeF = Free/PF  #fat
    FreeK = Free/PK  #kidney
    FreeSk = Free/PSk #skin
    FreeR = Free/PR  #rest of tissues
    FreeG = Free/PG  #gut
    
    return(list('QC'=QC, 'QCP'=QCP, 'QL'=QL, 'QF'=QF, 'QK'=QK, 
                'Qfil'=Qfil, 'QG'=QG, 'QSk'=QSk, 'QR'=QR,
                'VPlas'=VPlas, 'VL'=VL, 'VF'=VF, 'VK'=VK, 
                'Vfil'=Vfil, 'VG'=VG, 'VSk'=VSk, 'VR'=VR,
                'PL'=PL, 'PF'=PF, 'PK'=PK, 'PSk'=PSk,
                'PG'=PG, 'PR'=PR,
                'Tm'=Tm, 'Kt'=Kt, 'Free'=Free,
                'FreeL'=FreeL, 'FreeF'=FreeF, 'FreeK'=FreeK,
                'FreeSk'=FreeSk, 'FreeR'=FreeR, 'FreeG'=FreeG, 
                'kurine'=kurine,
                "admin.type" = admin.type,
                "admin.dose" = admin.dose, 
                "admin.time" = admin.time,
                'Drinking_rate'=Drinking_rate
    ))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    APlas<-0; AG<-0; AL<-0; AK<-0; AF<-0; Afil<-0; AStore<-0; AUrine<-0; ASk<-0;
    AR<-0; C_water<-0; 
    return(c('APlas'=APlas, 'AG'=AG, 'AL'=AL, 'AF'=AF, 'AK'=AK, 
             'Afil'=Afil, 'AStore'=AStore, 'AUrine'=AUrine, 'ASk'=ASk, 
             'AR'=AR, 'C_water'=C_water))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{

    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if(admin.type=='oral'){
        events <- data.frame(var = c(rep('C_water', ltimes)), 
                             time = admin.time,
                             value = admin.dose,
                             method = rep('rep',ltimes))
      }else if(admin.type=='iv'){
        events <- data.frame(var = c(rep('APlas', ltimes)), 
                             time = admin.time,
                             value = admin.dose,
                             method = rep('add',ltimes))
      }
    }

    #events <- events[order(events$time),] 
    return(list(data=events))
  })
}

custom.func <- function(){
  return()
}

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits, params)),{
    # Concentrations
    CPlas <- APlas/VPlas # Concentration in plasma 
    CG <- AG/VG # Concentration in gut 
    CL <- AL/VL # Concentration in liver
    CF <- AF/VF # Concentration in fat
    CK <- AK/VK # Concentration in kidney
    Cfil <- Afil/Vfil # Concentration in filtrate
    CSk <- ASk/VSk # Concentration in skin
    CR <- AR/VR
    
    # Plasma Compartment
    dAPlas <- QF*CF*FreeF + (QL+QG)*CL*FreeL + QR*CR*FreeR + QSk*CSk*FreeSk + 
      QK*CK*FreeK - QCP*CPlas*Free
    
    # Gut compartment
    dAG <- QG*(CPlas*Free - CG*FreeG) + C_water*Drinking_rate
    
    # Liver compartment
    dAL <- QL*CPlas*Free + QG*CG*FreeG - (QL+QG)*CL*FreeL 
    
    # Fat compartment
    dAF <- QF*(CPlas*Free - CF*FreeF)
    
    # Kidney compartment
    dAK <- QK*(CPlas*Free - CK*FreeK) + (Tm*Cfil)/(Kt+Cfil)
    
    # Filtrate compartment
    dAfil = Qfil*(CPlas*Free - Cfil) - (Tm*Cfil)/(Kt+Cfil)
    
    # Storage compartment
    dAStore <- Qfil*Cfil - kurine*AStore
    
    # Urine
    dAUrine <- kurine*AStore
    
    # Skin compartment
    dASk <- QSk*(CPlas*Free-CSk*FreeSk) #+ input4*DoseOn
    
    # Rest of the body
    dAR <- QR*(CPlas*Free - CR*FreeR)
    
    # Water concentration
    dC_water <- 0
    
    return(list(c('dAPlas'=dAPlas, 'dAG'=dAG, 'dAL'=dAL, 'dAF'=dAF,
                  'dAK'=dAK, 'dAfil'=dAfil, 'dAStore'=dAStore,
                  'dAUrine'=dAUrine, 'dASk'=dASk, 'dAR'=dAR, 'dC_water'=dC_water),
                'CPlas'=CPlas, 'CG'=CG, 'CL'=CL, 'CF'=CF,
                'CK'=CK, 'Cfil'=Cfil, 'CSk'=CSk, 'CR'=CR))
  })
}

################################################################################




# reproduce results in figure 8 and figure 9
BW = 70 # human body weight in kg
substance <- 'PFOA' # select substance: PFOA/PFOS
admin.type <- 'oral' # administration type values: oral/iv
admin.dose <- c(3.55,0) # administered dose in mg
admin.time <- c(0, 30*24*360) # time when doses are administered, in hours
Drinking_rate <- BW*11/1000/24 # L/h

user_input <- list('BW'=BW, 
                   'substance'=substance,
                   "admin.type" = admin.type,
                   "admin.dose" = admin.dose, 
                   "admin.time" = admin.time,
                   'Drinking_rate'=Drinking_rate)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,50*360*24,24)

solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           events = events,
                           method="lsodes",rtol = 1e-05, atol = 1e-05))
print(tail(solution))

#====================
#7. Upload on Jaqpot 
#===================
# Subset of features to be displayed on the user interface
predicted.feats <- c('APlas', 'AG', 'AL', 'AF',
                     'AK', 'Afil', 'AStore',
                     'AUrine', 'ASk', 'AR', 'C_water',
                     'CPlas', 'CG', 'CL', 'CF',
                     'CK', 'Cfil', 'CSk', 'CR')

# Log in Jaqpot server
jaqpotr::login.api()
# Deploy the model on the Jaqpot server to create a web service
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func, 
                     method = "bdf",url = "https://api.jaqpot.org/jaqpot/")
