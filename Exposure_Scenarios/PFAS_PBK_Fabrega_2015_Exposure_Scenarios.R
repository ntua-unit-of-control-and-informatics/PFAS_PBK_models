library(deSolve)

create.params <- function(user_input){
  with(as.list(user_input),{
    # Physiological parameters (from Brown, et al)
    # fractional blood flows to tissues
    QCC = 12.5		  # Cardiac blood output (L/h/kg^0.75)
    QFC = 0.052		  # Fraction cardiac output going to fat
    QLC = 0.19 #0.25 from Loccisano 2011 		  # Fraction cardiac output going to liver
    QKC = 0.175		  # Fraction cardiac output going to kiney
    #QFilC = 0.035		# Fraction cardiac output to the filtrate compartment (20% of kiney blood flow)
    QGC = 0.181		  # Fraction cardiac output going to gut
    QLuC = 0.034    # Fraction cardiac output going to lungs
    QBC = 0.117    # Fraction cardiac output going to brain
    QBmC = 0.05    # Fraction cardiac output going to bone marrow
    
    #fractional tissue volumes
    VLC = 0.026		  # Fraction liver volume
    VFC = 0.02 #0.214		  # Fraction fat volume
    VKC = 0.004		  # Fraction kidney volume
    VFilC = 0.0004	# Fraction filtrate compartment volume (10% of kidney volume)
    VGC = 0.0171		# Fraction gut volume
    VPlasC = 0.0428	# Fraction plasma volume
    VLuC = 0.014	# Fraction lungs volume
    VBC = 0.021	# Fraction lungs volume
    VBmC = 0.1429	# Fraction lungs volume
    
    # Scaling parameters
    QC = QCC*BW**0.75	#Cardiac output (L/h)
    Htc = 0#0.44        #hematocrit
    QCP = QC*(1-Htc)	# Plasma flow
    QL = QLC*QCP			# Plasma flow to liver (L/h)
    QF = QFC*QCP			# Plasma flow to fat (L/h)
    QK = QKC*QCP	    # Plasma flow to kiney (L/h)
    QFil = 0.2*QK   	# Plasma flow to filtrate compartment (L/h)# 20% of QK
    QG = QGC*QCP	    # Plasma flow to gut (L/h)
    QLu = QLuC*QCP	  # Plasma flow to lungs (L/h)
    QB = QBC*QCP	    # Plasma flow to brain (L/h)
    QBm = QBmC*QCP	  # Plasma flow to bone marrow (L/h)
    QR = QCP - QL - QF - QK - QFil - QG - QLu - QB - QBm	 # Plasma flow to rest of the boy (L/h)
    
    Qbal = QCP - (QL+QF+QK+QFil+QG+QR+QLu+QB+QBm)        # balance check 
    
    VL = VLC*BW			    # Liver volume (L)
    VF = VFC*BW			    # Fat volume (L)
    VK = VKC*BW			    # Kiney volume (L)
    VFil = VFilC*BW	    # Fitrate compartment volume (L)
    VG = VGC*BW			    # Gut volume (L)
    VPlas = VPlasC*BW		# Plasma volume (L)
    VLu = VLuC*BW			  # Lungs volume (L)
    VB = VBC*BW			    # Brain volume (L)
    VBm = VBmC*BW			    # Bone marrow volume (L)
    VR = 1*BW - VL - VF - VK - VFil - VG -VLu -VB - VBm - VPlas		# Rest of the boy volume (L) # Loccisano 2011 uses 0.84*BW
    
    Vbal = (1*BW)-(VL+VF+VK+VFil+VG+VPlas+VR+VLu+VB+VBm)         # Balance check
    
    loaded_params <- data.frame(matrix(data = c(1.38,	4.91,	0.74,	27.31,	12.27,	2456.3,	30,		0.001,	
                                                0.72,	0.28,	0.5,	1.27,	3.27,	2456.3,	3,		0.005,
                                                7.48,	0.11,	0.36,	2.13,	5.54,	86,		0.018,	0.03,
                                                0.001,	3.28,	0.58,	5.99,	9.46,	2456.3,	15,		0.001,
                                                128.8,	39.87,	201.6,	56.11,	6.27,	6.1,	5,		0.001,
                                                47.82,	110.71,	1.94,	24.98,	11.2,	368.4,	7,		0.015,
                                                4.23,	18.73,	0.37,	9.08,	0.62,	147.4,	0.116,	0.03,
                                                1.65,	2.66,	37.8,	19.47,	36.02,	491.3,	4.5,	0.01,
                                                0.001,	0.28,	43.68,	31.92,	11.57,	245.6,	0.6,	0.01,
                                                0.002,	0.25,	14.82,	4.61,	11.69,	2456.3,	9,		0.03,
                                                0.001,	0.001,	63.13,	24.95,	15.78,	24.6,	0.07,	0.02),
                                       nrow = 11, byrow = T))
    
    colnames(loaded_params) <- c('Liver',	'Bone_marrow',	'Brain',
                                 'Lung',	'Kidney',	'Tm',	'Kt', 'Free')
    rownames(loaded_params) <- c('PFBS', 'PFHxS', 'PFOS', 'PFDS', 'PFHxA', 'PFHpA',
                                 'PFOA', 'PFNA', 'PFDA', 'PFUnDA', 'PFTeDA')
    
    keep_params <- loaded_params[which(rownames(loaded_params)==substance),1:dim(loaded_params)[2]]
    
    Tm = keep_params$Tm #ug/h Fabrega (2015) Table 2
    Kt = keep_params$Kt #ug/L Fabrega (2015) Table 2
    Free = keep_params$Free #unitless Fabrega (2015) Table 2
    PL = keep_params$Liver # Partition Coefficient for Liver: Fabrega (2015) Table 2
    PF = 0.467 # Partition Coefficient for Fat: Fabrega (2014) Table 1
    PB = keep_params$Brain # Partition Coefficient for Brain: Fabrega (2015) Table 2
    PBm = keep_params$Bone_marrow # Partition Coefficient for Bone marrow: Fabrega (2015) Table 2
    PLu = keep_params$Lung # Partition Coefficient for Lung: Fabrega (2015) Table 2
    PK = keep_params$Kidney # Partition Coefficient for Kidney: Fabrega (2015) Table 2
    PG = 0.05 # Partition Coefficient for Gut: Loccisano (2011) Table 1
    PR = 0.12 # Partition Coefficient for Rest of body: Loccisano (2011) Table 1
    kurinec = 3e-04	#urinary elimination rate constant  (/h/kg^-0.25); 
    kurine = kurinec*BW**(-0.25) # Elimination rate (1/h)
    
    return(list('QC'=QC, 'QCP'=QCP, 'QL'=QL, 'QF'=QF, 'QK'=QK, 
                'QFil'=QFil, 'QG'=QG, 'QLu'=QLu, 'QB'=QB, 'QBm'=QBm, 'QR'=QR,
                'VPlas'=VPlas, 'VL'=VL, 'VF'=VF, 'VK'=VK, 
                'VFil'=VFil, 'VG'=VG, 'VLu'=VLu, 'VB'=VB, 'VBm'=VBm, 'VR'=VR,
                'PL'=PL, 'PF'=PF, 'PB'=PB, 'PBm'=PBm, 'PLu'=PLu, 'PK'=PK, 
                'PG'=PG, 'PR'=PR,
                'Tm'=Tm, 'Kt'=Kt, 'Free'=Free, 'kurine'=kurine,
                'admin.dose'=admin.dose, 'admin.time'=admin.time,
                'f_unabs'=f_unabs))
  })
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    
    APlas<-0; AG<-0; AL<-0; AF<-0; ALu<-0; AB<-0; ABm<-0
    AK<-0; AFil<-0; AStor<-0; AUrine<-0
    AR<-0
    
    return(c("APlas"=APlas, "AG"=AG, "AL"=AL, "AF"=AF, "ALu"=ALu, "AB"=AB, "ABm"=ABm,
             "AK"=AK, "AFil"=AFil, "AStor"=AStor, "AUrine"=AUrine,
             "AR"=AR))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    ltimes <- length(admin.time)
    lexposure <- length(admin.dose)
    
    events <- data.frame(var = rep('AG', ltimes),
                         time = admin.time,
                         value = (1-f_unabs)*admin.dose, method = rep('add',ltimes))
    return(list(data=events))
  })
}

custom.func <- function(){
  return()
}

ode.func <- function(time, inits, params, custom.func){
  with(as.list(c(inits,params)),{
    
    # Units
    # C_i: Concentration in tissue i (ug/kg) #cosidering density = 1g/ml
    # Intake: hourly ingestion of PFASs (ug/h)
    
    # Concentrations (ng/g)
    CPlas <- APlas/VPlas #total concentration of chemical in plasma
    CG <- AG/VG # Concentration in Gut
    CL <- AL/VL # Concentration in Liver
    CF <- AF/VF # Concentration in Fat
    CLu <- ALu/VLu # Concentration in Lungs
    CB <- AB/VB # Concentration in Brain
    CBm <- ABm/VBm # Concentration in Bone marrow
    CK <- AK/VK # Concentration in Kidney
    CFil <- AFil/VFil # Concentration in Filtration compartment
    CR <- AR/VR  # Concentration in Rest of the body compartment
    
    # Plasma compartment (ng)
    dAPlas <- QF*Free*CF/PF + (QL+QG)*Free*CL/PL + QR*Free*CR/PR + QK*Free*CK/PK + 
      QB*Free*CB/PB + QBm*Free*CBm/PBm + QLu*Free*CLu/PLu - QCP*Free*CPlas 
    
    # Gut compartment
    dAG <- QG*Free*(CPlas - CG/PG)
    
    # Liver compartment
    dAL <- QL*Free*CPlas + QG*Free*CG/PG - (QL+QG)*Free*CL/PL
    
    # Fat compartment
    dAF <- QF*Free*(CPlas - CF/PF)
    
    # Lungs Compartment
    dALu <- QLu*Free*(CPlas - CLu/PLu)
    
    # Brain Compartment
    dAB <- QB*Free*(CPlas - CB/PB)
    
    # Bone Marrow Compartment
    dABm <- QBm*Free*(CPlas - CBm/PBm)
    
    # Kidney compartment
    dAK <- QK*Free*(CPlas - CK/PK) + (Tm*CFil)/(Kt+CFil)
    
    # Filtrate compartment
    dAFil = QFil*(Free*CPlas - CFil) - (Tm*CFil)/(Kt+CFil)
    
    # Storage compartment for urine
    dAStor <- QFil*CFil - kurine*AStor
    
    # Amount excrete via urine
    dAUrine <- kurine*AStor
    
    # Rest of the body
    dAR <- QR*Free*(CPlas - CR/PR)
    
    
    return(list(c("dAPlas"=dAPlas, "dAG"=dAG, "dAL"=dAL, "dAF"=dAF, "dALu"=dALu, 
                  "dAB"=dAB, "dABm"=dABm, "dAK"=dAK, "dAFil"=dAFil, "dAStor"=dAStor,
                  "dAUrine"=dAUrine, "dAR"=dAR),
                "CPlas"=CPlas, "CG"=CG, "CL"=CL, "CF"=CF, "CLu"=CLu, "CB"=CB, "CBm"=CBm,
                "CK"=CK, "CFil"=CFil, "CR"=CR))
  })
}

########################################
# PFBS - Total intake
daily_intake <- (0.33+0.01+0.004)*70 # ng of PFAS per day
BW <- 70 # kg
substance <- 'PFBS'
f_unabs <- 0
admin.dose <- rep(daily_intake, 40*365 ) # administered dose in ug
admin.time <- seq(0, 40*365*24-1, 24) # time when doses are administered, in hours
user_input <- list('BW'=BW,
                   'substance'=substance,
                   "f_unabs"=f_unabs,
                   "admin.dose"=admin.dose,
                   "admin.time"= admin.time)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

sample_time <- seq(0,40*365*24-1,24)
start.time <- Sys.time()
total_solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                           events = events, 
                           method="lsodes",rtol = 1e-05, atol = 1e-05)) 
end.time <- Sys.time()
print(end.time-start.time)

########################################
# PFBS - Dietary intake
daily_intake <- (0.33)*70 # ng of PFAS per day
BW <- 70 # kg
substance <- 'PFBS'
f_unabs <- 0
admin.dose <- rep(daily_intake, 40*365 ) # administered dose in ug
admin.time <- seq(0, 40*365*24-1, 24) # time when doses are administered, in hours
user_input <- list('BW'=BW,
                   'substance'=substance,
                   "f_unabs"=f_unabs,
                   "admin.dose"=admin.dose,
                   "admin.time"= admin.time)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

sample_time <- seq(0,40*365*24-1,24)
start.time <- Sys.time()
dietary_solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                 events = events, 
                                 method="lsodes",rtol = 1e-05, atol = 1e-05)) 
end.time <- Sys.time()
print(end.time-start.time)



########################################
# PFBS - dwater intake
daily_intake <- (0.01)*70 # ng of PFAS per day
BW <- 70 # kg
substance <- 'PFBS'
f_unabs <- 0
admin.dose <- rep(daily_intake, 40*365 ) # administered dose in ug
admin.time <- seq(0, 40*365*24-1, 24) # time when doses are administered, in hours
user_input <- list('BW'=BW,
                   'substance'=substance,
                   "f_unabs"=f_unabs,
                   "admin.dose"=admin.dose,
                   "admin.time"= admin.time)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

sample_time <- seq(0,40*365*24-1,24)
start.time <- Sys.time()
dwater_solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                   events = events, 
                                   method="lsodes",rtol = 1e-05, atol = 1e-05)) 
end.time <- Sys.time()
print(end.time-start.time)


########################################
# PFOA - Dust intake
daily_intake <- (0.004)*70 # ng of PFAS per day
BW <- 70 # kg
substance <- 'PFBS'
f_unabs <- 0
admin.dose <- rep(daily_intake, 40*365 ) # administered dose in ug
admin.time <- seq(0, 40*365*24-1, 24) # time when doses are administered, in hours
user_input <- list('BW'=BW,
                   'substance'=substance,
                   "f_unabs"=f_unabs,
                   "admin.dose"=admin.dose,
                   "admin.time"= admin.time)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

sample_time <- seq(0,40*365*24-1,24)
start.time <- Sys.time()
dust_solution <- data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                                   events = events, 
                                   method="lsodes",rtol = 1e-05, atol = 1e-05)) 
end.time <- Sys.time()
print(end.time-start.time)

#----------
# PLOTS
#----------

color_codes <- scales::hue_pal()(4)
names(color_codes) <-  c("Total", "Water", "Food", "Dust")

library(ggplot2)
PFBS_plot <- ggplot()+
  geom_line(data = dietary_solution, aes(x = time/24/365, y = CL, color='Food'), size=1.3)+
  geom_line(data = dwater_solution, aes(x = time/24/365, y = CL, color='Water'), size=1.3)+
  geom_line(data = dust_solution, aes(x = time/24/365, y = CL, color='Dust'), size=1.3)+
  geom_line(data = total_solution, aes(x = time/24/365, y = CL, color='Total'), size=1.3)+
  
  labs(title = paste0('Concenrtation of PFBS in liver for different exposures'),
       y = 'Liver Concentration (ng/g)' , x = "Time (years)")+
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.y=element_text(size=22),
        axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
        axis.text.x=element_text(size=22),
        legend.title=element_text(hjust = 0.5,size=25), 
        legend.text=element_text(size=22),
        panel.border = element_rect(colour = "black", fill=NA, size=1.0)) + 
  
  scale_y_log10()+
  scale_x_continuous(limits=c(0, 40))+ 
  scale_color_manual("Exposure", values=color_codes)+
  theme(legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        axis.text = element_text(size = 14))
PFBS_plot

setwd('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/Exposure_Scenarios')
jpeg(file = paste0(getwd(), "/PFBS_plot.jpeg"),
     width = 16, height = 12, units = 'in', res = 150)
plot(PFBS_plot)
dev.off()