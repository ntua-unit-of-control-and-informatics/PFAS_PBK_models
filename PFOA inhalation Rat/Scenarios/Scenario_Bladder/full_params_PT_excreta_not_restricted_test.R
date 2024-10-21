
setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Scenarios/Scenario_Bladder/Training/AAFE/full_params_PT_excreta_not_restricted_test")
load("full_params_PT_excreta_not_restricted.RData")

change_parms <- function(params){
  with(as.list(params),{
    
  RAFOatp_k <- estimated_params[2]
  RAFOat1 <- estimated_params[4]
  RAFbile_transp <- estimated_params[6]
  
  RAFOatp_l <- estimated_params[7]
  RAFUrat <- RAFOatp_k
  RAFOat3 <- RAFOat1
  RAFOatp2_l <- RAFOatp_l
  RAFOatp_lu_ap <- estimated_params[8]
  RAFOatp_lu_bas <- RAFOatp_lu_ap
  RAFNtcp <- RAFOatp_l
  RAFOatp2_Int <- estimated_params[9]
  
  RAF_papp <- estimated_params[10]
  
  f_fabp_avail <- estimated_params[11]
  f_alb_avail <- estimated_params[12]
  
  koff_alb <-0.1
  koff_fabp <-  koff_alb
  koff_a2u <- koff_alb
  
  VmK_api <- 0
  VmK_baso <- 0
  KmK_baso <- 1e20
  KmK_api <-   1e20
  KLfabp <- (1.2e5+4e4+1.9e4)  #[L/mol]*1e-3 , value from Cheng et al. (2017)
  Ka <- 5.8e05 #  from Rue et
  
  return(list('VB'=VB, 'Vplasma'=Vplasma, 'VK'=VK, 'VKB'=VKB, 
              'VKF'=VKF, 'VKT'=VKT, 'VFil_rest'=VFil_rest, 'VPT' = VPT, 'VBladder' = VBladder,
              'VL'=VL, 'VLB'=VLB, 'VLF'=VLF, 'VLT'=VLT, 'VLbile'=VLbile,
              'VM'=VM, 'VMB'=VMB, 'VMF'=VMF, 'VMT'=VMT, 'VA'=VA, 'VAB'=VAB, 
              'VAF'=VAF, 'VAT'=VAT, 'VR'=VR, 'VRB'=VRB, 
              'VRF'=VRF, 'VRT'=VRT, 'VVen' = VVen,
              'VArt' = VArt, 'VLu'=VLu, 'VLuB'=VLuB, 'VLuF'=VLuF,
              'VLuAF'=VLuAF, 'VLuT'=VLuT,
              'VSP'=VSP, 'VSPB'=VSPB, 'VSPF'=VSPF, 'VSPT'=VSPT,
              'VH'=VH, 'VHB'=VHB, 'VHF'=VHF, 'VHT'=VHT,
              'VBr'=VBr, 'VBrB'=VBrB, 'VBrF'=VBrF, 'VBrT'=VBrT,
              'VGo'=VGo, 'VGoB'=VGoB, 'VGoF'=VGoF, 'VGoT'=VGoT,
              'VIN'=VIN, 'VINB'=VINB, 'VINF'=VINF, 'VINT'=VINT,
              'VST'=VST, 'VSTB'=VSTB, 'VSTF'=VSTF, 'VSTT'=VSTT,
              'VSTL'=VSTL, 'VINL'=VINL,
              'VSK'=VSK,'VSKB'=VSKB, 'VSKF'=VSKF, 'VSKT'=VSKT,
              'VBo'=VBo,'VBoB'=VBoB, 'VBoF'=VBoF, 'VBoT'=VBoT,
              
              "VPT" = VPT, "VTDL"= VTDL, "VThinAL"=VThinAL, "VThickAL"=VThickAL,
              "VDT" = VDT, "Vduct"=Vduct,
              
              'AK'=AK, 'AKG'=AKG, 'AL'=AL, 
              'AM'=AM, 'AA'=AA, 'AR'=AR, 'ALu'= ALu, 
              'ASP'=ASP, 'AH'=AH, 'ABr'=ABr, 'AST'= AST,
              'AIN'=AIN, 'AGo'=AGo,
              'ASK'= ASK, 'ABo'=ABo,
              
              "SKi" = SKi,"SLi" = SLi,"SSt" = SSt,"SIn" = SIn,
              "SMu" = SMu,"SAd" = SAd,"SRe" = SRe,"SLu" = SLu,
              "SSp" = SSp,"SHt" = SHt,"SBr" = SBr,"SGo" = SGo,
              "SSk" = SSk,"SBo" = SBo,
              
              
              'PeffK'=PeffK, 'PeffL'=PeffL, 
              'PeffA'=PeffA, 'PeffM'=PeffM, 'PeffR'=PeffR, 'PeffLu' = PeffLu,
              'PeffSP'=PeffSP, 'PeffH'=PeffH, 'PeffBr'=PeffBr, 'PeffST' = PeffST,
              'PeffIN'=PeffIN, 'PeffGo'=PeffGo,
              'PeffSK' = PeffSK,  'PeffBo' = PeffBo,  
              
              'Qcardiac'=Qcardiac, 'QBK'=QBK, 
              'QBL'=QBL, 'QBLtot'=QBLtot,
              'QBM'=QBM, 'QBA'=QBA,
              'QBR'=QBR, 'QBLu'=QBLu, 'Qfeces'=Qfeces, 'feces_density'=feces_density,
              'Qbile'=Qbile, 'QGFR'=QGFR,'Qurine'=Qurine, 
              'QBSP'=QBSP, 'QBH'=QBH, 'QBBr'=QBBr, 'QBST'=QBST,
              'QBIN'=QBIN, 'QGE'=QGE,
              'QBGo'=QBGo,
              'QBSK'=QBSK, 'QBBo'=QBBo, 'Hct' = Hct,
              
              "QparaKi" = QparaKi,"QparaLi" = QparaLi,"QparaSt" = QparaSt,"QparaIn" = QparaIn,
              "QparaMu" = QparaMu,"QparaAd" = QparaAd,"QparaRe" = QparaRe,"QparaLu" = QparaLu,
              "QparaSp" = QparaSp,"QparaHt" = QparaHt,"QparaBr" = QparaBr,"QparaGo" = QparaGo,
              "QparaSk" = QparaSk,"QparaBo" = QparaBo,
              
              "QPT" = QPT, 'QTDL' = QTDL,
              # "QTDL"= QTDL, "QTAL"= QTAL, "QThinAL"=QThinAL, "QThickAL"=QThickAL,
              # "QDT" = QDT, "Qduct"=Qduct,
              
              'CalbB_init'= CalbB_init, 'CalbKF_init'=CalbKF_init, 'CalbLF_init'=CalbLF_init,
              'CalbMF_init'=CalbMF_init, 'CalbAF_init'=CalbAF_init, 'CalbRF_init'=CalbRF_init,
              'CalbBoF_init'=CalbBoF_init, 'CalbLuF_init' =CalbLuF_init,
              'CalbLuAF_init'=CalbLuAF_init, 'CalbSPF_init' =CalbSPF_init,
              'CalbGoF_init' =CalbGoF_init, 'CalbHF_init' =CalbHF_init,
              'CalbBrF_init' =CalbBrF_init, 'CalbSTF_init' =CalbSTF_init,
              'CalbINF_init' =CalbINF_init, 'CalbSKF_init' =CalbSKF_init, 
              
              'CalbKB_init'=CalbKB_init,'CalbLB_init'=CalbLB_init,'CalbSTB_init'=CalbSTB_init,
              'CalbINB_init'=CalbINB_init, 'CalbMB_init'=CalbMB_init,'CalbAB_init'=CalbAB_init,
              'CalbRB_init'=CalbRB_init,'CalbBoB_init'=CalbBoB_init,
              'CalbLuB_init'=CalbLuB_init, 'CalbSPB_init'=CalbSPB_init,'CalbGoB_init'=CalbGoB_init,
              'CalbHB_init'=CalbHB_init,'CalbBrB_init'=CalbBrB_init,'CalbSKB_init'=CalbSKB_init,
              
              'Ca2uKT_init'=Ca2uKT_init,'CFabpKT_init'=CFabpKT_init,'CFabpLT_init'=CFabpLT_init, 
              
              'Ka'=Ka, 'Ka2u'=Ka2u, 'KLfabp'=KLfabp,
              
              "koff_alb" = koff_alb, "koff_a2u" = koff_a2u, "koff_fabp" = koff_fabp,
              "kon_alb" = kon_alb, "kon_a2u" = kon_a2u, "kon_fabp" = kon_fabp,
              
              'VmL_Oatp'=VmL_Oatp, 'KmL_Oatp'= KmL_Oatp, 'VmL_Ntcp'= VmL_Ntcp,
              'VmL_Oatp2'=VmL_Oatp2, 'KmL_Oatp2'= KmL_Oatp2, 
              'VmIn_Oatp2'=VmIn_Oatp2, 'KmIn_Oatp2'= KmIn_Oatp2,
              'KmL_Ntcp'= KmL_Ntcp,'VmK_Oatp'= VmK_Oatp, 
              'Vmbile_transp' = Vmbile_transp,    'Kmbile_transp' = Kmbile_transp,
              'KmK_Oatp'=KmK_Oatp, 'VmLu_Oatp_ap'= VmLu_Oatp_ap,
              'KmLu_Oatp_ap'=KmLu_Oatp_ap, 'VmLu_Oatp_bas'= VmLu_Oatp_bas,
              'KmLu_Oatp_bas'=KmLu_Oatp_bas,
              
              
              'VmK_Oat1'=VmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 'KmK_Oat1'=KmK_Oat1, 
              'VmK_Oat3'=VmK_Oat3, 'KmK_Oat3'=KmK_Oat3, 
              'VmK_Urat'=VmK_Urat, 'KmK_Urat'=KmK_Urat, 
              
              'KmK_baso' = KmK_baso, 'KmK_api' = KmK_api,
              'VmK_baso' = VmK_baso,'VmK_api' = VmK_api,
              
              'Papp' = Papp, 'P_passive' = P_passive,
              'kKFKT'=kKFKT, 'kFilKT'=kFilKT, 'kPTKT' = kPTKT,  
              'kLFLT'=kLFLT, 'kLTLbile'=kLTLbile,  'kAFAT'=kAFAT, 
              'kRFRT'=kRFRT,
              'kabST'=kabST, 
              'kMFMT'=kMFMT, 'kLuTLuF' =kLuTLuF, 'kLuTLuAF'=kLuTLuAF, 'kSPFSPT' =kSPFSPT,
              'kSTFSTT' =kSTFSTT, 'kINFINT' =kINFINT, 'kHFHT' =kHFHT,
              'kBrFBrT' =kBrFBrT, 'kGoFGoT' =kGoFGoT,
              'kSKFSKT' =kSKFSKT, 'kBoFBoT'=kBoFBoT,
              
              "admin.time" = admin.time, "admin.dose" = admin.dose,
              "admin.type" = admin.type, "MW"=MW
  ))
  
  })
}

setwd("C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat")

MW <- 414.07 #g/mol

# Read data
kudo_high_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_high_kudo_2007.xlsx")
kudo_low_dose <- openxlsx::read.xlsx("Data/IV_male_rats_tissues_low_kudo_2007.xlsx")
kim_IV_Mtissues <- openxlsx::read.xlsx("Data/PFOA_male_tissues_IV_kim_2016.xlsx")
kim_OR_Mtissues <- openxlsx::read.xlsx("Data/PFOA_male_tissues_ORAL_kim_2016.xlsx")
kim_IV_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_IV_kim_2016.xlsx")
kim_OR_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_ORAL_kim_2016.xlsx")
dzi_OR_Mtissues <- openxlsx::read.xlsx("Data/Dzierlenga_tissue_male_ORAL_2021.xlsx")
dzi_OR_Mtissues$Concentration_microM <- dzi_OR_Mtissues$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Ftissues <- openxlsx::read.xlsx("Data/Dzierlenga_tissue_female_ORAL_2021.xlsx")
dzi_OR_Ftissues$Concentration_microM <- dzi_OR_Ftissues$Concentration_microM* MW/1000 #convert from uM to ug/g
kim_OR_Mblood <- openxlsx::read.xlsx("Data/PFOA_male_blood_ORAL_kim_2016.xlsx")
kim_IV_Mblood <- openxlsx::read.xlsx("Data/PFOA_male_blood_IV_kim_2016.xlsx")
Lup_OR_Ftissues <- openxlsx::read.xlsx("Data/PFOA_female_tissues_Lupton_2020.xlsx")
Kemp_OR_Ffeces <- openxlsx::read.xlsx("Data/PFOA_Feces_female_oral_25_mg_per_kg-Loc.xlsx")
Kemp_OR_Mfeces <- openxlsx::read.xlsx("Data/PFOA_Feces_male_oral_25_mg_per_kg-Loc.xlsx")
Kemp_OR_Furine_low <- openxlsx::read.xlsx("Data/PFOA_Urine_female_oral_1_mg_per_kg.xlsx")
Kemp_OR_Furine_med <- openxlsx::read.xlsx("Data/PFOA_Urine_female_oral_5_mg_per_kg.xlsx")
Kemp_OR_Furine_high <- openxlsx::read.xlsx("Data/PFOA_Urine_female_oral_25_mg_per_kg.xlsx")
Kemp_OR_Murine_low <- openxlsx::read.xlsx("Data/PFOA_Urine_male_oral_1_mg_per_kg.xlsx")
Kemp_OR_Murine_med <- openxlsx::read.xlsx("Data/PFOA_Urine_male_oral_5_mg_per_kg.xlsx")
Kemp_OR_Murine_high <- openxlsx::read.xlsx("Data/PFOA_Urine_male_oral_25_mg_per_kg.xlsx")
dzi_IV_Mserum <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_IV_2021.xlsx")
dzi_IV_Mserum$Concentration_microM <- dzi_IV_Mserum$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_OR_Mserum_low <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_ORAL_low_2021.xlsx")
dzi_OR_Mserum_low$Concentration_microM <- dzi_OR_Mserum_low$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Mserum_medium <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_ORAL_medium_2021.xlsx")
dzi_OR_Mserum_medium$Concentration_microM <- dzi_OR_Mserum_medium$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_OR_Mserum_high <- openxlsx::read.xlsx("Data/Dzierlenga_serum_male_ORAL_high_2021.xlsx")
dzi_OR_Mserum_high$Concentration_microM <- dzi_OR_Mserum_high$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_IV_Fserum <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_IV_2021.xlsx")
dzi_IV_Fserum$Concentration_microM <- dzi_IV_Fserum$Concentration_microM* MW/1000 #convert from uM to ug/g
dzi_OR_Fserum_low <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_ORAL_low_2021.xlsx")
dzi_OR_Fserum_low$Concentration_microM <- dzi_OR_Fserum_low$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Fserum_medium <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_ORAL_medium_2021.xlsx")
dzi_OR_Fserum_medium$Concentration_microM <- dzi_OR_Fserum_medium$Concentration_microM * MW/1000 #convert from uM to ug/g
dzi_OR_Fserum_high <- openxlsx::read.xlsx("Data/Dzierlenga_serum_female_ORAL_high_2021.xlsx")
dzi_OR_Fserum_high$Concentration_microM <- dzi_OR_Fserum_high$Concentration_microM * MW/1000 #convert from uM to ug/g
kim_OR_Fblood <- openxlsx::read.xlsx("Data/PFOA_female_blood_ORAL_kim_2016.xlsx")
kim_IV_Fblood <- openxlsx::read.xlsx("Data/PFOA_female_blood_IV_kim_2016.xlsx")
gus_OR_Mblood <- openxlsx::read.xlsx("Data/Gustafsson 2022_PFOA_Plasma Male rats_Oral.xlsx")
gus_OR_Mtissues <- openxlsx::read.xlsx("Data/Gustafsson 2022_PFOA_Tissues Male rats_Oral.xlsx")

dataset <- list("df1" = kudo_high_dose, "df2" = kudo_low_dose, "df3" = kim_IV_Mtissues, "df4" = kim_OR_Mtissues,
                "df5" = kim_IV_Ftissues, "df6" = kim_OR_Ftissues, "df7" = dzi_OR_Mtissues, "df8" = dzi_OR_Ftissues,
                "df9" = kim_OR_Mblood, "df10" = kim_IV_Mblood, "df11" = Lup_OR_Ftissues, "df12" = Kemp_OR_Ffeces,
                "df13" = Kemp_OR_Mfeces, "df14" = Kemp_OR_Furine_low, "df15" = Kemp_OR_Furine_med, "df16" = Kemp_OR_Furine_high,
                "df17" = Kemp_OR_Murine_low, "df18" = Kemp_OR_Murine_med, "df19" = Kemp_OR_Murine_high, 
                "df20" = dzi_IV_Mserum, "df21" = dzi_OR_Mserum_low, "df22" = dzi_OR_Mserum_medium,
                "df23" = dzi_OR_Mserum_high, "df24" = dzi_IV_Fserum, "df25" = dzi_OR_Fserum_low, "df26" = dzi_OR_Fserum_medium,
                "df27" = dzi_OR_Fserum_high, "df28" = kim_OR_Fblood, "df29" = kim_IV_Fblood, "df30" = gus_OR_Mblood,
   "df31" = gus_OR_Mtissues)


# Set up simulations for the 5th case, i.e. kim (2016) IV female tissues
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)

params <- change_parms(params)

sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_IV_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]

# Set up simulations for the 6th case, i.e. kim (2016) ORAL female tissues
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
params <- change_parms(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Clungs", 
                                       "Cspleen", "Cheart")]


# Set up simulations for the 8th case, i.e. Dzierlenga (2021) ORAL female tissues
BW <- 0.2  # body weight (kg) not reported
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
params <- change_parms(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Ftissues <-  solution[, c("time", "Cliver","Ckidney", "Cbrain")]




##########################
#-------------------------
# Kemper 2003 (Loccisano)
#-------------------------
##########################

# Set up simulations for the 12th case, i.e.Kemper 2003 (Loccisano) ORAL female feces

sex <- "F"
BW <- 0.2 #kg
sample_time <- seq(0,672,1)
admin.type <-"oral"
admin.dose <- 25 * BW*1000 #ug
admin.time <- 0

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
params <- change_parms(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=c(0, 0.25, 0.5, 1, 1.5, 2, seq(4,100,2), seq(104,672,4))

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_Kemp_OR_Ffeces <-  solution[, c("time", "Mfeces")]


##########################
#-------------------------
# Kemper 2003 (Worley)
#-------------------------
##########################

# Set up simulations for the 14th case, i.e.Kemper 2003 (Worley) ORAL female urine LOW

sex <- "F"
BW <- 0.2 #kg
sample_time <- c(0, 0.25, 0.5, 1, 1.5, 2, seq(4,168,1), 169)
admin.type <-"oral"
admin.time <- 0

#Female, oral 1mg/kg dose
admin.dose <- 1 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
parameters <- change_parms(parameters)

events <- create.events(parameters)
inits <- create.inits (parameters)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))


preds_Kemp_OR_Furine_low <-  solution[, c("time", "Murine")]


# Set up simulations for the 15th case, i.e.Kemper 2003 (Worley) ORAL female urine MEDIUM

admin.dose <- 5 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
events <- create.events(parameters)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))


preds_Kemp_OR_Furine_med <-  solution[, c("time", "Murine")]


# Set up simulations for the 16th case, i.e.Kemper 2003 (Worley) ORAL female urine HIGH

admin.dose <- 25 * BW * 1000 #ug
parameters <-   create.params(list('BW'=BW,
                                   "admin.dose"= admin.dose,
                                   "admin.time" = admin.time, 
                                   "admin.type" = admin.type,
                                   "estimated_params" = estimated_params,
                                   "sex" = sex))
parameters <- change_parms(parameters)

events <- create.events(parameters)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = parameters, events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))


preds_Kemp_OR_Furine_high <-  solution[, c("time", "Murine")]



# Set up simulations for the 24th case, i.e. Dzierlenga 2021, IV female serum 
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
params <- change_parms(params)

inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, seq(4, 192, 2))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_IV_Fserum <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 25th case, i.e. Dzierlenga 2021, ORAL female serum low
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 40 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
params <- change_parms(params)

inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0, 96, 1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_low <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 26th case, i.e. Dzierlenga 2021, ORAL female serum medium
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
params <- change_parms(params)

inits <- create.inits(params)
events <- create.events(params)


sample_time <- c(0, 0.25, 0.5, 1, 1.5,2, seq(4, 192, 2))

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_medium <-  solution[, c("time", "Cplasma")]



# Set up simulations for the 27th case, i.e. Dzierlenga 2021, ORAL female serum high
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 320 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"  


user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
params <- change_parms(params)

inits <- create.inits(params)
events <- create.events(params)


sample_time <- seq(0, 96, 1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_dzi_OR_Fserum_high <-  solution[, c("time", "Cplasma")]



#################################################################################
#--------------------------------------------------------------------------------
#                                Kim 2016 female
#-------------------------------------------------------------------------------
#################################################################################

# Set up simulations for the 28th case, i.e. Kim (2016) ORAL male blood
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
params <- change_parms(params)

inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_OR_Fblood <-  solution[, c("time", "Cplasma")]


# Set up simulations for the 29th case, i.e. Kim (2016) IV male blood
BW <- 0.25  # body weight (kg) not reported
admin.dose_per_g <- 1 # administered dose in mg PFOA/kg BW 
admin.dose <- admin.dose_per_g*BW *1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "iv"
sex <- "F"

user_input <- list('BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "estimated_params" = estimated_params,
                   "sex" = sex)


params <- create.params(user_input)
params <- change_parms(params)

inits <- create.inits(params)
events <- create.events(params)


sample_time=seq(0,24,1)

solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,events = events,
                                    method="lsodes",rtol = 1e-03, atol = 1e-03))

preds_kim_IV_Fblood <-  solution[, c("time", "Cplasma")]






#convert ug/L, which is returned by the ODE function, to ug/g, which are the units in all the datasets

preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] <- preds_kim_IV_Ftissues[,2:dim(preds_kim_IV_Ftissues)[2]] /1000
preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] <- preds_kim_OR_Ftissues[,2:dim(preds_kim_OR_Ftissues)[2]] /1000
preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] <- preds_dzi_OR_Ftissues[,2:dim(preds_dzi_OR_Ftissues)[2]] /1000
preds_Kemp_OR_Ffeces[,2:dim(preds_Kemp_OR_Ffeces)[2]] <- preds_Kemp_OR_Ffeces[,2:dim(preds_Kemp_OR_Ffeces)[2]] /1000
preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] <- preds_Kemp_OR_Furine_low[,2:dim(preds_Kemp_OR_Furine_low)[2]] /1000
preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] <- preds_Kemp_OR_Furine_med[,2:dim(preds_Kemp_OR_Furine_med)[2]] /1000
preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] <- preds_Kemp_OR_Furine_high[,2:dim(preds_Kemp_OR_Furine_high)[2]] /1000
preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] <- preds_dzi_IV_Fserum[,2:dim(preds_dzi_IV_Fserum)[2]] /1000
preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] <- preds_dzi_OR_Fserum_low[,2:dim(preds_dzi_OR_Fserum_low)[2]] /1000
preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] <- preds_dzi_OR_Fserum_medium[,2:dim(preds_dzi_OR_Fserum_medium)[2]] /1000
preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] <- preds_dzi_OR_Fserum_high[,2:dim(preds_dzi_OR_Fserum_high)[2]] /1000
preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] <- preds_kim_OR_Fblood[,2:dim(preds_kim_OR_Fblood)[2]] /1000
preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] <- preds_kim_IV_Fblood[,2:dim(preds_kim_IV_Fblood)[2]] /1000


# ######################################################################################
#Plot the predictions against the observations
library(ggplot2) 

# Function that creates a plot given a compartment name and the respective predictions and observations
create.plots <- function(predictions, observations, compartment){  
  #Colours of observations and predictions
  cls <-  c("predictions" = "#56B4E9",  "Observations" = "#D55E00")
  
  ggplot(data = predictions)+
    geom_line( aes_string(x= "Time", y= rlang::expr(!!compartment), 
                          color = '"predictions"'),  size=1.5,alpha = 0.7) +
    geom_point(data=observations, aes_string(x="Time", y= rlang::expr(!!compartment), 
                                             color='"Observations"'), size=4)+
    labs(title = rlang::expr(!!compartment), 
         y = expression("PFOA concentration (" * mu* "g/g tissue)" ),
         x = "Time (hours)")+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual("", values=cls,
                       guide = guide_legend(override.aes =
                                              list(shape = c(16,NA),
                                                   linetype = c(0,1))))+
    theme_light() + 
    theme(legend.position=c(1,1), 
          legend.justification=c(0, 1), 
          legend.key.size = unit(1.5, 'cm'),  
          legend.title = element_text(size=14),
          axis.title=element_text(size=14),
          legend.text = element_text(size=14)
    )
  
}



# Convert Kim IV female tissues from long to wide format using reshape
experiment5 <- reshape(kim_IV_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment5) <- c("Time",kim_IV_Ftissues$Tissue )

# Convert Kim ORAL female tissues from long to wide format using reshape
experiment6 <- reshape(kim_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment6) <- c("Time",kim_OR_Ftissues$Tissue )


# Convert Dzierlenga ORAL female tissues from long to wide format using reshape
experiment8 <- reshape(dzi_OR_Ftissues[c("Tissue" ,"Time_hours", 
                                         "Concentration_microM")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment8) <- c("Time",unique(dzi_OR_Ftissues$Tissue))
# Convert Kemper ORAL female feces from long to wide format using reshape
experiment12 <- reshape(Kemp_OR_Ffeces[c("Tissue" ,"Time_h", 
                                         "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment12) <- c("Time",unique(Kemp_OR_Ffeces$Tissue))
experiment12$Feces = (experiment12$Feces/100)*0.2*25



# Convert Kemper ORAL female urine low from long to wide format using reshape
experiment14 <- reshape(Kemp_OR_Furine_low [c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment14) <- c("Time",unique(Kemp_OR_Furine_low$Tissue))
experiment14$Urine = (experiment14$Urine/100)*0.2*1


# Convert Kemper ORAL female urine med from long to wide format using reshape
experiment15 <- reshape(Kemp_OR_Furine_med[c("Tissue" ,"Time_h", 
                                             "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment15) <- c("Time",unique(Kemp_OR_Furine_med$Tissue))
experiment15$Urine = (experiment15$Urine/100)*0.2*5

# Convert Kemper ORAL female urine high from long to wide format using reshape
experiment16 <- reshape(Kemp_OR_Furine_high[c("Tissue" ,"Time_h", 
                                              "Cum_dose_%")], 
                        idvar = "Time_h", timevar = "Tissue", direction = "wide")
colnames(experiment16) <- c("Time",unique(Kemp_OR_Furine_high$Tissue))
experiment16$Urine = (experiment16$Urine/100)*0.2*25

#Convert Dzierlenga 2021, IV female serum from long to wide format using reshape
experiment24 <- reshape(dzi_IV_Fserum[c("Tissue" ,"Time_hours", 
                                        "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment24) <- c("Time",unique(dzi_IV_Fserum$Tissue))

#Convert Dzierlenga 2021, ORAL female serum low from long to wide format using reshape
experiment25 <- reshape(dzi_OR_Fserum_low[c("Tissue" ,"Time_hours", 
                                            "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment25) <- c("Time",unique(dzi_OR_Fserum_low$Tissue))

#Convert Dzierlenga 2021, ORAL female serum medium from long to wide format using reshape
experiment26 <- reshape(dzi_OR_Fserum_medium[c("Tissue" ,"Time_hours", 
                                               "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment26) <- c("Time",unique(dzi_OR_Fserum_medium$Tissue))

#Convert Dzierlenga 2021, ORAL female serum high from long to wide format using reshape
experiment27 <- reshape(dzi_OR_Fserum_high[c("Tissue" ,"Time_hours", 
                                             "Concentration_microM")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment27) <- c("Time",unique(dzi_OR_Fserum_high$Tissue))

#Convert Kim 2016, ORAL female serum long to wide format using reshape
experiment28 <- reshape(kim_OR_Fblood[c("Tissue" ,"Time_hours", 
                                        "Concentration_microg_per_g_organ")], 
                        idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment28) <- c("Time",unique(kim_OR_Fblood$Tissue))


#Convert Kim 2016, IV female serum long to wide format using reshape
experiment29<- reshape(kim_IV_Fblood[c("Tissue" ,"Time_hours", 
                                       "Concentration_microg_per_g_organ")], 
                       idvar = "Time_hours", timevar = "Tissue", direction = "wide")
colnames(experiment29) <- c("Time",unique(kim_IV_Fblood$Tissue))


# Put the experiments in a list
experiments <- list(
                    experiment5 = experiment5, experiment6 = experiment6,  experiment8 = experiment8,
                     experiment12 = experiment12,
                   experiment14 = experiment14, experiment15 = experiment15, experiment16 = experiment16,
                    experiment24 = experiment24,
                    experiment25 = experiment25, experiment26 = experiment26, experiment27 = experiment27,
                    experiment28 = experiment28, experiment29 = experiment29)



colnames(preds_kim_IV_Ftissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")
colnames(preds_kim_OR_Ftissues) <- c( "Time", "Liver",  "Kidney", "Lung",
                                      "Spleen", "Heart")

colnames(preds_dzi_OR_Ftissues) <- c("Time","Liver","Kidney","Brain")

colnames(preds_Kemp_OR_Ffeces) <- c ("Time", "Feces")

colnames(preds_Kemp_OR_Furine_low) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_med) <- c ("Time", "Urine")
colnames(preds_Kemp_OR_Furine_high) <- c ("Time", "Urine")

colnames(preds_dzi_IV_Fserum) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_low) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_medium) <- c ("Time", "Serum")
colnames(preds_dzi_OR_Fserum_high) <- c ("Time", "Serum")

colnames(preds_kim_IV_Fblood) <- c ("Time", "Plasma")
colnames(preds_kim_OR_Fblood) <- c ("Time", "Plasma")



# Create a list containing the corresponding predictions
simulations <- list( predictions5 = preds_kim_IV_Ftissues, predictions6 = preds_kim_OR_Ftissues,
                   predictions8 = preds_dzi_OR_Ftissues,  predictions12 = preds_Kemp_OR_Ffeces,
       predictions14 = preds_Kemp_OR_Furine_low, predictions15 = preds_Kemp_OR_Furine_med,
                    predictions16 =preds_Kemp_OR_Furine_high, 
                    predictions24 =preds_dzi_IV_Fserum, predictions24 =preds_dzi_OR_Fserum_low, predictions26 =preds_dzi_OR_Fserum_medium,
                    predictions27 =preds_dzi_OR_Fserum_high, predictions28 = preds_kim_OR_Fblood, 
                    predictions29 = preds_kim_IV_Fblood)


# Iterate over all existing experiments and create the accompanying plots
for(i in 1:length(experiments)){
  # Retrieve the corresponding observations and simulations
  observations <- experiments[[i]]
  predictions <- simulations[[i]]
  # Extract the compartment names
  compartments <- names(predictions)[2:length(predictions)]
  
  # Use lapply to iterate over the column names and create plots
  plots <- lapply(compartments, function(compartment) {
    create.plots(predictions, observations, compartment )
  })
  if(length(compartments) == 1){
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 1, nrow = 1,
                                               common.legend = TRUE, legend = "right"))
    
  }else{
    final_plot <- do.call(ggpubr::ggarrange, c(plots, ncol = 3, nrow = ceiling(length(plots) / 3),
                                               common.legend = TRUE, legend = "right"))
  }
  
  
  plot.margin=unit(c(0,0,0,0), "pt")
  
  
  # Save the plot with dynamically adjusted dimensions
  ggsave(paste0("experiment", i,".png"), plot = final_plot,
         device = 'png', dpi = 300,
         width = 13,
         height = 10,
         units = "in")
}



