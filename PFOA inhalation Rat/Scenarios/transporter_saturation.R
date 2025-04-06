library(ggplot2)

BW <- 0.2  # body weight (kg) not reported, # body weight (kg) not reported, based on 8 week female rats from https://animal.ncku.edu.tw/p/412-1130-16363.php?Lang=en
admin.dose_per_g <- 80 # administered dose in mg PFOA/kg BW http://127.0.0.1:20627/graphics/plot_zoom_png?width=1216&height=614
admin.dose <- admin.dose_per_g*BW*1e03 #ug PFOA
admin.time <- 0 #time when doses are administered, in hours
admin.type <- "oral"
sex <- "F"
alpha = 1#0.01 #kPtcTu*(CPTCf - alpha*CPT) 
kappa = 1#5 #RAFOat3 <- estimated_params[3]*kappa
lambda = 1 #Ka <-5.8e5*lambda 
lambda_1 = 3#3 #f_fabp_avail <- estimated_params[9]*lambda_1
lambda_2 = 0.1  #estimated_params[10]*lambda_2
eta = 1 #Papp <- estimated_params[8]*eta
eta_2 = 1 # k_gut_in = ( (Papp_gut/100) *AINL)*1000*eta_2
gamma = 1 #CL_int <- estimated_params[4]*gamma 
gamma2 = 1 #RAFOatp_l <- estimated_params[5]*gamma2
gamma3 = 1#0.1 # f_fabp_avail <- estimated_params[9]*gamma3
gamma4 = 1 #KLfabp <-  gamma4*(1.2e5+4e4+1.9e4)
user_input <- list("alpha" = alpha,
                   "kappa" = kappa,
                   "lambda" = lambda, 
                   "lambda_1" = lambda_1, 
                   "lambda_2" = lambda_2, 
                   "eta" = eta, 
                   "eta_2" = eta_2, 
                   "gamma" = gamma,
                   "gamma2" = gamma2,
                   "gamma3" = gamma3,
                   "gamma4" = gamma4,
                   "t_lag" = 2.5,
                    'BW'=BW,
                   "admin.dose"= admin.dose,
                   "admin.time" = admin.time, 
                   "admin.type" = admin.type,
                   "sex" = sex,
                   "estimated_params" = estimated_params)

params <-create_params(user_input)
inits <- create.inits(params)
events <- create.events(params)

# sample_time: a vector of time points to solve the ODEs
sample_time=seq(0,192,0.1)

# ode(): The solver of the ODEs
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params,
                                    events = events,
                                    method="lsodes",rtol = 1e-04, atol = 1e-04))

CPT <- solution$CPT 
CKBf <- solution$CKBf
CPTCf <- solution$CPTCf
CPTCb <- solution$CPTCb
CKBb <- solution$CKBb
CKB = solution$CKB
KmK_Oatp <- params$KmK_Oatp
KmK_Urat<- params$KmK_Urat
KmK_Oat1<- params$KmK_Oat1
KmK_Oat3<- params$KmK_Oat3

oatp_loading <- (CPT/(KmK_Oatp+CPT))

urat_loading <- (CPT/(KmK_Urat+CPT))
  
  
oat1_loading <- (CKBf/(KmK_Oat1+CKBf))

oat3_loading <- (CKBf/(KmK_Oat3+CKBf))

df <- data.frame(CPTCf = solution$CPTCf ,CPT = CPT, CKBf =CKBf ,
                 CPTCb = CPTCb, CKBb = CKBb,
                 CKB = CKB, Ckidney = solution$Ckidney, CKTrestf=solution$CKTrestf,
                 CDALCf=solution$CDALCf, CDTCf=solution$CDTCf,  CCDCf=solution$CCDCf,
                 CKF = solution$CKF, Cplasma = solution$Cplasma,
                 CDAL =  solution$CDAL, CDT = solution$CDT,  CCD = solution$CCD,
                 MDAL = solution$MDAL, Cliver = solution$Cliver )/1000
  df$time <- solution$time
  
  
  ggplot(data = df[-1,])+
    geom_line(aes(x = time, y =Cplasma, colour = "Cplasma"))+
    geom_point(data= dzi_OR_Fserum_medium, aes(x = Time_hours , y =Concentration_microM ))+
    scale_y_continuous(trans='log10')
  
  
ggplot(data = df[2:241,])+
  geom_line(aes(x = time, y =Ckidney, colour = "Ckidney"))+
  geom_point(data= dzi_OR_Ftissues[dzi_OR_Ftissues$Tissue == "Kidney",], 
             aes(x = Time_hours , y =Concentration_microM ))+
  scale_y_continuous(trans='log10')



ggplot(data = df[2:241,])+
  geom_line(aes(x = time, y =Cliver, colour = "Cliver"))+
  geom_point(data= dzi_OR_Ftissues[dzi_OR_Ftissues$Tissue == "Liver",], 
             aes(x = Time_hours , y =Concentration_microM ))+
  scale_y_continuous(trans='log10')

tail(solution$CLT/solution$CLB)
tail(solution$CLT/solution$CLF)
tail(solution$CKT/solution$CKB)
tail(solution$CKT/solution$CKF)
  # geom_line(aes(x = time, y =CPT, colour = "CPT"))+
  # geom_line(aes(x = time, y =CDAL, colour = "CDAL"))+
  # geom_line(aes(x = time, y =CDT, colour = "CDT"))+
  # geom_line(aes(x = time, y =CCD, colour = "CCD"))
  # 
  # geom_line(aes(x = time, y =CPTCf, colour = "CPTCf"))+
  # geom_line(aes(x = time, y =CKB, colour = "CKB"))+
  #geom_line(aes(x = time, y =CPT, colour = "CPT"))+
  # geom_line(aes(x = time, y =CKF, colour = "CKF"))+
  # geom_line(aes(x = time, y =Ckidney, colour = "Ckidney"))+
  # geom_line(aes(x = time, y =CKTrestf, colour = "CKTrestf"))+
  # geom_line(aes(x = time, y =CDALCf, colour = "CDALCf"))+
  # geom_line(aes(x = time, y =CDTCf, colour = "CDTCf"))+
  # geom_line(aes(x = time, y =CCDCf, colour = "CCDCf"))+
    






# dMKBf = QBK*CArtf - QBK*CKBf   - PeffK*A_peritubular_PTC*(CKBf-CPTCf) -
#   PeffK*A_peritubular_DTC*(CKBf-CDTCf) - QparaKi*(1-SKi)*CKBf -
#   (VmK_Oat1*CKBf/(KmK_Oat1+CKBf)) - (VmK_Oat3*CKBf/(KmK_Oat3+CKBf))+ (VmK_baso*CPTCf/(KmK_baso+CPTCf))+
#   koff_alb*CKBb*VKB - kon_alb*CalbKBf*CKBf*VKB



# 
# #proximal tubule  cells subcompartment
# dMPTCf =  QparaKi*(1-SKi)*CKBf + PeffK*A_peritubular_PTC*(CKBf-CPTCf) + 
#   kPtcF*(CKFf-CPTCf) - kPtcTu*(CPTCf - CPT)  +
#   (VmK_Oatp*CPT/(KmK_Oatp+CPT)) + (VmK_Urat*CPT/(KmK_Urat+CPT))+
#   (VmK_Oat1*CKBf/(KmK_Oat1+CKBf)) + (VmK_Oat3*CKBf/(KmK_Oat3+CKBf)) - 
#   (VmK_baso*CPTCf/(KmK_baso+CPTCf)) -(VmK_api*CPTCf/(KmK_api+CPTCf))-
#   (kon_a2u*Ca2uPTCf*CPTCf*VPTC + kon_fabp*CFabpPTCf*CPTCf*VPTC -
#      koff_fabp*CPTCb*VPTC - koff_a2u*CPTCb*VPTC) 




# #Proximal convoluted tubule
# dMPT =  QGFR*CArtf + kPtcTu*(CPTCf - CPT) - (VmK_Oatp*CPT/(KmK_Oatp+CPT)) - 
#   (VmK_Urat*CPT/(KmK_Urat+CPT)) + (VmK_api*CPTCf/(KmK_api+CPTCf))- QTDL*CPT


term_kb_1 <- params$QBK * solution$CArtf
term_kb_2 <- params$QBK * CKBf
term_kb_3 <- -params$PeffK * params$A_peritubular_PTC* (CKBf-CPTCf)
term_kb_4 <- params$QparaKi*(1-params$SKi)*CKBf 
term_kb_5 <-  ( params$VmK_Oat1*CKBf/( params$KmK_Oat1+CKBf)) +  ( params$VmK_Oat3*CKBf/( params$KmK_Oat3+CKBf))
term_ptc_6<-  params$kPtcTu*(CPTCf - CPT)
term_ptc_7<-   (params$VmK_Oatp*CPT/(params$KmK_Oatp+CPT)) + (params$VmK_Urat*CPT/(params$KmK_Urat+CPT))




df2 <- data.frame(term_kb_1 = term_kb_1 ,term_kb_2 = term_kb_2, time = solution$time, 
                  term_kb_3 =term_kb_3 , term_kb_4 = term_kb_4, term_kb_5 = term_kb_5,
                  term_ptc_6 = term_ptc_6, term_ptc_7 = term_ptc_7)
library(ggplot2)
ggplot(df2)+
  geom_line(aes(x = time, y =term_kb_1, colour = "term_kb_1"))+
  geom_line(aes(x = time, y =term_kb_2, colour = "term_kb_2"))+
  geom_line(aes(x = time, y =term_kb_3, colour = "term_kb_3"))+
  geom_line(aes(x = time, y =term_kb_4, colour = "term_kb_4"))+
  geom_line(aes(x = time, y =term_kb_5, colour = "term_kb_5"))+
  geom_line(aes(x = time, y =term_ptc_6, colour = "term_ptc_6"))+
  geom_line(aes(x = time, y =term_ptc_7, colour = "term_ptc_7"))

Ccell <- 104381.8#ug/L
Cbulk <- 89080.57#ug/L

D <- (5.46*1e-6)*1e-4*3600 #[cm2/s]* *1e-4*3600 --> [m2/h] 
L <- 2.29e-05 #m, diameter of tube
kinematic_viscosity <- 0.6959*1e-6*3600 # [mm2/s ]*1e-6*3600 --> [m2/h], at 37oC from https://wiki.anton-paar.com/en/water/
Sc <- kinematic_viscosity/D
cross_sectional_area <- (pi*2.29e-05^2)/4 #m^2 
Q <- ((0.0599184+0.03896499)/2)*1e-3/76000 #[L/h]*1e-3 -->[m3/h]
velocity <- Q/cross_sectional_area #m/h
Re <- velocity*L/kinematic_viscosity
Pe <- velocity*L/D
Sh <- 2+0.552* Re^(1/2)*Sc^(1/3)
k_conv <- Sh*D/L 
#k_conv <- 0.01
Ccell <- CPTCf  #ug/L
Cbulk <- CPT  #ug/L
Papp <-0.3354964 #M/H
C_wall <- (k_conv*params$APT*Cbulk+ Papp*params$APT)/(k_conv*params$APT+ Papp*params$APT)



df3 <- data.frame(Ccell = Ccell ,Cbulk = Cbulk, time = solution$time, 
                  C_wall =C_wall)
library(ggplot2)
ggplot(df3)+
  geom_line(aes(x = time, y =Ccell, colour = "Ccell"))+
  geom_line(aes(x = time, y =Cbulk, colour = "Cbulk"))+
  geom_line(aes(x = time, y =C_wall, colour = "C_wall"))

Q <- params$QPT
P <-0.2127459
k <- ((1/Q)+ (1/P))^-1



  
