# The code replicates the Kim et al.(2019) PBK model
# and simulates IV or Oral PFDA distribution for female or male rats of constant or varying
# body weight. The Bernstein et al. (2021) code was used for replicating the
# model and the parameter set was expanded to also describe male rats
# (the original code only described female rats). 
# We did not include simultaneous IV and oral because there are no such experiments


library(deSolve)

#=========================
#1. Parameters of the model
#=========================

create.params  <- function(user_input){
  with( as.list(user_input),{
    
    #------------------------------------------------------------------------------
    # PARAMETERS for the model (which are independent of time).
    # Default values are zero
    #   unless the parameter is used as a divisor
    #   (in which case the default value is 1)
    
    # Experimental parameters.
    gut2liver <- 0;       # Option to have gut lumen connect directly to liver
    #   (assumes no gi blood flow; 1.0 for yes, connected)
    
    
    # Biochemical parameters
    F_unabs <-  F_unabs #0.650;        # Fraction of unabsorbed dose
    Bmax <- 0.0;           # Saturable binding in liver, maximum binding capacity (mg/h)
    K_b <- 1.0;            # Saturable binding in liver, affinity constant (mg/L)
    k_off <- 0.0;          # binding in liver, dissociation rate constant (/h)
    K_f <- 0;              # Rate of excretion from GI tissue to fecal storage (/h)
    K_unabs <- 0.0;        # Rate unabsorbed fraction goes to fecal storage (/h)
    K_t <- 30.051;            # Transporter affinity constant (mg/L)
    
    if(sex == "M"){
      T_mc <- 83953/0.25^0.75/1e06;           # Transporter maximum (mg/h)
      K_abs <- 3.84;          # Oral absorption rate (/h)
      K_uc <- 0.391/(0.25^(-0.25));           # Urinary elimination rate (/h)
      K_ustc <- 0.148/(0.25^(-0.25));         # Rate constant to urinary storage compartment (/h)
      K_fstc <- 0.382/(0.25^(-0.25));         # Rate constant from fecal storage to fecal elim. (/h)
      F_free <- 0.00118;       # Free fraction of chemical in blood
      K_bilec <- 0.000995/(0.25^(-0.25)); # Biliary excretion rate (/h) adjusted by BW
      
      #Bernstein et al. (2021) did not include male rats in their analysis. Thus,
      # we scaled the female values using the ratio of male to female partition coefficients
      # from digitizing the corresponding plot. Since the plot did not include any 
      # rob compartment, we left the latter value unchanged.
      # Partition coefficients.
      K_li <- 0.607 * 0.745/0.609;           # Liver/plasma partition coefficient
      K_ki <- 0.2328 * 0.170/0.234;           # Kidney/blood partition coefficient
      K_gi <- 0.0102 * 0.007/0.012;           # GI Tract/blood partition coefficient
      K_lu <-  0.1002 * 0.087/0.101;          # Lungs/blood partition coefficient
      K_ht <- 0.0436 * 0.0399/0.423;          # Heart/blood partition coefficient
      K_br <- 1.0;          # Brain/blood partition coefficient
      K_ad <- 1.0;          # Adipose/blood partition coefficient
      K_bm <- 1.0;          # Bone marrow/blood partition coefficient
      K_rb <-  0.0289 ;           # Rest of Body/blood partition coefficient
      
    }else if(sex == "F"){
      T_mc <- 0.11384419;           # Transporter maximum (mg/h)
      K_abs <- 1.440;          # Oral absorption rate (/h)
      K_uc <- 0.48153972;           # Urinary elimination rate (/h)
      K_ustc <-  0.53598694;         # Rate constant to urinary storage compartment (/h)
      K_fstc <- 0.27647875;         # Rate constant from fecal storage to fecal elim. (/h)
      F_free <- 0.000122;       # Free fraction of chemical in blood
      K_bilec <- 0.00119/(0.25^(-0.25)); # Biliary excretion rate (/h) adjusted by BW
      
      # Partition coefficients.
      K_li <- 0.607;           # Liver/plasma partition coefficient
      K_ki <- 0.2328;           # Kidney/blood partition coefficient
      K_gi <- 0.0102;           # GI Tract/blood partition coefficient
      K_lu <-  0.1002;          # Lungs/blood partition coefficient
      K_ht <- 0.0436;          # Heart/blood partition coefficient
      K_br <- 1.0;          # Brain/blood partition coefficient
      K_ad <- 1.0;          # Adipose/blood partition coefficient
      K_bm <- 1.0;          # Bone marrow/blood partition coefficient
      K_rb <-  0.0289;           # Rest of Body/blood partition coefficient
    }
    
    # Blood flow rates to compartments as fraction of cardiac output
    Q_cardiacc <- 7.29734198;
    Q_lic <- 0.32093023;
    Q_kic <- 0.20155039 ;
    Q_filc <- 0.10077519;
    Q_gic <- 0.17441860;
    Q_luc <- 1.0;
    Q_htc <- 0.09069767 ;
    Q_brc <- 0.0;
    Q_adc <- 0.0;
    Q_bmc <- 0.0;
    Q_rbc <- 0.20465116;
    
    # Compartment volumes as fraction of total body volume
    V_blc <- 0.0816;
    V_kic <- 0.0072;
    V_filc <- 0.00072;
    V_lic <- 0.0320;
    V_gic <- 0.040;
    V_luc <- 0.004;
    V_htc <- 0.0032;
    V_brc <- 1.0;
    V_adc <- 1.0;
    V_bmc <- 1.0;
    V_rbc <- 0.536 ;
    
    return(list("sex" = sex, "admin.type" = admin.type,
                "admin.dose" = admin.dose, 
                "admin.time" = admin.time,
                "BW" =  BW, "BW.times" = BW.times,
                "gut2liver" = gut2liver,
                "K_li" = K_li, "K_ki" = K_ki, "K_gi" = K_gi, "K_lu" = K_lu,
                "K_ht" = K_ht, "K_br" = K_br, "K_ad" = K_ad, "K_bm" = K_bm,
                "K_rb" = K_rb,
                "T_mc" = T_mc, "K_t" = K_t, "K_uc" = K_uc, "K_abs" = K_abs,
                "F_unabs" = F_unabs, "K_unabs" = K_unabs, "K_ustc" = K_ustc,
                "K_fstc" = K_fstc, "Bmax" = Bmax, "K_b" = K_b, "k_off" = k_off,
                "F_free" = F_free, "K_f" = K_f, "K_bilec" = K_bilec,
                "Q_cardiacc" = Q_cardiacc, "Q_lic" = Q_lic, "Q_kic" = Q_kic,
                "Q_filc" = Q_filc, "Q_gic" = Q_gic, "Q_luc" = Q_luc,
                "Q_htc" = Q_htc, "Q_brc" = Q_brc,"Q_adc" = Q_adc,
                "Q_bmc" = Q_bmc, "Q_rbc" = Q_rbc,
                "V_blc" = V_blc, "V_kic" = V_kic,
                "V_filc" = V_filc, "V_lic" = V_lic,"V_gic" = V_gic,
                "V_luc" = V_luc, "V_htc" = V_htc, "V_brc" = V_brc,
                "V_adc" = V_adc, "V_bmc" = V_bmc, "V_rbc" = V_rbc))
    
  })
}

#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters){
  with( as.list(parameters),{
    
    A_bl <- 0 #iv_dose * BW;
    A_glumen <- 0 #oral_dose_init * BW * (1 - F_unabs);
    A_gi <- 0.0;
    A_li <- 0.0;
    A_lib <- 0.0;
    A_ki <- 0.0;
    A_fil <- 0.0;
    A_lu <- 0.0;
    A_ht <- 0.0;
    A_br <- 0.0;
    A_ad <- 0.0;
    A_bm <- 0.0;
    A_rb <- 0.0;
    A_fst <- 0 #oral_dose_init * BW * F_unabs;
    A_ust <- 0.0;
    A_fecal <- 0.0;
    A_urine <- 0.0;
    iv_dose_cont <- 0.0;
    A_in <- 0 #oral_dose_init * BW + iv_dose * BW;
    
    
    return(c( "A_bl" = A_bl,"A_glumen" = A_glumen,"A_gi" = A_gi,
              "A_li" = A_li,"A_lib" = A_lib,"A_ki" = A_ki,"A_fil" = A_fil,
              "A_lu" = A_lu,"A_ht" = A_ht,"A_br" = A_br,
              "A_ad" = A_ad,"A_bm" = A_bm,"A_rb" = A_rb,
              "A_fst" = A_fst, "A_ust" = A_ust, "A_fecal" = A_fecal,
              "A_urine" = A_urine,"iv_dose_cont" = iv_dose_cont,"A_in" = A_in))
  })
}

#===================
#3. Events function
#===================

create.events <- function(parameters){
  with(as.list(parameters), {
    
    # Calculate number of administrated doses and corresponding administration time
    ldose <- length(admin.dose)
    ltimes <- length(admin.time)
    # If not equal, then stop 
    if (ltimes != ldose){
      stop("The times of administration should be equal in number to the doses")
    }else{
      if (admin.type == "iv"){
        vector_of_doses <- c()
        for (dose in admin.dose) {
          added_doses <- rep(dose, 2)
          vector_of_doses <- c(vector_of_doses, added_doses)
        }
        events <- list(data = rbind(data.frame(var = c("A_bl", "A_in"),  time = rep(admin.time, each = 2), 
                                               value = vector_of_doses, method = c("add")) ))
      }else if (admin.type == "oral"){
        vector_of_doses <- c()
        for (dose in admin.dose) {
          added_doses <- c(dose*(1 - F_unabs), dose*F_unabs, dose)
          vector_of_doses <- c(vector_of_doses, added_doses)
        }
        
        events <- list(data = rbind(data.frame(var = c("A_glumen", "A_fst", "A_in"),  time = rep(admin.time, each = 3), 
                                               value = vector_of_doses, method = c("add")) ))
      }
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
    
    # Body weight (kg)
    if(length(BW)>1){
      # Body weight (kg)
      BW_fnc <- approxfun(BW.times, BW, method = "linear", rule = 2)
      BW_out <- BW_fnc(time)
    }else{
      BW_out <- BW
    }
    
    
    # Compute rate constants that scale with BW
    # Resorption maximum (mg/h) adjusted by BW
    T_m = T_mc * BW_out ^ 0.75;
    # Biliary excretion rate (/h) adjusted by BW
    K_bile = K_bilec *  BW_out^ (-0.25);
    # Urinary elimination rate (/h) adjusted by BW
    K_u = K_uc *  BW_out ^ (-0.25);
    # Fecal elimination rate (/h) adjusted by BW
    K_fst = K_fstc *  BW_out ^ (-0.25);
    # Rate constant to urinary storage (/h) adjusted by BW
    K_ust = K_ustc * BW_out^ (-0.25);
    
    #Computed parameters based on BW
    # Cardiac output and blood flows to compartments (L/h)
    Q_cardiac = Q_cardiacc * BW_out^0.75;
    Q_li = Q_lic*Q_cardiac;
    Q_ki = Q_kic*Q_cardiac;
    Q_fil = Q_filc*Q_cardiac;
    Q_gi = Q_gic*Q_cardiac;
    Q_lu = Q_luc*Q_cardiac;
    Q_ht = Q_htc*Q_cardiac;
    Q_br = Q_brc*Q_cardiac;
    Q_ad = Q_adc*Q_cardiac;
    Q_bm = Q_bmc*Q_cardiac;
    Q_rb = Q_rbc*Q_cardiac;
    
    # Compartment volumes (L; assumed approx. equal to kg)
    V_bl = V_blc*BW_out;
    V_ki = V_kic*BW_out;
    V_fil = V_filc*BW_out;
    V_li = V_lic*BW_out;
    V_gi = V_gic*BW_out;
    V_lu = V_luc*BW_out;
    V_ht = V_htc*BW_out;
    V_br = V_brc*BW_out;
    V_ad = V_adc*BW_out;
    V_bm = V_bmc*BW_out;
    V_rb = V_rbc*BW_out;
    
    # Total amount (currently in system + metabolized) (mg)
    A_totl = A_glumen + A_bl + A_gi + A_fecal + A_li + A_ki + A_fil + A_ust +
      A_fst + A_urine + A_lu + A_ht + A_br + A_ad + A_bm + A_rb + A_lib;
    # Concentrations (mg/L)
    C_bl = A_bl / V_bl;
    C_gi = A_gi / V_gi;
    C_li = A_li / V_li;
    C_ki = A_ki / V_ki;
    C_fil = A_fil / V_fil;
    C_lu = A_lu / V_lu;
    C_ht = A_ht / V_ht;
    C_br = A_br / V_br;
    C_ad = A_ad / V_ad;
    C_bm = A_bm / V_bm;
    C_rb = A_rb / V_rb;
    
    #=========================================
    # Time rates of change of amounts (ODEs)
    #========================================
    
    # Amount in blood (plasma) (mg)
    dA_bl <- F_free*((Q_li+Q_gi)*C_li/K_li + Q_ki*C_ki/K_ki + 
                       Q_lu*C_lu/K_lu + Q_ht*C_ht/K_ht + Q_br*C_br/K_br+ 
                       Q_ad*C_ad/K_ad + Q_bm*C_bm/K_bm + Q_rb*C_rb/K_rb) -
      F_free*C_bl*(Q_li + Q_ki + Q_gi + Q_fil + Q_lu + Q_ht + Q_br +
                     Q_ad + Q_bm + Q_rb) + iv_dose_cont / 24 * BW_out; 
    
    # Amount gut lumen (mg)
    dA_glumen <- -K_abs*A_glumen - K_unabs*A_glumen;
    
    # Amount gut tissue (mg)
    dA_gi <- F_free*(Q_gi*C_bl - Q_gi*C_gi/K_gi) - K_f*A_gi +
      (1-gut2liver)*K_abs*A_glumen;
    
    # Amount liver (mg)
    dA_li <- F_free*(Q_li*C_bl + Q_gi*C_gi/K_gi - (Q_li+Q_gi)*C_li/K_li) -
      K_bile*A_li + gut2liver*K_abs*A_glumen -
      (Bmax*C_li*F_free/K_li)/(K_b + C_li*F_free/K_li) + k_off*A_lib;
    
    # Amount bound in liver (mg)
    dA_lib <- (Bmax*C_li*F_free/K_li)/(K_b + C_li*F_free/K_li) - k_off*A_lib;
    
    # Amount kidney (mg)
    dA_ki <- F_free*(Q_ki*C_bl - Q_ki*C_ki/K_ki) + T_m*C_fil/(K_t + C_fil);
    
    # Amount filtrate (mg)
    dA_fil <- F_free*Q_fil*C_bl - T_m*C_fil/(K_t + C_fil) - K_ust*V_fil*C_fil;
    
    # Amount lungs (mg)
    dA_lu <- F_free*(Q_lu*C_bl - Q_lu*C_lu/K_lu);
    
    # Amount heart (mg)
    dA_ht <- F_free*(Q_ht*C_bl - Q_ht*C_ht/K_ht);
    
    # Amount brain (mg)
    dA_br <- F_free*(Q_br*C_bl - Q_br*C_br/K_br);
    
    # Amount adipose (mg)
    dA_ad <- F_free*(Q_ad*C_bl - Q_ad*C_ad/K_ad);
    
    # Amount bone marrow (mg)
    dA_bm <- F_free*(Q_bm*C_bl - Q_bm*C_bm/K_bm);
    
    dA_rb <- F_free*(Q_rb*C_bl - Q_rb*C_rb/K_rb);
    
    # Amount fecal storage (mg)
    dA_fst <- K_bile*A_li + K_unabs*A_glumen + K_f*A_gi - K_fst*A_fst;
    
    # Amount urinary storage (mg)
    dA_ust <- K_ust*C_fil*V_fil - K_u*A_ust;
    
    # Amount excreted via fecal elimination (mg)
    dA_fecal <- K_fst*A_fst;
    
    # Amount excreted via urinary elimination (mg)
    dA_urine <- K_u*A_ust;
    
    # Time rate of change of continuous intravenous dose.
    div_dose_cont <- 0;
    
    # Time rate of change of cumulative amount that has entered organism.
    dA_in <- iv_dose_cont / 24 * BW_out;
    
    # Mass balance.
    A_bal = A_in - A_totl;
    
    list(c("dA_bl" = dA_bl, "dA_glumen" = dA_glumen, "dA_gi" = dA_gi,
           "dA_li" = dA_li, "dA_lib" = dA_lib, "dA_ki" = dA_ki,"dA_fil" = dA_fil,
           "dA_lu" = dA_lu, "dA_ht" = dA_ht, "dA_br" = dA_br,
           "dA_ad" = dA_ad, "dA_bm" = dA_bm, "dA_rb" = dA_rb,
           "dA_fst" = dA_fst, "dA_ust" = dA_ust, "dA_fecal" = dA_fecal,
           "dA_urine" = dA_urine, "div_dose_cont" = div_dose_cont, "dA_in" = dA_in), 
         "A_totl" = A_totl,  
         "A_bal" = A_bal, "C_li" = C_li, "C_gi" = C_gi, "C_ki" = C_ki,
         "C_fil" = C_fil, "C_lu" = C_lu, "C_ht" = C_ht, "C_br" = C_br,
         "C_ad" = C_ad, "C_bm" = C_bm, "C_rb" = C_rb, "C_bl" = C_bl,
         "BW_out" = BW_out, "F_free" = F_free, "Q_cardiac" = Q_cardiac, "Q_li" = Q_li,
         "Q_ki" = Q_ki, "Q_gi" = Q_gi, "Q_lu" = Q_lu, "Q_ht" = Q_ht,
         "Q_br" = Q_br, "Q_ad" = Q_ad, "Q_bm" = Q_bm, "Q_rb" = Q_rb,
         "Q_fil" = Q_fil, "T_m" = T_m, "K_f" = K_f, "K_bile" = K_bile,
         "K_u" = K_u, "K_fst" = K_fst, "K_ust" = K_ust, "V_bl" = V_bl,
         "V_ki" = V_ki, "V_fil" = V_fil, "V_li" = V_li, "V_gi" = V_gi,
         "V_lu" = V_lu, "V_ht" = V_ht, "V_br" = V_br, "V_ad" = V_ad,
         "V_bm" = V_bm, "V_rb" = V_rb)
    
  })
}

#=============
#6. User input 
#=============
# Parameters for reproducing example of sheet "FKimRecreateBW" in "PFAS_template_parameters_PFDA.xlsx" of Bernstein et al.2021
sex <- "F" # rat sex, values: M/F
BW <- 0.25# rat body weight in kg
BW.times <-0 # Times corresponding to BW vector. If BW is constant, then type 0
admin.type <-  "oral" # administration type values: iv/oral
admin.dose <- 1*BW  # administered dose in mg
admin.time <- 0 # time when doses are administered, in hours
F_unabs <-   0.650 # Fraction of unabsorbed dose


user_input <- list( "admin.type" = admin.type,
                    "admin.dose" = admin.dose, 
                    "admin.time" = admin.time,
                    "BW"=BW, "BW.times" = BW.times,
                    "F_unabs" = F_unabs, "sex" = sex)


params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,150*24,0.1)  #hours

solution <-  data.frame(ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                 events = events, method="bdf",rtol = 1e-05, atol = 1e-05))
print(tail(solution))

# Male, multiple oral doses
sex <- "M" # rat sex, values: M/F
BW <- 0.25# rat body weight in kg
BW.times <-0 # Times corresponding to BW vector. If BW is constant, then type 0
admin.type <-  "oral" # administration type values: iv/oral
admin.dose <- c(1 * BW, 1 * BW, 10 * BW)  # administered dose in mg
admin.time <- c(0,5,10) # time when doses are administered, in hours
F_unabs <-   0.650 # Fraction of unabsorbed dose

user_input <- list( "admin.type" = admin.type,
                    "admin.dose" = admin.dose, 
                    "admin.time" = admin.time,
                    "BW"=BW, "BW.times" = BW.times,
                    "F_unabs" = F_unabs, "sex" = sex)

params <- create.params(user_input)
inits <- create.inits(params)
events <- create.events(params)
sample_time <- seq(0,15,0.1)  #hours

solution <-  ode(times = sample_time,  func = ode.func, y = inits, parms = params,
                 events = events, method="bdf",rtol = 1e-05, atol = 1e-05)
plot(solution[,1],solution[,3], type = "l")

#====================
#7. Upload on Jaqpot 
#===================
# Subset of features to be displayed on the user interface
predicted.feats <- c("A_li", "A_gi", "A_ki", "A_fil", "A_rb", "A_bl", "A_lu", "A_ht",
                     "A_fecal", "A_urine",  "A_fst",  "A_ust", "A_glumen" = "A_glumen",
                     "C_li", "C_gi", "C_ki", "C_fil", "C_rb", "C_bl","C_lu", "C_ht",
                     "BW_out")

# Deploy the model on the Jaqpot server to create a web service
jaqpotr::deploy.pbpk(user.input = user_input,out.vars = predicted.feats,
                     create.params = create.params,  create.inits = create.inits,
                     create.events = create.events, custom.func = custom.func,
                     ode.fun = ode.fun, envFile = "")