  library(deSolve)
  library(tidyverse)
  # library(jaqpotr)  # only needed for Jaqpot upload, not by the Shiny app
  #=========================
  # 1. Parameters of the model
  #=========================

#'  
#' @export
#'  
create.params <- function(user_input) {
  with( as.list(user_input),{

    # Replaced %>% select(all_of(chemical)) with base R subsetting [, chemical, drop = FALSE]
    # Contents of estimated_parameters_PopGen.csv
    variables <- data.frame(
      PFHpA   = c(364.06, 0.0045, 3300, 263, 4500, 60, 1136.21342943039, 0.000142164477918033, 0.00631107236193145, 3.8, 2.78584860940939, 1.512, 10),
      PFOA    = c(414.07, 0.00245, 2500, 137.5, 4500, 47, 1, 0.00045, 0.00235, 3.96, 10, 1.512, 50.88),
      PFNA    = c(464.08, 0.00187, 1000, 42, 8500, 58, 0.4045538, 7e-04, 0.0021553774486375, 4.02, 3180.077, 1.512, 50.88),
      PFDA    = c(514.09, 0.000618, 2500, 137.5, 6000, 39, 7.78143006004393, 6.74263619631535e-05, 0.00261140549164202, 4.01, 15486.9201267857, 1.512, 50.88),
      PFBS    = c(300.1, 0.0365, 2500, 137.5, 4500, 47, 12.3656911519399, 0.000386421164690321, 0, 18.42, 190.6249146875, 1.512, 50.88),
      PFHxS   = c(403.11, 0.000695, 5800, 778, 7300, 92, 7.34410734176636, 1e-04, 0, 3.76, 0.630791833381653, 1.512, 50.88),
      PFOS    = c(500.13, 0.0045, 2500, 137.5, 2200, 48, 1, 0.0014, 0, 3.78, 225, 1.512, 50.88),
      DONA    = c(378.07, 0.0193, 2500, 137.5, 4500, 39, 0.243517057527131, 6.50582714720544e-05, 0.0257317294556275, 3.71, 35.9865797071135, 1.512, 50.88),
      HFPO_DA = c(330.0489, 0.0491, 2500, 137.5, 4500, 47, 1305.45789452388, 8.10939128423146e-06, 0.0873472097851557, 19.58, 43.8342106997064, 1.512, 50.88),
      PFBA    = c(214.04, 0.229, 2500, 137.5, 4500, 47, 1, 0.0001178515625, 0, 4.01, 200, 1.512, 50.88),
      PFHxA   = c(314.05, 0.038, 2500, 137.5, 4500, 47, 0.029798235492761, 2.10509023692685e-06, 0.097821027749457, 3.99, 71.6921170325873, 1.512, 50.88)
    )[, chemical, drop = FALSE]

    # Set the row names to match the 'Parameters' column
    rownames(variables) <- c("MW", "Free", "Vmax_baso_invitro", "Km_baso", "Vmax_apical_invitro", "Km_apical", "RAFbaso", "RAFapi", "kbilec", "admin_dose_bolus", "keffluxc", "kurinec", "kabsc")

    # Contents of shiny_app_population/PCs_human_replaced.csv
    PC <- data.frame(
      DONA    = c(0.298730574951549, 0.292474891124872, 0.111622255587086, 0.116569258573112, 0.106399963739007, 0.141692100982492, 0.1941693387432, 0.113133210469352, 0.380709559520465, 0.142250732230917, 1, 1.64143748626212, 0.0850006573309938, 0.153685559151805),
      HFPO_DA = c(0.471005008459102, 0.514779965369464, 0.184642687522042, 0.197206093001642, 0.174190691566295, 0.547360264077368, 0.250986777510723, 0.186244768659733, 0.604424013633604, 0.220753777215081, 1, 1.59131067236414, 0.119002839282567, 0.213331314975323),
      PFBA    = c(0.396679756602288, 0.469464282259187, 0.205080014985862, 0.211181085376696, 0.183480477763335, 0.222643167101742, 0.263486323098763, 0.195617533859863, 0.512701355956668, 0.224412267845041, 1, 1.61281348758351, 0.0944081001506644, 0.235717159939647),
      PFBS    = c(0.4302271, 0.959157701947092, 0.254968312875282, 0.282684418139592, 0.217717605023404, 0.29795547911487, 0.273045463634937, 0.240324652072097, 0.989215861903838, 0.281968134883447, 1, 1.60336619815795, 0.11414314946686, 0.271123743346798),
      PFDA    = c(0.23, 0.5, 0.520427078733946, 0.36, 0.441815944289152, 3.327, 0.434435367210863, 0.23, 1.1, 0.564202015421272, 1, 1.47554616251317, 0.217429049819238, 0.483429060391693),
      PFHpA   = c(0.202879340315232, 0.180877682681212, 0.0896933317854623, 0.0907720077120878, 0.0900993291907227, 0.120228982178295, 0.183347431492427, 0.0942338124843852, 0.280099126014041, 0.12131561210594, 1, 1.64513756657515, 0.0815440034781734, 0.135653273194568),
      PFHxA   = c(0.206107554198993, 0.191377826099034, 0.09774633994465, 0.0986375340050109, 0.0967067829815432, 0.134201744794219, 0.189558379342699, 0.101239374916713, 0.286229401054081, 0.12810636587219, 1, 1.64291039800941, 0.0818348883352994, 0.142818229163531),
      PFHxS   = c(0.455194155149325, 0.231078431372549, 0.138724632000269, 0.149657895794101, 0.1221470381763, 0.4837, 1.09803921568627, 0.133159668041511, 0.533026074061412, 0.588235294117647, 1, 1.64899062138046, 0.0809898715850802, 0.176017483447799),
      PFNA    = c(0.2380952, 0.199821428571429, 0.385375096691658, 0.439168238146343, 0.319629662224778, 2.42033, 1.52678571428571, 0.1790368, 1.64211491882327, 0.732142857142857, 1, 1.56374610578011, 0.150157789934354, 0.375860643763351),
      PFOA    = c(0.23333333, 0.115099349092155, 0.633333333333333, 0.158278588395683, 0.133811934378827, 0.56, 0.79, 0.144735437423378, 0.547613753046703, 0.56, 1, 1.63115429770854, 0.0936242052382557, 0.433333333333333),
      PFOS    = c(0.08, 0.209916607288945, 0.666666666666667, 0.719197187017213, 0.518657820130813, 2.89171262002202, 1.79436175546935, 0.196078431372549, 1.2, 0.936742934051144, 1, 1.45858823589366, 0.233222746137853, 0.686274509803922),
      PFUnDA  = c(1.20160041562323, 1.39093833943251, 0.477928551939649, 0.534574341412838, 0.503998549654399, 0.68327594360335, 0.533718816720171, 0.533448959646928, 1.70088000194443, 0.607115817727957, 1, 1.20759092980041, 0.404841651445938, 0.437281917824021)
    )[, chemical, drop = FALSE]

    # Set Organs as row names
    rownames(PC) <- c("Adipose", "Brain", "Gonads", "Gut", "Heart", "Liver", "Lung", "Muscle", "Skin", "Spleen", "Blood", "Plasma", "Bone", "Pancreas")


    if  (missing(time_scale) || is.na(time_scale)){
      time_scale <- 1/24
    }else if (time_scale == "minutes"){
      time_scale <- 60
    }else if (time_scale == "hours"){
      time_scale <- 1
    }else if (time_scale == "days"){
      time_scale <- 1/24
    }else if (time_scale == "weeks"){
      time_scale <- (1/24)/7
    }else if (time_scale == "months"){
      time_scale <- ((1/24)/30)
    }else if (time_scale == "years"){
      time_scale <- ((1/24)/365)
    }
    inv_time_scale <- 1/time_scale

    # --- Cardiac Output and Blood Flow (as fraction of cardiac output) ---
    QCC <- 12.5 * inv_time_scale  # cardiac output in L/day/kg^0.75; Brown 1997
    QLC <- 0.065  # fraction blood flow to liver; Brown 1997
    QKC <- 0.175  # fraction blood flow to kidney; Brown 1997
    QAdiC <-0.05  # fraction blood flow to adipose tissue; Brown 1997
    QBraC <-0.12  # fraction blood flow to brain; Brown 1997
    QGonC <-0.0005  # fraction blood flow to gonads; ICRP 2002
    QHeaC <-0.04  # fraction blood flow to heart; Brown 1997
    QLunC <-0.025  # fraction blood flow to lung; Brown 1997
    QMusC <-0.17  # fraction blood flow to muscle; Brown 1997
    QSkiC <-0.05  # fraction blood flow to skin; Brown 1997
    QSplC <-0.03  # fraction blood flow to spleen; ICRP 2002
    QPanC <-0.01  # fraction blood flow to pancreas; ICRP 2002
    QGIC<-0.15  # fraction blood flow to GI tract; ICRP 2002

    Htc <- 0.42  # hematocrit; 

    # --- Tissue Volumes ---
    VplasC <- 0.0428  # fraction vol. of plasma (L/kg BW); Davies 1993
    VLC <- 0.026  # fraction vol. of liver (L/kg BW); Brown 1997
    VKC <- 0.004  # fraction vol. of kidney (L/kg BW); Brown 1997
    VAdiC <- 0.214  # fraction vol. of adipose tissue (L/kg BW); Brown 1997
    VBraC <- 0.02  # fraction vol. of brain (L/kg BW); Brown 1997
    VGonC <- 0.0005 # fraction vol. of gonads (L/kg BW); 
    VHeaC <- 0.005  # fraction vol. of heart (L/kg BW); Brown 1997
    VLunC <- 0.008  # fraction vol. of lung (L/kg BW); Brown 1997
    VMusC <- 0.4  # fraction vol. of muscle (L/kg BW); Brown 1997
    VSkiC <- 0.037  # fraction vol. of skin (L/kg BW); Brown 1997
    VSplC <- 0.002  # fraction vol. of spleen (L/kg BW); ICRP 2002
    VPanC <- 0.002 # fraction vol. of pancreas (L/kg BW); ICRP 2002
    VGIC <- 0.014  # fraction vol. of GI tract (L/kg BW); Brown 1997

    VfilC <- 4e-4  # fraction vol. of filtrate (L/kg BW)
    VPTCC <- 1.35e-4  # vol. of proximal tubule cells (L/g kidney)

    # --- Chemical Specific Parameters ---
    MW <- variables["MW",]  # PFAS molecular mass (g/mol)
    Free <- variables["Free",]  # free fraction in plasma (Smeltz 2023); fitted to model

    # --- Kidney Transport Parameters ---
    Vmax_baso_invitro <- variables["Vmax_baso_invitro", ]  # Vmax of basolateral transporter (pmol/mg protein/min)
    Km_baso <- variables["Km_baso", ] * variables["MW", ]  # Km of basolateral transporter (ug/L)
    Vmax_apical_invitro <- variables["Vmax_apical_invitro", ]  # Vmax of apical transporter (pmol/mg protein/min)
    Km_apical <- variables["Km_apical", ] * variables["MW", ]  # Km of apical transporter (ug/L)
    RAFbaso <- variables["RAFbaso", ]  # relative activity factor, basolateral transporters (male)
    RAFapi <- variables["RAFapi", ] # relative activity factor, apical transporters (male); fitted to model
    protein <- 2.0e-6  # amount of protein in proximal tubule cells (mg protein/cell)
    GFRC <- 24.19 * inv_time_scale  # glomerular filtration rate (L/day/kg kidney); Corley 2005

    # --- Partition Coefficients (from Allendorf 2021) ---
    PR   <- (1-Htc) * PC["Bone",] # rest of body:blood#; fitted to model
    PAdi <- (1-Htc) * PC["Adipose",] # adipose tissue:blood;
    PBra <- (1-Htc) * PC["Brain",]  # brain:blood;
    PGon <- (1-Htc) * PC["Gonads",]  # gonads:blood;
    PGI  <- (1-Htc) * PC["Gut",]# GI tract:blood;
    PHea <- (1-Htc) * PC["Heart",]  # heart:blood;
    PL   <- (1-Htc) * PC["Liver",]  # liver:blood;
    PLun <- (1-Htc) * PC["Lung",]  # lung:blood;
    PMus <- (1-Htc) * PC["Muscle",] # muscle:blood;
    PSki <- (1-Htc) * PC["Skin",]  # skin:blood;
    PSpl <- (1-Htc) * PC["Spleen",]  # spleen:blood;
    PPan <- (1-Htc) * PC["Pancreas",]  # pancreas:blood; estimated as average of GI tract and spleen

    # --- Rate Constants ---
    kdif <- 0.001 * inv_time_scale  # diffusion rate from proximal tubule cells (L/day)
    kabsc <-  2.12 * inv_time_scale  # rate of absorption from small intestine (1/(day*BW^-0.25))
    kunabsc <- 7.06e-5 * inv_time_scale  # rate of unabsorbed dose to feces (1/(day*BW^-0.25)); fitted to model 
    # estimated_parameters.csv stores keffluxc and kbilec per day, hence /24 (as in shiny_app_population)
    keffluxc <-variables["keffluxc", ] * inv_time_scale / 24  # rate of efflux from PTC to blood (1/(day*BW^-0.25))
    kbilec <- variables["kbilec", ] * inv_time_scale / 24 # biliary elimination rate (1/(day*BW^-0.25)); fitted to model
    kurinec <- 0.063 * inv_time_scale  # urinary elimination rate (1/(day*BW^-0.25))

    # --- Water Consumption ---
    water_consumption <- 1.36  # L/day

    # === SCALED PARAMETERS (calculated from above) ===

    # Cardiac output and blood flows
    QC <- QCC * (BW^0.75) * (1 - Htc)  # cardiac output in L/day; adjusted for plasma
    QK <- (QKC * QC)  # plasma flow to kidney (L/day)
    QL <- (QLC * QC)  # plasma flow to liver (L/day)
    QAdi <- (QAdiC * QC)  # plasma flow to adipose tissue (L/day)
    QBra <- (QBraC * QC)  # plasma flow to brain (L/day)
    QGon <- (QGonC * QC)  # plasma flow to gonads (L/day)
    QHea <- (QHeaC * QC)  # plasma flow to heart (L/day)
    QLun <- (QLunC * QC)  # plasma flow to lung (L/day)
    QMus <- (QMusC * QC)  # plasma flow to muscle (L/day)
    QSki <- (QSkiC * QC)  # plasma flow to skin (L/day)
    QSpl <- (QSplC * QC)  # plasma flow to spleen (L/day)
    QPan <- (QPanC * QC)  # plasma flow to pancreas (L/day)
    QGI <- (QGIC * QC)  # plasma flow to GI tract (L/day)
    QR <- QC - QK - QL - QAdi - QBra - QGon - QHea - QLun - QMus - QSki - QSpl - QPan - QGI  # plasma flow to rest of body (L/day)

    QBal <- QC - (QK + QL + QR + QAdi + QBra + QGon + QHea + QLun + QMus + QSki + QSpl + QPan + QGI)  # Balance check; should equal zero

    # Tissue Volumes
    VPlas <- VplasC * BW  # volume of plasma (L)
    VK <- VKC * BW  # volume of kidney (L)
    MK <- VK * 1.0 * 1000  # mass of the kidney (g)
    VKb <- VK * 0.16  # volume of blood in the kidney (L); Brown 1997
    Vfil <- VfilC * BW  # volume of filtrate (L)
    VL <- VLC * BW  # volume of liver (L)
    ML <- VL * 1.05 * 1000  # mass of the liver (g)
    VAdi <-VAdiC * BW  # volume of adipose tissue (L)
    VBra <- VBraC * BW  # volume of brain (L)
    VGon <- VGonC * BW  # volume of gonads (L)
    VHea <- VHeaC * BW  # volume of heart (L) 
    VLun <- VLunC * BW  # volume of lung (L)
    VMus <- VMusC * BW  # volume of muscle (L)
    VSki <- VSkiC * BW  # volume of skin (L)
    VSpl <- VSplC * BW  # volume of spleen (L)
    VPan <- VPanC * BW  # volume of pancreas (L)
    VGI <- VGIC * BW  # volume of GI tract (L)


    # Kidney Parameters
    PTC <- VKC * 1000 * 6e7  # number of PTC (cells/kg BW)
    VPTC <- VK * 1000 * VPTCC  # volume of proximal tubule cells (L)
    MPTC <- VPTC * 1000  # mass of the proximal tubule cells (g)
    VR <- (0.93 * BW) - VPlas- VPTC - Vfil - VL -VAdi - VBra - VGon - VHea - VLun - VMus - VSki - VSpl - VPan - VGI  # volume of rest of body (L)
    VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas + VAdi + VBra + VGon + VHea + VLun + VMus + VSki + VSpl + VPan + VGI)  # Balance check; should equal zero

    Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6)  * inv_time_scale 
    Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * inv_time_scale 
    Vmax_baso <- Vmax_basoC * BW^0.75  # (ug/day)
    Vmax_apical <- Vmax_apicalC * BW^0.75  # (ug/day)
    kbile <- kbilec * BW^(-0.25)  # biliary elimination; liver to feces storage (/day)
    kurine <- kurinec * BW^(-0.25)  # urinary elimination, from filtrate (/day)
    kefflux <- keffluxc * BW^(-0.25)  # efflux clearance rate, from PTC to blood (/day)
    GFR <- GFRC*VK # glomerular filtration rate, (L/day)Abraham's personal value

    # GI Tract Parameters
    kabs <- kabsc * BW^(-0.25)  # rate of absorption from small intestine (/day)
    kunabs <- kunabsc * BW^(-0.25)  # rate of unabsorbed dose to feces (/day)
    

    return(list(
      "BW"=BW, "QBAL" = QBal, "VBAL" = VBal,
      "Free" = Free, "QC" = QC, "QK" = QK, "QL" = QL, "QR" = QR,
      "QAdi" = QAdi, "QBra" = QBra, "QGon" = QGon, "QHea" = QHea,
      "QLun" = QLun, "QMus" = QMus, "QSki" = QSki, "QSpl" = QSpl,
      "QPan" = QPan, "QGI" = QGI,
      "VPlas" = VPlas, "VKb" = VKb, "Vfil" = Vfil, "VL" = VL,
      "VR" = VR, "ML" = ML,
      "VAdi" = VAdi, "VBra" = VBra, "VGon" = VGon, "VHea" = VHea,
      "VLun" = VLun, "VMus" = VMus, "VSki" = VSki, "VSpl" = VSpl, "VPan" = VPan,
      "VGI" = VGI, "MK" = MK,
      "VPTC" = VPTC, "Vmax_baso" = Vmax_baso, "Vmax_apical" = Vmax_apical,
      "kdif" = kdif, "Km_baso" = Km_baso, "Km_apical" = Km_apical,
      "kbile" = kbile, "kurine" = kurine, "kefflux" = kefflux,
      "GFR" = GFR, "kabs" = kabs, "kunabs" = kunabs, 
      "PR" = PR, "PAdi" = PAdi, "PBra" = PBra, "PGon" = PGon, "PGI" = PGI,
      "PHea" = PHea, "PL" = PL, "PLun" = PLun, "PMus" = PMus, "PSki" = PSki,
      "PSpl" = PSpl, "PPan" = PPan,
      "ingestion" = ingestion,
      "ingestion_time" = ingestion_time ,
      "time_scale" = time_scale, "exp_type" = exp_type))

  })
}
#===============================================
#2. Function to create initial values for ODEs
#===============================================
#'  
#' @export
#'  
create.inits <- function(parameters){
  with( as.list(parameters),{
    "AR" = 0; "AAdi"=0; "ABra"=0; "AGon"=0;
    "AHea"=0; "ALun"=0; "AMus"=0; "ASki"=0; "ASpl"=0; "APan"=0;
    "Adif" = 0; "A_baso" = 0; "AKb" = 0;
    "ACl" = 0; "Aefflux" = 0;
    "A_apical" = 0; "APTC" = 0; "Afil" = 0;
    "Aurine" = 0; "ALumen" = 0; "AGI" = 0;
    "AabsLumen" = 0; "Afeces" = 0;
    "AL" = 0; "Abile" = 0; "Aplas_free" = 0;
    "ingestion" = 0; 

    return(c("AR" = AR, "AAdi"=AAdi, "ABra"=ABra, "AGon"=AGon,
             "AHea"=AHea, "ALun"=ALun, "AMus"=AMus, "ASki"=ASki,
             "ASpl"=ASpl, "APan"=APan,
             "Adif" = Adif, "A_baso" = A_baso, "AKb" = AKb,
             "ACl" = ACl, "Aefflux" = Aefflux,
             "A_apical" = A_apical, "APTC" = APTC, "Afil" = Afil,
             "Aurine" = Aurine, "ALumen" = ALumen, "AGI" = AGI,
             "AabsLumen" = AabsLumen, "Afeces" = Afeces, 
             "AL" = AL, "Abile" = Abile, "Aplas_free" = Aplas_free,
             "ingestion" = ingestion))
  })
}

#===================
# 3. Events function
#===================
#'  
#' @export
#'  
create.events <- function(parameters){
  with(as.list(parameters), {
 
      lingest <- length(ingestion)
      lingesttimes <- length(ingestion_time)
      if (lingest != lingesttimes) {
        #stop("The times of ingestion rate change should be equal in vector of ingestion")
      }
      if(exp_type == "pharmacokinetics"){
        # For pharmacokinetic studies, add the dose directly to the stomach compartment
        events <- list(data = data.frame(var = c("ALumen"),  time = ingestion_time, 
                                         value = ingestion, method = c("add")))
      }else if(exp_type == "continuous"){
                 # For continuous exposure scenarios, set the ingestion rate in the model  
        events <- list(data = data.frame(var = c("ingestion"),  time = ingestion_time, 
                                         value = ingestion, method = c("rep")))
      }else{
        stop("admin_type should be either 'bolus', or 'oral'")
      }
      return(events)

    }
  )
}

#==================
# 4. Custom functions
#==================
#'  
#' @export
#'  
custom.func <- function(){
  return()
}

#==============
# 5. ODEs System
#==============

#'  
#' @export
#'  
 
  ode.func <- function(time, inits, params, custom.func) {
    with(as.list(c(inits, params)), {

      # Concentrations in various compartments
      
      CR <- AR / VR  # concentration in rest of body (ug/L)
      CVR <- CR / PR  # concentration in venous blood leaving rest of body (ug/L)

      CAdi <- AAdi / VAdi  # concentration in adipose tissue (ug/L)
      CVAdi <- CAdi / PAdi  # concentration in venous blood leaving adipose
      
      CBra <- ABra / VBra  # concentration in brain tissue (ug/L)
      CVBra <- CBra / PBra  # concentration in venous blood leaving brain
      
      CGon <- AGon / VGon  # concentration in gonads tissue (ug/L)
      CVGon <- CGon / PGon  # concentration in venous blood leaving gonads
     
      
      CHea <- AHea / VHea  # concentration in heart tissue (ug/L)
      CVHea <- CHea / PHea  # concentration in venous blood
    
      
      CLun <- ALun / VLun  # concentration in lung tissue (ug/L)
      CVLun <- CLun / PLun  # concentration in venous blood leaving lung
      
      
      CMus <- AMus / VMus  # concentration in muscle tissue (ug/L)
      CVMus <- CMus / PMus  # concentration in venous blood leaving muscle
      

      CSki <- ASki / VSki  # concentration in skin tissue (ug/L)
      CVSki <- CSki / PSki  # concentration in venous blood leaving skin
      
      CSpl <- ASpl / VSpl  # concentration in spleen tissue (ug/L)
      CVSpl <- CSpl / PSpl  # concentration in venous blood leaving spleen
     

      Cpan <- APan / VPan  # concentration in pancreas tissue (ug/L)
      CVPan <- Cpan / PPan  # concentration in venous blood leaving pancreas
      
      CGI <- AGI / VGI  # concentration in GI tract (ug/L)
      CVGI <- CGI / PGI  # concentration in venous blood leaving GI tract
     
      
      CKb <- AKb / VKb  # concentration in kidney blood (ug/L)
      CVK <- CKb  # concentration in venous blood leaving kidney (ug/L)
      CPTC <- APTC / VPTC  # concentration in PTC (ug/L)
      Cfil <- Afil / Vfil  # concentration in filtrate (ug/L)
      
      CL <- AL / VL  # concentration in the liver (ug/L)
      CLiver <- AL / ML  # concentration in the liver (ug/g)
      CVL <- CL / PL  # concentration in the venous blood leaving the liver (ug/L)
      
      
      CA_free <- Aplas_free / VPlas  # free concentration in plasma (ug/L)
      CA <- CA_free / Free  # concentration of total PFAS in plasma (ug/L)
      
      # Rest of Body (Tis)
      dAR <- QR * (CA - CVR) * Free

      # Adipose Tissue (Adi)
      dAAdi <- QAdi * (CA - CVAdi) * Free

      # Brain Tissue (Bra)
      dABra <- QBra * (CA - CVBra) * Free

      # Gonads Tissue (Gon)
      dAGon <- QGon * (CA - CVGon) * Free

      # Heart Tissue (Hea)
      dAHea <- QHea * (CA - CVHea) * Free

      # Lung Tissue (Lun)
      dALun <- QLun * (CA - CVLun) * Free
      
      # Muscle Tissue (Mus)
      dAMus <- QMus * (CA - CVMus) * Free

      # Skin Tissue (Ski)
      dASki <- QSki * (CA - CVSki) * Free

      # Spleen Tissue (Spl)
      dASpl <- QSpl * (CA - CVSpl) * Free

      # Pancreas Tissue (Pan)
      dAPan <- QPan * (CA - CVPan) * Free
      
      # Kidney
      # Kidney Blood (Kb)
      dAdif <- kdif * (CKb - CPTC)
      dA_baso <- (Vmax_baso * CKb) / (Km_baso + CKb)
      dAKb <- QK * (CA - CVK) * Free - CA * GFR * Free - dAdif - dA_baso
      dACl <- CA * GFR * Free

      # Proximal Tubule Cells (PTC)
      dAefflux <- kefflux * APTC
      dA_apical <- (Vmax_apical * Cfil) / (Km_apical + Cfil)
      dAPTC <- dAdif + dA_apical + dA_baso - dAefflux

      # Filtrate (Fil)
      dAfil <- CA * GFR * Free - dA_apical - Afil * kurine

      # Urinary elimination
      dAurine <- kurine * Afil

      # GI Tract (Absorption site of oral dose)
      # Stomach_lumen
      dALumen <- ingestion + - kabs * ALumen- kunabs * ALumen 


      # GI Tract - Stomach,Small Intestine,Large Intestine 
      dAGI <- kabs * ALumen + QGI*(CA - CVGI) * Free  
      dAabsLumen <- kabs * ALumen


      # Feces compartment
      dAfeces <-  kunabs * ALumen + kbile * AL 

      # Liver
      dAL <- QL * (CA - CVL) * Free - kbile * AL + QGI*(CVGI - CVL) * Free +
            QSpl*(CVSpl-CVL)*Free +QPan*(CVPan-CVL)*Free
      
      dAbile <- kbile * AL
      amount_per_gram_liver <- CLiver

      # Plasma compartment
      dAplas_free <- (QR * CVR * Free) + (QK * CVK * Free) + (QL * CVL * Free) +
      (QAdi*CVAdi*Free) + (QBra*CVBra*Free) + (QGon*CVGon*Free) + (QHea*CVHea*Free) +
      (QLun*CVLun*Free) +(QMus*CVMus*Free) +(QSki*CVSki*Free)+(QSpl*CVL*Free)+(QPan*CVL*Free)+
      (QGI*CVL*Free)- (QC * CA * Free) + dAefflux

      
      dingestion <- 0

      # Mass Balance Check
      Atissue <- Aplas_free + AR + AAdi + ABra + AGon + AHea+ ALun +
      AMus + ASki + ASpl+ APan + AKb + Afil + ALumen+ AGI + APTC + AL 
      Aloss <- Aurine + Afeces
      Atotal <- Atissue + Aloss

      list(c(
          "dAR" = dAR, "dAAdi" = dAAdi, "dABra" = dABra, "dAGon" = dAGon,
          "dAHea" = dAHea, "dALun" = dALun, "dAMus" = dAMus, "dASki" = dASki,
          "dASpl" = dASpl, "dAPan"=dAPan,
          "dAdif" = dAdif, "dA_baso" = dA_baso, "dAKb" = dAKb,
          "dACl" = dACl, "dAefflux" = dAefflux,
          "dA_apical" = dA_apical, "dAPTC" = dAPTC, "dAfil" = dAfil,
          "dAurine" = dAurine, "dALumen" = dALumen, "dAGI" = dAGI,
          "dAabsLumen" = dAabsLumen, "dAfeces" = dAfeces,
          "dAL" = dAL, "dAbile" = dAbile, "dAplas_free" = dAplas_free,
          "dingestion" = dingestion
        ),
        "amount_per_gram_liver" = amount_per_gram_liver,
        "Atissue" = Atissue, "Aloss" = Aloss, "Atotal" = Atotal,
        "CR" = CR,  "CVR" = CVR, 
        "CAdi"=CAdi,  "CVAdi"=CVAdi,
        "CBra"=CBra,  "CVBra"=CVBra,
        "CGon"=CGon,  "CVGon"=CVGon,
        "CHea"=CHea, "CVHea"=CVHea, 
        "CLun"=CLun, "CVLun"=CVLun, 
        "CMus"=CMus, "CVMus"=CVMus, 
        "CSki"=CSki, "CVSki"=CVSki, 
        "CSpl"=CSpl, "CVSpl"=CVSpl, 
        "Cpan"=Cpan, "CVPan"=CVPan, 
        "CKb" = CKb,  "CVK" = CVK, "CPTC" = CPTC,
        "Cfil" = Cfil, 
        "CL" = CL, "CVL" = CVL, 
        "CGI" = CGI, "CVGI" = CVGI, 
        "CA_free" = CA_free, "CA" = CA)
  })
}
