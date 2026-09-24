  library(deSolve)
  library(tidyverse)

  #=========================
  # 1. Parameters of the model
  #=========================

#'  
#' @export
#'  
  .create.params <- function(user_input) {
  with( as.list(user_input),{

variables <-data.frame(
  PFHpA = c(364.06, 4.50E-03, 3300, 263, 4500, 60, 2.63E+01,
            1.79E-04, 7.76E-04, 3.8, 3.72E+03),
  PFOA = c(414.07, 2.45E-03, 2500, 137.5, 4500, 47, 0.007757613, 
           0.005284911, 2.40E-03, 3.96, 2.4),
  PFNA = c(464.08, 1.87E-03, 1000, 42, 8500, 58, 2.161081781, 
           0.000943089, 0.002219551, 4.02, 1.73968375),
  PFDA = c(514.09, 6.18E-04, 2500, 137.5, 6000, 39, 78.37437317, 
           7.29E-05, 0.002076146, 4.01, 0.156992915),
  PFBS = c(300.1, 3.65E-02, 2500, 137.5, 4500, 47, 0.00E+00, 
           0.00054072, 0.00E+00, 18.42, 5.525079871),
  PFHxS = c(403.11, 6.95E-04, 5800, 778, 7300, 92, 44.76000058, 
            0.000165342, 2.40E-03, 3.76, 0.959152914),
  PFOS = c(500.13, 0.000753, 2500, 137.5, 2200, 48, 0.000424071, 
           0.700484202, 0, 3.78, 0.069918193),
  DONA = c(378.07, 0.0268, 2500, 137.5, 4500, 39, 5.68E-05, 
           8.11E-05, 0.015495299, 3.71, 98.46893334),
  HFPO_DA = c(330.0489, 0.00385, 2500, 137.5, 4500, 47, 
          4.39E-02, 4.75E-06, 6.12E-02, 19.58, 1.41),
  PFBA = c(214.04, 2.29E-01, 2500, 137.5, 4500, 47, 0, 
           0.000224546, 0.00E+00, 4.01, 50.4078006),
  PFHxA = c(314.05, 3.80E-02, 2500, 137.5, 4500, 47, 0.011660009, 
            1.00E-07, 6.73E-02, 3.99, 1.22444103)
) %>% select(all_of(chemical)) 

# Set the row names to match the 'Parameters' column
rownames(variables) <- c('MW', 'Free', 'Vmax_baso_invitro', 'Km_baso', 'Vmax_apical_invitro', 
                            'Km_apical', 'RAFbaso', 'RAFapi', 'kbilec', 'admin_dose_bolus', 'keffluxc')

 PC<- data.frame(
  DONA = c(0.124947389, 0.093417353, 0.453895406, 0.39404558, 0.516525519, 0.260245376, 0.802893425, 0.245729715, 0.333327203, 0.554155414),
  HFPO_DA = c(0.132911327, 0.167241461, 0.335472967, 0.300497404, 0.363279275, 0.993256185, 0.642857005, 0.350626131, 0.297788612, 0.405961151),
  PFBA = c(0.129165456, 0.127429659, 0.278016846, 0.126434927, 0.325021589, 0.361492503, 0.603411764, 0.156031705, 0.304498056, 0.378143683),
  PFBS = c(0.197947634, 0.350683638, 0.429465629, 0.408625829, 0.450282602, 0.458664426, 0.673479649, 0.220412593, 0.409674109, 0.478156735),
  PFDA = c(0.050113807, 0.497299996, 0.89483373, 0.018004527, 0.073924732, 3.409634891, 0.040734602, 0.02085036, 0.094637375, 0.840556154),
  PFHpA = c(0.098067143, 0.054941827, 0.527979623, 0.231463023, 0.583097047, 1.168406124, 0.841703023, 0.279995413, 0.39715958, 0.623833077),
  PFHxA = c(0.145771589, 0.054376066, 0.485393186, 0.213538328, 0.543916041, 0.263938131, 0.814340373, 0.260989886, 0.390107977, 0.59141182),
  PFHxS = c(0.254116384, 0.053076379, 0.430516787, 0.378514537, 0.504146317, 0.244750669, 0.079558592, 0.118922026, 0.318553055, 0.529997757),
  PFNA = c(0.136107702, 0.166395692, 0.662236244, 0.655397113, 0.064056921, 2.566865548, 1.487738844, 0.016102782, 0.666595171, 0.02198743),
  PFOA = c(0.063812225, 0.158920559, 0.406905162, 0.11934646, 0.460418226, 1.261933405, 0.188237566, 0.109461608, 0.172308549, 0.985281684),
  PFOS = c(0.311979033, 0.311341927, 0.560484488, 1.131641191, 1.02340325, 1.071556994, 0.239982717, 0.104992399, 1.18549103, 1.05033605)
) %>% select(all_of(chemical)) 

# Set Organs as row names
rownames(PC) <- c("Adipose", "Brain", "Gonads", "Gut", "Heart", "Liver",
"Lung", "Muscle", "Skin", "Spleen")

                       
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
    QCC <- 12.5 * inv_time_scale  # cardiac output in L/time scale/kg^0.75; Brown 1997
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
    GFRC <- 24.19 * inv_time_scale  # glomerular filtration rate (L/time scale/kg kidney); Corley 2005

    # --- Partition Coefficients (from Allendorf 2021) ---
    PR <- 0.1 # rest of body:blood#; fitted to model
    PAdi <- (1-Htc) * PC["Adipose",] # adipose tissue:blood;
    PBra <- (1-Htc) * PC["Brain",]  # brain:blood;
    PGon <- (1-Htc) * PC["Gonads",]  # gonads:blood;
    PGI <- (1-Htc) * PC["Gut",]# GI tract:blood;
    PHea <- (1-Htc) * PC["Heart",]  # heart:blood;
    PL<- (1-Htc) * PC["Liver",]  # liver:blood;
    PLun <- (1-Htc) * PC["Lung",]  # lung:blood;
    PMus <- (1-Htc) * PC["Muscle",] # muscle:blood;
    PSki <- (1-Htc) * PC["Skin",]  # skin:blood;
    PSpl <- (1-Htc) * PC["Spleen",]  # spleen:blood;
    PPan <- (PGI+PSpl)/2  # pancreas:blood; estimated as average of GI tract and spleen

    # --- Rate Constants ---
    kdif <- 0.001 * inv_time_scale  # diffusion rate from proximal tubule cells (L/time scale)
    kabsc <-  2.12 * inv_time_scale  # rate of absorption from small intestine (1/(day*BW^-0.25))
    kunabsc <- 7.06e-5 * inv_time_scale  # rate of unabsorbed dose to feces (1/(day*BW^-0.25)); fitted to model 
    keffluxc <-variables["keffluxc", ] * inv_time_scale  # rate of efflux from PTC to blood (1/(day*BW^-0.25))
    kbilec <- variables["kbilec", ] * inv_time_scale # biliary elimination rate (1/(day*BW^-0.25)); fitted to model 
    kurinec <- 0.063 * inv_time_scale  # urinary elimination rate (1/(day*BW^-0.25))
    
    # --- Water Consumption ---
    water_consumption <- 1.36  # L/time scale

    # === SCALED PARAMETERS (calculated from above) ===

    # Cardiac output and blood flows
    QC <- QCC * (BW^0.75) * (1 - Htc) # cardiac output in L/time scale; adjusted for plasma
    QK <- (QKC * QC)  # plasma flow to kidney (L/time scale)
    QL <- (QLC * QC)  # plasma flow to liver (L/time scale)
    QAdi <- (QAdiC * QC)  # plasma flow to adipose tissue (L/time scale)
    QBra <- (QBraC * QC)  # plasma flow to brain (L/time scale)
    QGon <- (QGonC * QC)  # plasma flow to gonads (L/time scale)
    QHea <- (QHeaC * QC)  # plasma flow to heart (L/time scale)
    QLun <- (QLunC * QC)  # plasma flow to lung (L/time scale)
    QMus <- (QMusC * QC)  # plasma flow to muscle (L/time scale)
    QSki <- (QSkiC * QC)  # plasma flow to skin (L/time scale)
    QSpl <- (QSplC * QC)  # plasma flow to spleen (L/time scale)
    QPan <- (QPanC * QC)  # plasma flow to pancreas (L/time scale)
    QGI <- (QGIC * QC)  # plasma flow to GI tract (L/time scale)
    QR <- QC - QK - QL - QAdi - QBra - QGon - QHea - QLun - QMus - QSki - QSpl - QPan - QGI  # plasma flow to rest of body (L/time scale)

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
    VBal <- (0.93 * BW) - (VR + VL + VPTC + Vfil + VPlas + VAdi 
    + VBra + VGon + VHea + VLun + VMus + VSki + VSpl + VPan + VGI)  # Balance check; should equal zero

    Vmax_basoC <- (Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24 * inv_time_scale 
    Vmax_apicalC <- (Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1e6) * 24 * inv_time_scale
    Vmax_baso <- Vmax_basoC * BW^0.75  # (ug/time scale)
    Vmax_apical <- Vmax_apicalC * BW^0.75  # (ug/time scale)
    kbile <- kbilec * BW^(-0.25)  # biliary elimination; liver to feces storage (/time scale)
    kurine <- kurinec * BW^(-0.25)  # urinary elimination, from filtrate (/time scale)
    kefflux <- keffluxc * BW^(-0.25)  # efflux clearance rate, from PTC to blood (/time scale)
    GFR <- GFRC*VK # glomerular filtration rate, (L/time scale)
    # GI Tract Parameters
    kabs <- kabsc * BW^(-0.25)  # rate of absorption from small intestine (/time scale)
    kunabs <- kunabsc * BW^(-0.25)  # rate of unabsorbed dose to feces (/time scale)


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
      "PSpl" = PSpl, "PPan" = PPan,"water_consumption" = water_consumption,
      "admin_type" = admin_type, "admin_dose" = admin_dose,
      "admin_time" = admin_time, "ingestion" = ingestion,
      "ingestion_time" = ingestion_time ,  "duration" = duration,
      "time_scale" = time_scale, "exp_type" = exp_type))
    
   })
}
  #===============================================
  #2. Function to create initial values for ODEs 
  #===============================================
#'  
#' @export
#'  
  .create.inits <- function(parameters){
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
  .create.events <- function(parameters){
  with(as.list(parameters), {
    if (tolower(admin_type) == "iv") {
      ldose <- length(admin_dose)
      ltimes <- length(admin_time)
      if (ltimes != ldose) {
        stop("The times of administration should be equal in number to the doses")
      }
      events <- list(data = data.frame(
        var = Aplas_free, time = admin_time,
        value = admin_dose, method = "add"
      ))

    } else if (admin_type == "bolus") {
      ldose <- length(admin_dose)
      ltimes <- length(admin_time)
      if (ltimes != ldose) {
        stop("The times of administration should be equal in number to the doses")
      }
      events <- list(data = data.frame(
        var = "ALumen", time = admin_time,
        value = admin_dose, method = "add"
      ))   
    } else if (admin_type == "oral") {
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
      stop("admin_type should be either 'iv', 'bolus', or 'oral'")
    }
      return(events)
  
  }}
  )
  }

  #==================
  # 4. Custom functions
  #==================
  #'  
  #' @export
  #'  
  .custom.func <- function(){
  return()
}

  #==============
  # 5. ODEs System
  #==============

#'  
#' @export
#'  
  .ode.func <- function(time, inits, params, custom.func) {
    with(as.list(c(inits, params)), {

      # Concentrations in various compartments
      
      CR <- AR / VR  # concentration in rest of body (ug/L)
      CVR <- CR / PR  # concentration in venous blood leaving rest of body (ug/L)
      CR_free <- CR * Free  # free concentration in rest of body (ug/L)

      CAdi <- AAdi / VAdi  # concentration in adipose tissue (ug/L)
      CVAdi <- CAdi / PAdi  # concentration in venous blood leaving adipose
      CAdi_free <- CAdi * Free  # free concentration in adipose tissue (ug/L)
      
      CBra <- ABra / VBra  # concentration in brain tissue (ug/L)
      CVBra <- CBra / PBra  # concentration in venous blood leaving brain
      CBra_free <- CBra * Free  # free concentration in brain tissue (ug/L)
      
      CGon <- AGon / VGon  # concentration in gonads tissue (ug/L)
      CVGon <- CGon / PGon  # concentration in venous blood leaving gonads
      CGon_free <- CGon * Free  # free concentration in gonads tissue (ug/L)
      
      CHea <- AHea / VHea  # concentration in heart tissue (ug/L)
      CVHea <- CHea / PHea  # concentration in venous blood
      CHea_free <- CHea * Free  # free concentration in heart tissue (ug/L)
      
      CLun <- ALun / VLun  # concentration in lung tissue (ug/L)
      CVLun <- CLun / PLun  # concentration in venous blood leaving lung
      CLun_free <- CLun * Free  # free concentration in lung tissue (ug/L)
      
      CMus <- AMus / VMus  # concentration in muscle tissue (ug/L)
      CVMus <- CMus / PMus  # concentration in venous blood leaving muscle
      CMus_free <- CMus * Free  # free concentration in muscle tissue (ug/L)

      CSki <- ASki / VSki  # concentration in skin tissue (ug/L)
      CVSki <- CSki / PSki  # concentration in venous blood leaving skin
      CSki_free <- CSki * Free  # free concentration in skin tissue (ug/L)
      
      CSpl <- ASpl / VSpl  # concentration in spleen tissue (ug/L)
      CVSpl <- CSpl / PSpl  # concentration in venous blood leaving spleen
      CSpl_free <- CSpl * Free  # free concentration in spleen tissue (ug/L)

      Cpan <- APan / VPan  # concentration in pancreas tissue (ug/L)
      CVPan <- Cpan / PPan  # concentration in venous blood leaving pancreas
      CPan_free <- Cpan * Free  # free concentration in pancreas tissue (ug/L)

      CGI <- AGI / VGI  # concentration in GI tract (ug/L)
      CVGI <- CGI / PGI  # concentration in venous blood leaving GI tract
      CGI_free <- CGI * Free  # free concentration in GI tract (ug/L)
      
      CKb <- AKb / VKb  # concentration in kidney blood (ug/L)
      CVK <- CKb  # concentration in venous blood leaving kidney (ug/L)
      CPTC <- APTC / VPTC  # concentration in PTC (ug/L)
      Cfil <- Afil / Vfil  # concentration in filtrate (ug/L)
      
      CL <- AL / VL  # concentration in the liver (ug/L)
      CLiver <- AL / ML  # concentration in the liver (ug/g)
      CVL <- CL / PL  # concentration in the venous blood leaving the liver (ug/L)
      CL_free <- CL * Free  # free concentration in the liver (ug/L)
      
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
        "CR" = CR, "CR_free" = CR_free, "CVR" = CVR, 
        "CAdi"=CAdi, "CAdi_free"=CAdi_free, "CVAdi"=CVAdi,
        "CBra"=CBra, "CBra_free"=CBra_free, "CVBra"=CVBra,
        "CGon"=CGon, "CGon_free"=CGon_free, "CVGon"=CVGon,
        "CHea"=CHea, "CVHea"=CVHea, "CHea_free"=CHea_free,
        "CLun"=CLun, "CVLun"=CVLun, "CLun_free"=CLun_free,
        "CMus"=CMus, "CVMus"=CVMus, "CMus_free"=CMus_free,
        "CSki"=CSki, "CVSki"=CVSki, "CSki_free"=CSki_free,
        "CSpl"=CSpl, "CVSpl"=CVSpl, "CSpl_free"=CSpl_free,
        "Cpan"=Cpan, "CVPan"=CVPan, "Cpan_free"=CPan_free,
        "CKb" = CKb,  "CVK" = CVK, "CPTC" = CPTC,
        "Cfil" = Cfil, 
        "CL" = CL, "CVL" = CVL, "CL_free" = CL_free,
        "CGI" = CGI, "CVGI" = CVGI, "CGI_free" = CGI_free,
        "CA_free" = CA_free, "CA" = CA)
  })
}
