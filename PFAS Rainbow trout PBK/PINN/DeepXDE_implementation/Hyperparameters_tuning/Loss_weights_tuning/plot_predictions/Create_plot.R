library(ggplot2)
library(gridExtra)

create.params <- function(user.input){
  with(as.list(user.input),{
    # Transform input temperature into Kelvin scale
    Texp <- 273 + Texp # K
    
    Tref <- 273 + c(6,12,18) # Reference Temperature K - Grech et al.2018
    keep_ref_value <- which.min(abs(Tref - Texp))
    
    # Cardiac output reference value at T = 6 C (Barron et al. 1987, Table II)
    F_card_ref_6 <- 1.188 # ml/h/g 
    # Cardiac output reference value at T = 12 C (Barron et al. 1987, Table II)
    F_card_ref_12 <- 2.322 # ml/h/g 
    # Cardiac output reference value at T = 18 C (Barron et al. 1987, Table II)
    F_card_ref_18 <- 3.75 # ml/h/g 
    F_card_ref_values <- c(F_card_ref_6,F_card_ref_12,F_card_ref_18)
    F_card_ref <- F_card_ref_values[keep_ref_value]
    
    # Body weight reference value at T = 6 C (Barron et al. 1987, Table II)
    BW_ref_6 <- 270.1 # g 
    # Body weight reference value at T = 12 C (Barron et al. 1987, Table II)
    BW_ref_12 <- 296.4 # g 
    # Body weight reference value at T = 18 C (Barron et al. 1987, Table II)
    BW_ref_18 <- 414.5 # g
    BW_ref_values <- c(BW_ref_6,BW_ref_12,BW_ref_18)
    BW_ref <- BW_ref_values[keep_ref_value]
    
    # Arrhenius Temperature function 
    TA <- 6930 # Arrhenius Temperature K - Grech et al.2018
    Tr <- Tref[which.min(abs(Tref - Texp))]
    KT <- exp(TA/Tr - TA/Texp)
    
    # Load the xlsx file with the physiological params pf rainbow trout
    phys.params <- openxlsx::read.xlsx('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 1)
    
    # Keep only the physiological parameters from the paper of Vidal et al. 2019
    # fw are the fractions of tissue_weight/total_weight (unitless)
    fw <- phys.params[phys.params$Source=='Vidal et al. 2019', c('Liver', 'Blood', 'Skin',
                                                                 'Muscle', 'Gills', 'Kidney', 'Viscera')]
    fw_Liver <- fw$Liver
    fw_Blood <- fw$Blood
    fw_Skin <- fw$Skin
    fw_Muscle <- fw$Muscle
    fw_Gills <- fw$Gills
    fw_Kidney <- fw$Kidney
    fw_Viscera <- fw$Viscera
    fw_lumen <- 0.012
    
    # Load the xlsx file with the physiological params pf rainbow trout
    phys.params <- openxlsx::read.xlsx('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 2)
    
    # Keep only the physiological parameters from the paper of Vidal et al. 2019
    # fb are the fractions of blood flow of each tissue (unitless)
    fb <- phys.params[phys.params$Source=='Vidal et al. 2019', c('Liver', 'Skin', 'Muscle',
                                                                 'Gills', 'Kidney', 'Viscera')]
    
    fb_Liver <- fb$Liver
    fb_Skin <- fb$Skin
    fb_Muscle <- fb$Muscle
    fb_Gills <- fb$Gills
    fb_Kidney <- fb$Kidney
    fb_Viscera <- fb$Viscera
    
    #Ku <- 0.13 # 1/h
    
    # Reabsorption coefficients from bile to intestine
    # estimated by Cao et al., 2022
    # K_urine = Cl_urine/f_reab_urine estimated by Ng et al., 2013 (unitless)
    if(substance=='PFOA'){
      a <- 0.138 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.30 # Cao et al., 2022
      K_urine <- 2.08 
      Cl_urine <- 0.029*3600 # 1/h (Sun et al., 2022)
    }else if(substance=='PFNA'){
      a <- 0.522 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.34 # Cao et al., 2022
      K_urine <- 1.35
      Cl_urine <- 0.050*3600 # 1/h (Sun et al., 2022)
    }else if(substance=='PFBS'){
      a <- 0.0598 # Goeritz et al.2013
      f_reab_hep <- 0.23 # Cao et al., 2022
      K_urine <- 5.88
      Cl_urine <- 0.023*3600 # 1/h (Sun et al., 2022) # Assumed equal to PFHxS
    }else if(substance=='PFHxS'){
      a <- 0.558 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.30 # Cao et al., 2022
      K_urine <- 5.88
      Cl_urine <- 0.023*3600 # 1/h (Sun et al., 2022)
    }else if(substance=='PFOS'){
      a <- 0.721 # Sun et al., 2022, Goeritz et al.2013
      f_reab_hep <- 0.42 # Cao et al., 2022
      K_urine <- 1.35
      Cl_urine <- 0.050*3600 # 1/h (Sun et al., 2022)
    }
    
    # Bile flow coefficient
    Q_bile_coef <- 7.5e-05 # ml/g BW/h Grosell et al., 2000
    Q_urine_coef <- 2.755e-03 # ml/h/g of BW Urinary flow rate (Curtis et al., 1981)
    V_urine_coef <- 2.2e-03 # ml/g of BW Urine volume inside urinary bladder (Curtis et al., 1981)
    
    a_skin <- 0.9 # 90% of venous blood of skin was assumed to flow directly to kidney (Nichols et al.1996)
    a_muscle <- 0.6 # 60% of venous blood of muscle was assumed to flow directly to kidney (Nichols et al.1996)
    
    plasma <- 0.7
    
    return(list('F_card_ref'=F_card_ref, 'BW_ref'=BW_ref, 'KT'=KT, 
                
                'fw_Liver'=fw_Liver, 'fw_Blood'=fw_Blood, 'fw_Skin'=fw_Skin, 
                'fw_Muscle'=fw_Muscle, 'fw_Gills'=fw_Gills, 'fw_Kidney'=fw_Kidney,
                'fw_Viscera'=fw_Viscera, 'fw_lumen'=fw_lumen, 
                
                'fb_Liver'=fb_Liver, 'fb_Skin'=fb_Skin, 'fb_Muscle'=fb_Muscle,
                'fb_Gills'=fb_Gills, 'fb_Kidney'=fb_Kidney, 'fb_Viscera'=fb_Viscera,
                
                'a_skin'=a_skin, 'a_muscle'=a_muscle,
                'Q_bile_coef'=Q_bile_coef,
                'Q_urine_coef'=Q_urine_coef, 'V_urine_coef'=V_urine_coef,
                'K_urine'=K_urine, 'Cl_urine'=Cl_urine,
                'f_reab_hep'=f_reab_hep, 'plasma'=plasma,"a"=a))
  })
}
fish_weight <- function(time){
  x <- c(0,28,56)*24 
  y <- c(314, 655, 808)
  
  if(time <= x[1]){
    w = y[1]
  }else if(time >= x[3]){
    w = y[3]
  }else if(time >= x[1] & time < x[2]){
    w = approx(x=x[1:2],y=y[1:2], xout = time)$y
  }else if(time >= x[2] & time < x[3]){
    w = approx(x=x[2:3],y=y[2:3], xout = time)$y
  }
  return(w)
}


user.input <- list('substance'='PFOS',
                   'Texp'=15)

time_specific_params <- function(time, params){
  with(as.list(params),{
    BW = fish_weight(time)
    # Calculate the mass of each tissue - g
    w_blood <- fw_Blood*BW*plasma     # Blood mass - g
    w_liver <- fw_Liver*BW     # Liver mass - g
    w_skin <- fw_Skin*BW       # Skin weight - g
    w_muscle <- fw_Muscle*BW   # Muscle weight - g
    w_gills <- fw_Gills*BW     # Gills weight - g
    w_kidney <- fw_Kidney*BW   # Kidney weight - g
    w_viscera <- fw_Viscera*BW # Viscera weight - g
    w_lumen <- fw_lumen*BW
    w_art <- 1/3*w_blood
    w_venous <- 2/3*w_blood
    w_carcass <- BW - (w_blood/plasma + w_liver + w_skin + w_muscle +
                         w_gills + w_kidney + w_viscera + w_lumen)
    
    return(list('w_blood'=w_blood, 'w_liver'=w_liver, 'w_skin'=w_skin,
                'w_muscle'=w_muscle, 'w_gills'=w_gills, 'w_kidney'=w_kidney,
                'w_viscera'=w_viscera, 'w_lumen'=w_lumen, 'w_art'=w_art,
                'w_venous'=w_venous, 'w_carcass'=w_carcass))
  })
}

# Load experimental data
# Directory of folder with saved data files
data_dir <- '/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015'

# Load PFOS data
#---------------
PFOS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFOS.xlsx'))
PFOS_data$Time <- PFOS_data$Time*24

params <- create.params(user.input)


create_plot <- function(file_name, legend_position, title_name){
  # Load the predictions of the PINN on experimental time points
  pinn_preds = read.csv(file = paste0('/Users/vassilis/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PINN/DeepXDE_implementation/Hyperparameters_tuning/Loss_weights_tuning/plot_predictions/',file_name) )
  
  # Keep only the time points and columns of interest
  pinn_preds <- pinn_preds[,c('Time', 'M_liver', 'M_blood', 'M_skin', 'M_muscle', 'M_gills', 
                              'M_kidney', 'M_carcass')]
  
  # Tranform the predicted variables from mass into concentrations (ug -> ug PFOS/kg tissue)
  pinn_preds_conc <- pinn_preds
  for (i in 1:dim(pinn_preds_conc)[1]) {
    current_time_params <- time_specific_params(i, params)
    pinn_preds_conc[i, 2:8] = pinn_preds[i, 2:8]/c(current_time_params$w_liver, current_time_params$w_blood,
                                                   current_time_params$w_skin, current_time_params$w_muscle,
                                                   current_time_params$w_gills, current_time_params$w_kidney,
                                                   current_time_params$w_carcass)*1000
  }
  
  names(pinn_preds_conc)[-1] <- c('C_Liver', 'C_Blood', 'C_Skin', 'C_Muscle', 'C_Gills', 'C_Kidney', 'C_Carcass')
  
  color_codes <- scales::hue_pal()(7)
  names(color_codes) <-  c('Liver', 'Blood', 'Skin', 'Muscle', 'Gills', 'Kidney',
                           'Carcass')
  
  plot <- ggplot()+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Liver, color='Liver'), size=0.8)+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Blood, color='Blood'), size=0.8)+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Skin, color='Skin'), size=0.8)+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Muscle, color='Muscle'), size=0.8)+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Gills, color='Gills'), size=0.8)+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Kidney, color='Kidney'), size=0.8)+
    geom_line(data = pinn_preds_conc, aes(x = Time, y = C_Carcass, color='Carcass'), size=0.8)+
    
    geom_point(data = PFOS_data, aes(x = Time, y = Liver, color = 'Liver'), size=5)+
    geom_point(data = PFOS_data, aes(x = Time, y = Blood, color = 'Blood'), size=5)+
    geom_point(data = PFOS_data, aes(x = Time, y = Skin, color = 'Skin'), size=5)+
    geom_point(data = PFOS_data, aes(x = Time, y = Muscle, color = 'Muscle'), size=5)+
    geom_point(data = PFOS_data, aes(x = Time, y = Gills, color = 'Gills'), size=5)+
    geom_point(data = PFOS_data, aes(x = Time, y = Kidney, color = 'Kidney'), size=5)+
    geom_point(data = PFOS_data, aes(x = Time, y = Carcass, color = 'Carcass'), size=5)+
    #scale_y_log10(limits = c(1, 600))+
    #ylim(c(1, 600))+
    geom_vline(xintercept = 28*24, size=1.0)+ # end of uptake
    
    labs(title = bquote(omega* ~'='~ .(title_name)),
         y = expression('Concentration ('*mu*'g/kg)') , x = "Time (h)")+
    theme(plot.title = element_text(hjust = 0.5,size=20),
          axis.title.y =element_text(hjust = 0.5,size=15),#,face="bold"),
          axis.text.y=element_text(size=13),
          axis.title.x =element_text(hjust = 0.5,size=15),#,face="bold"),
          axis.text.x=element_text(size=13),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0)) +
    
    scale_color_manual("Tissues", values=color_codes)+
    theme(legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          legend.position = legend_position,
          legend.background = element_rect(fill="lightgrey", 
                                           size=0.5, linetype="solid",
                                           colour ="black"))
  return(plot)
}


files <- c('0.01_predictions.csv', '0.1_predictions.csv', '1.0_predictions.csv',
           '2.0_predictions.csv', '5.0_predictions.csv', '10.0_predictions.csv',
           '100.0_predictions.csv')

plots_list = list()
counter <-1
plots_titles <- c('0.01', '0.1', '1.0', '2.0', '5.0', '10.0', '100.0')
for (i in files){
  plot <- create_plot(i, 'none', plots_titles[counter])
  plots_list[[i]] <- plot
  assign(paste0('plot_', as.character(counter)), plot)
  counter = counter+1
}

plot_with_legend <- create_plot(files[1], 'bottom', plots_titles[1])

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legendxa
shared_legend <- extract_legend(plot_with_legend)
# Draw plots with shared legend
final_plot <- grid.arrange(arrangeGrob(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7, ncol = 3),
                           shared_legend, nrow = 2, heights = c(10, 1))
