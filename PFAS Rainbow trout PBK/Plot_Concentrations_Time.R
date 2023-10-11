library(ggplot2)
library(gridExtra)

create_plot <- function(substance, legend_postition){  
  # Keep the experimental data of current subsstance
  data_df <- data_list[[substance]]
  # Generate the sd values
  CV <- 50/100
  errors_df <- data.frame(matrix(NA, nrow = nrow(data_df), ncol = ncol(data_df)))
  for (k in 1:nrow(data_df)) {
    for (j in 2:ncol(data_df)) {
      set.seed(100)
      errors_df[k,j] <- abs(rnorm(1, data_df[k,j]*CV, 1))
    }
  }
  errors_df[,1] <- data_df[,1]
  colnames(errors_df) <- colnames(data_df) 
  
  current_params <- substance_specific_params[substance,]
  
  user.input <- list('substance'=substance,
                     'Texp'=15,
                     'admin.dose_dietary'=admin.dose_dietary,
                     'admin.time_dietary'=admin.time_dietary)
  params <- create.params(user.input)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- seq(0,56*24,2)
  
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, 
                                      parms = c(current_params, common_params, params),
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  # Keep the predictions only for the time points at which there are available data
  predictions_df <- solution[,c('time' , 'C_Liver', 'C_Blood',
                                'C_Skin', 'C_Muscle', 'C_Gills',
                                'C_Kidney', 'C_Carcass')]
  predictions_df[,-1] <- predictions_df[,-1]*1000
  exp_data <- data_list[[substance]]
  
  compartments <- colnames(exp_data)[2:8]
  color_codes <- scales::hue_pal()(length(compartments))
  names(color_codes) <-  colnames(exp_data)[2:8]
  
  plot <- ggplot()+
    geom_line(data = predictions_df, aes(x = time, y = C_Liver, color='Liver'), size=0.8)+
    geom_line(data = predictions_df, aes(x = time, y = C_Blood, color='Blood'), size=0.8)+
    geom_line(data = predictions_df, aes(x = time, y = C_Skin, color='Skin'), size=0.8)+
    geom_line(data = predictions_df, aes(x = time, y = C_Muscle, color='Muscle'), size=0.8)+
    geom_line(data = predictions_df, aes(x = time, y = C_Gills, color='Gills'), size=0.8)+
    geom_line(data = predictions_df, aes(x = time, y = C_Kidney, color='Kidney'), size=0.8)+
    geom_line(data = predictions_df, aes(x = time, y = C_Carcass, color='Carcass'), size=0.8)+
    
    geom_point(data = exp_data, aes(x = Time, y = Liver, color = 'Liver'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Blood, color = 'Blood'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Skin, color = 'Skin'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Muscle, color = 'Muscle'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Gills, color = 'Gills'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Kidney, color = 'Kidney'), size=5)+
    geom_point(data = exp_data, aes(x = Time, y = Carcass, color = 'Carcass'), size=5)+
    #scale_y_log10(limits = c(1, 600))+
    #ylim(c(1, 600))+
    geom_vline(xintercept = 28*24, size=1.0)+ # end of uptake
    
    
    
    labs(title = paste0(substance),
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
          legend.position = legend_postition,
          legend.background = element_rect(fill="lightgrey", 
                                           size=0.5, linetype="solid",
                                           colour ="black"))
  
  #print(plot)
  return(plot)}

plots_list = list()
for (i in 1:length(substances)) {
  substance <- substances[i]
  plot <- create_plot(substance, 'none')
  plots_list[[i]] <- plot
  assign(paste0('plot_', as.character(i)), plot)
}

plot_with_legend <- create_plot('PFOS', 'bottom')

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(plot_with_legend)
# Draw plots with shared legend
final_plot <- grid.arrange(arrangeGrob(plot_1, plot_2, plot_3, plot_4, plot_5, ncol = 2),
                            shared_legend, nrow = 2, heights = c(10, 1))

