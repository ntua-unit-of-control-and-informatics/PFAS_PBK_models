single_plot <- function(i, legend_postition, global_optimum, thetas,
                        analysis_results = test, alpha = 0.95,
                        df = 1){
  output <- analysis_results$Likelihood_profiles
  data_to_plot <- output[[i]]$plik
  current_param <- names(data_to_plot)[1]
  names(data_to_plot)[1] <- "Parameter"
  optimal_value <- data.frame(thetas[i], global_optimum)
  names(optimal_value) <- c("Parameter", "Likelihood")
  data_to_plot[,1] <- data_to_plot[,1]
  cols <- c('Optimal value'='blue', 'Threshold'='red', 'Likelihood'='black')
  create_plot <- ggplot()+
    geom_hline(aes(yintercept=global_optimum + qchisq(alpha,df), color = "Threshold"), linetype="dashed", size=1)+
    geom_hline(aes(yintercept=global_optimum , color = "Optimal value"), linetype="dashed", size=1)+
    geom_line(data = data_to_plot,  aes(x=Parameter, y=Likelihood, color = 'Likelihood'), size=1.3)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=7)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, colour="pink", size=6.5)+
    geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=3.5)+
    
    scale_colour_manual(name='Lines', values = cols)+
    
    labs(#title = paste0( current_param), #"Profile Likelihood of ",
      y = expression(paste(chi^2, "(", theta, ")")) , x = x_axis_names[i])+
    theme(plot.title = element_text(hjust = 0.5,size=15), 
          axis.title.y =element_text(hjust = 0.5,size=15),
          axis.text.y=element_text(size=13),
          axis.title.x =element_text(hjust = 0.5,size=15),
          axis.text.x=element_text(size=13),
          panel.border = element_rect(colour = "black", fill=NA, size=1.0),
          
          legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          legend.position = legend_postition,
          legend.background = element_rect(fill="lightgrey",
                                           size=0.5, linetype="solid",
                                           colour ="black"))
  return(create_plot)
}


x_axis_names <- c(TeX("$P_{liver}$"), TeX("$P_{muscle}$"), TeX("$P_{kidney}$"), 
                  TeX("$P_{skin}$"), TeX("$P_{gills}$"),TeX("$P_{carcass}$"),
                  TeX("$P_{viscera}$"))
plot_list <- list()
for (i in 1:length(test$Likelihood_profiles)) {
  current_plot <- single_plot(i,'none', global_optimum, thetas)
  plot_list[[i]] <- plot
  assign(paste0('plot_', as.character(i)), current_plot)
}

plot_with_legend <- single_plot(1, c(0.9,5.3), global_optimum, thetas)


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
final_plot <- grid.arrange(arrangeGrob(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6, plot_7, nrow = 2),
                           shared_legend, ncol = 1, heights = c(10, 1))
final_plot
