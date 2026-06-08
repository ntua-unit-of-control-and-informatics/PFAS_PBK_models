#Reverse dosimetry runs for exposure reconstruction






library(tidyverse)
library(ggplot2)
library(deSolve)
library(nloptr)
library(patchwork)

source("Extended_Model.r")
source("Exposure.r")

# Dataframe containing the 5th, 50th, and 95th 
# percentiles of body weight for USA population
# Anthropometric Reference Data for
# Children and Adults: United States,
# August 2021–August 2023


#
# Data declaration 
#
# weight_data <- data.frame(
#   sex = c(
#     "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male"#,
#     #"Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female"
#   ),
#   age_group = c(
#     "20 and older", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80 and older"#,
#     #"20 and older", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80 and older"
#   ),
#   Percentile_5th = c(
#     60.9, 57.0, 63.6, 66.0, 61.0, 62.6, 60.8, 58.0#,
#     #50.7, 47.6, 51.4, 52.8, 53.1, 51.6, 49.8, 46.8
#   ),
#   Percentile_50th = c(
#     86.6, 83.7, 87.1, 88.4, 87.8, 86.9, 86.6, 79.9#,
#     #73.5, 68.9, 75.0, 74.8, 78.5, 74.3, 71.7, 65.7
#   ),
#   Percentile_95th = c(
#     131.8, 128.3, 132.5, 141.9, 138.7, 123.9, 121.3, 108.5#,
#     #119.2, 121.3, 125.9, 123.7, 121.4, 109.3, 103.5, 94.2
#   )
# )

Olsen<- read.csv("C:/Users/fotis/Documents/GitHub/PFAS_PBK_models/Extended PFAS PBK model/Olsen_biomonitoring.csv")

average_pfas_concentration <- function(data, chemical, year) {
  # Filter for the target chemical and year
  subset <- data %>% filter(Chemical == chemical, time == year)
  
  if (nrow(subset) == 0) {
    message(sprintf("No data found for %s in %s.", chemical, year))
    return(NA)
  }
  
  # Calculate & return the mean (ignores NA values if any)
  mean(subset$concentration, na.rm = TRUE)
}


#Set up functions

.Olsens_run <- function(chemical, BW, Olsen_bio_mon, time_scale = "years",
                        optim_tol = 1e-3, solver = "lsodes", rtol = 1e-4, atol = 1e-4){
  # Filter the biomonitoring data for the specified chemical and age group
  bio_mon <- Olsen_bio_mon[Olsen_bio_mon$Chemical == chemical , ]  
  bio_mon$time <- bio_mon$time - 1995
  bio_mon<- subset(bio_mon, select = c("time", "avg_concentration"))
  colnames(bio_mon) <- c("bio_time", "bio_con") 
  
  # Run reverse dosimetry to estimate ingestion rates
  results <- .reverse_dosimetry_exposure(
    chemical = chemical,
    BW = BW,
    duration = max(bio_mon$bio_time),  # duration based on the maximum time in biomonitoring data
    time_scale = time_scale,
    bio_mon = bio_mon,
    optim_tol = optim_tol,
    solver = solver,
    rtol = rtol,
    atol = atol
  )

return(results)
}

all_results <- list() 

Olsen_bio_mon<- Olsen %>%
  group_by(Chemical, time) %>%
  summarise(avg_concentration = mean(concentration, na.rm = TRUE))

print(Olsen_bio_mon)


for (chemical in unique(Olsen_bio_mon$Chemical)) {
   cat("Running reverse dosimetry for chemical:", chemical, "\n")

    results<- .Olsens_run(chemical = chemical, BW = 86.6, Olsen_bio_mon = Olsen_bio_mon)

     all_results[[chemical]] <- results

    print(results)
    
    }


saveRDS(all_results, file = "Olsen_reverse_dosimetry_results.rds")

all_results <- readRDS("Olsen_reverse_dosimetry_results.rds")
# Extract biomonitoring time points per chemical/age_group

results_df <- map_dfr(
  names(all_results),
  function(chem) {
    map_dfr(
      names(all_results[[chem]]),
      function(age) {
        res <- all_results[[chem]]
        df <- NULL
        
        # 1. If result is already a dataframe
        if (is.data.frame(res)) df <- res
        # 2. If result is a list, find the first dataframe inside it
        else if (is.list(res)) {
          dfs <- keep(res, is.data.frame)
          if (length(dfs) > 0) df <- dfs[[1]]
          # 3. Fallback: convert named vectors of equal length to dataframe
          else if (all(sapply(res, is.atomic))) {
            df <- tryCatch(as.data.frame(res), error = function(e) NULL)
          }
        }
        
        # Skip empty/failed extractions
        if (is.null(df) || nrow(df) == 0) return(tibble())
        
        # Add identifiers and move them to the front
        df %>%
          mutate(chemical = chem) %>%
          relocate(chemical, .before = 1)
      }
    )
  }
)

theme_elsevier <- function(base_size = 10, base_family = "Times") {
  (theme_minimal(base_size = base_size, base_family = base_family)
    + theme(
      # Panel & background
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title = element_text(face = "bold", size = rel(1.1)),
      axis.text = element_text(colour = "black", size = rel(0.95)),
      
      # Legend
      legend.position = "right",
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_blank(),
      legend.title = element_text(face = "bold", size = rel(1)),
      legend.text = element_text(size = rel(0.95)),
      legend.box = "vertical",
      
      # Titles & captions
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(size = rel(1), colour = "grey30", margin = margin(b = 12)),
      plot.caption = element_text(size = rel(0.85), colour = "grey50", hjust = 1),
      
      # Margins & layout
      plot.margin = margin(20, 20, 15, 15, "pt"),
      panel.spacing = unit(1, "lines"),
      
      # Ensure black-and-white print compatibility
      strip.background = element_rect(fill = "grey95", colour = "black"),
      strip.text = element_text(face = "bold", size = rel(0.95))
    ))
}
# Create a list to store individual plots
chemical_plots <- list()

# Loop through each chemical

  
  p <- ggplot(results_df, aes(x = bio_time+1995, y = exposure*1000/(365*86),colour = chemical)) +
    geom_point(size = 5, alpha = 0.9 , shape="diamond")+
        labs(
      title = paste("Exposure vs. Time"),
      x = "Time (years)",
      y = "log10(Exposure (ng/kg/day))"
    ) +
     scale_color_viridis_d(option = "D") +
    scale_y_log10() +
    theme_elsevier(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "right"
    )
  p

  
  # Save individual file
  ggsave(filename = paste0("combined_plot.png"), 
         plot = p, width = 8, height = 6)




combined_plot <- 
  chemical_plots[["PFOS"]] + 
  chemical_plots[["PFHxS"]] + 
  chemical_plots[["PFOA"]] + 
  chemical_plots[["PFNA"]] + 
  chemical_plots[["PFDA"]] +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) &  
  theme(
    legend.position = "bottom",
    #legend.justification = "right",
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.margin = margin(10, 10, 10, 10),
    legend.title = element_text(face = "bold")
  )

print(combined_plot)

ggsave(filename = "combined_plot.png", plot = combined_plot, width = 12, height = 17)
