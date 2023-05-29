setwd('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK')

# Experimental data from Falk et al.2015
#---------------------------------------
# The concentrations in the data are given in ug PFAS/kg tissue units. 
# The time is given in days and will be transformed in hours, to be compatible
# with the model

# Directory of folder with saved data files
data_dir <- 'C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015'

# Load PFOS data
#---------------
PFOS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFOS.xlsx'))
PFOS_data$Time <- PFOS_data$Time*24

#-------------------------------------------------------------------------------

# Suppose that the input is a dataframe whose the 1st column is the timepoints (x)
# and each column contains the concentrations of a tissue (y).

df <- PFOS_data

Weighting.func <- function(df, a_first=2, a_last=2, a_tp=2){
  x_values <- df[,1] # Take the x values 
  y_values <- df[,-1] # Take th y values
  
  # Number of x values
  Nx <- dim(df)[1] 
  # Number of y columns (subtract one column because it's the x values)
  Ny <- dim(y_values)[2]
  
  # Transform the data df to a list of dataframes. Each dataframe will have 
  # the x values and y values of a tissue.
  sub_list <- list()
  
  # Create a list to save the weights
  weights_list <- list()
  for (i in 1:Ny) { # Loop over the y outputs
    
    # Create a 2-column dataframe with x and y values
    sub_list[[i]] <- cbind(x_values, y_values[,i])
    names(sub_list)[i] <- colnames(y_values)[i]
    colnames(sub_list[[i]]) <- colnames(df)[c(1,i+1)]
    
    weights_list <- append(weights_list, NA)
    # Select a df from the sub_list to calulcate its weights
    sub_df <- sub_list[[i]]
    N <- dim(sub_df)[1] # Number of x values of sub_df
    
    # Estimate the weight of the 1st point
    weights_list[[i]][1] <- a_first*(abs(sub_df[1,1] - sub_df[2,1]) + abs(sub_df[1,2] - sub_df[2,2]))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
    # Estimate the weight of the Last point
    weights_list[[i]][N] <- a_last*(abs(sub_df[N,1] - sub_df[N-1,1]) + abs(sub_df[N,2] - sub_df[N-1,2]))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
    
    # next_point indicates if j+1 point is a turning point
    next_point <- 0
    prev_point<- 0
    for (j in 2:(N-2)) {
      if(next_point){ # j is defined as turning point in previous loop
        a <- a_tp
        weights_list[[i]][j] <- (abs(sub_df[j,1] - sub_df[j-1,1]) + abs(sub_df[j,1] - sub_df[j+1,1]) +
                                   a*(abs(sub_df[j,2] - sub_df[j-1,2]) + abs(sub_df[j,2] - sub_df[j+1,2])))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
        
        # Check if the j+1 point is turning point to define next_point value for the next loop
        if((sub_df[j+1,2]<sub_df[j,2] & sub_df[j+1,2]<sub_df[j+2,2]) | (sub_df[j+1,2]>sub_df[j,2] & sub_df[j+1,2]>sub_df[j+2,2])){
          prev_point <- T
          next_point <- T
        }else{
          # prev_point is used to identify if the point of the next loop (which is not a turning point) is after
          # a turning point
          prev_point <- T
          next_point <- F
        }
      }else{ # j is not defined as turning point in previous loop
        # Check if the j+1 point is turning point to define next_point value for the next loop
        # or check if the previous point was a turning point
        if((sub_df[j+1,2]<sub_df[j,2] & sub_df[j+1,2]<sub_df[j+2,2]) | (sub_df[j+1,2]>sub_df[j,2] & sub_df[j+1,2]>sub_df[j+2,2])){
          a <- a_tp
          prev_point <- F
          next_point <- T
        }else{
          a <- ifelse(prev_point, a_tp, 1)
          prev_point <- F
          next_point <- F
        }
        weights_list[[i]][j] <- (abs(sub_df[j,1] - sub_df[j-1,1]) + abs(sub_df[j,1] - sub_df[j+1,1]) +
                                   a*(abs(sub_df[j,2] - sub_df[j-1,2]) + abs(sub_df[j,2] - sub_df[j+1,2])))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
      }
    }
    
    a <- ifelse(next_point | prev_point, a_tp, 1)
    j <- N-1
    weights_list[[i]][j] <- (abs(sub_df[j,1] - sub_df[j-1,1]) + abs(sub_df[j,1] - sub_df[j+1,1]) +
                               a*(abs(sub_df[j,2] - sub_df[j-1,2]) + abs(sub_df[j,2] - sub_df[j+1,2])))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
  }
  names(weights_list) <- colnames(y_values)
  return(weights_list)
}

weight_values <- Weighting.func(df)
weight_values



