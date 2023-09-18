library(parallel)
library(deSolve)
library(nloptr)
library(ggplot2)
library(gridExtra)
#=====================================#
#  Weighted Sum of Squared Residuals  #
#=====================================#

WSSR <- function(observed, predicted, weights, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted) || !is.list(weights)){
    stop(" The observations, predictions and weights must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted) || length(observed) != length(weights)){
    stop(" The observations, predictions and weights must have the same compartments")
  }
  
  # Define the number of observed outputs
  N_outputs <- length(predicted)
  # Define the number of observations per output
  N_obs <- rep(NA, N_outputs)
  
  # A vector to store the values of the weighted squared sum residuals of each compartment
  outputs_res <- c()
  for (i in 1:N_outputs) { # loop over the observed outputs
    N_obs[i] <- length(observed[[i]])
    
    # Check that all observed, predicted and weights vectors have the same length
    if(N_obs[i] != length(predicted[[i]]) || N_obs[i] != length(weights[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    # The number of observations for output i
    N <- N_obs[i]
    
    # Initiate a variable to estimate the sum of squared residuals for output j
    sq_weighted_res_j <- 0
    for (j in 1:N) { #loop over the experimental points i compartment i
      sq_weighted_res_j <- sq_weighted_res_j + ((observed[[i]][j] - predicted[[i]][j]) / weights[[i]][j])^2   
    }
    outputs_res[i] <- sq_weighted_res_j
  }
  
  WSSR_results <- sum(outputs_res)
  
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(WSSR_results) <- comp.names
  }else if (!is.null(names(observed))){
    names(WSSR_results) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(WSSR_results) <- names(predicted)
  } else if (!is.null(names(weights)) && is.null(comp.names) ){
    names(WSSR_results) <- names(weights)
  }
  
  return(WSSR_results)
}

#======================#
#  Profile Likelihood  #
#======================#

profile_likelihood <- function(obj_f, 
                               i,
                               thetas,
                               thetas_names, 
                               constant_params = NULL,
                               data_df,
                               errors_df,
                               lb, ub, N_samples,
                               alpha, df, q, global_optimum, 
                               min_step_coef, max_step_coef,
                               break_at_bounds = FALSE,
                               # nlopt settings for the main optimization problem
                               opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                           "xtol_rel" = 1e-06, 
                                           "ftol_rel" = 1e-06,
                                           "ftol_abs" = 0.0,
                                           "xtol_abs" = 0.0 ,
                                           "maxeval" = 300,
                                           "print_level" = 1),
                               # nlopt settings for the estimation of theta_step
                               opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                      "xtol_rel" = 1e-05, 
                                                      "ftol_rel" = 1e-05,
                                                      "ftol_abs" = 0.0,
                                                      "xtol_abs" = 0.0 ,
                                                      "maxeval" = 50,
                                                      "print_level" = 0),
                               create_txt = FALSE){
  
  # INPUT VARIABLES:
  # obj_f             function that returns the value of the objective function
  # i                 integer defining the position of the parameter to estimate the profile-likelihood
  # thetas            vector containing the optimal values of parameters 
  # thetas_names      vector of characters with the names of the parameters
  # constant_params   vector with any extra constant parameters of the model and their names
  # data_df           dataframe containing the data used in the obj_f 
  # errors_df         dataframe with the meassured experimental errors (or SD values) 
  # lb                vector with the lower bounds of the "thetas" parameters
  # ub                vector with the upper bounds of the "thetas" parameters
  # N_samples         integer defining the number of samples taken around the theta optimal
  #                   value (N samples will be taken from each side of the theta)
  # alpha             probability of chi-squared to estimate the quantile (it is
  #                   the "p" variable of qchisq() function)
  # df                degrees of freedom of qchisq()
  # global_optimum    scalar defining the global minimum of the objective function
  # q                 a variable used in the estimation of the adaptive step (see Raue et al., 2009)
  # max_step_coef     coefficient defining the maximum permitted step  
  # min_step_coef     coefficient defining the minimum permitted step
  # break_at_bounds   logical; if TRUE the the sampling of the parameter stops because the bounds were exceeded
  # opts              list with the options selected for the minimization of the objective function
  #                   (check the nloptr package for more details)
  # opts_theta_step   list with the options selected for the estimation of the adaptive step
  #                   (check the nloptr package for more details)
  # create_txt        logical; if TRUE a txt file will be created at the current working directory 
  #                   to save the samples of the parameters and the corresponding values of the objective function
  
  # Estimate the Delta_alpha parameter
  Delta_alpha <- qchisq(alpha, df)
  
  # Take the name of the i-th parameter
  theta_i_name <- names(thetas)[i] # the name of the theta_i parameter
  
  # Function to estimate the theta_step by solving the equation
  # chi^2(theta) - chi^2(theta_hat) - q*Delta_alpha = 0 
  theta_step_estimation <- function(theta_step, theta_last, obj_f, constant_params, index, current_score, q, Delta_alpha){
    i <- index
    x <- theta_last
    x[i] <- x[i] + theta_step
    chi2_last <- current_score
    
    chi2 <- obj_f(x = x, constant_theta = NULL, constant_theta_name = NULL, params_names = names(x),
                  constant_params = constant_params,
                  data_df = data_df, errors_df = errors_df)
    
    return(abs(chi2 - chi2_last - q*Delta_alpha))
  }
  
  # Set the threshold. The threshold is estimated as global_optimum + Delta_alpha 
  threshold =  global_optimum + Delta_alpha
  
  if(create_txt) sink(paste0("check_progress_",thetas_names[i], ".txt")) #keep a txt to check the progress while running
  
  # Forward search
  cat("Forward search begins...", "\n")
  # Initiate a counter for the number of iterations
  iter_counter <- 0
  
  # Inititate a parameter to check if the threshold is exceeded in each iteration
  current_score <- global_optimum
  
  # Create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  forward_results_df <- data.frame(matrix(NA, nrow=N_samples, ncol = 2))
  colnames(forward_results_df) <- c(theta_i_name, "Likelihood")
  
  # take as constant_theta the parameter that will be fixed to a constant value
  constant_theta <- thetas[i]
  
  # Initiate theta_last and set it equal with the theta_hat
  theta_last <- thetas
  
  # set the current score equal global optimum to estimate the first step
  current_score <- global_optimum 
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr::nloptr(x0 = 0.25*thetas[[i]],
                                  eval_f = theta_step_estimation,
                                  lb	= min_step_coef*thetas[[i]],
                                  ub = max_step_coef*thetas[[i]],
                                  opts = opts_theta_step,
                                  theta_last = theta_last, index = i,
                                  current_score = current_score, q = q, Delta_alpha=Delta_alpha,
                                  obj_f=obj_f, 
                                  constant_params=constant_params)$solution
    
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(thetas[[i]])){
      theta_step <- max_step_coef*abs(thetas[[i]])
    }else if(theta_step < min_step_coef*abs(thetas[[i]])){
      theta_step <- min_step_coef*abs(thetas[[i]])
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta + theta_step
    
    #Check if the constant_theta exceeded the corresponding upper boundary
    if(constant_theta > ub[i]){
      f_exit = 3
      break
    } 
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    # define the lower and upper bounds of the parameters
    set.seed(12312)
    optimization<- nloptr::nloptr(x0 = x0,
                                  eval_f = obj_f,
                                  lb	= lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  constant_params=constant_params,
                                  data_df = data_df,
                                  errors_df = errors_df)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step , constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
    forward_results_df[iter_counter,] <- c(constant_theta, optimization$objective)
    
    # update current score
    current_score <- optimization$objective
    
    #update theta_last vector with the new optimized values
    for (k in 1:length(theta_last)) {
      if(k==i){
        theta_last[k] <- constant_theta
      }else if(k>i){
        theta_last[k] <- optimization$solution[k-1]   
      }else{
        theta_last[k] <- optimization$solution[k]} 
    }
    
  }
  if(iter_counter >= N_samples){
    f_exit <- 1
  }else if(current_score > threshold){
    f_exit <- 2
  }
  cat("Forward search ended after ", iter_counter, "iterations.", "\n")
  
  
  
  # Backward search
  cat("Backward search begins...", "\n")
  # Initiate a counter for the number of iterations
  iter_counter <- 0
  # Inititate a parameter to check if the threshold is exceeded in each iteration
  current_score <- global_optimum
  # Create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  backward_results_df <- data.frame(matrix(NA, nrow=N_samples, ncol = 2))
  colnames(backward_results_df) <- c(theta_i_name, "Likelihood")
  # take as constant_theta the parameter that will be fixed to a constant value
  constant_theta <- thetas[i]
  # Initiate theta_last and set it equal with the theta_hat
  theta_last <- thetas
  # set the current score equal global optimum to estimate the first step
  current_score <- global_optimum 
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr::nloptr(x0 = 0.25*thetas[[i]],
                                  eval_f = theta_step_estimation,
                                  lb	= min_step_coef*thetas[[i]],
                                  ub = max_step_coef*thetas[[i]],
                                  opts = opts_theta_step,
                                  theta_last = theta_last, index = i,
                                  current_score = current_score, q = q, Delta_alpha=Delta_alpha,
                                  obj_f=obj_f, 
                                  constant_params=constant_params)$solution
    
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(thetas[[i]])){
      theta_step <- max_step_coef*abs(thetas[[i]])
    }else if(theta_step < min_step_coef*abs(thetas[[i]])){
      theta_step <- min_step_coef*abs(thetas[[i]])
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta - theta_step
    #Check if the constant_theta exceeded the corresponding lower boundary
    if(constant_theta < lb[i]){
      b_exit <- 3
      break
    } 
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    set.seed(12312)
    optimization<- nloptr::nloptr(x0 = x0,
                                  eval_f = obj_f,
                                  lb	= lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  constant_params=constant_params,
                                  data_df = data_df,
                                  errors_df = errors_df)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step , constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
    backward_results_df[iter_counter,] <- c(constant_theta, optimization$objective)
    
    # update current score
    current_score <- optimization$objective
    
    #update theta_last vector with the new optimized values
    for (k in 1:length(theta_last)) {
      if(k==i){
        theta_last[k] <- constant_theta
      }else if(k>i){
        theta_last[k] <- optimization$solution[k-1]   
      }else{
        theta_last[k] <- optimization$solution[k]} 
    }
  }
  if(iter_counter >= N_samples){
    b_exit <- 1
  }else if(current_score > threshold){
    b_exit <- 2
  }
  cat("Backward search ended after ", iter_counter, "iterations.", "\n")
  
  results_df <- rbind(backward_results_df, forward_results_df, c(thetas[i], global_optimum))
  results_df <- results_df[order(results_df[,1]),]
  results_df <- results_df[complete.cases(results_df),]
  if(create_txt) sink()
  return(list("plik"=results_df, "b_exit"=b_exit, "f_exit"=f_exit))
}

#============================#
#  Identifiability Analysis  #
#============================#

Identifiability_analysis <- function(obj_f, thetas, thetas_names, data_df, errors_df,
                                     lb, ub,
                                     N_samples = 50,
                                     alpha = 0.95, df = 1,
                                     q = 0.5,
                                     global_optimum, 
                                     min_step_coef = 1e-03, max_step_coef = 0.2,
                                     N_cores,
                                     constant_params = NULL,
                                     exported_to_cluster = NULL,
                                     break_at_bounds = FALSE,
                                     # nlopt settings for the main optimization problem
                                     opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                                 "xtol_rel" = 1e-06, 
                                                 "ftol_rel" = 1e-06,
                                                 "ftol_abs" = 0.0,
                                                 "xtol_abs" = 0.0 ,
                                                 "maxeval" = 300,
                                                 "print_level" = 0),
                                     # nlopt settings for the estimation of theta_step
                                     opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                            "xtol_rel_step" = 1e-05, 
                                                            "ftol_rel_step" = 1e-05,
                                                            "ftol_abs_step" = 0.0,
                                                            "xtol_abs_step" = 0.0 ,
                                                            "maxeval_step" = 50,
                                                            "print_level_step" = 0),
                                     create_txt = TRUE){
  
  # INPUT VARIABLES:
  # obj_f             function that returns the value of the objective function
  # thetas            vector containing the optimal values of parameters 
  # thetas_names      vector of characters with the names of the parameters
  # data_df           dataframe containing the data used in the obj_f 
  # errors_df         dataframe with the meassured experimental errors (or SD values) 
  # lb                vector with the lower bounds of the "thetas" parameters
  # ub                vector with the upper bounds of the "thetas" parameters
  # N_samples         integer defining the number of samples taken around the theta optimal
  #                   value (N samples will be taken from each side of the theta)
  # alpha             probability of chi-squared to estimate the quantile (it is
  #                   the "p" variable of qchisq() function)
  # df                degrees of freedom of qchisq()
  # q                 a variable used in the estimation of the adaptive step (see Raue et al., 2009)
  # global_optimum    scalar defining the global minimum of the objective function
  # min_step_coef     coefficient defining the minimum permitted step
  # max_step_coef     coefficient defining the maximum permitted step  
  # N_cores           integer defining the number of cores for the parallel processing
  # constant_params   vector with any extra constant parameters of the model and their names
  # break_at_bounds   logical; if TRUE the the sampling of the parameter stops because the bounds were exceeded
  # opts              list with the options selected for the minimization of the objective function
  #                   (check the nloptr package for more details)
  # opts_theta_step   list with the options selected for the estimation of the adaptive step
  #                   (check the nloptr package for more details)
  # create_txt        logical; if TRUE a txt file will be created at the current working directory 
  #                   to save the samples of the parameters and the corresponding values of the objective function
  
  
  # Number of parameters tested in the identifiability analysis
  N_parameters <- length(thetas)
  # prepare the input for parallel processing
  X <- vector("list", N_parameters)
  for(i in 1:N_parameters){
    X[[i]] <- list(index=i, thetas=thetas, thetas_names=thetas_names,
                   constant_params = constant_params,
                   data_df = data_df,
                   errors_df = errors_df,
                   lb=lb, ub=ub, N_samples=N_samples, 
                   alpha=alpha, df=df,
                   global_optimum = global_optimum,
                   q = q,
                   max_step_coef = max_step_coef,
                   min_step_coef = min_step_coef,
                   break_at_bounds = break_at_bounds,
                   opts = opts,
                   opts_theta_step=opts_theta_step,
                   create_txt = create_txt)
  }
  
  start.time <- Sys.time()
  # Set up the cluster.
  cluster <- makeCluster(N_cores)
  
  # parallel_func is a wrapper of the profile_likelihood function in order to 
  # take all the input variable of the profile_likelihood as a list and un-list
  # the to provide them to profile_likelihood. It just returns the output
  # of the profile_likelihood function.
  parallel_func <- function(X){
    with(as.list(X),{
      profile_likelihood(obj_f, 
                         i=index,
                         thetas,
                         thetas_names, 
                         constant_params,
                         data_df,
                         errors_df,
                         lb, ub, 
                         N_samples,
                         alpha,
                         df,
                         q, 
                         global_optimum, 
                         min_step_coef,
                         max_step_coef,
                         break_at_bounds,
                         # nlopt settings for the main optimization problem
                         opts,
                         # nlopt settings for the estimation of theta_step
                         opts_theta_step,
                         create_txt)
    })
  }
  # Export to the cluster any function or parameter that the obj_f needs to work.
  clusterExport(cl=cluster, c(names(exported_to_cluster),"obj_f", "profile_likelihood"))
  output <- parLapply(cluster, X, parallel_func)
  # Terminate the cluster.
  stopCluster(cluster)
  total.duration <- Sys.time() - start.time
  
  # Estimation of the confidence intervals
  ci_estimation <- function(pl_results, alpha=0.95, df=1, global_optimum){
    # pl_results should be a dataframe with 2 columns containing the theta values 
    # and the corresponding chi^2 values
    
    # As confidence interevals of each theta are considered the borders of the
    # following refion: 
    # {theta | chi^2(theta) - chi^2(theta_hat) < Delta_alpha} where
    Delta_alpha = qchisq(alpha, df)
    threshold = Delta_alpha + global_optimum
    # Estimate the lower confidence interval
    # Check if the threshold was exceeded at the backwards search. If true, then
    # the lower bound is the last theta under the threshold, else it is considered 
    # -Inf.
    if(pl_results[1,2] > threshold){
      # Use those 2 theta_i values that yield to the minimum distance from the threshold
      # and interpoalte between them in order to estimate the lower bound 
      # of the parameter. 
      
      # This is the last value of the parameter that yields chi^2 lower than the 
      # threshold
      theta_under_threshold <- pl_results[2,1]
      # This is the chi^2 value that this parameters has
      chi2_under_threshold <- pl_results[2,2]
      # This is the last value of the parameter that yields chi^2 grater than the 
      # threshold
      theta_over_threshold <- pl_results[1,1]
      # This is the chi^2 value that this parameters has
      chi2_over_threshold <- pl_results[1,2]
      
      slope <- (chi2_over_threshold - chi2_under_threshold)/(theta_over_threshold - theta_under_threshold)
      dy <- threshold - chi2_under_threshold
      dx <- dy/slope
      lower_bound <- theta_under_threshold + dx
      
    }else{
      lower_bound <- -Inf
    }
    
    # Estimate the upper confidence interval
    # Check if the threshold was exceeded at the forward search. If true, then
    # the lower bound is the last theta under the threshold, else it is considered 
    # +Inf.
    if(pl_results[dim(pl_results)[1],2] > threshold){
      # Use those 2 theta_i values that yield to the minimum distance from the threshold
      # and interpoalte between them in order to estimate the lower bound 
      # of the parameter. 
      
      # This is the last value of the parameter that yields chi^2 lower than the 
      # threshold
      theta_under_threshold <- pl_results[(dim(pl_results)[1]-1),1]
      # This is the chi^2 value that this parameters has
      chi2_under_threshold <- pl_results[(dim(pl_results)[1]-1),2]
      # This is the last value of the parameter that yields chi^2 grater than the 
      # threshold
      theta_over_threshold <- pl_results[dim(pl_results)[1],1]
      # This is the chi^2 value that this parameters has
      chi2_over_threshold <- pl_results[dim(pl_results)[1],2]
      
      slope <- (chi2_over_threshold - chi2_under_threshold)/(theta_over_threshold - theta_under_threshold)
      dy <- threshold - chi2_under_threshold
      dx <- dy/slope
      upper_bound <- theta_under_threshold + dx
    }else{
      upper_bound <- +Inf
    }
    
    return(list("Lower_bound" = lower_bound,
                "Upper_bound" = upper_bound)) 
  }
  
  # Collect all the results of interest and present them in a dataframe.
  results_df <- data.frame(matrix(NA, ncol = 6, nrow = length(thetas)))
  rownames(results_df) <- thetas_names
  colnames(results_df) <- c("Non-Identifiability" , "Optimal", "Lower_Bound", "Upper_Bound", "Exit_code_B", "Exit_code_F")
  for (i in 1:length(thetas)) {
    pl_results <- output[[i]]$'plik'
    confidence_intervals <- ci_estimation(pl_results,  alpha = alpha, df=df, global_optimum)
    results_df$Lower_Bound[i] <- confidence_intervals$Lower_bound
    results_df$Upper_Bound[i] <- confidence_intervals$Upper_bound
    results_df$Optimal[i] <- thetas[[i]]
    results_df$Exit_code_B[i] <- output[[i]]$b_exit
    results_df$Exit_code_F[i] <- output[[i]]$f_exit
    if(confidence_intervals$Lower_bound == -Inf & confidence_intervals$Upper_bound == Inf){
      results_df$'Non-Identifiability'[i] <- "Structural"
    }else if(confidence_intervals$Lower_bound == -Inf & confidence_intervals$Upper_bound != Inf |
             confidence_intervals$Lower_bound != -Inf & confidence_intervals$Upper_bound == Inf){
      results_df$'Non-Identifiability'[i] <- "Practical"
    }else{
      results_df$'Non-Identifiability'[i] <- "Identifiable"
    }
  }
  
  return(list("Likelihood_profiles" = output, "Identiafiability_Analysis" = results_df, "Total_duration" = total.duration))
}

################################################################################

#=======================#
#  PBK Model Functions  #
#=======================#

create.params <- function(user.input){
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
  phys.params <- openxlsx::read.xlsx('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 1)
  
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
  phys.params <- openxlsx::read.xlsx('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Rainbow trout Physiological parameters/Rainbow trout Physiological parameters.xlsx', sheet = 2)
  
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
  #Free = 3.2e-02
  
  # Reabsorption coefficients from bile to intestine
  # estimated by Cao et al., 2022
  # K_urine = Cl_urine/f_reab_urine estimated by Ng et al., 2013 (unitless)
  if(user.input$substance=='PFOA'){
    a <- 0.138
    f_reab_hep <- 0.30
    K_urine <- 2.08 
    Free <- 0.385
  }else if(user.input$substance=='PFNA'){
    a <- 0.522
    f_reab_hep <- 0.34
    K_urine <- 1.35
    Free <- 0.622
  }else if(user.input$substance=='PFBS'){
    a <- 0.3
    f_reab_hep <- 0.23
    K_urine <- 10.41
    Free <- 0.1 # assumed
  }else if(user.input$substance=='PFHxS'){
    a <- 0.558
    f_reab_hep <- 0.30
    K_urine <- 5.88
    Free <- 0.217
  }else if(user.input$substance=='PFOS'){
    a <- 0.721
    f_reab_hep <- 0.42
    K_urine <- 1.35
    Free <- 0.819
  }
  
  # Bile flow coefficient
  Q_bile_coef <- 7.5e-05 # ml/g BW/h Grosell et al., 2000
  Q_urine_coef <- 2.755e-03 # ml/h/g of BW Urinary flow rate
  V_urine_coef <- 2.2e-03 # ml/g of BW Urine volume inside urinary bladder
  
  a_skin <- 0.9 # 90% of venous blood of skin was assumed to flow directly to kidney (Nichols et al.1996)
  a_muscle <- 0.6 # 60% of venous blood of muscle was assumed to flow directly to kidney (Nichols et al.1996)
  
  plasma <- 0.7
  
  return(list('F_card_ref'=F_card_ref, 'BW_ref'=BW_ref, 'KT'=KT, 
              'admin.dose_dietary'=admin.dose_dietary, 
              'admin.time_dietary'=admin.time_dietary,
              
              'fw_Liver'=fw_Liver, 'fw_Blood'=fw_Blood, 'fw_Skin'=fw_Skin, 
              'fw_Muscle'=fw_Muscle, 'fw_Gills'=fw_Gills, 'fw_Kidney'=fw_Kidney,
              'fw_Viscera'=fw_Viscera, 'fw_lumen'=fw_lumen, 
              
              'fb_Liver'=fb_Liver, 'fb_Skin'=fb_Skin, 'fb_Muscle'=fb_Muscle,
              'fb_Gills'=fb_Gills, 'fb_Kidney'=fb_Kidney, 'fb_Viscera'=fb_Viscera,
              
              'a_skin'=a_skin, 'a_muscle'=a_muscle,
              'Q_bile_coef'=Q_bile_coef,
              'Q_urine_coef'=Q_urine_coef, 'V_urine_coef'=V_urine_coef,
              'K_urine'=K_urine,
              'f_reab_hep'=f_reab_hep, 'plasma'=plasma, "Free"=Free, "a"=a))
}

create.inits <- function(parameters){
  with(as.list(parameters),{
    M_art<-0; M_venous<-0;
    M_gills<-0; M_lumen=0; M_lumen_2=0; M_viscera<-0; M_liver<-0; M_kidney<-0;
    M_muscle<-0; M_skin<-0; M_carcass<-0; M_storage<-0; M_urine<-0; M_feces<-0; M_input<-0 
    
    return(c('M_art'=M_art, 'M_venous'=M_venous, 'M_gills'=M_gills, 'M_lumen'=M_lumen,
             'M_lumen_2'=M_lumen_2, 'M_viscera'=M_viscera, 'M_liver'=M_liver, 'M_kidney'=M_kidney, 
             'M_muscle'=M_muscle, 'M_skin'=M_skin, 'M_carcass'=M_carcass,
             'M_storage'=M_storage,
             'M_urine'=M_urine, 'M_feces'=M_feces, 'M_input'=M_input))
  })
}

create.events <- function(parameters){
  with(as.list(parameters),{
    
    # Calculate number of administrated doses and corresponding administration time
    ldose_dietary <- length(admin.dose_dietary)
    ltimes_dietary <- length(admin.time_dietary)
    
    # If not equal, then stop 
    if (ltimes_dietary != ldose_dietary){
      stop("The times of administration should be equal in number to the doses")
    }else{
      events <- data.frame(var = c(rep(c('M_lumen', 'M_input'), ltimes_dietary)), 
                           time = sort(rep(admin.time_dietary,2)),
                           value = rep(admin.dose_dietary,each=2),
                           method = 'add')
    }
    return(list(data=events))
  })
}



ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
    # Body weight - g
    #BW <- 1000*a_bio*L^b_bio
    BW <- fish_weight(time)
    
    # Total cardiac output ml/h
    Q_total <- F_card_ref*KT*(BW/BW_ref)^(-0.1)*BW*plasma  
    
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
    
    # Calculate the regional blood flows - ml/h
    Q_liver <- fb_Liver*Q_total     # Liver blood flow - ml/h
    Q_skin <- fb_Skin*Q_total      # Skin blood flow - ml/h
    Q_muscle <- fb_Muscle*Q_total   # Muscle blood flow - ml/h
    Q_gills <- Q_total #fb_Gills*BW     # Gills blood flow - ml/h
    Q_kidney <- fb_Kidney*Q_total   # Kidney blood flow - ml/h
    Q_viscera <- fb_Viscera*Q_total # Viscera blood flow - ml/h
    Q_carcass <- Q_total - (Q_liver + Q_skin + Q_muscle + 
                              Q_kidney + Q_viscera)
    
    # Calculate the absolute bile flow rate - ml/h
    Q_bile <- Q_bile_coef*BW
    # Calculate Urinary flow rate - ml/h
    Q_urine <- Q_urine_coef*BW
    
    # Calculate urine volume  - ml 
    v_urine <- V_urine_coef*BW
    
    # Calculate f_reab_urine based on Cl_urine and K_urine - 1/h
    f_reab_urine <- Cl_urine/K_urine
    
    # Tissue concentrations ug PFAS/g tissue
    C_gills <- M_gills/w_gills
    C_viscera <- M_viscera/w_viscera
    C_liver <- M_liver/w_liver
    C_kidney <- M_kidney/w_kidney
    C_muscle <- M_muscle/w_muscle 
    C_skin <- M_skin/w_skin
    C_carcass <- M_carcass/w_carcass
    C_lumen <- (M_lumen+M_lumen_2)/w_lumen
    C_art <- M_art/w_art
    C_venous <- M_venous/w_venous
    C_blood <- (M_art + M_venous)/w_blood
    C_storage <- M_storage/v_urine
    
    # Arterial Blood
    dM_art <- Free*Q_gills*C_gills/P_gills - 
      (Q_viscera + Q_liver + Q_kidney +
         Q_muscle + Q_skin + Q_carcass)*Free*C_art
    
    dM_venous <- - Free*Q_total*C_venous + 
      ((Q_liver + Q_viscera)*C_liver/P_liver +
         (Q_kidney + a_muscle*Q_muscle + a_skin*Q_skin)*C_kidney/P_kidney +
         (1-a_muscle)*Q_muscle*C_muscle/P_muscle +
         (1-a_skin)*Q_skin*C_skin/P_skin + Q_carcass*C_carcass/P_carcass)*Free
    
    # Gills
    dM_gills <- Q_gills*Free*(C_venous - C_gills/P_gills)
    
    dM_input=0
    
    # Viscera lumen - Available PFAS for absorption and elimination
    dM_lumen = f_reab_hep*Q_bile*C_liver - Ku*a*M_lumen - Cl_feces*(1-a)*M_lumen 
    
    # Viscera lumen_2- Unavailable PFAS for absorption. Can be only eliminated.
    dM_lumen_2 = (1-f_reab_hep)*Q_bile*C_liver - Cl_feces*M_lumen_2 
    
    # Viscera tissue
    dM_viscera <- Q_viscera*Free*(C_art - C_viscera/P_viscera) + Ku*a*M_lumen 
    
    # Liver
    dM_Liver <- Q_liver*Free*C_art + Q_viscera*Free*C_viscera/P_viscera - 
      (Q_liver + Q_viscera)*Free*C_liver/P_liver - Q_bile*C_liver
    
    # Kidney
    dM_kidney <- Q_kidney*Free*C_art - 
      (Q_kidney + a_muscle*Q_muscle + a_skin*Q_skin)*Free*C_kidney/P_kidney + 
      a_muscle*Q_muscle*Free*C_muscle/P_muscle + 
      a_skin*Q_skin*Free*C_skin/P_skin - Cl_urine*M_kidney + f_reab_urine*M_storage
    
    # Muscle
    dM_muscle <- Q_muscle*Free*(C_art - C_muscle/P_muscle)
    
    # Skin
    dM_skin <- Q_skin*Free*(C_art - C_skin/P_skin)
    
    # Carcass 
    dM_carcass <- Q_carcass*Free*(C_art - C_carcass/P_carcass)
    
    # Urine storage
    dM_storage <- Cl_urine*M_kidney - f_reab_urine*M_storage - Q_urine*C_storage
    
    # Urine
    dM_urine <- Q_urine*C_storage
    
    # Feces
    dM_feces <- Cl_feces*((1-a)*M_lumen + M_lumen_2)
    
    Mass_balance <- M_input - (M_art + M_venous + M_gills + M_lumen + M_lumen_2 + 
                                 M_viscera + M_liver + M_kidney + M_muscle + 
                                 M_skin + M_carcass + M_storage + M_urine + M_feces)
    
    return(list(c('dM_art'=dM_art, 'dM_venous'=dM_venous, 
                  'dM_gills'=dM_gills, 'dM_lumen'=dM_lumen, 'dM_lumen_2'=dM_lumen_2,
                  'dM_viscera'=dM_viscera, 'dM_Liver'=dM_Liver, 
                  'dM_kidney'=dM_kidney, 'dM_muscle'=dM_muscle,
                  'dM_skin'=dM_skin, 'dM_carcass'=dM_carcass, 'dM_storage'=dM_storage,
                  'dM_urine'=dM_urine, 'dM_feces'=dM_feces, 'dM_input'=dM_input),
                'C_Gills'=C_gills, 'C_Viscera'=C_viscera,
                'C_Liver'=C_liver, 'C_Kidney'=C_kidney, 'C_Muscle'=C_muscle,
                'C_Skin'=C_skin, 'C_Carcass'=C_carcass, 'C_Lumen'=C_lumen,
                'C_Blood'=C_blood*plasma,
                'Mass_balance'=Mass_balance, 'BW'=BW))
  })
}

#======================#
#  Objective Function  #
#======================#
obj_f <- function(x, constant_theta, constant_theta_name, params_names, constant_params=NULL,
                          data_df, errors_df){
  
  if(!is.null(constant_theta)){
    if(length(constant_theta_name) != length(constant_theta)){
      stop("The constant_theta_name vector must be of equal length with the constant_theta vector")
    }
    for (j in 1:length(constant_theta)){
      assign(constant_theta_name[j], constant_theta[[j]])
    }  
  }
  # Assign the values of the x vector to the corresponding parameters
  if(length(x) != length(params_names)){
    stop("The params_names must be of equal length with the x vector")
  }
  for (k in 1:length(x)) {
    assign(params_names[k], x[k])
  }
  
  if(!is.null(constant_params)){
    for (k in 1:length(constant_params)) {
      assign(names(constant_params)[k], constant_params[[k]])
    }
  }
  
  params <- create.params(user.input)
  inits <- create.inits(params)
  events <- create.events(params)
  sample_time <- seq(0,56*24,2)
  # Time of measurement of selected PFAS 
  exp_time <- data_df[,1]
  
  pbk_params <- c('P_liver'=P_liver[[1]], 'P_muscle'=P_muscle[[1]], 'P_kidney'=P_kidney[[1]], 
                  'P_skin'=P_skin[[1]], 'P_gills'=P_gills[[1]], 'P_carcass'=P_carcass[[1]],
                  'P_viscera'=P_viscera[[1]], 'Cl_feces'=Cl_feces[[1]], 'Cl_urine'=Cl_urine[[1]],
                  'Ku'=Ku[[1]])
  
  solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = c(pbk_params,params),
                                      events = events,
                                      method="lsodes",rtol = 1e-03, atol = 1e-03))
  
  
  
  output_variables <- c('C_Liver', 'C_Blood', 
                        'C_Skin', 'C_Muscle', 'C_Gills',
                        'C_Kidney', 'C_Carcass')
  
  if(sum(solution$time %in% exp_time) == length(exp_time)){
    results <- solution[which(solution$time %in% exp_time), output_variables]*1000
  }else{
    stop(print("Length of predictions is not equal to the length of data"))
  }
  score <- c()
  # Estimate WSSR for every output variable with available data
  for (i in 1:length(output_variables)) {
    score[i] <- WSSR(list(data_df[1,i+1]), list(results[1,i]), list(errors_df[1,i+1]))
  }
  average_score <- mean(score)
  return(average_score)
}

profile_likelihood_plots <- function(analysis_results, thetas, global_optimum, alpha = 0.95,
                                     df = 1){
  output <- analysis_results$Likelihood_profiles
  plot_list <- list()
  for (i in 1:length(output)) {
    data_to_plot <- output[[i]]$plik
    current_param <- names(data_to_plot)[1]
    names(data_to_plot)[1] <- "Parameter"
    optimal_value <- data.frame(thetas[i], global_optimum)
    names(optimal_value) <- c("Parameter", "Likelihood")
    
    plot <- ggplot()+
      geom_hline(yintercept=global_optimum + qchisq(alpha,df), linetype="dashed", color = "red", size=1)+
      #geom_hline(yintercept=global_optimization$objective + qchisq(0.95,1), linetype="dashed", color = "green", size=1)+
      geom_hline(yintercept=global_optimum , linetype="dashed", color = "blue", size=1)+
      geom_line(data = data_to_plot,  aes(x=Parameter, y=Likelihood), color = 'black', size=2)+
      #geom_smooth(data = data_to_plot,  aes(x=Parameter, y=Likelihood), method = "loess", span = 0.5, se =0, color = 'black', size=2)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=11)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, colour="pink", size=10)+
      geom_point(data = optimal_value,  aes(x=Parameter, y=Likelihood), shape=18, size=5)+
      
      #scale_y_log10()+
      #ylim(c(6,NA))+
      
      labs(title = paste0( current_param), #"Profile Likelihood of ",
           y = expression(paste(chi^2, "(", theta, ")")) , x = current_param)+
      theme(plot.title = element_text(hjust = 0.5,size=30), 
            axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.y=element_text(size=22),
            axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
            axis.text.x=element_text(size=22),
            legend.title=element_text(hjust = 0.5,size=25), 
            legend.text=element_text(size=22),
            panel.border = element_rect(colour = "black", fill=NA, size=1.0))
    
    #print(plot)
    plot_list[[i]] <- plot
  }
  # Arrange and print the plots in a grid
  grid.arrange(grobs = plot_list, nrow = 3)
}

################################################################################
################################################################################
################################################################################
#                           END OF FUNCTIONS                                   #
################################################################################
################################################################################
################################################################################


###########################
#   GLOBAL OPTIMIZATION   #
###########################

# # initial values of parameters
# x0 <- c(1.696437, 0.116814, 0.5202888, 0.2893776, 0.2529077, 0.1154041,
#         6.845947, 0.007199194, 0.1281913, 0.6346485)
# # x0 <- c(1.6, 0.1, 0.5, 0.26, 0.2, 0.1,
# #         6.8, 0.007, 0.1, 0.6)
# # x0 <- rep(1,10)
# names(x0) <- c('P_liver', 'P_muscle', 'P_kidney', 
#                'P_skin', 'P_gills', 'P_carcass', 'P_viscera',
#                'Cl_feces', 'Cl_urine', 'Ku')
# params_names <- names(x0)
# # Total body weight of fish
# Texp <- 15 # C
# 
# # Fish feeding
# # The fish were fed once per day. The daily food intake was 2.6% of the mean body weight.
# # The amount of food was  consumed immediately and completely. The uptake phase
# # lasted 28 days. Nominal concentration of PFAS in spiked food is 500 ug/kg.
# 
# # We estimate the mean body weight of the fish from day 0 till day 28 in order to 
# # estimate the added mount of food everyday.
# 
# fish_weight <- function(time){
#   x <- c(0,28,56)*24 
#   y <- c(314, 655, 808)
#   
#   if(time <= x[1]){
#     w = y[1]
#   }else if(time >= x[3]){
#     w = y[3]
#   }else if(time >= x[1] & time < x[2]){
#     w = approx(x=x[1:2],y=y[1:2], xout = time)$y
#   }else if(time >= x[2] & time < x[3]){
#     w = approx(x=x[2:3],y=y[2:3], xout = time)$y
#   }
#   return(w)
# }
# 
# # Time points of added food
# admin.time_dietary <- seq(0,27*24,24)
# # Calculate fish weight over time (g)
# fish_weights <- unlist(lapply(admin.time_dietary, fish_weight))
# # Multiply fish_weights * g daily_food_intake/g of BW * Concentration (ug/g of food)
# admin.dose_dietary <- fish_weights*2.6/100*500/1000
# 
# user.input <- list('substance'='PFOS',
#                    'Texp'=Texp,
#                    'admin.dose_dietary'=admin.dose_dietary,
#                    'admin.time_dietary'=admin.time_dietary)
# 
# params <- create.params(user.input)
# inits <- create.inits(params)
# events <- create.events(params)
# sample_time <- seq(0,56*24,1)

# # Experimental data from Falk et al.2015
# #---------------------------------------
# # The concentrations in the data are given in ug PFAS/kg tissue units. 
# # The time is given in days and will be transformed in hours, to be compatible
# # with the model
# 
# # Directory of folder with saved data files
# data_dir <- 'C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015'
# # Load PFOA data
# #---------------
# PFOS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFOS.xlsx'))
# PFOS_data$Time <- PFOS_data$Time*24
# data_df <- PFOS_data
# 
# # Consider a Coefficient of Variation of the data points (CV = sd/mean)
# CV <- 20/100
# errors_df <- data.frame(matrix(NA, nrow = nrow(data_df), ncol = ncol(data_df)))
# for (i in 1:nrow(data_df)) {
#   for (j in 2:ncol(data_df)) {
#     set.seed(100)
#     errors_df[i,j] <- abs(rnorm(1, data_df[i,j]*CV, 1))
#   }
# }
# errors_df[,1] <- data_df[,1]
# colnames(errors_df) <- colnames(data_df) 

# lb <- rep(1e-05,length(x0))
# ub <- rep(20, length(x0))
# global_optimization <- nloptr::nloptr(x0 = x0,
#                                       eval_f = obj_f,
#                                       lb	= lb,
#                                       ub = ub,
#                                       opts = list("algorithm" = "NLOPT_LN_SBPLX", #"NLOPT_LN_NELDERMEAD" "NLOPT_LN_SBPLX" 
#                                                   "xtol_rel" = 1e-06, 
#                                                   "ftol_rel" = 1e-06,
#                                                   "ftol_abs" = 0.0,
#                                                   "xtol_abs" = 0.0 ,
#                                                   "maxeval" = 1000,
#                                                   "print_level" = 1),
#                                       constant_theta = NULL,
#                                       constant_theta_name = NULL,
#                                       params_names = params_names,
#                                       constant_params=NULL,
#                                       data_df = data_df,
#                                       errors_df = errors_df)




################################
#   Identifiability Analysis   #
################################

# Experimental data from Falk et al.2015
#---------------------------------------
# The concentrations in the data are given in ug PFAS/kg tissue units. 
# The time is given in days and will be transformed in hours, to be compatible
# with the model

# Directory of folder with saved data files
data_dir <- 'C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/PFAS_Data/Falk et al.2015'
# Load PFOA data
#---------------
PFOS_data <- openxlsx::read.xlsx(paste0(data_dir,'/','PFOS.xlsx'))
PFOS_data$Time <- PFOS_data$Time*24
data_df <- PFOS_data

# Consider a Coefficient of Variation of the data points (CV = sd/mean)
CV <- 20/100
errors_df <- data.frame(matrix(NA, nrow = nrow(data_df), ncol = ncol(data_df)))
for (i in 1:nrow(data_df)) {
  for (j in 2:ncol(data_df)) {
    set.seed(100)
    errors_df[i,j] <- abs(rnorm(1, data_df[i,j]*CV, 1))
  }
}
errors_df[,1] <- data_df[,1]
colnames(errors_df) <- colnames(data_df) 

setwd('C:/Users/vassi/Documents/GitHub/PFAS_PBK_models/PFAS Rainbow trout PBK/Identifiability_Analysis')
global_optimum <- 0.04261145 #global_optimization$objective
x_opt <- c(1.2270596458, 0.1961113100, 0.4291278784, 0.2532434467, 0.2536738621,
           0.1013013492, 6.3254528747, 0.0007452958, 0.0489045306, 0.4031006573)

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

# Time points of added food
admin.time_dietary <- seq(0,27*24,24)
# Calculate fish weight over time (g)
fish_weights <- unlist(lapply(admin.time_dietary, fish_weight))
# Multiply fish_weights * g daily_food_intake/g of BW * Concentration (ug/g of food)
admin.dose_dietary <- fish_weights*2.6/100*500/1000
user.input <- list('substance'='PFOS',
                   'Texp'=15,
                   'admin.dose_dietary'=admin.dose_dietary,
                   'admin.time_dietary'=admin.time_dietary)

thetas <- (x_opt)
thetas_names <- c('P_liver', 'P_muscle', 'P_kidney', 
                  'P_skin', 'P_gills', 'P_carcass', 'P_viscera',
                  'Cl_feces', 'Cl_urine', 'Ku')
names(thetas) <- thetas_names

# lower bounds of parameters
lb <- (rep(1e-07,length(thetas)))
# upper bounds of parameters
ub <- (rep(100, length(thetas)))

names(thetas) <- thetas_names

constant_params <- NULL
substance <- "PFOS" 
Texp <- 15
exported_to_cluster = list("fish_weight"=fish_weight,
                           "create.inits"=create.inits,
                           "create.params"=create.params,
                           "create.events"=create.events,
                           "ode.func"=ode.func,
                           "WSSR"=WSSR,
                           "substance"=substance, 
                           "Texp"=Texp,
                           "admin.dose_dietary"=admin.dose_dietary,
                           "admin.time_dietary"=admin.time_dietary,
                           "user.input"=user.input)



test <- Identifiability_analysis(obj_f = obj_f,
                                 thetas=thetas,
                                 thetas_names=thetas_names,
                                 data_df=data_df ,
                                 errors_df=errors_df,
                                 lb=lb ,
                                 ub=ub,
                                 N_samples = 3,
                                 alpha = 0.95 ,
                                 df = 1,
                                 q = 0.5,
                                 global_optimum = global_optimum ,
                                 min_step_coef = 0.01 ,
                                 max_step_coef = 0.5,
                                 N_cores = 8,
                                 constant_params = constant_params,
                                 exported_to_cluster = exported_to_cluster,
                                 break_at_bounds = TRUE,
                                 # nlopt settings for the main optimization problem
                                 opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                             "xtol_rel" = 1e-03, 
                                             "ftol_rel" = 1e-03,
                                             "ftol_abs" = 0.0,
                                             "xtol_abs" = 0.0 ,
                                             "maxeval" = 75,
                                             "print_level" = 0),
                                 # nlopt settings for the estimation of theta_step
                                 opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                        "xtol_rel" = 1e-03,
                                                        "ftol_rel" = 1e-03,
                                                        "ftol_abs" = 0.0,
                                                        "xtol_abs" = 0.0 ,
                                                        "maxeval" = 40,
                                                        "print_level" = 0),
                                 create_txt = TRUE)




profile_likelihood_plots(analysis_results = test, thetas, global_optimum, alpha = 0.95,
                         df = 1)



#### plot
# x_opt <- global_optimization$solution
# params <- c(create.params(user.input), x_opt)
# inits <- create.inits(params)
# events <- create.events(params)
# 
# solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func, y = inits, parms = c(x_opt, params),
#                                     events = events,
#                                     method="lsodes",rtol = 1e-05, atol = 1e-05))
# 
# # Keep the predictions only for the time points at which there are available data
# predictions_df <- solution[,c('time' , 'C_Liver', 'C_Blood',
#                               'C_Skin', 'C_Muscle', 'C_Gills',
#                               'C_Kidney', 'C_Carcass')]
# 
# predictions_df[,-1] <- predictions_df[,-1]*1000
# exp_data <- data_df
# compartments <- colnames(exp_data)[2:8]
# color_codes <- scales::hue_pal()(length(compartments))
# names(color_codes) <-  colnames(exp_data)[2:8]
# 
# plot <- ggplot()+
#   geom_line(data = predictions_df, aes(x = time, y = C_Liver, color='Liver'), size=1.3)+
#   geom_line(data = predictions_df, aes(x = time, y = C_Blood, color='Blood'), size=1.3)+
#   geom_line(data = predictions_df, aes(x = time, y = C_Skin, color='Skin'), size=1.3)+
#   geom_line(data = predictions_df, aes(x = time, y = C_Muscle, color='Muscle'), size=1.3)+
#   geom_line(data = predictions_df, aes(x = time, y = C_Gills, color='Gills'), size=1.3)+
#   geom_line(data = predictions_df, aes(x = time, y = C_Kidney, color='Kidney'), size=1.3)+
#   geom_line(data = predictions_df, aes(x = time, y = C_Carcass, color='Carcass'), size=1.3)+
#   
#   geom_point(data = exp_data, aes(x = Time, y = Liver, color = 'Liver'), size=5)+
#   geom_point(data = exp_data, aes(x = Time, y = Blood, color = 'Blood'), size=5)+
#   geom_point(data = exp_data, aes(x = Time, y = Skin, color = 'Skin'), size=5)+
#   geom_point(data = exp_data, aes(x = Time, y = Muscle, color = 'Muscle'), size=5)+
#   geom_point(data = exp_data, aes(x = Time, y = Gills, color = 'Gills'), size=5)+
#   geom_point(data = exp_data, aes(x = Time, y = Kidney, color = 'Kidney'), size=5)+
#   geom_point(data = exp_data, aes(x = Time, y = Carcass, color = 'Carcass'), size=5)+
#   #scale_y_log10(limits = c(1, 600))+
#   #ylim(c(1, 600))+
#   geom_vline(xintercept = 28*24, size=1.0)+ # end of uptake
#   
#   
#   
#   labs(title = paste0("Tissues Predicted vs Observed Concentrations of ", "PFOS"),
#        y = 'Concentration (ug/kg)' , x = "Time (hours)")+
#   theme(plot.title = element_text(hjust = 0.5,size=30), 
#         axis.title.y =element_text(hjust = 0.5,size=25,face="bold"),
#         axis.text.y=element_text(size=22),
#         axis.title.x =element_text(hjust = 0.5,size=25,face="bold"),
#         axis.text.x=element_text(size=22),
#         legend.title=element_text(hjust = 0.5,size=25), 
#         legend.text=element_text(size=22),
#         panel.border = element_rect(colour = "black", fill=NA, size=1.0)) + 
#   
#   scale_color_manual("Tissues", values=color_codes)+
#   theme(legend.key.size = unit(1.5, 'cm'),  
#         legend.title = element_text(size=14),
#         legend.text = element_text(size=14),
#         axis.text = element_text(size = 14))
#   plot