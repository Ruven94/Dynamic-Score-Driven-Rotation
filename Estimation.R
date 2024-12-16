dynamic_score_driven_rot_estimation <- function(et,
                                                df1_range = seq(3, 12),
                                                df2_range = seq(3, 12)){
  ###############################################################################
  # Estimation Setup
  ###############################################################################
  library(parallel)
  
  # Detect number of cores for parallelization
  n_cores <- detectCores()
  
  # Initialize a result container
  result_matrix <- matrix(0, nrow = 1, ncol = 14)
  result_df <- data.frame(result_matrix)
  colnames(result_df) <- c("df1", "df2", "omega", "beta", "kappa", "likelihood", 
                           "score", "df1", "df2", "omega", "beta", "kappa", 
                           "likelihood", "score")
  
  # Define parameter grids
  omega_vec <- seq(0, 0.9 * pi, 0.1 * pi)
  beta_vec <- seq(0.70, 1 - 0.025, 0.025) 
  kappa_vec <- seq(0, 0.3, 0.025)
  
  ###############################################################################
  # run_estimation Function
  #
  # This function:
  # 1. Given df1, df2, and data 'et', it:
  #    - Evaluates lambda_est_max over a grid of (omega, beta, kappa) to find 
  #      a good initialization.
  #    - Sorts these initial evaluations by their likelihood.
  #    - Uses the top results as starting points for the maxBHHH optimization.
  # 2. Returns the final estimates and related values.
  #
  # Arguments:
  #   df1, df2 : degrees of freedom parameters
  #   et       : data matrix (innovations)
  #   length   : number of top initializations from the grid search to refine 
  #              (default: 2)
  #
  # Returns:
  #   A vector of results including parameter estimates, likelihood and score.
  ###############################################################################
  run_estimation <- function(df1, df2, et, length = 2) {
    df <- c(df1, df2)
    
    # Compute all initial likelihood values over parameter grid
    N <- length(kappa_vec) * length(omega_vec) * length(beta_vec)
    init_matrix <- matrix(0, nrow = N, ncol = 4)
    
    i <- 1
    for (omega in omega_vec) {
      for (beta in beta_vec) {
        for (kappa in kappa_vec) {
          theta_test <- c(omega, beta, kappa)
          init_matrix[i, 1:3] <- theta_test
          init_matrix[i, 4] <- lambda_est_max(theta_test, et, df)
          i <- i + 1
        }
      }
    }
    
    # Sort initial results by likelihood (descending)
    sorted_init <- init_matrix[order(-init_matrix[,4]),]
    
    # First candidate from the grid search
    theta_afterinit <- sorted_init[1, 1:3]
    
    # Run the optimization using maxBHHH with the first starting point
    result <- maxBHHH(fn = lambda_est_max, grad = grad_function, 
                      start = theta_afterinit, et = et, df = df,
                      control = list(tol = 1e-15, reltol = 1e-30, gradtol = 1e-15, 
                                     steptol = 1e-15, iterlim = 1000))
    
    likelihood <- lambda_est_max(theta = result$estimate, et = et, df = df)
    score_value <- colSums(grad_function(result$estimate, et = et, df = df))
    
    # If we only need the best result
    if (length == 1) {
      return(c(df1, df2, result$estimate[1], result$estimate[2], result$estimate[3], 
               likelihood, score_value))
    }
    
    # If we want two results, then also run optimization using the second best start
    if (length == 2) {
      theta_afterinit2 <- sorted_init[2, 1:3]
      
      result2 <- maxBHHH(fn = lambda_est_max, grad = grad_function, 
                         start = theta_afterinit2, et = et, df = df,
                         control = list(tol = 1e-15, reltol = 1e-30, gradtol = 1e-15, 
                                        steptol = 1e-15, iterlim = 1000))
      likelihood_2 <- lambda_est_max(theta = result2$estimate, et = et, df = df)
      score_value_2 <- colSums(grad_function(result2$estimate, et = et, df = df))
      
      return(c(df1, df2, 
               result$estimate[1], result$estimate[2], result$estimate[3], likelihood, score_value, 
               df1, df2, 
               result2$estimate[1], result2$estimate[2], result2$estimate[3], likelihood_2, score_value_2))
    }
  }
  
  ###############################################################################
  # Parallelized Estimation Over Grids of df1 and df2
  ###############################################################################
  
  # Define combinations of df1 and df2
  df_combinations <- expand.grid(df1 = df1_range, df2 = df2_range)
  
  # Start time for tracking performance
  start <- Sys.time()
  
  # Parallel run over all df_combinations
  results <- mclapply(1:nrow(df_combinations), function(i) {
    run_estimation(df_combinations$df1[i], df_combinations$df2[i], et = et, length = 2)
  }, mc.cores = n_cores)
  
  # Combine results into a data frame
  result_df <- do.call(rbind, results)
  
  # Assign column names appropriately
  colnames(result_df) <- c("df1", "df2", "omega", "beta", "kappa", "likelihood", 
                           "score omega","score beta", "score gamma", 
                           "df1", "df2", "omega", "beta", "kappa", "likelihood", 
                           "score omega","score beta", "score gamma")
  
  end <- Sys.time()
  print(end - start)
  
  # Reformat results by stacking rows appropriately (adjusting indices as needed)
  result_df <- rbind(result_df[,1:9], result_df[,10:18])
  
  # Sort results by likelihood
  sorted_result <- result_df[order(-result_df[, 6]), 1:6]
  return(sorted_result)
}
