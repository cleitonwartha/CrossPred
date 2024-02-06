#' Generalized bootstrapper
#' 
#' @description Computes bootstrap resamples of a statistic
#' 
#' @param x A numeric vector.
#' @param y NULL (default) or a numeric vector for statistics that require a second
#' numeric vector
#' @param fun A character vector of the desired statistic function. For instance,
#' to calculate the variance of values in the numeric vector, pass "var".
#' @param boot.reps The number of bootstrapping replications.
#' @param alpha The significance level for calculating a confidence interval.
#' 
#' @import boot
#' 
#' @export
#' 
bootstrap <- function(x, y = NULL, fun, boot.reps = 1000, alpha = 0.05) {
  
  # Error handling
  boot.reps <- as.integer(boot.reps)
  
  # Error if the function does not exist
  stopifnot(exists(fun))
  # Otherwise get the function
  fun_torun <- get(fun)
  
  # Prob must be between 0 and 1
  alpha_check <- alpha > 0 | alpha < 1
  
  if (!alpha_check)
    stop("'alpha' must be between 0 and 1.")
  
  # Combine the data
  mat <- cbind(x, y)
  # How many columns
  n_col <- ncol(mat)
  
  
  ## Define a function for bootstrapping a single vector
  boot_fun_vec <- function(data, i) {
    # Use the index i to sample the data
    data_i <- data[i,,drop = FALSE]
    # Execute the function
    fun_torun(data_i[,1])
    
  }
  
  ## Define a function for bootstrapping a two vectors
  boot_fun_mat <- function(data, i) {
    # Use the index i to sample the data
    data_i <- data[i,]
    # Execute the function
    fun_torun(data_i[,1], data_i[,2])
    
  }
  
  # Calculate the base statistic
  base_stat <- if (n_col > 1) fun_torun(mat[,1], mat[,2]) else fun_torun(mat[,1])
  
  # If the base is not NA, proceed
  if (!is.na(base_stat)) {
    
    # Perform the bootstrapping
    if (n_col > 1) {
      boot_results <- boot(data = mat, statistic = boot_fun_mat, R = boot.reps)
      
    } else {
      boot_results <- boot(data = mat, statistic = boot_fun_vec, R = boot.reps)
      
    }
    
    # Standard error
    se <- sd(boot_results$t)
    # Bias
    bias <- mean(boot_results$t) - base_stat
    
    
    # Confidence interval
    ci_upper <- quantile(boot_results$t, 1 - (alpha / 2))
    ci_lower <- quantile(boot_results$t, (alpha / 2))
    
  } else {
    
    se <- bias <- ci_lower <- ci_upper <- NA
    
  }
  
  # Assemble list and return
  data.frame(statistic = fun, base = base_stat, se = se, bias = bias,
             ci_lower = ci_lower, ci_upper = ci_upper, row.names = NULL)
}