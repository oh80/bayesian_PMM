get_estimated_values <- function(missing_rate, missing_type){
  # read sample and data
  file_name <- paste0(missing_type , "_", missing_rate, ".obj")
  path <- here::here("Liner_regression", "03_analyze","output", "bayesian_PMM", file_name)
  sample <- readRDS(path)
  
  path <- here::here("Liner_regression", "02_build","data", file_name)
  data <- readRDS(path)
  
  code_path <- here::here("Liner_Regression","03_analyze","code","Bayesian_PMM","compute_complete_data_params.R")
  source(code_path)
  
  #get estimated values
  weight <- get_weight(data)
  burn_in <- 100
  beta_mean <- get_sample_mean(sample, burn_in)
  beta_0 <- get_OLS_extimater(data)
  
  estimated_values <- compute_weighted_mean(weight, beta_mean, beta_0)
  
  return(estimated_values)
}


get_weight <- function(data){
  N <- length(data$Y)
  output <- (data |> dplyr::count(R))/N
  return(output)
}


get_sample_mean <- function(sample, burn_in){
  col_names <- paste0("beta", seq(1,10,by=1))
  mean_matrix <- matrix(ncol = 10, nrow = 10)
  for (j in 1:10) {
    sample_j <- sample[[j]]
    beta <- sample_j$beta_sample |> t()
    N <- length(beta[,1])
    for (k in 1:10) {
      mean_matrix[j,k] <- mean(beta[burn_in:N,k])
    }
  }
  output <- mean_matrix |> as.data.frame() |> 
    magrittr::set_colnames(col_names)
  
  return(output)
}


compute_weighted_mean <- function(weight, beta_mean, beta_0){
  weight <- weight$n
  col_names <- paste0("beta", seq(1,10,by=1))
  mean_matrix <- matrix(ncol = 10, nrow = 1)
  for (j in 1:10) {
    mean_matrix[1,j] <- sum(beta_mean[,j] * weight[2:11]) +  beta_0[j] * weight[1]
  }
  output <- mean_matrix |> as.data.frame() |> 
    magrittr::set_colnames(col_names)
  return(output)
}


get_standard_error <- function(missing_rate, missing_type){
  # read sample and data
  file_name <- paste0(missing_type , "_", missing_rate, ".obj")
  path <- here::here("Liner_regression", "03_analyze","output", "bayesian_PMM", file_name)
  sample <- readRDS(path)
  
  
  path <- here::here("Liner_regression", "02_build","data", file_name)
  data <- readRDS(path)
  
  code_path <- here::here("Liner_Regression","03_analyze","code","Bayesian_PMM","compute_complete_data_params.R")
  source(code_path)
  
  #get standard error
  weight <- get_weight(data)
  beta_var <- get_sample_var(sample)
  beta_var_0 <- get_ols_var(data)
  
  standard_error <- compute_standard_error(weight, beta_var, beta_var_0)
  
  #save 
  
  return(standard_error)
}


get_sample_var <- function(sample){
  col_names <- paste0("var_",seq(1, 10, by = 1))
  se_matrix <- matrix(nrow = 10, ncol = 10)
  for (j in 1:10) {
    beta <- sample[[j]]$beta_sample |> t()
    for (k in 1:10) {
      se_matrix[j,k] <- var(beta[,k])
    }
  }
  output <- se_matrix |> as.data.frame() |> 
    magrittr::set_colnames(col_names)
  return(output)
}


compute_standard_error <- function(weight, beta_var, beta_var_0){
  weight <- weight$n
  col_names <- paste0("se ", seq(1,10,by=1))
  se_matrix <- matrix(ncol = 10, nrow = 1)
  for (j in 1:10) {
    se_matrix[1,j] <- sqrt(sum(beta_var[,j] + beta_var_0[j]) )
  }
  output <- se_matrix |> as.data.frame() |> 
    magrittr::set_colnames(col_names)
  return(output)
}



get_estimated_values(0.1, "MCAR")
