gibbs_sampler <- function(data, j, sample_size, beta_0, Sigma_0, sigma){
  
  complete_data <- data |> dplyr::filter(R == 0)
  missing_data <- data |> dplyr::filter(R == j) |> dplyr::select(-R)
  n <- length(missing_data$Y)
  missing_col <- paste0("x",j)
  X <- missing_data |> dplyr::select(-Y)
  X_obs_matrix <- X |> dplyr::select(-all_of(missing_col)) |> as.matrix()
  Y <- missing_data |> dplyr::select(Y)|> as.matrix()
  
  #initial value set
  initial_value <- mean(complete_data[,j])
  
  missing_value_sample <- data.frame()
  missing_value_sample[1:n,1] <- rep(initial_value, n)
  
  beta_sample <- data.frame()
  
  #beta prior params
  prior_precision <- solve(Sigma_0)
  prior_mean <- prior_precision  %*% beta_0
  
  #missing_value prior params
  prior_var <- var(complete_data[,j])
  prior_mean_X <- rep(mean(complete_data[,j]), n)
  
  #gibbs sammpling
  for (i in 1:sample_size) {
    
    X[, j] <- missing_value_sample[1:n,i] 
    X <- X |> as.matrix()
    X_square <- t(X) %*% X

    #beta posterior params
    posterior_precision <- prior_precision + X_square/sigma
    posterior_mean <- solve(posterior_precision) %*% (prior_mean + (t(X) %*% Y)/sigma)

    #sampling
    beta_sample[1:10, i] <- mvtnorm::rmvnorm(n = 1, mean =  posterior_mean, sigma = solve(posterior_precision)) |> t()

    #missing value posterior params
    posterior_var <- 1/(1/prior_var +  beta_sample[j, i]/sigma)
    beta_obs <- beta_sample[-j,i] |> as.matrix()
    beta_mis <- beta_sample[j,i]
    coef_first_order_term <- prior_mean_X/prior_var +  rep(beta_mis/sigma, n) * (Y + X_obs_matrix%*%beta_obs)
    posterior_mean_X <- rep(posterior_var, n) * coef_first_order_term
    
    #sampling
    missing_value_sample[1:n, i+1]<- mvtnorm::rmvnorm(n = 1, mean =  posterior_mean_X, sigma = rep(posterior_var,n) *diag(n)) |> t() 
  }
  output <- list()
  output$missing_value_sample <- missing_value_sample
  output$beta_sample <- beta_sample
  return(output)
}




