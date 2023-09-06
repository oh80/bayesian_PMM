gibbs_sampler <- function(data, j, sample_size, beta_0, Sigma_0, sigma){
  
  complete_data <- data |> dplyr::filter(R == 0)
  missing_data <- data |> dplyr::filter(R == j) |> dplyr::select(-R)
  n <- length(missing_data$Y)
  missing_col <- paste0("x",j)
  X_obs <- missing_data |> dplyr::select(-Y)
  X_obs_matrix <- X_obs |> dplyr::select(-all_of(missing_col)) |> as.matrix()
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
    X <- X_obs |> dplyr::mutate( "{missing_col}":= missing_value_sample[,j]) |> as.matrix()
    X_square <- t(X) %*% X

    #neta posterior params
    posterior_precision <- prior_precision + X_square/sigma
    posterior_mean <- solve(posterior_precision) %*% (prior_mean + (t(X) %*% Y)/sigma)

    #sampling
    beta_sample[1:10, i] <- mvtnorm::rmvnorm(n = 1, mean =  posterior_mean, sigma = solve(posterior_precision)) |> t()

    #missing value posterior params
    posterior_var <- 1/prior_var +  beta_sample[j, i]/sigma
    posterior_mean_X <- rep(posterior_var,n) * (prior_mean_X/prior_var + (rep(beta_sample[j,i],n) * as.vector(Y) + rep(beta_sample[j,i],n) * X_obs_matrix %*% beta_sample[-j,i])/sigma)

    #sampling
    missing_value_sample[1:n, i+1] <- mvtnorm::rmvnorm(n = 1, mean =  posterior_mean_X, sigma = posterior_var) |> t()
  }
  return(missing_value_sample)
}



data <- readRDS("Liner_Regression/02_build/data/MCAR_0.1.obj")
beta_0 <- get_OLS_extimater(data)
sigma <- var(data$Y)
Sigma_0 <- get_Sigma0(data, g = 10, sigma)

missing_value_sample<- gibbs_sampler(data, j = 1, sample_size = 10, beta_0, Sigma_0, sigma)
