main <- function(){
  #read data
  missing_type <- "MCAR"
  missing_rate <- 0.3
  
  file_name <- paste0(missing_type,"_", missing_rate, ".", "obj")
  path <- here::here("RCT","02_build","data", missing_type ,file_name)
  data <- readRDS(path)
  
  #read functions
  code_path <- here::here("RCT","03_analyze","code","bayesian_PMM","gibbs_sampler.R")
  source(code_path)

  #complete data giggs sampler
  R3_params_sample <- R3_gibbs_sampler(data, sample_size = 10)
  #R3_estimated_values <- get_estimated_values(R3_params_sample)
  
  #R=2 data gibbs sampler
  R2_prior_parameter <- get_R2_prior_parameter(R3_params_sample)
  R2_sample <- R2_gibbs_sampler(data, sample_size = 10, R2_prior_parameter)
  
  
  return(R2_sample)
}


get_estimated_values <- function(R3_params_sample){
  n <- length(R3_params_sample$beta1)
  beta1 <- R3_params_sample$beta1[[1]]/n
  for (i in 2:n){
    add <- R3_params_sample$beta1[[i]]/n
    beta1 <- beta1+add
  }
  
  n <- length(R3_params_sample$beta2)
  beta2 <- R3_params_sample$beta2[[1]]/n
  for (i in 2:n){
    add <- R3_params_sample$beta2[[i]]/n
    beta2 <- beta2+add
  }
  
  n <- length(R3_params_sample$beta3)
  beta3 <- R3_params_sample$beta3[[1]]/n
  for (i in 2:n){
    add <- R3_params_sample$beta3[[i]]/n
    beta3 <- beta3+add
  }
  
  output <- list(beta1, beta2, beta3)
  return(output)
}


get_R2_prior_parameter <- function(sample){
  #beta3
  n <- length(sample$beta3)
  mu_30 <- sample$beta3[[1]]/n
  for (i in 2:n){
    add <- sample$beta3[[i]] / n
    mu_30 <- mu_30 + add
  }
  
  Sigma_30 <- t(sample$beta3[[1]] - mu_30)  %*% (sample$beta3[[1]] - mu_30) / n
  for (i in 2:n){
    add <- t(sample$beta3[[i]] - mu_30)  %*% (sample$beta3[[i]] - mu_30) / n
    Sigma_30 <- Sigma_30 + add
  }
  
  #H3
  H3_sample <- c(1/sample$H3[[1]][1,1])
  for (i in 2:n){
    add <- c(1/sample$H3[[i]][1,1])
    H3_sample <- c(H3_sample, add)
  }
  H3_params <- MASS::fitdistr(H3_sample, "gamma")
  v_H30 <- H3_params$estimate[1]
  sigma_H30 <- H3_params$estimate[2]
  
  #other params
  #beta1
  Sigma_10 <- 100*solve(monomvn::rwish(7, diag(7)))
  mu_10 <- rep(0,7)
  
  #beta2
  Sigma_20 <- 100*solve(monomvn::rwish(7, diag(7)))
  mu_20 <- rep(0,7)
  
  #T
  sigma_t0 <- 100
  
  #G
  v_G0 <- 7+2
  S_G0 <- 0.1*diag(7)
  
  #H1 H2 
  v_H0 <- 7
  sigma_H0 <- 100
  
  output <- list(mu_10, Sigma_10, mu_20, Sigma_20, t(mu_30), Sigma_30, sigma_t0,
                 v_G0, S_G0, v_H0, sigma_H0, v_H0, sigma_H0, v_H30, sigma_H30)
  return(output)
}


R2_sample <- main()


