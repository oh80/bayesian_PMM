#### R2_prior_parameter####

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

#### R1_prior_parameter####

get_R1_prior_parameter <- function(R2_sample, R3_sample){
  #beta2
  n <- length(R2_sample$beta2)
  mu_20 <- R2_sample$beta2[[1]]/n
  for (i in 2:n){
    add <- R2_sample$beta2[[i]] / n
    mu_20 <- mu_20 + add
  }
  
  Sigma_20 <- t(R2_sample$beta2[[1]] - mu_20)  %*% (R2_sample$beta2[[1]] - mu_20) / n
  for (i in 2:n){
    add <- t(R2_sample$beta2[[i]] - mu_20)  %*% (R2_sample$beta2[[i]] - mu_20) / n
    Sigma_20 <- Sigma_20 + add
  }
  
  #beta3
  n <- length(R3_sample$beta3)
  mu_30 <- R3_sample$beta3[[1]]/n
  for (i in 2:n){
    add <- R3_sample$beta3[[i]] / n
    mu_30 <- mu_30 + add
  }
  
  Sigma_30 <- t(R3_sample$beta3[[1]] - mu_30)  %*% (R3_sample$beta3[[1]] - mu_30) / n
  for (i in 2:n){
    add <- t(R3_sample$beta3[[i]] - mu_30)  %*% (R3_sample$beta3[[i]] - mu_30) / n
    Sigma_30 <- Sigma_30 + add
  }
  
  #H2
  H2_sample <- c(1/R2_sample$H2[[1]][1,1])
  for (i in 2:n){
    add <- c(1/R2_sample$H2[[i]][1,1])
    H2_sample <- c(H2_sample, add)
  }
  H2_params <- MASS::fitdistr(H2_sample, "gamma")
  v_H20 <- H2_params$estimate[1]
  sigma_H20 <- H2_params$estimate[2]
  
  #H3
  H3_sample <- c(1/R3_sample$H3[[1]][1,1])
  for (i in 2:n){
    add <- c(1/R3_sample$H3[[i]][1,1])
    H3_sample <- c(H3_sample, add)
  }
  H3_params <- MASS::fitdistr(H3_sample, "gamma")
  v_H30 <- H3_params$estimate[1]
  sigma_H30 <- H3_params$estimate[2]
  
  #other params
  #beta1
  Sigma_10 <- 100*solve(monomvn::rwish(7, diag(7)))
  mu_10 <- rep(0,7)
  
  #T
  sigma_t0 <- 100
  
  #G
  v_G0 <- 7+2
  S_G0 <- 0.1*diag(7)
  
  #H1 
  v_H0 <- 7
  sigma_H0 <- 100
  
  output <- list(mu_10, Sigma_10, t(mu_20), Sigma_20, t(mu_30), Sigma_30, sigma_t0,
                 v_G0, S_G0, v_H0, sigma_H0, v_H20, sigma_H20, v_H30, sigma_H30)
  return(output)
}

#### R0_prior_parameter####

get_R0_prior_parameter <- function(R1_sample, R2_sample, R3_sample){
  #beta1
  n <- length(R1_sample$beta1)
  mu_10 <- R1_sample$beta1[[1]]/n
  for (i in 2:n){
    add <- R1_sample$beta1[[i]] / n
    mu_10 <- mu_10 + add
  }
  
  Sigma_10 <- t(R1_sample$beta1[[1]] - mu_10)  %*% (R1_sample$beta1[[1]] - mu_10) / n
  for (i in 2:n){
    add <- t(R1_sample$beta1[[i]] - mu_10)  %*% (R1_sample$beta1[[i]] - mu_10) / n
    Sigma_10 <- Sigma_10 + add
  }
  
  #beta2
  n <- length(R2_sample$beta2)
  mu_20 <- R2_sample$beta2[[1]]/n
  for (i in 2:n){
    add <- R2_sample$beta2[[i]] / n
    mu_20 <- mu_20 + add
  }
  Sigma_20 <- t(R2_sample$beta2[[1]] - mu_20)  %*% (R2_sample$beta2[[1]] - mu_20) / n
  for (i in 2:n){
    add <- t(R2_sample$beta2[[i]] - mu_20)  %*% (R2_sample$beta2[[i]] - mu_20) / n
    Sigma_20 <- Sigma_20 + add
  }
  
  #beta3
  n <- length(R3_sample$beta3)
  mu_30 <- R3_sample$beta3[[1]]/n
  for (i in 2:n){
    add <- R3_sample$beta3[[i]] / n
    mu_30 <- mu_30 + add
  }
  
  Sigma_30 <- t(R3_sample$beta3[[1]] - mu_30)  %*% (R3_sample$beta3[[1]] - mu_30) / n
  for (i in 2:n){
    add <- t(R3_sample$beta3[[i]] - mu_30)  %*% (R3_sample$beta3[[i]] - mu_30) / n
    Sigma_30 <- Sigma_30 + add
  }
  
  #H1
  H1_sample <- c(1/R1_sample$H1[[1]][1,1])
  for (i in 2:n){
    add <- c(1/R1_sample$H1[[i]][1,1])
    H1_sample <- c(H1_sample, add)
  }
  H1_params <- MASS::fitdistr(H1_sample, "gamma")
  v_H10 <- H1_params$estimate[1]
  sigma_H10 <- H1_params$estimate[2]
  
  #H2
  H2_sample <- c(1/R2_sample$H2[[1]][1,1])
  for (i in 2:n){
    add <- c(1/R2_sample$H2[[i]][1,1])
    H2_sample <- c(H2_sample, add)
  }
  H2_params <- MASS::fitdistr(H2_sample, "gamma")
  v_H20 <- H2_params$estimate[1]
  sigma_H20 <- H2_params$estimate[2]
  
  #H3
  H3_sample <- c(1/R3_sample$H3[[1]][1,1])
  for (i in 2:n){
    add <- c(1/R3_sample$H3[[i]][1,1])
    H3_sample <- c(H3_sample, add)
  }
  H3_params <- MASS::fitdistr(H3_sample, "gamma")
  v_H30 <- H3_params$estimate[1]
  sigma_H30 <- H3_params$estimate[2]
  
  #other params
  #T
  sigma_t0 <- 100
  
  #G
  v_G0 <- 7+2
  S_G0 <- 0.1*diag(7)
  
  output <- list(t(mu_10), Sigma_10, t(mu_20), Sigma_20, t(mu_30), Sigma_30, sigma_t0,
                 v_G0, S_G0, v_H10, sigma_H10, v_H20, sigma_H20, v_H30, sigma_H30)
  return(output)
}
