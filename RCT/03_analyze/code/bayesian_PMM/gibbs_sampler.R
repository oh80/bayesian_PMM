R3_gibbs_sampler <- function(data, sample_size){
  set.seed(428)
  
  data <- data |> dplyr::filter(R == 3)
  n <- length(data$X1)
  X <- data |> dplyr::select(Treatment, X1:X5) |> dplyr::mutate("constant" = rep(1, n)) |> 
    dplyr::select(constant, Treatment, X1:X5) |> as.matrix()
  Y1 <- data |> dplyr::select(Y1)|> as.matrix()
  Y2 <- data |> dplyr::select(Y2)|> as.matrix()
  Y3 <- data |> dplyr::select(Y3)|> as.matrix()
  
  #initial values set
  beta2 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  beta3 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  H1 <- 10*diag(n) |> as.matrix()
  H2 <- 10*diag(n) |> as.matrix()
  H3 <- 10*diag(n) |> as.matrix()
  t <- diag(7) |> as.matrix()
  G <- solve(monomvn::rwish(7, diag(7))) |> as.matrix()
  
  initial_value <- list(beta2, beta3, H1, H2, H3, t, G)
  
  #prior params set
  #beta1
  Sigma_10 <- 100*solve(monomvn::rwish(7, diag(7)))
  
  #beta2
  Sigma_20 <- 100*solve(monomvn::rwish(7, diag(7)))
  
  #beta3
  Sigma_30 <- 100*solve(monomvn::rwish(7, diag(7)))
  
  #T
  sigma_t0 <- 100
  
  #G
  v_G0 <- 7
  S_G0 <- 100*solve(monomvn::rwish(7, diag(7)))
  
  #H 
  v_H0 <- 7
  sigma_H0 <- 100
  
  prior_params <- list(Sigma_10, Sigma_20, Sigma_30, sigma_t0,  v_G0, S_G0, v_H0, sigma_H0)
  
  #gibbs sample
  output <- parmater_sampler(X, Y1, Y2, Y3, initial_value, prior_params, sample_size)
  return(output)
}


parmater_sampler <- function(X, Y1, Y2, Y3, initial_value, prior_params, sample_size){
  beta1_sample <- list()
  beta2_sample <- list() |> append(initial_value[1])
  beta3_sample <- list() |> append(initial_value[2])
  H1_sample <- list() |> append(initial_value[3])
  H2_sample <- list() |> append(initial_value[4])
  H3_sample <- list() |> append(initial_value[5])
  T_sample <- list() |> append(initial_value[6])
  G_sample <- list() |> append(initial_value[7])
  
  for(i in 1:sample_size){
    #beta1 
    Sigma_1 <- solve(t(X) %*% solve(H1_sample[[i]]) %*% X + solve(prior_params[[1]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_1 <- Sigma_1 %*% (t(X) %*% solve(H1_sample[[i]]) %*% Y1
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta2_sample[[i]]))
    
    beta1_sample[[i]] <- mvtnorm::rmvnorm(n = 1, mean = mu_1, sigma = Sigma_1) 
    
    #beta2
    Sigma_2 <- solve(t(X) %*% solve(H2_sample[[i]]) %*% X + solve(prior_params[[2]]) + solve(G_sample[[i]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_2 <- Sigma_2 %*% (t(X) %*% solve(H2_sample[[i]]) %*% Y2 +
                           solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta1_sample[[i]])
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta3_sample[[i]]))
    
    beta2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_2, sigma = Sigma_2) 
    
    #beta3
    Sigma_3 <- solve(t(X) %*% solve(H3_sample[[i]]) %*% X + solve(prior_params[[3]]) + solve(G_sample[[i]]))
    mu_3 <- Sigma_3 %*% (t(X) %*% solve(H3_sample[[i]]) %*% Y3 +
                           solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta2_sample[[i+1]]))
    
    beta3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_3, sigma = Sigma_3)
    
    #T
    t_sample <-  matrix(0, 7, 7)
    for (t in 1:7) {
      sigma_t <- 1/(1/prior_params[[4]] + 
        ((beta1_sample[[i]][t]))^2/G_sample[[i]][t,t] +
        (beta2_sample[[i+1]][t])^2/G_sample[[i]][t,t])
      mu_t <- sigma_t * (beta1_sample[[i]][t]*beta2_sample[[i+1]][t]/G_sample[[i]][t,t] +
                           beta2_sample[[i+1]][t]*beta3_sample[[i+1]][t]/G_sample[[i]][t,t])
      t_sample[t,t] <- rnorm(n = 1, mu_t, sigma_t) |> as.matrix()
    }
    T_sample[[i+1]] <- t_sample
    
    #G
    v_G <- prior_params[[5]] + 20
    S_G <- prior_params[[6]] + (t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]]))%*%
      t(t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]])) + 
      (t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]]))%*%
      t(t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]])) 
    
    G_sample[[i+1]] <- solve(monomvn::rwish(v_G, S_G)) |> as.matrix()
    
    #H1
    v_H1 <- (prior_params[[7]] + length(Y1))/2
    sigma_H1 <- prior_params[[7]] * prior_params[[8]] +
      (Y1- X %*% t(beta1_sample[[i]])) %*% t(Y1- X %*% t(beta1_sample[[i]]))
    sigma_1 <- 1/rgamma(n = 1 ,v_H1, sigma_H1)
    H1_sample[[i+1]] <- sigma_1 * diag(length(Y1)) 
    #H2
    v_H2 <- (prior_params[[7]] + length(Y2))/2
    sigma_H2 <- prior_params[[7]] * prior_params[[8]] +
      (Y2- X %*% t(beta2_sample[[i+1]])) %*% t(Y2- X %*% t(beta2_sample[[i+1]]))
    sigma_2 <- 1/rgamma(n = 1 ,v_H2, sigma_H2)
    H2_sample[[i+1]] <- sigma_2 * diag(length(Y2))
    
    #H3
    v_H3 <- (prior_params[[7]] + length(Y3))/2
    sigma_H3 <- prior_params[[7]] * prior_params[[8]] +
      (Y3- X %*% t(beta3_sample[[i+1]])) %*% t(Y3- X %*% t(beta3_sample[[i+1]]))
    sigma_3 <- 1/rgamma(n = 1 ,v_H3, sigma_H3)
    H3_sample[[i+1]] <- sigma_3 * diag(length(Y3))
  }
  output <- list()
  output$beta1 <- beta1_sample
  output$beta2 <- beta2_sample
  output$beta3 <- beta3_sample
  output$t <- T_sample
  output$G <- G_sample
  output$H1 <- H1_sample
  output$H2 <- H2_sample
  output$H3 <- H3_sample
  
  return(output)
}




