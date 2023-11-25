####R3_gibbs_sampler####
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
  mu_10 <- rep(0,7)
  
  #beta2
  Sigma_20 <- 100*solve(monomvn::rwish(7, diag(7)))
  mu_20 <- rep(0,7)
  
  #beta3
  Sigma_30 <- 100*solve(monomvn::rwish(7, diag(7)))
  mu_30 <- rep(0,7)
  
  #T
  sigma_t0 <- 100
  
  #G
  v_G0 <- 7
  S_G0 <- 0.1*diag(7)
  
  #H1 H2 H3
  v_H0 <- 7
  sigma_H0 <- 100
  
  
  prior_params <- list(mu_10, Sigma_10, mu_20, Sigma_20, mu_30, Sigma_30, sigma_t0,
                       v_G0, S_G0, v_H0, sigma_H0, v_H0, sigma_H0, v_H0, sigma_H0)
  
  #gibbs sample
  output <- R3_sampler(X, Y1, Y2, Y3, initial_value, prior_params, sample_size)
  return(output)
}


R3_sampler <- function(X, Y1, Y2, Y3, initial_value, prior_params, sample_size){
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
    Sigma_1 <- solve(t(X) %*% solve(H1_sample[[i]]) %*% X + solve(prior_params[[2]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_1 <- Sigma_1 %*% (t(X) %*% solve(H1_sample[[i]]) %*% Y1 
                         + solve(prior_params[[2]]) %*% prior_params[[1]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta2_sample[[i]]))
    
    beta1_sample[[i]] <- mvtnorm::rmvnorm(n = 1, mean = mu_1, sigma = Sigma_1) 
    
    #beta2
    Sigma_2 <- solve(t(X) %*% solve(H2_sample[[i]]) %*% X + solve(prior_params[[4]]) + solve(G_sample[[i]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_2 <- Sigma_2 %*% (t(X) %*% solve(H2_sample[[i]]) %*% Y2 +
                           solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta1_sample[[i]])
                         + solve(prior_params[[4]]) %*% prior_params[[3]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta3_sample[[i]]))
    
    beta2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_2, sigma = Sigma_2) 
    
    #beta3
    Sigma_3 <- solve(t(X) %*% solve(H3_sample[[i]]) %*% X + solve(prior_params[[6]]) + solve(G_sample[[i]]))
    mu_3 <- Sigma_3 %*% (t(X) %*% solve(H3_sample[[i]]) %*% Y3 +
                           + solve(prior_params[[6]]) %*% prior_params[[5]]
                           + solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta2_sample[[i+1]]))
    
    beta3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_3, sigma = Sigma_3)
    
    #T
    t_sample <-  matrix(0, 7, 7)
    for (t in 1:7) {
      sigma_t <- 1/(1/prior_params[[7]] + 
        ((beta1_sample[[i]][t]))^2/G_sample[[i]][t,t] +
        (beta2_sample[[i+1]][t])^2/G_sample[[i]][t,t])
      mu_t <- sigma_t * (beta1_sample[[i]][t]*beta2_sample[[i+1]][t]/G_sample[[i]][t,t] +
                           beta2_sample[[i+1]][t]*beta3_sample[[i+1]][t]/G_sample[[i]][t,t])
      t_sample[t,t] <- rnorm(n = 1, mu_t, sigma_t) |> as.matrix()
    }
    T_sample[[i+1]] <- diag(7)
    
    #G
    v_G <- prior_params[[8]] + 2
    S_G <- prior_params[[9]] + (t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]]))%*%
      t(t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]])) + 
      (t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]]))%*%
      t(t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]])) 
    
    G_sample[[i+1]] <- solve(monomvn::rwish(v_G, solve(S_G))) |> as.matrix()
    
    #H1
    v_H1 <- (prior_params[[10]] + length(Y1))
    sigma_H1 <- prior_params[[10]] * prior_params[[11]] +
      (Y1- X %*% t(beta1_sample[[i]])) %*% t(Y1- X %*% t(beta1_sample[[i]]))
    sigma_1 <- 1/rgamma(n = 1 ,v_H1/2, sigma_H1/2)
    H1_sample[[i+1]] <- sigma_1 * diag(length(Y1)) 
    #H2
    v_H2 <- (prior_params[[12]] + length(Y2))
    sigma_H2 <- prior_params[[12]] * prior_params[[13]] +
      (Y2- X %*% t(beta2_sample[[i+1]])) %*% t(Y2- X %*% t(beta2_sample[[i+1]]))
    sigma_2 <- 1/rgamma(n = 1 ,v_H2/2, sigma_H2/2)
    H2_sample[[i+1]] <- sigma_2 * diag(length(Y2))
    
    #H3
    v_H3 <- (prior_params[[14]] + length(Y3))
    sigma_H3 <- prior_params[[14]] * prior_params[[15]] +
      (Y3- X %*% t(beta3_sample[[i+1]])) %*% t(Y3- X %*% t(beta3_sample[[i+1]]))
    sigma_3 <- 1/rgamma(n = 1 ,v_H3/2, sigma_H3/2)
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

####R2_gibbs_sampler####

R2_gibbs_sampler <- function(data, sample_size, prior_params){
  set.seed(428)
  R3_data <- data|> dplyr::filter(R == 3)
  
  data <- data |> dplyr::filter(R == 2)
  n <- length(data$X1)
  X <- data |> dplyr::select(Treatment, X1:X5) |> dplyr::mutate("constant" = rep(1, n)) |> 
    dplyr::select(constant, Treatment, X1:X5) |> as.matrix()
  Y1 <- data |> dplyr::select(Y1)|> as.matrix()
  Y2 <- data |> dplyr::select(Y2)|> as.matrix()
  
  #initial values set
  beta2 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  beta3 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  H1 <- 10*diag(n) |> as.matrix()
  H2 <- 10*diag(n) |> as.matrix()
  H3 <- 10*diag(n) |> as.matrix()
  t <- diag(7) |> as.matrix()
  G <- solve(monomvn::rwish(7, diag(7))) |> as.matrix()
  Y3 <- mvtnorm::rmvnorm(n = 1, mean = X %*% t(beta3), sigma = H3) |> t()
  
  initial_value <- list(beta2, beta3, H1, H2, H3, t, G, Y3)
  
  #gibbs sample
  output <- R2_sampler(X, Y1, Y2, initial_value, prior_params, sample_size)
  return(output)
}


R2_sampler <- function(X, Y1, Y2, initial_value, prior_params, sample_size){
  beta1_sample <- list()
  beta2_sample <- list() |> append(initial_value[1])
  beta3_sample <- list() |> append(initial_value[2])
  H1_sample <- list() |> append(initial_value[3])
  H2_sample <- list() |> append(initial_value[4])
  H3_sample <- list() |> append(initial_value[5])
  T_sample <- list() |> append(initial_value[6])
  G_sample <- list() |> append(initial_value[7])
  Y3_sample <- list()|> append(initial_value[8])
  
  for(i in 1:sample_size){
    #beta1 
    Sigma_1 <- solve(t(X) %*% solve(H1_sample[[i]]) %*% X + solve(prior_params[[2]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_1 <- Sigma_1 %*% (t(X) %*% solve(H1_sample[[i]]) %*% Y1 
                         + solve(prior_params[[2]]) %*% prior_params[[1]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta2_sample[[i]]))
    
    beta1_sample[[i]] <- mvtnorm::rmvnorm(n = 1, mean = mu_1, sigma = Sigma_1) 
    
    #beta2
    Sigma_2 <- solve(t(X) %*% solve(H2_sample[[i]]) %*% X + solve(prior_params[[4]]) + solve(G_sample[[i]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_2 <- Sigma_2 %*% (t(X) %*% solve(H2_sample[[i]]) %*% Y2 +
                           solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta1_sample[[i]])
                         + solve(prior_params[[4]]) %*% prior_params[[3]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta3_sample[[i]]))
    
    beta2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_2, sigma = Sigma_2) 
    
    #beta3
    Sigma_3 <- solve(t(X) %*% solve(H3_sample[[i]]) %*% X + solve(prior_params[[6]]) + solve(G_sample[[i]]))
    mu_3 <- Sigma_3 %*% (t(X) %*% solve(H3_sample[[i]]) %*% Y3_sample[[i]] +
                           + solve(prior_params[[6]]) %*% prior_params[[5]]
                         + solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta2_sample[[i+1]]))
    
    beta3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_3, sigma = Sigma_3)
    
    #T
    t_sample <-  matrix(0, 7, 7)
    for (t in 1:7) {
      sigma_t <- 1/(1/prior_params[[7]] + 
                      ((beta1_sample[[i]][t]))^2/G_sample[[i]][t,t] +
                      (beta2_sample[[i+1]][t])^2/G_sample[[i]][t,t])
      mu_t <- sigma_t * (beta1_sample[[i]][t]*beta2_sample[[i+1]][t]/G_sample[[i]][t,t] +
                           beta2_sample[[i+1]][t]*beta3_sample[[i+1]][t]/G_sample[[i]][t,t])
      t_sample[t,t] <- rnorm(n = 1, mu_t, sigma_t) |> as.matrix()
    }
    T_sample[[i+1]] <- diag(7)
    
    #G
    v_G <- prior_params[[8]] + 2
    S_G <- prior_params[[9]] + (t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]]))%*%
      t(t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]])) + 
      (t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]]))%*%
      t(t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]])) 
    
    G_sample[[i+1]] <- solve(monomvn::rwish(v_G, solve(S_G))) |> as.matrix()
    
    #H1
    v_H1 <- (prior_params[[10]] + length(Y1))
    sigma_H1 <- prior_params[[10]] * prior_params[[11]] +
      (Y1- X %*% t(beta1_sample[[i]])) %*% t(Y1- X %*% t(beta1_sample[[i]]))
    sigma_1 <- 1/rgamma(n = 1 ,v_H1/2, sigma_H1/2)
    H1_sample[[i+1]] <- sigma_1 * diag(length(Y1)) 
    #H2
    v_H2 <- (prior_params[[12]] + length(Y2))
    sigma_H2 <- prior_params[[12]] * prior_params[[13]] +
      (Y2- X %*% t(beta2_sample[[i+1]])) %*% t(Y2- X %*% t(beta2_sample[[i+1]]))
    sigma_2 <- 1/rgamma(n = 1 ,v_H2/2, sigma_H2/2)
    H2_sample[[i+1]] <- sigma_2 * diag(length(Y2))
    
    #H3
    v_H3 <- (prior_params[[14]] + length(Y3_sample[[i]]))
    sigma_H3 <- prior_params[[14]] * prior_params[[15]] +
      (Y3_sample[[i]] - X %*% t(beta3_sample[[i+1]])) %*% t(Y3_sample[[i]] - X %*% t(beta3_sample[[i+1]]))
    sigma_3 <- 1/rgamma(n = 1 ,v_H3/2, sigma_H3/2)
    H3_sample[[i+1]] <- sigma_3 * diag(length(Y3_sample[[i]]))
    
    #Y3 
    Y3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean =  X %*% t(beta3_sample[[i+1]]),
                                         sigma = H3_sample[[i+1]]) |>  t()
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
  output$Y3 <- Y3_sample
  
  return(output)
}

####R1_gibbs_sampler####

R1_gibbs_sampler <- function(data, sample_size, prior_params){
  set.seed(428)
  R3_data <- data|> dplyr::filter(R == 3)
  R2_data <- data|> dplyr::filter(R == 2)
  
  data <- data |> dplyr::filter(R == 1)
  n <- length(data$X1)
  X <- data |> dplyr::select(Treatment, X1:X5) |> dplyr::mutate("constant" = rep(1, n)) |> 
    dplyr::select(constant, Treatment, X1:X5) |> as.matrix()
  Y1 <- data |> dplyr::select(Y1)|> as.matrix()

  
  #initial values set
  beta2 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  beta3 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  H1 <- 10*diag(n) |> as.matrix()
  H2 <- 10*diag(n) |> as.matrix()
  H3 <- 10*diag(n) |> as.matrix()
  t <- diag(7) |> as.matrix()
  G <- solve(monomvn::rwish(7, diag(7))) |> as.matrix()
  Y2 <- mvtnorm::rmvnorm(n = 1, mean = X %*% t(beta2), sigma = H2) |> t()
  Y3 <- mvtnorm::rmvnorm(n = 1, mean = X %*% t(beta3), sigma = H3) |> t()

  
  initial_value <- list(beta2, beta3, H1, H2, H3, t, G, Y2, Y3)
  
  #gibbs sample
  output <- R1_sampler(X, Y1, initial_value, prior_params, sample_size)
  return(output)
}

R1_sampler <- function(X, Y1, initial_value, prior_params, sample_size){
  beta1_sample <- list()
  beta2_sample <- list() |> append(initial_value[1])
  beta3_sample <- list() |> append(initial_value[2])
  H1_sample <- list() |> append(initial_value[3])
  H2_sample <- list() |> append(initial_value[4])
  H3_sample <- list() |> append(initial_value[5])
  T_sample <- list() |> append(initial_value[6])
  G_sample <- list() |> append(initial_value[7])
  Y2_sample <- list()|> append(initial_value[8])
  Y3_sample <- list()|> append(initial_value[9])
  
  for(i in 1:sample_size){
    #beta1 
    Sigma_1 <- solve(t(X) %*% solve(H1_sample[[i]]) %*% X + solve(prior_params[[2]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_1 <- Sigma_1 %*% (t(X) %*% solve(H1_sample[[i]]) %*% Y1 
                         + solve(prior_params[[2]]) %*% prior_params[[1]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta2_sample[[i]]))
    
    beta1_sample[[i]] <- mvtnorm::rmvnorm(n = 1, mean = mu_1, sigma = Sigma_1) 
    
    #beta2
    Sigma_2 <- solve(t(X) %*% solve(H2_sample[[i]]) %*% X + solve(prior_params[[4]]) + solve(G_sample[[i]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_2 <- Sigma_2 %*% (t(X) %*% solve(H2_sample[[i]]) %*% Y2_sample[[i]] +
                           solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta1_sample[[i]])
                         + solve(prior_params[[4]]) %*% prior_params[[3]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta3_sample[[i]]))
    
    beta2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_2, sigma = Sigma_2) 
    
    #beta3
    Sigma_3 <- solve(t(X) %*% solve(H3_sample[[i]]) %*% X + solve(prior_params[[6]]) + solve(G_sample[[i]]))
    mu_3 <- Sigma_3 %*% (t(X) %*% solve(H3_sample[[i]]) %*% Y3_sample[[i]] +
                           + solve(prior_params[[6]]) %*% prior_params[[5]]
                         + solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta2_sample[[i+1]]))
    
    beta3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_3, sigma = Sigma_3)
    
    #T
    t_sample <-  matrix(0, 7, 7)
    for (t in 1:7) {
      sigma_t <- 1/(1/prior_params[[7]] + 
                      ((beta1_sample[[i]][t]))^2/G_sample[[i]][t,t] +
                      (beta2_sample[[i+1]][t])^2/G_sample[[i]][t,t])
      mu_t <- sigma_t * (beta1_sample[[i]][t]*beta2_sample[[i+1]][t]/G_sample[[i]][t,t] +
                           beta2_sample[[i+1]][t]*beta3_sample[[i+1]][t]/G_sample[[i]][t,t])
      t_sample[t,t] <- rnorm(n = 1, mu_t, sigma_t) |> as.matrix()
    }
    T_sample[[i+1]] <- diag(7)
    
    #G
    v_G <- prior_params[[8]] + 2
    S_G <- prior_params[[9]] + (t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]]))%*%
      t(t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]])) + 
      (t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]]))%*%
      t(t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]])) 
    
    G_sample[[i+1]] <- solve(monomvn::rwish(v_G, solve(S_G))) |> as.matrix()
    
    #H1
    v_H1 <- (prior_params[[10]] + length(Y1))
    sigma_H1 <- prior_params[[10]] * prior_params[[11]] +
      (Y1- X %*% t(beta1_sample[[i]])) %*% t(Y1- X %*% t(beta1_sample[[i]]))
    sigma_1 <- 1/rgamma(n = 1 ,v_H1/2, sigma_H1/2)
    H1_sample[[i+1]] <- sigma_1 * diag(length(Y1)) 
    
    #H2
    v_H2 <- (prior_params[[12]] + length(Y2_sample[[i]]))
    sigma_H2 <- prior_params[[12]] * prior_params[[13]] +
      (Y2_sample[[i]]- X %*% t(beta2_sample[[i+1]])) %*% t(Y2_sample[[i]]- X %*% t(beta2_sample[[i+1]]))
    sigma_2 <- 1/rgamma(n = 1 ,v_H2/2, sigma_H2/2)
    H2_sample[[i+1]] <- sigma_2 * diag(length(Y2_sample[[i]]))
    
    #H3
    v_H3 <- (prior_params[[14]] + length(Y3_sample[[i]]))
    sigma_H3 <- prior_params[[14]] * prior_params[[15]] +
      (Y3_sample[[i]] - X %*% t(beta3_sample[[i+1]])) %*% t(Y3_sample[[i]] - X %*% t(beta3_sample[[i+1]]))
    sigma_3 <- 1/rgamma(n = 1 ,v_H3/2, sigma_H3/2)
    H3_sample[[i+1]] <- sigma_3 * diag(length(Y3_sample[[i]]))
    
    #Y2
    Y2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean =  X %*% t(beta2_sample[[i+1]]),
                                         sigma = H2_sample[[i+1]]) |>  t()
    
    #Y3 
    Y3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean =  X %*% t(beta3_sample[[i+1]]),
                                         sigma = H3_sample[[i+1]]) |>  t()
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
  output$Y2 <- Y3_sample
  output$Y3 <- Y3_sample
  
  return(output)
}

####R0_gibbs_sampler####

R0_gibbs_sampler <- function(data, sample_size, prior_params){
  set.seed(428)
  R3_data <- data|> dplyr::filter(R == 3)
  R2_data <- data|> dplyr::filter(R == 2)
  R1_data <- data|> dplyr::filter(R == 1)
  
  data <- data |> dplyr::filter(R == 1)
  n <- length(data$X1)
  X <- data |> dplyr::select(Treatment, X1:X5) |> dplyr::mutate("constant" = rep(1, n)) |> 
    dplyr::select(constant, Treatment, X1:X5) |> as.matrix()
  
  #initial values set
  beta2 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  beta3 <- mvtnorm::rmvnorm(n = 1, mean = rep(0, 7), sigma = 5*diag(7)) |> as.matrix() 
  H1 <- 10*diag(n) |> as.matrix()
  H2 <- 10*diag(n) |> as.matrix()
  H3 <- 10*diag(n) |> as.matrix()
  t <- diag(7) |> as.matrix()
  G <- solve(monomvn::rwish(7, diag(7))) |> as.matrix()
  Y1 <- rnorm(n = n, mean = mean(R1_data$Y1), sd = sqrt(var(R1_data$Y1)))
  Y2 <- mvtnorm::rmvnorm(n = 1, mean = X %*% t(beta2), sigma = H2) |> t()
  Y3 <- mvtnorm::rmvnorm(n = 1, mean = X %*% t(beta3), sigma = H3) |> t()
  
  
  initial_value <- list(beta2, beta3, H1, H2, H3, t, G, Y1, Y2, Y3)
  
  #gibbs sample
  output <- R0_sampler(X, initial_value, prior_params, sample_size)
  return(output)
}

R0_sampler <- function(X, initial_value, prior_params, sample_size){
  beta1_sample <- list()
  beta2_sample <- list() |> append(initial_value[1])
  beta3_sample <- list() |> append(initial_value[2])
  H1_sample <- list() |> append(initial_value[3])
  H2_sample <- list() |> append(initial_value[4])
  H3_sample <- list() |> append(initial_value[5])
  T_sample <- list() |> append(initial_value[6])
  G_sample <- list() |> append(initial_value[7])
  Y1_sample <- list()|> append(initial_value[8])
  Y2_sample <- list()|> append(initial_value[9])
  Y3_sample <- list()|> append(initial_value[10])
  
  for(i in 1:sample_size){
    #beta1 
    Sigma_1 <- solve(t(X) %*% solve(H1_sample[[i]]) %*% X + solve(prior_params[[2]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_1 <- Sigma_1 %*% (t(X) %*% solve(H1_sample[[i]]) %*% Y1_sample[[i]]
                         + solve(prior_params[[2]]) %*% prior_params[[1]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta2_sample[[i]]))
    
    beta1_sample[[i]] <- mvtnorm::rmvnorm(n = 1, mean = mu_1, sigma = Sigma_1) 
    
    #beta2
    Sigma_2 <- solve(t(X) %*% solve(H2_sample[[i]]) %*% X + solve(prior_params[[4]]) + solve(G_sample[[i]])
                     + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% T_sample[[i]])
    mu_2 <- Sigma_2 %*% (t(X) %*% solve(H2_sample[[i]]) %*% Y2_sample[[i]] +
                           solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta1_sample[[i]])
                         + solve(prior_params[[4]]) %*% prior_params[[3]]
                         + t(T_sample[[i]]) %*% solve(G_sample[[i]])%*% t(beta3_sample[[i]]))
    
    beta2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_2, sigma = Sigma_2) 
    
    #beta3
    Sigma_3 <- solve(t(X) %*% solve(H3_sample[[i]]) %*% X + solve(prior_params[[6]]) + solve(G_sample[[i]]))
    mu_3 <- Sigma_3 %*% (t(X) %*% solve(H3_sample[[i]]) %*% Y3_sample[[i]] +
                           + solve(prior_params[[6]]) %*% prior_params[[5]]
                         + solve(G_sample[[i]]) %*% T_sample[[i]] %*% t(beta2_sample[[i+1]]))
    
    beta3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean = mu_3, sigma = Sigma_3)
    
    #T
    t_sample <-  matrix(0, 7, 7)
    for (t in 1:7) {
      sigma_t <- 1/(1/prior_params[[7]] + 
                      ((beta1_sample[[i]][t]))^2/G_sample[[i]][t,t] +
                      (beta2_sample[[i+1]][t])^2/G_sample[[i]][t,t])
      mu_t <- sigma_t * (beta1_sample[[i]][t]*beta2_sample[[i+1]][t]/G_sample[[i]][t,t] +
                           beta2_sample[[i+1]][t]*beta3_sample[[i+1]][t]/G_sample[[i]][t,t])
      t_sample[t,t] <- rnorm(n = 1, mu_t, sigma_t) |> as.matrix()
    }
    T_sample[[i+1]] <- diag(7)
    
    #G
    v_G <- prior_params[[8]] + 2
    S_G <- prior_params[[9]] + (t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]]))%*%
      t(t(beta2_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta1_sample[[i]])) + 
      (t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]]))%*%
      t(t(beta3_sample[[i+1]]) - T_sample[[i+1]] %*% t(beta2_sample[[i+1]])) 
    
    G_sample[[i+1]] <- solve(monomvn::rwish(v_G, solve(S_G))) |> as.matrix()
    
    #H1
    v_H1 <- (prior_params[[10]] + length(Y1_sample[[i]]))
    sigma_H1 <- prior_params[[10]] * prior_params[[11]] +
      (Y1_sample[[i]]- X %*% t(beta1_sample[[i]])) %*% t(Y1_sample[[i]]- X %*% t(beta1_sample[[i]]))
    sigma_1 <- 1/rgamma(n = 1 ,v_H1/2, sigma_H1/2)
    H1_sample[[i+1]] <- sigma_1 * diag(length(Y1_sample[[i]])) 
    
    #H2
    v_H2 <- (prior_params[[12]] + length(Y2_sample[[i]]))
    sigma_H2 <- prior_params[[12]] * prior_params[[13]] +
      (Y2_sample[[i]]- X %*% t(beta2_sample[[i+1]])) %*% t(Y2_sample[[i]]- X %*% t(beta2_sample[[i+1]]))
    sigma_2 <- 1/rgamma(n = 1 ,v_H2/2, sigma_H2/2)
    H2_sample[[i+1]] <- sigma_2 * diag(length(Y2_sample[[i]]))
    
    #H3
    v_H3 <- (prior_params[[14]] + length(Y3_sample[[i]]))
    sigma_H3 <- prior_params[[14]] * prior_params[[15]] +
      (Y3_sample[[i]] - X %*% t(beta3_sample[[i+1]])) %*% t(Y3_sample[[i]] - X %*% t(beta3_sample[[i+1]]))
    sigma_3 <- 1/rgamma(n = 1 ,v_H3/2, sigma_H3/2)
    H3_sample[[i+1]] <- sigma_3 * diag(length(Y3_sample[[i]]))
    
    #Y1
    Y1_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean =  X %*% t(beta1_sample[[i]]),
                                         sigma = H1_sample[[i+1]]) |>  t()
    
    #Y2
    Y2_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean =  X %*% t(beta2_sample[[i+1]]),
                                         sigma = H2_sample[[i+1]]) |>  t()
    
    #Y3 
    Y3_sample[[i+1]] <- mvtnorm::rmvnorm(n = 1, mean =  X %*% t(beta3_sample[[i+1]]),
                                         sigma = H3_sample[[i+1]]) |>  t()
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
  output$Y1 <- Y1_sample
  output$Y2 <- Y2_sample
  output$Y3 <- Y3_sample
  
  return(output)
}
