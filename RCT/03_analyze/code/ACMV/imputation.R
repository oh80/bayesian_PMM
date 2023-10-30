decison_disutribution <- function(data, imputed_pattern){
  R3_data <- data |> dplyr::filter(R == 3)
  R2_data <- data |> dplyr::filter(R == 2)
  R1_data <- data |> dplyr::filter(R == 1)
  R0_data <- data |> dplyr::filter(R == 0)
  
  if(imputed_pattern == 1){
    missing_data <- R1_data
    n <- length(R1_data$X1)
    
    alpha3 <- length(R3_data$X1)
    alpha2 <- length(R2_data$X1)
    
    mean3 <- c(mean(R3_data$Y0), mean(R3_data$Y1))
    mean2 <- c(mean(R2_data$Y0), mean(R2_data$Y1))
    
    cov3 <- matrix(c(var(R3_data$Y0),cov(R3_data$Y1, R3_data$Y0),
                    cov(R3_data$Y0, R3_data$Y1), var(R3_data$Y1)),
                  nrow = 2)
    cov2 <- matrix(c(var(R2_data$Y0),cov(R2_data$Y0, R2_data$Y1),
                    cov(R2_data$Y0, R2_data$Y1), var(R2_data$Y1)),
                  nrow = 2)
    
    output <- rep(0, n)
    for (i in 1:n) {
      set.seed(i)
      vec <- c(missing_data[i,]$Y0, missing_data[i,]$Y1)
      omega <- (alpha3 * mvtnorm::dmvnorm(vec, mean = mean3, sigma = cov3)) /
        ((alpha3 * mvtnorm::dmvnorm(vec, mean = mean3, sigma = cov3)) + 
        (alpha2 * mvtnorm::dmvnorm(vec, mean = mean2, sigma = cov2)))
      
      u <- runif(n = 1, min = 0, max = 1)
      if(u <= omega){
        output[i] <- 3
      }else{
        output[i] <- 2
      }}
    return(output)}
  
  if(imputed_pattern == 0){
    missing_data <- R0_data
    n <- length(R0_data$X1)
    
    alpha3 <- length(R3_data$X1)
    alpha2 <- length(R2_data$X1)
    alpha1 <- length(R1_data$X1)
    
    mean3 <- mean(R3_data$Y0)
    mean2 <- mean(R2_data$Y0)
    mean1 <- mean(R1_data$Y0)
    
    sd3 <- sd(R3_data$Y0)
    sd2 <- sd(R2_data$Y0)
    sd1 <- sd(R1_data$Y0)
    
    output <- rep(0, n)
    for (i in 1:n) {
      set.seed(i)
      prob3 <- dnorm(missing_data[i,]$Y0, mean3, sd3)
      prob2 <- dnorm(missing_data[i,]$Y0, mean2, sd2)
      prob1 <- dnorm(missing_data[i,]$Y0, mean1, sd1)
      
      omega3 <- (alpha3 * prob3) / (alpha3 * prob3 + alpha2 * prob2 + alpha1 * prob1)
      omega2 <- (alpha2 * prob2) / (alpha3 * prob3 + alpha2 * prob2 + alpha1 * prob1)
      omega1 <- (alpha1 * prob1) / (alpha3 * prob3 + alpha2 * prob2 + alpha1 * prob1)
      u <- runif(n = 1, min = 0, max = 1)
      if(u <= omega1){
        output[i] <- 1
      }
      if(omega1 <  u &  u < (omega2 + omega1)){
        output[i] <- 2
      }
      else{
        output[i] <- 3
      }}
    return(output)}
  
}


impute_data <- function(data, imputed_pattern, delta){
  R3_data <- data |> dplyr::filter(R == 3)
  R2_data <- data |> dplyr::filter(R == 2)
  R1_data <- data |> dplyr::filter(R == 1)
  R0_data <- data |> dplyr::filter(R == 0)
  
  if(imputed_pattern == 2){
    missing_data <- R2_data
    n <- length(missing_data$X1)
    formula <- "Y3 ~ X1 + X2 + X3 + X4 + X5 + Y0 + Y1 + Y2 + Treatment "
    model  <- lm(data = R3_data, formula = formula)
    
    for(i in 1:n){
      coefs <- model$coefficients |> as.vector()
      covariables <- missing_data[i,] |> dplyr::select(X1,X2,X3,X4,X5, Y0, Y1, Y2, Treatment)
      sd <- sqrt(var(R3_data$Y3))
      missing_data[i,]$Y3 <- rnorm(n = 1, mean = (coefs[1]  + sum(coefs[2:10]* covariables)) + delta, sd = sd)
    }
    return(missing_data)
    }
    
    if(imputed_pattern == 1){
      missing_data <- R1_data
      n <- length(missing_data$X1)
      formula <- "Y2 ~ X1 + X2 + X3 + X4 + X5 + Y0 + Y1 + Treatment "
      model3  <- lm(data = R2_data, formula = formula)
      model2  <- lm(data = R2_data, formula = formula)
      
      use_model <- decison_disutribution(data, imputed_pattern = 1)
      
      for(i in 1:n){
        if(use_model[i] == 3){
          coefs <- model3$coefficients |> as.vector()
          sd <- sqrt(var(R3_data$Y2))
        }else{
          coefs <- model2$coefficients |> as.vector()
          sd <- sqrt(var(R2_data$Y2))
        }
        covariables <- missing_data[i,] |> dplyr::select(X1,X2,X3,X4,X5, Y0, Y1, Treatment)
        mean <- (coefs[1]  + sum(coefs[2:9]* covariables))
        missing_data[i,]$Y2 <- rnorm(n = 1, mean = mean + delta, sd = sd)
      }
      return(missing_data)
    }
    
    if(imputed_pattern == 0){
      missing_data <- R0_data
      n <- length(missing_data$X1)
      formula <- "Y1 ~ X1 + X2 + X3 + X4 + X5 + Y0  + Treatment "
      model3  <- lm(data = R2_data, formula = formula)
      model2  <- lm(data = R2_data, formula = formula)
      model1  <- lm(data = R1_data, formula = formula)
      
      use_model <- decison_disutribution(data, imputed_pattern = 0)
      
      for(i in 1:n){
        if(use_model[i] == 3){
          coefs <- model3$coefficients |> as.vector()
          sd <- sqrt(var(R3_data$Y1))
        }
        if(use_model[i] == 2){
          coefs <- model2$coefficients |> as.vector()
          sd <- sqrt(var(R2_data$Y1))
        }
        else{
          coefs <- model1$coefficients |> as.vector()
          sd <- sqrt(var(R1_data$Y1))
        }
        covariables <- missing_data[i,] |> dplyr::select(X1,X2,X3,X4,X5, Y0, Treatment)
        mean <- (coefs[1]  + sum(coefs[2:8]* covariables))
        missing_data[i,]$Y1 <- rnorm(n = 1, mean = mean + delta, sd = sd)
      }
      return(missing_data)
    }
  }


i <- impute_data(data, 0, 0)
