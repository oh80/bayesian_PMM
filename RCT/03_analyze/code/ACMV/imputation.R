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

a <- decison_disutribution(data, imputed_pattern = 0)

