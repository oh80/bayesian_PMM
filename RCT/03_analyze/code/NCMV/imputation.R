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
    model  <- lm(data = R2_data, formula = formula)
    
    for(i in 1:n){
      coefs <- model$coefficients |> as.vector()
      covariables <- missing_data[i,] |> dplyr::select(X1,X2,X3,X4,X5, Y0, Y1, Treatment)
      sd <- sqrt(var(R2_data$Y2))
      mean <- (coefs[1]  + sum(coefs[2:9]* covariables))
      missing_data[i,]$Y2 <- rnorm(n = 1, mean = mean + delta, sd = sd)
    }
    return(missing_data)
  }
  
  if(imputed_pattern == 0){
    missing_data <- R0_data
    n <- length(missing_data$X1)
    formula <- "Y1 ~ X1 + X2 + X3 + X4 + X5 + Y0 + Treatment "
    model  <- lm(data = R1_data, formula = formula)
    
    for(i in 1:n){
      coefs <- model$coefficients |> as.vector()
      covariables <- missing_data[i,] |> dplyr::select(X1,X2,X3,X4,X5,Y0, Treatment)
      sd <- sqrt(var(R1_data$Y1))
      mean <- (coefs[1]  + sum(coefs[2:8]* covariables))
      missing_data[i,]$Y1 <- rnorm(n = 1, mean = mean + delta, sd = sd)
    }
    return(missing_data)
  }
}


single_imputation <- function(data, delta){
  pattern <- data$R
  n <- length(data$X1)
  
  complete_data <- data |> dplyr::filter(R == 3)
  R2_data <- data |> dplyr::filter(R == 2)
  
  #one time imputation
  R3_imputed_data <- impute_data(data, imputed_pattern = 2, delta = delta)
  R2_impute_data <- impute_data(data, imputed_pattern = 1, delta = delta) |> 
    dplyr::mutate(R = R + 1)
  R1_impute_data <- impute_data(data, imputed_pattern = 0, delta = delta) |> 
    dplyr::mutate(R = R + 1)
  
  #two times imputation
  R3R2_imputed_data <- impute_data(data = dplyr::bind_rows(complete_data, R2_impute_data),
                                   imputed_pattern = 2, delta = 0) |> 
    dplyr::mutate(R = R - 1)
  
  R1R2_imputed_data <- impute_data(data = dplyr::bind_rows(complete_data, R2_data, R1_impute_data),
                                   imputed_pattern = 1, delta = 0) |> 
    dplyr::mutate(R = R + 1)
  
  #three times imputation
  R1R2R3_imputed_data <- impute_data(data = dplyr::bind_rows(complete_data, R1R2_imputed_data),
                                   imputed_pattern = 2, delta = 0) |> 
    dplyr::mutate(R = R - 2)
  
  output <- dplyr::bind_rows(complete_data, R3_imputed_data, R3R2_imputed_data, R1R2R3_imputed_data)
  
  
  return(output)
}




