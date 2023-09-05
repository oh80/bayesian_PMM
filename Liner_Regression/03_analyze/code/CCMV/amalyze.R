main <- function(){
  #read_data
  missing_rate <- 0.1
  missing_type <- "MCAR"
  file_name <- paste0(missing_type , "_", missing_rate, ".obj")
  path <- here::here("Liner_regression", "02_build","data", file_name)
  data <- readRDS(path)
  
  #param set
  set.seed(428)
  D = 10
  
  #get imputation model
  coefficients <- get_coeffcients(data)
  sigma <- get_sigma(data)
  
  #imputation 
  pseudo_complete_data <- multiple_imputation(data, coefficients, sigma, D)
  
  #estimation
  estimated_values <- esimate_each_data(pseudo_complete_data)
  
  #combine each estimated values
  estimate_result <- combine(estimated_values)
  
  return(estimate_result)
}


get_coeffcients <- function(data){
  complete_data <- data |> dplyr::filter(R == 0) |> dplyr::select(-R)
  
  param_list <- list()
  for (j in 1:10) {
    imputed_col <- paste0("x",j)
    formula <- as.formula(paste0(imputed_col, "~ ."))
    lm <- lm(formula = formula, complete_data)
    param_list[[j]] <- coef(lm)
  }
  output <- param_list
  return(output)
}


get_sigma <- function(data){
  complete_data <- data |> dplyr::filter(R == 0) |> dplyr::select(-R)
  sd_list <- list()
  for (j in 1:10) {
    sd_list[j] <- sd(complete_data[,j])
  }
  output <- sd_list
  return(output)
}


get_single_imputation <- function(data, coefficients, sigma){
  output <- data |> dplyr::filter(R == 0)
  for(j in 1:10){
    missing_col <- paste0("x",j)
    missing_data <- data |> dplyr::filter(R == j) |> dplyr::select(-all_of(missing_col))
    n <- length(missing_data$Y)
    coeffcient <- coefficients[[j]]
    imputation <- rep(0,n)
      for (i in 1:n) {
        predict <- sum(coeffcient[1] + coeffcient[2:11] * missing_data[i, 2:11])
        imputation[i] <- rnorm(n = 1, mean =predict, sd = unlist(sigma[j]) )
      }
    missing_data <- missing_data |> dplyr::mutate("{missing_col}" := scale(imputation))
    output <- output |> dplyr::bind_rows(missing_data)
  }
  return(output)
}


multiple_imputation <- function(data, coefficients, sigma, D){
  output <- list()
  for (d in 1:D) {
    set.seed(d)
    output[[d]] <- get_single_imputation(data, coefficients, sigma)
  }
  return(output)
}


esimate_each_data <- function(data_list){
  D <- length(data_list)
  output <- list()
  for (d in 1:D) {
    data <- data_list[[d]] |> dplyr::select(-R)
    output[[d]] <- estimatr::lm_robust(formula = Y ~ ., data = data,
                                       se_type = "stata")
  }
  return(output)
}


combine <- function(estimate_list){
  D = length(estimate_list)
  
  combined_estimation <- rep(0,11)
  for (d in 1:D) {
    combined_estimation <- combined_estimation + estimate_list[[d]]$coefficients/D
  }
  
  var <- rep(0,11)
  for (d in 1:D) {
    var <- var + 1/D * (estimate_list[[d]]$std.error)^2 + 
      (1+1/D)*(1/(D-1))*(combined_estimation-estimate_list[[d]]$coefficients)^2
  }
  
  output <- list()
  output$estimation <- combined_estimation
  output$var <- var
  return(output)
}


results <- main()

