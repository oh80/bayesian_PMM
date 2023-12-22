main <- function(){
  #read data
  path <- here::here("Liner_Regression", "01_data", "data", "row_data.obj")
  data <- readRDS(path)
  
  #get loglikelihood
  missing_type <- "NMAR"
  
  bayesian_PMM_loglikelihood <- get_bayesian_PMM_loglikelihood(data, missing_type)
  CCMV_loglikelihood <- get_CCMV_loglikelihood(data, missing_type)
  
  #get AIC and BIC
  bayesian_PMM_AIC_BIC <- get_criteria(data, bayesian_PMM_loglikelihood)
  CCMV_AIC_BIC <- get_criteria(data, CCMV_loglikelihood)
  
  #get table
  criteria_df <- cobine_result(bayesian_PMM_loglikelihood, CCMV_loglikelihood,
                               bayesian_PMM_AIC_BIC,  CCMV_AIC_BIC)
  table <- get_table(criteria_df)
  
  #save
  save(table,  missing_type)
}


get_bayesian_PMM_loglikelihood <- function(data, missing_type){
  code_path1 <- here::here("Liner_Regression","03_analyze","code","Bayesian_PMM","compute_complete_data_params.R")
  code_path2 <-  here::here("Liner_Regression","04_report","code","bayesian_PMM","get_estimated_values.R")
  source(code_path1)
  source(code_path2)
  
  data_path1 <- here::here("Liner_Regression","02_build","data", paste0(missing_type, "_0.1.obj"))
  data_path2 <- here::here("Liner_Regression","02_build","data", paste0(missing_type, "_0.25.obj"))
  data_path3 <- here::here("Liner_Regression","02_build","data", paste0(missing_type, "_0.5.obj"))
  data1 <- readRDS(data_path1)
  data2 <- readRDS(data_path2)
  data3 <- readRDS(data_path3)
  
  N <- length(data$Y)
  
  estimated_values1 <- get_estimated_values(missing_rate = 0.1, missing_type)
  mean1 <- rep(0, N)
  for(i in 1:N){
    mean1[i] <- sum(estimated_values1 * data[i,1:10])
  }
  sigma1 <- get_sigma(data1)
  log_likelihood1 <- 0
  for (i in 1:N) {
    log_likelihood1 <- log_likelihood1 + log(dnorm(data$Y[i], mean1[i], sigma1))
  }
  
  estimated_values2 <- get_estimated_values(missing_rate = 0.25, missing_type)
  mean2 <- rep(0, N)
  for(i in 1:N){
    mean2[i] <- sum(estimated_values2 * data[i,1:10])
  }
  sigma2 <- get_sigma(data2)
  log_likelihood2 <- 0
  for (i in 1:N) {
    log_likelihood2 <- log_likelihood2 + log(dnorm(data$Y[i], mean2[i], sigma2))
  }
  
  estimated_values3 <- get_estimated_values(missing_rate = 0.5, missing_type)
  mean3 <- rep(0, N)
  for(i in 1:N){
    mean3[i] <- sum(estimated_values3 * data[i,1:10])
  }
  sigma3 <- get_sigma(data3)
  log_likelihood3 <- 0
  for (i in 1:N) {
    log_likelihood3 <- log_likelihood3 + log(dnorm(data$Y[i], mean3[i], sigma3))
  }
  
  output <- c(log_likelihood1, log_likelihood2, log_likelihood3)
  return(output)
}


get_CCMV_loglikelihood <- function(data, missing_type){
  result_path1 <- here::here("Liner_Regression","03_analyze","output","CCMV", paste0(missing_type, "_0.1.obj"))
  result_path2 <- here::here("Liner_Regression","03_analyze","output","CCMV", paste0(missing_type, "_0.25.obj"))
  result_path3 <- here::here("Liner_Regression","03_analyze","output","CCMV", paste0(missing_type, "_0.5.obj"))
  result1 <- readRDS(result_path1)
  result2 <- readRDS(result_path2)
  result3 <- readRDS(result_path3)
  
  data_path1 <- here::here("Liner_Regression","02_build","data", paste0(missing_type, "_0.1.obj"))
  data_path2 <- here::here("Liner_Regression","02_build","data", paste0(missing_type, "_0.25.obj"))
  data_path3 <- here::here("Liner_Regression","02_build","data", paste0(missing_type, "_0.5.obj"))
  data1 <- readRDS(data_path1)
  data2 <- readRDS(data_path2)
  data3 <- readRDS(data_path3)
  
  code_path <- here::here("Liner_Regression","03_analyze","code","Bayesian_PMM","compute_complete_data_params.R")
  source(code_path)
  
  N <- length(data$Y)
  
  estimated_values1 <- result1$estimated_values[2:11]
  mean1 <- rep(0, N)
  for(i in 1:N){
    mean1[i] <- sum(estimated_values1 * data[i,1:10])
  }
  sigma1 <- get_sigma(data1)
  log_likelihood1 <- 0
  for (i in 1:N) {
    log_likelihood1 <- log_likelihood1 + log(dnorm(data$Y[i], mean1[i], sigma1))
  }
  
  estimated_values2 <- result2$estimated_values[2:11]
  mean2 <- rep(0, N)
  for(i in 1:N){
    mean2[i] <- sum(estimated_values2 * data[i,1:10])
  }
  sigma2 <- get_sigma(data2)
  log_likelihood2 <- 0
  for (i in 1:N) {
    log_likelihood2 <- log_likelihood2 + log(dnorm(data$Y[i], mean2[i], sigma2))
  }
  
  estimated_values3 <- result3$estimated_values[2:11]
  mean3 <- rep(0, N)
  for(i in 1:N){
    mean3[i] <- sum(estimated_values3 * data[i,1:10])
  }
  sigma3 <- get_sigma(data3)
  log_likelihood3 <- 0
  for (i in 1:N) {
    log_likelihood3 <- log_likelihood3 + log(dnorm(data$Y[i], mean3[i], sigma3))
  }
  
  output <- c(log_likelihood1, log_likelihood2, log_likelihood3)
  return(output)
}


get_criteria <- function(data, loglikelihood){
  N <- length(data$Y)
  k <- 100 + 10
  
  AIC <- -2*loglikelihood + 2*k
  BIC <- -2*loglikelihood + k*log(N)
  
  output <- c(AIC, BIC)
  return(output)
}


cobine_result <- function(loglikelihood1, loglikelihood2, criteria1, criteria2){
  method <- c("proposed", " ", " ", "CCMV", " ", " ")
  missing_rate <- c(rep(c("0.1", "0.25", "0.5"), 2))
  log_likelihood <- c(loglikelihood1, loglikelihood2)
  AIC <- c(criteria1[1:3], criteria2[1:3])
  BIC <- c(criteria1[4:6], criteria2[4:6])
  
  output <- data.frame("method" = method,
                       "rate" = missing_rate,
                       "log_likelihood" = log_likelihood,
                       "AIC" = AIC,
                       "BIC" = BIC) |> 
    dplyr::mutate(log_likelihood = format(log_likelihood, digits = 2),
                  AIC = format(AIC, digits = 2), 
                  BIC = format(BIC, digits = 2))
  return(output)
}


get_table <- function(data_frame){
  output <- data_frame |> kableExtra::kable(booktabs = TRUE, format = "latex")
  
  return(output)
}


save <- function(table,  missing_type){
  file_name <- paste0(missing_type,"_criteria.txt")
  path <- here::here("Liner_Regression","04_report","output","common", file_name)
  writeLines(table, path)
}


main()